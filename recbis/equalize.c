/*
   This file is part of RECBIS:
   A library for parallel computation of medians and for
   parallel data distribution using recursive bisection.

   Written by Tomas Oppelstrup at the Lawrence Livermore
   National Laboratory.
*/

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <mpi.h>

#include "integer.h"
#include "genptr.h"

#include "exchange.h"
#include "txmap.h"

#include "rectimer.h"
#include "functimer.h"


#include <stdarg.h>
int dbprintf(char *s,...) {
  return 0;
  int pid,n = 0;
  MPI_Comm_rank(MPI_COMM_WORLD,&pid);
  if(pid >= 0) {
    va_list args;
    va_start(args,s);
    n = vprintf(s,args);
    va_end(args);
  }
  return n;
}

integer equalize_sched(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
		       integer nloc,txmap_struct *txmap) {
  STARTTIMER;
  const integer np = (pid1-pid0)/pstride;
  const integer g_pid = (pid-pid0)/pstride;
  const integer pcut = (np+1)/2;

  integer nleft,nright,ntot;
  integer ngleft,ngright;

  assert(np > 0);
  if(np == 1) { STOPTIMER; return nloc; }

  /* Recursive equalization over left and right half of processor set */
  if(g_pid < pcut) {
    integer i0,i1;
    nleft  = equalize_sched(pid0,pid0+pcut*pstride,pid,pstride,comm,
			    nloc,txmap);
    partition(nleft,pcut,g_pid,&i0,&i1);
    nloc = i1 - i0;
  } else {
    integer i0,i1;
    nright = equalize_sched(pid0+pcut*pstride,pid1,pid,pstride,comm,
			    nloc,txmap);
    partition(nright,np-pcut,g_pid-pcut,&i0,&i1);
    nloc = i1 - i0;
  }

  /* Receive atom count from other half of processor set */ {
    MPI_Request reqlist[4];
    int nreq = 0;

    if(g_pid < pcut) {
      integer p;
      MPI_Irecv(&nright,sizeof(integer),MPI_BYTE,
		pid0 + (pcut+g_pid%(np-pcut))*pstride,717,comm,reqlist + nreq++);
      for(p = pcut+g_pid; p<np; p+=pcut)
	MPI_Isend(&nleft,sizeof(integer),MPI_BYTE,
		  pid0+p*pstride,717,comm,reqlist + nreq++);
    } else {
      integer p;
      MPI_Irecv(&nleft,sizeof(integer),MPI_BYTE,
	       pid0 + ((g_pid-pcut)%pcut)*pstride,717,comm,reqlist + nreq++);
      for(p = g_pid-pcut; p<pcut; p+=(np-pcut))
	MPI_Isend(&nright,sizeof(integer),MPI_BYTE,
		 pid0+p*pstride,717,comm,reqlist + nreq++);
    }
    assert(nreq <= (int) (sizeof(reqlist)/sizeof(*reqlist)));
    MPI_Waitall(nreq,reqlist,MPI_STATUSES_IGNORE);
  }
  
  ntot = nleft + nright;

  /* Calculate left-right partition */ {
    integer dummy;
    partition(ntot,np,pcut,&ngleft,&dummy);
    ngright = ntot - ngleft;
  }
  
  /* Exchange items with other half of processor set */
  if(g_pid < pcut) {
    exchange_sched(pid0,pid0+pcut*pstride,pid,pstride,
		   pid0+pcut*pstride,pid1,pstride,comm,
		   nleft,ngleft,nright,ngright,
		   txmap);
  } else {
    exchange_sched(pid0+pcut*pstride,pid1,pid,pstride,
		   pid0,pid0+pcut*pstride,pstride,comm,
		   nright,ngright,nleft,ngleft,
		   txmap);
  }

#if 0
  /* Calculate number of local elements...*/ {
    integer i0,i1;
    partition(ntot,np,g_pid,&i0,&i1);
    /* This can be deduced rom ntot anyway */
    nloc_p[0] = i1 - i0;
  }
#endif
  STOPTIMER;
  return ntot;
}

extern int wrapcond[4];
int wrapcond[4] = {0,0,0,0};

integer equalize_new(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
		 integer nloc_p[1],integer nlocmax,integer sz,void *data) {
  STARTTIMER;
  const integer np = (pid1-pid0)/pstride;
  const integer g_pid = (pid-pid0)/pstride;
  txmap_struct txmap = { 0 , 0 , NULL };
  integer nloc = nloc_p[0],ntot;

  const int globbarr = 0;
  static thandle tcomm = NULL;


  if(globbarr) {
    static thandle tbarr = NULL;
    tbarr = rectimer_start(tbarr,"prebarrier");
    MPI_Barrier(comm);
    rectimer_stop(tbarr);
  }

  ntot = equalize_sched(pid0,pid1,pid,pstride,comm,
		 nloc,&txmap);

  if(globbarr) {
    static thandle tbarr = NULL;
    tbarr = rectimer_start(tbarr,"midbarrier");
    MPI_Barrier(comm);
    rectimer_stop(tbarr);
  }


  tcomm = rectimer_start(tcomm,"exchange");
  /* Execute communication schedule */ {
    integer n = txmap.nused,ns = 0,nr = 0;
    txmap_entry
      *sends = (txmap_entry *) malloc(sizeof(*sends) * n),
      *recvs = txmap.list;
    
    /* Separate communication items into send and receive lists */ {
      integer scount = 0,rcount = 0,i;
      for(i = 0; i<n; i++) {
	if(recvs[i].count < 0) {
	  integer c = -recvs[i].count;
	  recvs[nr  ].count = c;
	  recvs[nr++].pid   = recvs[i].pid;
	  rcount += c;
	} else if(recvs[i].count > 0) {
	  integer c = recvs[i].count;
	  sends[ns++] = recvs[i];
	  scount += c;
	}
      }
      assert(nloc + rcount - scount >= 0);
    }


    /* Submit send and receive calls */ {
      MPI_Request sreq[ns],rreq[nr];
      integer
	sptrlist[ns],rptrlist[nr];
      integer
	rptr = nloc,
	rptr_end = nlocmax,
	
	sptr = 0,
	sptr_end = nloc,
	
	wrap_flag = 0,
	recv_wrap = -1,
	send_wrap = -1,
	
	ir = 0,ir0 = 0,
	is = 0,is0 = 0;
      while(ir < nr || is < ns) {
	for(; ir<nr; ir++) {
	  
	  if((rptr + recvs[ir].count)/nlocmax > (rptr)/nlocmax) {
	    dbprintf("pid=%2lld: Receive wrap enter: nlocmax = %lld, ir=%d, nr=%d, rptr = %lld(%lld), count = %lld, "
		     "sw=%lld,rw=%lld,wf=%lld\n",
		     pid,nlocmax,ir,nr,
		     rptr,rptr%nlocmax,recvs[ir].count,send_wrap,recv_wrap,wrap_flag);
	    /* Wrap around. Advance ptr to first beginning of buffer */
	    assert(wrap_flag == 0);
	    assert(send_wrap < 0);
	    assert(recv_wrap < 0);
	    recv_wrap = rptr;
	    
	    /* Figure out where sends will wrap (if they do)... */ {
	      integer i,ipos;
	      send_wrap = recv_wrap;
	      ipos = sptr;
	      for(i = is; i<ns; i++) {
		if(ipos + sends[i].count > rptr) {
		  send_wrap = ipos;
		  break;
		} else {
		  ipos += sends[i].count;
		}
	      }
	    }

	    dbprintf("pid=%2lld: Receive wrap ok: nlocmax = %lld, ir=%d, nr=%d, rptr = %lld(%lld), count = %lld, "
		     "sw=%lld,rw=%lld,wf=%lld\n",
		     pid,nlocmax,ir,nr,
		     rptr,recvs[ir].count,send_wrap,recv_wrap,wrap_flag);


	    rptr += nlocmax - send_wrap%nlocmax;
	  }
	  if(rptr + recvs[ir].count > rptr_end) break;
	  MPI_Irecv(idx_ptr(data,rptr % nlocmax,sz),recvs[ir].count * sz,MPI_BYTE,
		    recvs[ir].pid,938,comm,rreq + ir);
	  dbprintf("pid=%2lld: Posting receive %lld of %lld in range [%lld,%lld)\n",
		   pid,ir,nr,rptr,rptr+recvs[ir].count);
	  rptr += recvs[ir].count;
	  rptrlist[ir] = rptr;
	}
	
	for(; is<ns; is++) {
	  if(sptr <= recv_wrap && (sptr + sends[is].count) > recv_wrap) {
	    /* Straddling last receive, need to copy end block */
	    assert(send_wrap == sptr);
	    sptr += nlocmax - (sptr % nlocmax);
	    wrap_flag = 1;
	  }
	  if(sptr + sends[is].count > sptr_end) break;
	  if(wrap_flag == 1) {
	    static thandle twrap = NULL;
	    twrap = rectimer_start(twrap,"wrapcopy");
	    memmove(data,idx_ptr(data,send_wrap,sz),(recv_wrap-send_wrap)*sz);
	    rectimer_stop(twrap);
	    dbprintf("pid=%2lld: Issuing send wrap move [%lld,%lld) = %lld\n",
		     pid,send_wrap,recv_wrap,recv_wrap-send_wrap);
	    wrap_flag = 0;
	    send_wrap = -1;
	    recv_wrap = -1;
	  }
	  MPI_Isend(idx_ptr(data,sptr % nlocmax,sz),sends[is].count * sz,MPI_BYTE,
		    sends[is].pid,938,comm,sreq + is);

	  dbprintf("pid=%2lld: Posting send %lld of %lld in range [%lld,%lld)\n",
		   pid,is,ns,sptr,sptr+sends[is].count);

	  sptr += sends[is].count;
	  sptrlist[is] = sptr;
	}
	
	/* Wait for a send or receive to complete */
	if(ir0 < ir && is0 < is) {
	  MPI_Request req_r_s[2];
	  int i,ndone,donelist[2];
	  req_r_s[0] = rreq[ir0];
	  req_r_s[1] = sreq[is0];
	  MPI_Waitsome(2,req_r_s,&ndone,donelist,MPI_STATUSES_IGNORE);
	  for(i = 0; i<ndone; i++) {
	    if(donelist[i] == 0) {
	      sptr_end = rptrlist[ir0++];
	      dbprintf("pid=%2lld: Receive %lld complete. sptr_end = %lld\n",
			pid,ir0-1,sptr_end);
	    } else {
	      rptr_end = sptrlist[is0++] + nlocmax;
	      dbprintf("pid=%2lld: Send %lld complete. sptr_end = %lld\n",
			pid,is0-1,rptr_end);	    }
	  }
	} else if(ir0 < ir) {
	  MPI_Wait(rreq + ir0,MPI_STATUS_IGNORE);
	  sptr_end = rptrlist[ir0++];
	  dbprintf("pid=%2lld: Receive %lld complete. sptr_end = %lld\n",
		    pid,ir0-1,sptr_end);
	} else if(is0 < is) {
	  MPI_Wait(sreq + is0,MPI_STATUS_IGNORE);
	  rptr_end = sptrlist[is0++] + nlocmax;
	  dbprintf("pid=%2lld: Send %lld complete. sptr_end = %lld\n",
		    pid,is0-1,rptr_end);
	}
      }

      /* Wait for pending communication to finish */
      MPI_Waitall(ns-is0,sreq + is0,MPI_STATUSES_IGNORE);
      MPI_Waitall(nr-ir0,rreq + ir0,MPI_STATUSES_IGNORE);
      
      /* Solve for nloc after communication */ {
	integer i0,i1;
	partition(ntot,np,g_pid,&i0,&i1);
	nloc = i1 - i0;
	nloc_p[0] = nloc;
      }


      /* Compact buffer... */ {
	/* This essentially copies the entire result buffer.
	   It could be optimized by only copying data from
	   the end of the buffer to the hole between the
	   receive and the send pointer. To do this one must
	   however calculate exactly where the receive just
	   before the last receive buffer wrap ended, since
	   there is no valid data between that point and the
	   end of the buffer. Actually, this can also be
	   calculated from knowing the final number of
	   elements, nloc, which we do know!
	*/
	static thandle tcompact = NULL;


	tcompact = rectimer_start(tcompact,"compact");
	if(0) {
	  if(rptr/nlocmax > sptr/nlocmax) {
	    printf("pid = %06d : buffer wrap.\n",(int) pid);
	    memmove(idx_ptr(data,rptr % nlocmax,sz),
		    idx_ptr(data,sptr % nlocmax,sz),
		    (nlocmax - sptr%nlocmax)*sz);
	  } else {
	    printf("pid = %06d : no buffer wrap.\n",(int) pid);
	    memmove(data,
		    idx_ptr(data,sptr % nlocmax,sz),
		    (rptr - sptr)*sz);
	  }
	} else {
	  if(rptr/nlocmax > sptr/nlocmax) {
	    const integer
	      rmod = rptr % nlocmax,
	      smod = sptr % nlocmax;
	    if(nloc < smod) {
	      wrapcond[0]++;
	      memcpy(idx_ptr(data,rmod,sz),
		     idx_ptr(data,smod,sz),
		     (nloc - rmod)*sz);
	    } else {
	      wrapcond[1]++;
	      memcpy(idx_ptr(data,rmod,sz),
		     idx_ptr(data,nloc,sz),
		     (smod - rmod)*sz);
	    }
	  } else {
	    const integer smod = sptr % nlocmax;
	    if(smod > rptr-sptr) {
	      wrapcond[2]++;
	      memcpy(data,idx_ptr(data,smod,sz),(rptr - sptr)*sz);
	    } else {
	      wrapcond[3]++;
	      memcpy(data,idx_ptr(data,rptr-sptr,sz),smod*sz);
	    }
	  }
	  
	}
	rectimer_stop(tcompact);
      }
    }
  }

  rectimer_stop(tcomm);

  if(globbarr) {
    static thandle tbarr = NULL;
    tbarr = rectimer_start(tbarr,"postbarrier");
    MPI_Barrier(comm);
    rectimer_stop(tbarr);
  }


  STOPTIMER;
  return ntot;
}


integer equalize_old(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
		     integer nloc_p[1],integer nlocmax,integer sz,void *data_p) {
  STARTTIMER;
  const integer np = (pid1-pid0)/pstride;
  const integer g_pid = (pid-pid0)/pstride;
  const integer pcut = (np+1)/2;

  integer nleft,nright,ntot;
  integer ngleft,ngright;

  if(np == 1) { STOPTIMER; return nloc_p[0]; }

  /* Recursive equalization over left and right half of processor set */
  if(g_pid < pcut)
    nleft  = equalize_old(pid0,pid0+pcut*pstride,pid,pstride,comm,
			  nloc_p,nlocmax,sz,data_p);
  else
    nright = equalize_old(pid0+pcut*pstride,pid1,pid,pstride,comm,
			  nloc_p,nlocmax,sz,data_p);

  /* Receive atom count from other half of processor set */ {
    MPI_Status status;
    integer p;

    if(g_pid < pcut) {
      for(p = pcut+g_pid; p<np; p+=pcut)
	MPI_Send(&nleft,sizeof(integer),MPI_BYTE,
		 pid0+p*pstride,717,comm);
      MPI_Recv(&nright,sizeof(integer),MPI_BYTE,
	       pid0 + (pcut+g_pid%(np-pcut))*pstride,717,comm,&status);    
    } else {
      MPI_Recv(&nleft,sizeof(integer),MPI_BYTE,
	       pid0 + ((g_pid-pcut)%pcut)*pstride,717,comm,&status);
      for(p = g_pid-pcut; p<pcut; p+=(np-pcut))
	MPI_Send(&nright,sizeof(integer),MPI_BYTE,
		 pid0+p*pstride,717,comm);
    }
  }
  
  ntot = nleft + nright;

  /* Calculate left-right partition */ {
    integer dummy;
    partition(ntot,np,pcut,&ngleft,&dummy);
    ngright = ntot - ngleft;
  }

  /* Exchange items with other half of processor set */
  if(g_pid < pcut) {
    exchange(pid0,pid0+pcut*pstride,pid,pstride,
	     pid0+pcut*pstride,pid1,pstride,comm,
	     nleft,ngleft,nright,ngright,
	     nlocmax,sz,data_p);
  } else {
    exchange(pid0+pcut*pstride,pid1,pid,pstride,
	     pid0,pid0+pcut*pstride,pstride,comm,
	     nright,ngright,nleft,ngleft,
	     nlocmax,sz,data_p);
  }

  /* Calculate number of local elements...*/ {
    integer i0,i1;
    partition(ntot,np,g_pid,&i0,&i1);
    nloc_p[0] = i1 - i0;
  }

  STOPTIMER;
  return ntot;
}

integer equalize(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
		 integer nloc_p[1],integer nlocmax,integer sz,void *data) {
  STARTTIMER;
  integer retval;
  retval = equalize_new(pid0,pid1,pid,pstride,comm,
			nloc_p,nlocmax,sz,data);
  STOPTIMER;
  return retval;
}
