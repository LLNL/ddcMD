#include "primes.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/resource.h>
#include <assert.h>

#define LONG64 unsigned long long

#define MIN(A,B) ((A) < (B) ? (A) : (B))

static int _initialized = 0;
static LONG64 _blockSize = 0;
static LONG64 _taskId;
static LONG64 _nTasks;
static LONG64 _upperBound = 0;
static LONG64 _prime = 1ll;  
static LONG64 _smallestPrime = (2ll<<30) + 1ll; // Find primes larger than this. 


static unsigned isPrime1(LONG64 n);
static void setSearchInterval(void);

/** Initialize prime generation. */
void prime_init(unsigned blockSize, unsigned taskId, unsigned nTasks)
{
   _blockSize = blockSize;
   _taskId = taskId;
   _nTasks = nTasks;
	_initialized = 1;
}

/** Return next prime.  Setup the next block if we exhaust the assigned
 *  range. */
LONG64 nextPrime(void)
{
	assert(_initialized !=0 );
   do
   {
      _prime +=2;
      if (_prime >= _upperBound)
			setSearchInterval();
   }
   while ( !isPrime1(_prime) ); 
   return _prime; 
}

/**
 *  Sets _upperBound and _prime for the next search interval for this
 *  task.  _prime is an odd number at the bottom of the interval (not
 *  necessarily prime) and _upperBound is an odd number at the top of
 *  the interval.
*/
static void setSearchInterval(void)
{
	static LONG64 iBlock=0;

	_upperBound = (iBlock*_nTasks+_taskId)*_blockSize + _smallestPrime;
	_prime = _upperBound - _blockSize;
	
	if (_upperBound%2 == 0)
		_upperBound -= 1;
	if (_prime%2 == 0)
		_prime += 1;
	++iBlock;
}

static LONG64 llog2(LONG64 a)
{
	LONG64 mask=0x8000000000000000ll;
	for (LONG64 s=64ll;s>0ll;s--)
	{
		if (a&mask) return s-1;
		mask = mask >> 1; 
	}
	return 0xffffffffffffffffll; 
}
static LONG64 MMUL(LONG64 a , LONG64 b, LONG64 N,unsigned n)
{
	LONG64 sum =0.0; 
	for (unsigned i=0;i<n;i++)
	{
		sum += (a & 1)*b; 
		sum +=(sum & 1)*N; 
		sum >>= 1; 
		a >>= 1; 
	}
	if (sum > N ) sum-=N; 
	return sum; 
}
static LONG64 two2n(unsigned n, LONG64 N)
{
/*
	LONG64 p=(1 << n); 
	if (p > N) p-=N; 
*/
	LONG64 p=(1ll<<n); 
	if (p > N) p-=N; 
	for (int s=0;s<(int)n;s++) 
	{
		p = (2*p);
		if (p > N) p-=N; 
	}
	return p; 

}
static LONG64 ipow1(LONG64 base, LONG64 exp,LONG64 N,unsigned n)
{
    LONG64 result = 1ll;
    while (exp)
    {
	
    	if (exp & 1ll) { result=MMUL(result,base,N,n); }
        exp >>= 1ll;
		base=MMUL(base,base,N,n); 
    }
    return result;
}

static unsigned isPrime1(LONG64 N)
{
	LONG64 primes[]={2,3,5,7,11,13,17};
	if (N%3ll==0) return 0; 
	if (N%5ll==0) return 0; 
	if (N%7ll==0) return 0; 
	if (N%11ll==0) return 0; 
	if (N%13ll==0) return 0; 
	unsigned isprime; 
	unsigned n = 1ll + llog2(N); 
	LONG64 s=N-1ll;
	LONG64 p2=  two2n(n,N);  
	unsigned r; 
	for (r=0;s%2==0;r++) s=s>>1;
	for (unsigned i=0;i<7 && primes[i]<(N-1);i++)
	{	
		LONG64 a = primes[i]; 
		isprime=0; 
		LONG64 x = ipow1(a,s,N,n); 
		if (x==1 || x == (N-1)) 
		{
			isprime=1;
			continue;
		}
		for (unsigned j=1;j<r;j++)
		{
			x = MMUL(x,x,N,n); 
			x=MMUL(x,p2,N,n); 
			if (x == 1) return 0;
			if (x ==( N-1)) { isprime=1; break; }
		}
		if (isprime == 0) return 0;
	}
	return 1;
}


#ifdef COMPILE_UNUSED
/** Find an interval starting at "x0" large enough to contain at least n
 *  primes. This is only a crude estimate. It mostly works but there may
 *  be cases that fail.  */
static LONG64  estimateUpperBound(unsigned n, LONG64 x0)
{
   double x1; 
   double l0 = log(x0*1.0); 
   double d = 1.0/l0 - 1.0/(l0*l0);  // Estimate of density of primes at x0; 
   if (n < 1000) n=1000;             // Choose n large enough so that interval is large enough that density makes sense
   x1 = x0 + n/d;                    // First estimate of upperBound to interval. 
   l0 = log(x1*1.0);                  
   d = 1.0/l0 - 1.0/(l0*l0);         // Density estimate at end of interval. 
   x1 = x0 + 1.2*n/d;                // Second estimate of upperBound. 1.2 is a fudge factor. 
   LONG64 upperBound = x1+1.0;  
   if (upperBound%2ll==0) upperBound++;  //Make sure upperBound is odd. 
   return upperBound ; 
}

static LONG64 NRMMS(LONG64 a , LONG64 b, LONG64 N,LONG64 s)
{
	LONG64 sum =0.0; 
	for (LONG64 i=0;i<s;i++)
	{
		sum+= (a & 1)*b; 
		sum+=(sum&1)*N; 
		sum = sum/2; 
		a = a >> 1; 
	}
	return sum; 
}
static LONG64 mult(LONG64 a , LONG64 b, LONG64 N)
{
	//a = a%N;
	//b = b%N;
	
	LONG64 i = a >> 32;
	LONG64 k = b >> 32;
	LONG64 j = a & 0xffffffffll;
	LONG64 l = b & 0xffffffffll;
	LONG64 ik = (i*k)%N;
	LONG64 il = (i*l)%N;
	LONG64 jk = (j*k)%N;
	LONG64 jl = (j*l)%N;
	LONG64 t = ( il + jk); 
	//for (int s=0;s<64;s+=16) ik = (65536*ik)%N;
	//for (int s=0;s<32;s+=16) t = (65536*t)%N;
	for (int s=0;s<64;s+=1) ik = (2*ik)%N;
	for (int s=0;s<32;s+=1) t = (2*t)%N;
	LONG64 r = (ik + t);
	r = (r + jl)%N;
//	printf("%15.15llu -(( %15.15llu * %15.15llu ) %% %15.15llu)\n",r,a,b,n); 
	return r; 
}
static LONG64 ipow(LONG64 base, LONG64 exp,LONG64 N)
{
    LONG64 result = 1ll;
    while (exp)
    {
        if (exp & 1ll) result = mult(result,base,N);
        exp >>= 1ll;
        base = mult(base,base,N);
    }

    return result;
}

static LONG64 ipow2(LONG64 base, LONG64 exp,LONG64 N,unsigned n)
{
    LONG64 cnt=0; 
    LONG64 p2 =1; 
    p2 = 1 << (n-2); 
    p2=(p2*2)%N; 
    p2=(p2*2)%N; 
    int ck =1; 
    LONG64 result = 1ll;
    while (exp)
    {
	
    	if (exp & 1ll) { result=MMUL(result,base,N,n); cnt+=ck; }
        exp >>= 1ll;
	base=MMUL(base,base,N,n); 
	ck = 2*ck; 
    }
    for (unsigned i=0;i<cnt*n;i++) result = (result * 2)%N;
    return result;
}
static unsigned isPrime(LONG64 N)
{
	LONG64 primes[]={2,3,5,7,11,13,17};
//	unsigned n = 1 + llog2(N); 
	
	unsigned i,j;
	LONG64 s,r;
	s = N-1;
	for (r=0;s%2==0;r++) s=s>>1;
	unsigned isprime; 
	for (i=0;i<7 && primes[i]<(N-1);i++)
	{	
		isprime=0; 
		LONG64 a = primes[i]; 
		LONG64 x = ipow(a,s,N); 
		if (x ==1 || x == (N-1)) 
		{
			isprime=1;
			continue;
		}
		for (j=1;j<r;j++)
		{
			x=mult(x,x,N); 
			if (x == 1) return 0;
			if (x==N-1) { isprime=1; break; }
		}
		if (isprime == 0) return 0;
	}
	return 1;
}



static LONG64 trial_divide(LONG64 N ,LONG64 max) 
{
  // Trial divides the positive integer N by the primes from 2 to max
  // Returns the first prime divisor found, or 0 if none found
  // Note: if N < max^2 is a prime, then N will be returned. 
  if (N%2ll == 0) return 0ll;
  if (N%3ll == 0) return 0ll;
  // No need to go past the square root of our number
  unsigned Stop = MIN(sqrt(1.0*N),max);
  // Okay, lets "wheel factor" alternately adding 2 and 4
  LONG64 di=2;
  for(LONG64 i=5; i<=Stop; i+=di) { di = 6 - di; if (N%i == 0) return 0ll; }
  if (N >= max*max) return 2ll; 
  return 1ll;
}

static LONG64 nextSafePrime(LONG64 prime)
// return first safe prime larger than prime. prime doesn't need to be a
// prime but does need to be such that (prime%12) == 1.
{
	while (1)
	{
		if (isPrime(prime+=12))
		{  
       		LONG64 p = (prime-1)/2; 
	   		if (isPrime1(p)) return prime; 
		}
	}
}

static double nthPrime(LONG64 n,double *lb, double *ub) // Estimate of the nth prime 
{
	double l = log(1.0*n); 
	double ll = log(log(1.0*n)); 
	double nP = l+ll-1.0+(ll-2)/l-(ll*ll-6*ll+11)/(2*l*l);
	nP *=n; 
	*lb =  n*(l + ll-1); 
	*ub =  n*(l + ll-1+1.8*ll/l); 
	return nP; 
	
}
#endif

#if 0
int main(int iargc, char *argv[])
{
//#include "4kSafePrimeBins.h"
	unsigned id=atol(argv[1]); 
	unsigned ntask=atol(argv[2]); 
	unsigned nlocal=atol(argv[3]); 
	unsigned kbin = id>>12;
	unsigned kpos = id % 4096; 
	unsigned isprime ;
	unsigned cntPrimes=0; 
	unsigned cntSafePrimes=0; 
	unsigned issafe=0; 
	LONG64 lpe; 
	unsigned n;
  	LONG64  p ; 
	LONG64 N; 
	LONG64 prime; 
	prime = (1ll << 62) - 1ll+(1ll<<32)*id;
	printf("%lf\n",log(prime*1.0)/log(2.0)); 
	
    unsigned cnt =0;
	double t1 = pcputime(); 
	for (unsigned i=0;i<nlocal;i++)
	{
		prime = nextPrime(); 
		printf("%u %llu\n",cnt,prime); 
		cnt++; 
	}
	double t2 = pcputime(); 
	double dt = t2 -t1; 
	printf("%llu %f\n",prime,dt); 
/*
	for ( n=safePrimeBins[kbin];cntSafePrimes <=kpos; n+=12)
	{
		N =n; 
		isprime = isPrime1(N);
		issafe=0; 
		if (isprime)
		{  
       		p = (N-1)/2; 
	   		issafe = isPrime1(p);
       		if (issafe) 
			{
				if (cntSafePrimes==kpos) break; 
				cntSafePrimes++; 
			}
			cntPrimes += 1; 
		}
	}
 	for (lpe=1;lpe<101;lpe+=2) 
	{ 
		LONG64 k2 = ipow(lpe, 2, N);
		LONG64 kp = ipow(lpe, p, N);
   		if ((k2 != 1) && (kp != 1) )break;
   	} 
	double t2 = pcputime(); 
	double dt = t2-t1; 
	printf("%u %u %u %llu %f\n",cntPrimes,cntSafePrimes,n,lpe,dt); 
	LONG64 aN=lpe; 
	unsigned cnt=0; 
	for (unsigned  k = 1; k<p;k++)
	{
		N = 2*k+1; 
		aN = (aN * lpe)%n;
		aN = (aN * lpe)%n;
		isprime = isPrime1(N);
		if (isprime) 
		{
			printf("%u %llu %llu\n",cnt,N,aN);
			cnt ++; 
			if (cnt == nlocal) break; 
		}
	}
*/
}
#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
