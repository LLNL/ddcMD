enum ExpandbufferEnum { REALLOC, FREE_MALLOC};
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <unistd.h>
#include <sys/resource.h>
#include <time.h>
#include <errno.h>
#include<sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <assert.h>
#ifndef __APPLE__
#include <malloc.h>
#endif
#include "ddcMalloc.h"
#include "pio.h"
#include "heap.h"
#include "mpiUtils.h"
#include "gid.h"

SIGNED64 simulate_getLoop(void *); 

static struct
{
	int currentLength, maxLength, ElementSize;
	void *StartPtr;
	void *StartPtrCopy;
	char *name;
}

*_bufferTable = NULL;
static int _tableSize = 0;
static int _memorysize = 0;
static int _mode = REALLOC; 
static int _tableCapacity = 0;
static int _tableCapacityIncrement = 128;
static int _lastLength=0; 

int getRank(int);

#ifndef __APPLE__
static struct mallinfo minfoStart[64] ;
static struct mallinfo minfoEnd[64] ;
#endif
void memorystart(char *label, int key)
{
#ifndef __APPLE__
   minfoStart[key] = mallinfo();
#endif   
}
void memoryend(char *label,int key)
{
#ifndef __APPLE__
   minfoEnd[key] = mallinfo();
	SIGNED64 loop=simulate_getLoop(NULL); 
   int delta = minfoEnd[key].uordblks-minfoStart[key].uordblks; 
	if (getRank(0) ==0 && delta != 0 ) printf("%"PRId64" %s %d %d %d\n", loop,label ,delta,minfoEnd[key].arena,minfoEnd[key].uordblks); 
#endif
}
void memoryinfo(FILE*file)
{
#ifdef __APPLE__
   fprintf(file, "In routine memoryinfo no mallinfo on OS X.  Sorry.\n");
#else
   struct mallinfo minfo;
   minfo = mallinfo();
   fprintf(file, "Malloc: heap size =%f MB used=%f MB unused =%f MB unused blocks =%d\n", minfo.arena/(1024*1024.), minfo.uordblks/(1024*1024.), minfo.fordblks/(1024*1024.), minfo.ordblks);
   heapSummary(file);
   fflush(file);
#endif
}
void memoryinfo_pio(PFILE*file)
{
#ifdef __APPLE__
   Pprintf(file, "No mallinfo on OS X.  Sorry.\n");
#else
   struct mallinfo minfo;
   minfo = mallinfo();
   Pprintf(file, "Malloc: heap size =%f MB used=%f MB unused =%f MB unused blocks =%d\n", minfo.arena/(1024*1024.), minfo.uordblks/(1024*1024.), minfo.fordblks/(1024*1024.), minfo.ordblks);
   heapSummary_pio(file);
#endif
}

void ExpandBuffersPrint(FILE*file)
{
	int i, currentLength, maxLength;
	
	fprintf(file, "ExpandBuffer memory allocation:  %d kbytes on proc: %d\n", _memorysize/1024, getRank(0));
	fprintf(file, "Index   Name    Address  Current Length (kB)  Max Length (kB)  Element Size(B)\n");;
	currentLength = maxLength = 0;
	for (i = 0; i < _tableSize; i++)
	{
	   fprintf(file, "%3d %16s %8p %10d %10d %10d\n", 
		i, _bufferTable[i].name, _bufferTable[i].StartPtr, _bufferTable[i].currentLength/1024, _bufferTable[i].maxLength/1024, _bufferTable[i].ElementSize);
		currentLength += _bufferTable[i].currentLength;
		maxLength += _bufferTable[i].maxLength;
	}
	fprintf(file, "%3d %16s %8p %10d %10d\n", i, "Totals", NULL, currentLength/1024, maxLength/1024);
	memoryinfo(file);
	fflush(file);
	return ;
}
void ExpandBuffersPrint_pio(PFILE*file)
{
	int i, currentLength, maxLength;
	
	Pprintf(file, "ExpandBuffer memory allocation:  %d kbytes on proc: %d\n", _memorysize/1024, getRank(0));
	Pprintf(file, "Index   Name    Address  Current Length (kB)  Max Length (kB)  Element Size(B)\n");;
	currentLength = maxLength = 0;
	for (i = 0; i < _tableSize; i++)
	{
		Pprintf(file, "%3d %16s %8p %10d %10d %10d\n", i, _bufferTable[i].name, _bufferTable[i].StartPtr, _bufferTable[i].currentLength/1024, _bufferTable[i].maxLength/1024,
			_bufferTable[i].ElementSize);
		currentLength += _bufferTable[i].currentLength;
		maxLength += _bufferTable[i].maxLength;
	}
	Pprintf(file, "%3d %16s %8p %10d %10d\n", i, "Totals", NULL, currentLength/1024, maxLength/1024);
	memoryinfo_pio(file);
}

void ExpandBuffersResize(void)
{
   return;
	int index, size, max;
	_memorysize = 0;
	for (index = 0; index < _tableSize; index++)
	{
	   max = _bufferTable[index].currentLength;
	   size = _bufferTable[index].ElementSize;
	   if (max==0) max = size; 
	   if (_bufferTable[index].StartPtr != _bufferTable[index].StartPtrCopy)
	   {
	      printf("This can't happen. (task %d)\n"
		     "StartPtr != StartPtrCopy in ExpandBufferResize\n"
		     "StartPtr=%p StartPtrCopy=%p index=%d\n",
		     getRank(0),
		     _bufferTable[index].StartPtr, _bufferTable[index].StartPtrCopy,
		     index);
	      exit(1);
	   }
	   _bufferTable[index].StartPtr = ddcRealloc(_bufferTable[index].StartPtr, max);
	   _bufferTable[index].StartPtrCopy = _bufferTable[index].StartPtr;
	   _bufferTable[index].maxLength = max;
	   _memorysize += max;
	}
}

int ExpandBuffersMaxLength(void *buffer)
{
	int index;
	if (buffer == NULL) return 0;
	for (index = 0; index < _tableSize; index++) if (_bufferTable[index].StartPtr == buffer) break;
	if (index == _tableSize) return 0;
	return _bufferTable[index].maxLength;
}
void ExpandBuffersMode(char *mode_string)
{
	_mode = REALLOC;
	if (strcmp(mode_string,"FREE_MALLOC")==0) _mode = FREE_MALLOC; 
}

int  ExpandBuffersLastLength(void)
{
	return _lastLength;
}
void ExpandBuffersFree(void *buffer)
{
	FILE *file = stderr;
	if (buffer == NULL) return; 
   int index = -1; 
   for (index = 0; index < _tableSize; index++)
   if (_bufferTable[index].StartPtr == buffer) break;
   if (index == _tableSize)
   {
      fprintf(file, "ExpandBuffers Error:\n"
	      "Request to free non-null pointer that is not in _bufferTable.\n   buffer=%p", buffer);
      fflush(file);
      exit(1);
   }
   _memorysize -=  _bufferTable[index].maxLength ;
   free(_bufferTable[index].name); 
  _bufferTable[index].name = NULL; 
  _bufferTable[index].currentLength = 0;
  _bufferTable[index].ElementSize =     0;
  _bufferTable[index].StartPtr = NULL;
  _bufferTable[index].StartPtrCopy = NULL;
  _bufferTable[index].maxLength = 0;
  _lastLength = 0; 
	if (index == _tableSize-1)
   {
		_tableSize--;
     	_lastLength = _bufferTable[_tableSize-1].maxLength ;
   }
}
void *ExpandBuffers(void *buffer, int size, int n, int increment, char *string, ...)
{
	va_list ap;
	char *name;
	int cnt, i;
	int max;
	FILE *file = stderr;
	int index = -1;
	assert(increment > 0); // we divide by increment
	if (buffer != NULL)
	{
	   for (index = 0; index < _tableSize; index++)
	      if (_bufferTable[index].StartPtr == buffer) break;
	   // It is an error to pass in a non-NULL pointer unless the
	   // pointer was previously allocated by ExpandBuffer (and can
	   // therefore be located in _bufferTable).
	   if (index == _tableSize)
	   {
	      va_start(ap, string);
	      name = va_arg(ap, char *);
	      va_end(ap);
	      fprintf(file, "ExpandBuffers Error on task %d\n"
		      "  Request to expand non-null pointer that is not in _bufferTable.\n"
		      "  buffer=%p  name=%s\n",
		      getRank(0), buffer, name);
	      fflush(file);
	      exit(1);
	   }
	}
	else // buffer == NULL
	{
		_tableSize++;
		if (_tableCapacity < _tableSize)
		{
		   _tableCapacity += _tableCapacityIncrement;
		   _bufferTable = ddcRealloc(_bufferTable, (_tableCapacity)*sizeof(*_bufferTable));
		}
		if (_bufferTable == NULL)
		{
			fprintf(file, "ID=%d _bufferTable =%p index=%d _tableSize=%d ", getRank(0), _bufferTable, index, _tableSize);
			perror("ExpandBuffer(a):");
			fflush(file);
			exit(1);
		}
		index = _tableSize - 1;
		_bufferTable[index].maxLength = 0;
		_bufferTable[index].StartPtr = NULL;
		_bufferTable[index].StartPtrCopy = NULL;
		va_start(ap, string);
		name = va_arg(ap, char *);
		va_end(ap);
		if (name > (char *)10) _bufferTable[index].name = strdup(name);
		else
		{
			_bufferTable[index].name = ddcMalloc(24);
			sprintf(_bufferTable[index].name,"_bufferTable%4.4d", index);
		}
		max = 0;
	}
	
	_bufferTable[index].currentLength = n*size;
	_bufferTable[index].ElementSize = size;
	if (n*size > _bufferTable[index].maxLength || _bufferTable[index].maxLength == 0)
	{
		max = MAX(_bufferTable[index].maxLength, size*(n/increment + 1)*increment);
		cnt = max - _bufferTable[index].maxLength;
		if ((_mode == FREE_MALLOC)  && (buffer != NULL)) {free(buffer); buffer=NULL;}
		_bufferTable[index].StartPtr = buffer = ddcRealloc(buffer, max);
		_bufferTable[index].StartPtrCopy = buffer;
		if (buffer == NULL )
		{
			char errfilename[1024]; 
			sprintf(errfilename,"Expandbuffer-%6.6d",getRank(0));
			file=fopen(errfilename,"w");
			fprintf(file, "Error in ExpandBuffer on task = %d\n",getRank(0));
			fprintf(file, "realloc call at %s returned NULL pointer\n", string);
			perror("System Error message is   ");
			fprintf(file, "Failed trying to expand buffer %d by %d bytes\n", index, cnt);
			fprintf(file, "Call requested buffer to be sized to %d elements of size %d elements\n", n, size);
			fprintf(file, "The requested expand of buffer %d would of required a total of %d bytes\n", index, _memorysize + cnt);
			//ExpandBuffersPrint(file);
			//ddcMemReport(file);
			fclose(file); 
			abortAll(1);
		}
		else 
		{
			_memorysize += cnt;
			for (i = _bufferTable[index].maxLength; i < max; i++)
				((char *)buffer)[i] = '\0';
			_bufferTable[index].maxLength = max;
		}
		_lastLength = _bufferTable[index].maxLength; 
	}
	return buffer;
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
