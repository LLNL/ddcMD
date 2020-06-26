#ifndef GET_REMOTE_DATA_H
#define GET_REMOTE_DATA_H

typedef void (*fcn_type) (void* commBuf, int nItems, unsigned* itemIndices,
			  int* ddcSendList, void* localData);

void getRemoteData(unsigned* iLocal, unsigned iLocalSize, void* localData,
		   fcn_type loaderFcn,
		   unsigned itemSize, void* remoteData,
		   int* nSend, int* nRecv);

#endif // #ifdef GET_REMOTE_DATA_H
