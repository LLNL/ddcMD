/* This acts as an interface to the central data broker. */
#ifndef CDBINTERFACE_H
#define CDBINTERFACE_H

#if USE_DB
#include <libdatabroker.h>
//#include "test_utils.h"
#endif

#include "energyInfo.h"
#include "system.h"
#include "gid.h"
//#include "simulate.h"

/* NOCDB: do not use DB
 * DB: Use DB for LJ Enhanced sampling
 * TRANS_CHARMM: use for transforming molecules with Charmm potential
 */
#define TEST_DB( function, expect ) ( (function)==(expect)? 0 : 1 )
#define TEST_NOT( function, expect ) ( 1 - TEST_DB( function, expect ) )

enum DATABROKER_CLASS{ NOCDB, CDB, TRANS_CHARMM, HYCOP_DB };


//#include "system.h"

/* Need:
 * Namespace for DB
 * ...

 */
typedef struct databroker_st
{
  //Type of data manager
  char *type;
  //Required for object struct
  char *value;
  char *name;
  char *objclass;
  void *parent; 
  //Name for the central data broker
  char *name_space;
  //Name for this instance of ddcMD
  char *mdName;
#if USE_DB
  //handle to the cdb
  DBR_Handle_t handle; 
#endif
  //Do I want to use the CDB?
  int useDb;
  //Do I start the CDB? 
  //(If false but we are using CDB, attach)
  int startDb;
  //Do I stop the CDB?
  int stopDb;
  //Update cdb values frequency
  //For CDB types that update within simulation loop 
  int updateRate;

} DATABROKER;

/* Initialize Data Broker struct
 * with information from object.data
 */
DATABROKER *dataBroker_init(void*, char *);

/* Either start the DB or attatch
 * to an existing DB name space
 */
int connectToDB(DATABROKER*);
/* Closes the DB with name DB->name_space*/
int shutdownDB(DATABROKER*);

/*Keep trying to connect o DB until 
 * successful.
 */
int waitUntilConnectToDB(DATABROKER*);

/*Synchronously read array of integers
 * accessed by a single key word.
 * Arguments:
 * db: data broker object
 * key: key to access key value pair
 * delim: deliminator to process value content
 * array: array to store final values in
 * size: size of integer array 
 * Returns:
 *  error value
 */
int readIntegerArray(DATABROKER*, char*, char*, int*, int);

/*Synchronously read array of integers
 * accessed by an ordered set of key words
 * Arguments:
 * db: data broker object
 * key_base: base name for keys
 * array: array to store final values in
 * size: size of integer array 
 * 
 * Returns: error value
 */
int readIntegerArray_orderedKeys(DATABROKER*, char*, int*, int);

/*Synchronously write array of integers to db
 * using a single key word.
 * Arguments:
 * db: data broker object
 * key: key to access key value pair
 * delim: deliminator to separate value content
 * array: input integer array
 * size: size of integer array 
 * 
 * Returns: error value
 */
int putIntegerArray(DATABROKER*, char*, char*, int*, int);

/* Asynchronously write array of integers to db
 * using a single key word.
 * Arguments:
 * db: data broker object
 * key: key to access key value pair
 * delim: deliminator to separate value content
 * array: input integer array
 * size: size of integer array 
 * 
 * Returns: error value
 */
int putIntegerArrayAsync(DATABROKER*, char*, char*, int*, int);

/* Synchronously write array of integers
 * accessed by an ordered set of key words
 * Arguments:
 * db: data broker object
 * key_base: base name for keys
 * array: input integer array
 * size: size of integer array 
 * 
 * Returns: error value
 */
int putIntegerArray_orderedKeys(DATABROKER*, char*, int*, int);

/* Asynchronously write array of integers
 * accessed by an ordered set of key words
 * Arguments:
 * db: data broker object
 * key_base: base name for keys
 * array: input integer array
 * size: size of integer array 
 * 
 * Returns: error value
 */
int putIntegerArrayAsync_orderedKeys(DATABROKER*, char*, int*, int);

/* ddcMD stores average energies
 * This function posts the average energies
 * for a given timestep to the DB with
 * name DB->name_space
 * ETYPE [=] struct to store average energy components
 *   in internal energy
 * int loop [=] MD simulation loop (not time)  
 */
int postAverageEnergies(DATABROKER*, ETYPE, SIGNED64, unsigned int);
//int postAverageEnergies(SIMULATE*);

/* ddcMD stores average energies
 * This function posts the total energies
 * for a given timestep to the DB with
 * name DB->name_space
 * ETYPE [=] struct to store average energy components
 *   in internal energy
 * loop [=] MD simulation loop  
 */
int postTotalEnergies(DATABROKER*, ETYPE, SIGNED64);

//convert everything to a character
//array and include units 
//use names ddcWorkDir-PotentialE etc. 
//int postEnergy(SYSTEM* system) {return -1;};

/* Read weights for species 1 and 2 from db */
//int readTransformWeights(DATABROKER*, double*, double*, int); 

/* Function to update parms specific to different system structs */
int dbUpdateInitialParms(DATABROKER*, SYSTEM*);

/*Update based on simulation loop and updateRate for the data broker
 * This is a general function with different calls depending 
 * on the DB type.
 * Input: databroker struct, sys struct, loop number, flag to
 * post total energies
 */
int updateDataBrokerContent(DATABROKER *, SYSTEM*, SIGNED64, int);

/* Puts a flag in the central databroker to indicate that
 * this instance of ddcMD is completing a simulation.
 * Expected input:
 * char * keyBase: base name for key
 * int loop: current MD simulation loop
 */
int postCompletion(DATABROKER *, char *,  SIGNED64);

#endif
