#ifndef TESTINGUTILS_H
#define TESTINGUTILS_H
#include "commandLineOptions.h"
#include "CuTest.h"
#include <mpi.h>
#include "object.h"
//#include "masterFactory.h"
#include "simulate.h"
#include "routineManager.h"
#include "primes.h"
#include "units.h"
#include "codata.h"
#include <stdlib.h>
#include <stdio.h>
#include "mpiUtils.h"
#include "masters.h"
#include "ddcenergy.h"
#include "parityHandler.h"
#include "utilities.h"
#include "hpmWrapper.h"
#include "auxNeighbor.h"
#include <string.h>
SIMULATE* setupSimulate( int i);
SIMULATE* setupSimulate2();
void setupEnergyCall(SIMULATE * simulate);
void setupFirstEnergyCall(SIMULATE * simulate);
void finishEnergyCall(SIMULATE * simulate);
#endif
