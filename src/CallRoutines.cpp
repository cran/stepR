#include <Rcpp.h>

#include "DynamicPrograms.h"

#include "Data.h"
#include "DataGauss.h"
#include "DataHsmuce.h"
#include "DataMDependentPS.h"

#include "IntervalSystem.h"
#include "IntervalSystemAll.h"
#include "IntervalSystemAllLengths.h"
#include "IntervalSystemDyaLen.h"
#include "IntervalSystemDyaLenLengths.h"
#include "IntervalSystemDyaPar.h"
#include "IntervalSystemDyaParLengths.h"

using namespace Rcpp;

// [[Rcpp::export(name = ".callRoutines")]]
RObject callRoutines(RObject observations,
                     int routineType, List argumentsListRoutine,
                     int dataType, List argumentsListData,
                     int intervalSystemType, List argumentsListIntervalSystem) {
  Data* data = NULL;
  
  switch (dataType) {
  case 0:
    DataGauss::setData(observations, argumentsListData);
    data = new DataGauss();
    break;
  case 10:
    DataMDependentPS::setData(observations, argumentsListData);
    data = new DataMDependentPS();
    break;
  case 20:
    DataHsmuce::setData(observations);
    data = new DataHsmuce();
    break;
  default:
    stop("dataType %d is not defined", dataType);
  }

  IntervalSystem* intervalSystem = NULL;
  switch(intervalSystemType) {
  case 0:
    intervalSystem = new IntervalSystemAll(data -> getN());
    break;
  case 1:
    intervalSystem = new IntervalSystemAllLengths(data -> getN(), argumentsListIntervalSystem);
    break;
  case 10:
    intervalSystem = new IntervalSystemDyaLen(data -> getN());
    break;
  case 11:
    intervalSystem = new IntervalSystemDyaLenLengths(data -> getN(), argumentsListIntervalSystem);
    break;
  case 20:
    intervalSystem = new IntervalSystemDyaPar(data -> getN());
    break;
  case 21:
    intervalSystem = new IntervalSystemDyaParLengths(data -> getN(), argumentsListIntervalSystem);
    break;
  default:
    data -> cleanUpStaticVariables();
    delete data;
    stop("intervalSystemType %d is not defined", intervalSystemType);
  }
  
  RObject ret;
  switch (routineType) {
  case 0:
    ret = intervalSystem -> computeMultiscaleStatisticNull(data);
    break;
  case 1:
    ret = intervalSystem -> computeMultiscaleStatistic(data, argumentsListRoutine);
    break;        
  case 2:
    Data::setCriticalValues(argumentsListRoutine);
    ret = intervalSystem -> computeBounds(data); 
    break;
  case 3:
    Data::setCriticalValues(argumentsListRoutine);
    ret = fitSimpleDynamicProgram(data, intervalSystem);
    break;
  case 4:
    Data::setCriticalValues(argumentsListRoutine);
    ret = fitIntervalDynamicProgram(data, intervalSystem);  
    break;
  case 5:
    Data::setCriticalValues(argumentsListRoutine);
    ret = fitBandDynamicProgram(data, intervalSystem);
    break;
  default:
    delete intervalSystem;
    data -> cleanUpStaticVariables();
    delete data;
    stop("routineType %d is not defined", routineType);
  }
  
  delete intervalSystem;
  data -> cleanUpStaticVariables();
  delete data;
  
  return ret;
}
