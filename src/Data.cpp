#include "Data.h"

/***************
* class Data
* abstract class maintaining the data
* Florian Pein, 2015
***************/
NumericVector Data::criticalValues_;

void Data::setCriticalValues(const List &input) {
  criticalValues_ = input["q"];
}

Data::~Data() {}
