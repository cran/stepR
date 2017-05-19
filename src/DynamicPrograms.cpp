#include "DynamicPrograms.h"

#include "SingleBounds.h"
#include "LocalOptimum.h"

#include <vector>
#include <stack>

/*************
* fitSimpleDynamicProgram
* fits the data with a simple dynamic program
****************/
List fitSimpleDynamicProgram(Data * const data, IntervalSystem * const intervalSystem) {
  unsigned int n = data -> getN();
  std::vector<LocalOptimum> optimum(n);           // optimum[i]: current local solution to reach i
  std::vector<SingleBounds> intersectedBounds(n); // intersectedBounds[i]:
                                                  // intersected bounds [i, current right index]
  
  unsigned int j = 0u;
  do {
    checkUserInterrupt();
    // find constant solution over [0, ..., j]
    data -> reset();
    data -> addRight(j);
    if (intervalSystem -> isInIntervalSystem(j, j)) {
      intersectedBounds[j].intersect(data -> computeSingleBounds());
    }
    unsigned int i = j;
    while (i > 0u) {
      --i;
      intersectedBounds[i].intersect(intersectedBounds[i + 1u]);
      data -> addLeft(i);
      if (intervalSystem -> isInIntervalSystem(i, j)) {
        intersectedBounds[i].intersect(data -> computeSingleBounds());
      }
    }    
    
    if (!intersectedBounds[0u].feasible()) {
      // no constant solution on [0, ..., j] possible
      break;
    }
    
    optimum[j] = data -> computeLocalOptimum(0, j, intersectedBounds[0u], NULL);
    ++j;
  } while (j < n);
  
  // dynamic program
  // find solution with K change-points
  unsigned int K = 0u;
  while (j < n) {
    ++K;
    unsigned int startPointRightIndex = j;
    
    // forward step
    do {
      checkUserInterrupt();
      // find constant solution over [startPointRightIndex, ..., j]
      data -> reset();
      data -> addRight(j);
      if (intervalSystem -> isInIntervalSystem(j, j)) {
        intersectedBounds[j].intersect(data -> computeSingleBounds());
      }
      unsigned int i = j;
      while (i > startPointRightIndex) {
        --i;
        intersectedBounds[i].intersect(intersectedBounds[i + 1u]);
        data -> addLeft(i);
        if (intervalSystem -> isInIntervalSystem(i, j)) {
          intersectedBounds[i].intersect(data -> computeSingleBounds());
        }
      }    
      
      if (!intersectedBounds[startPointRightIndex].feasible()) {
        // no constant solution on [startPointRightIndex, ..., j] possible
        break;
      }
      
      optimum[j] = data -> computeLocalOptimum(startPointRightIndex, j,
                                               intersectedBounds[startPointRightIndex],
                                               &optimum[startPointRightIndex - 1u]);
      
      // backward step
      while (i > 0u) {
        --i;
        // find constant solution over [i, ..., j]
        intersectedBounds[i].intersect(intersectedBounds[i + 1u]); 
        data -> addLeft(i);
        if (intervalSystem -> isInIntervalSystem(i, j)) {
          intersectedBounds[i].intersect(data -> computeSingleBounds());
        }
        
        if (!intersectedBounds[i].feasible()) {
          // no constant solution on [i, ..., j] possible
          break;
        }
        
        optimum[j].update(data -> computeLocalOptimum(i, j, intersectedBounds[i], &optimum[i - 1u]));
      }
      ++j;
    } while (j < n);
  }
  
  // return result
  IntegerVector leftIndex = IntegerVector(K + 1u);
  IntegerVector rightIndex = IntegerVector(K + 1u);
  NumericVector value = NumericVector(K + 1u);
  List ret = List::create(_["leftIndex"] = leftIndex, _["rightIndex"] = rightIndex, 
                          _["value"] = value);

  // return cost as attribute
  ret.attr("cost") = optimum[n - 1u].costs();
  
  // backtracking
  LocalOptimum localOptimum = optimum[n - 1u];
  unsigned int k = K;
  do {
    value[k] = localOptimum.estimatedValue();
    leftIndex[k] = localOptimum.leftIndex() + 1;
    rightIndex[k] = localOptimum.rightIndex() + 1; // + 1 for R-style index
    
    if (k > 0u) {
      localOptimum = *localOptimum.lastSegment();
    } else {
      break;
    }
    --k;
  } while (true);
  
  return ret;
}

/*************
* fitIntervalDynamicProgram
* fits the data with a dynamic program restricted to the confidence intervals
****************/
List fitIntervalDynamicProgram(Data * const data, IntervalSystem * const intervalSystem) {
  unsigned int n = data -> getN();
  std::vector<SingleBounds> intersectedBounds(n);
  
  // computation of left condidence intervals
  // intersectedBounds[i] will contain the intersected bounds from i to left.top() - 1
  std::stack<unsigned int> left;  // left confidence intervals starting with n, will be reordered
  left.push(n);
  
  unsigned int i = n;
  bool feasible = false;
  while (i > 0u) {
    --i;
    checkUserInterrupt();
    
    if (feasible) {
      intersectedBounds[i].intersect(intersectedBounds[i + 1u]);
    }
    
    data -> reset();
    for (unsigned int j = i; j < left.top(); ++j) {
      data -> addRight(j);
      
      if (intervalSystem -> isInIntervalSystem(i, j)) {
        intersectedBounds[i].intersect(data -> computeSingleBounds());
      }
    }    
    
    if (intersectedBounds[i].feasible()) {
      feasible = true;
    } else {
      feasible = false;
      left.push(i + 1u);
      intersectedBounds[i].reset();
      ++i;
    }
  }
  
  // reordering the left limits of the confidence intervals
  unsigned int K = left.size() - 1u;
  IntegerVector leftConfInt = IntegerVector(K);
  
  for (unsigned int k = 0u; k < K; ++k) {
    leftConfInt[k] = left.top();
    left.pop();
  }
  
  IntegerVector rightConfInt = IntegerVector(K);
  std::vector<LocalOptimum> optimum(n);    // optimum[i]: current local solution to reach i
  // intersectedBounds[i] contains the intersected bounds from i to leftConfInt[k] - 1u
  
  if (K > 0u) {
    // segment [0, leftConfInt[0] - 1]
    data -> reset();
    for (int i = 0; i < leftConfInt[0u]; ++i) {
      data -> addRight(i);
    }
    optimum[leftConfInt[0u] - 1] = data -> computeLocalOptimum(0, leftConfInt[0u] - 1,
                                                               intersectedBounds[0u], NULL);
    
    // first segment [0, j] until no constant solution is feasible
    for (unsigned int j = leftConfInt[0u]; j < n; ++j) {
      checkUserInterrupt();
      data -> reset();
      unsigned int i = j + 1u;
      while (i > 0u) {
        --i;
        data -> addLeft(i);
        if (intervalSystem -> isInIntervalSystem(i, j)) {
          intersectedBounds[0u].intersect(data -> computeSingleBounds());
        }
      }
      
      if (!intersectedBounds[0u].feasible()) {
        // j cannot reached without a change-point
        rightConfInt[0u] = j;
        break;
      }
      optimum[j] = data -> computeLocalOptimum(0, j, intersectedBounds[0u], NULL);
    }    
    
    // dynammic program
    for (unsigned int k = 1u; k < K; ++k) {
      data -> reset();
      
      for (unsigned int i = leftConfInt[k] - 1u; i > static_cast<unsigned int>(rightConfInt[k - 1u]); --i) {
        data -> addLeft(i);
      }
      
      // backward step for rightIndex = leftConfInt[k] - 1u
      for (unsigned int i = rightConfInt[k - 1u]; i >= static_cast<unsigned int>(leftConfInt[k - 1u]); --i) {
        data -> addLeft(i);
        optimum[leftConfInt[k] - 1u].update(data -> computeLocalOptimum(i, leftConfInt[k] - 1u,
                                                                        intersectedBounds[i],
                                                                        &optimum[i - 1u]));
      }
      
      // step forward
      for (unsigned int j = leftConfInt[k]; j < n; ++j) {
        checkUserInterrupt();
        data -> reset();
        
        for (unsigned int i = j; i >= static_cast<unsigned int>(rightConfInt[k - 1u]); --i) {
          data -> addLeft(i);
          if (intervalSystem -> isInIntervalSystem(i, j)) {
            intersectedBounds[rightConfInt[k - 1u]].intersect(data -> computeSingleBounds());
          }
        }        
        
        if (!intersectedBounds[rightConfInt[k - 1u]].feasible()) {
          // no constant solution on [startPointRightIndex, ..., j] possible
          rightConfInt[k] = j;
          break;
        }
        
        optimum[j] = data -> computeLocalOptimum(rightConfInt[k - 1u], j,
                                                 intersectedBounds[rightConfInt[k - 1u]],
                                                 &optimum[rightConfInt[k - 1u] - 1u]);
        
        // backward step
        for (unsigned int i = rightConfInt[k - 1u] - 1; i >= static_cast<unsigned int>(leftConfInt[k - 1u]);
          --i) {
          // find constant solution over [i, ..., j]
          intersectedBounds[i].intersect(intersectedBounds[i + 1u]); 
          data -> addLeft(i);
          if (intervalSystem -> isInIntervalSystem(i, j)) {
            intersectedBounds[i].intersect(data -> computeSingleBounds());
          }
          
          if (!intersectedBounds[i].feasible()) {
            // no constant solution on [i, ..., j] possible
            break;
          }
          
          optimum[j].update(data -> computeLocalOptimum(i, j, intersectedBounds[i], &optimum[i - 1u]));
        }
      }
    }
    
    // last segment
    data -> reset();
    for (unsigned int i = n - 1u; i >= static_cast<unsigned int>(rightConfInt[K - 1u]); --i) {
      data -> addLeft(i);
    }    
    optimum[n - 1u] = data -> computeLocalOptimum(rightConfInt[K - 1u], n - 1u, 
                                                  intersectedBounds[rightConfInt[K - 1u]],
                                                  &optimum[rightConfInt[K - 1u] - 1u]);
    
    for (unsigned int i = rightConfInt[K - 1u] - 1u; i >= static_cast<unsigned int>(leftConfInt[K - 1u]);
      --i) {
      checkUserInterrupt();
      data -> addLeft(i);
      optimum[n - 1u].update(data -> computeLocalOptimum(i, n - 1u, intersectedBounds[i], &optimum[i - 1u]));
    }
  }
  
  // return result
  IntegerVector leftIndex = IntegerVector(K + 1u);
  IntegerVector rightIndex = IntegerVector(K + 1u);
  NumericVector value = NumericVector(K + 1u);
  List ret = List::create(_["leftIndex"] = leftIndex, _["rightIndex"] = rightIndex, 
                          _["value"] = value,
                          _["leftConfInt"] = leftConfInt, _["rightConfInt"] = rightConfInt);
  
  for (unsigned int i = 0; i < K; ++i) {
    ++leftConfInt[i];   // C to R-style
    ++rightConfInt[i];  // C to R-style
  }
  
  // backtracking
  if (K > 0u) {
    LocalOptimum localOptimum = optimum[n - 1u];  
    unsigned int k = K;
    do {
      value[k] = localOptimum.estimatedValue();
      leftIndex[k] = localOptimum.leftIndex() + 1;
      rightIndex[k] = localOptimum.rightIndex() + 1; // + 1 for R-style index
      
      if (k > 0u) {
        localOptimum = *localOptimum.lastSegment();
      } else {
        break;
      }
      --k;
    } while (true);
    // return cost as attribute
    ret.attr("cost") = optimum[n - 1u].costs();
  } else {
    LocalOptimum optimum = data -> computeLocalOptimum(0u, n - 1u, intersectedBounds[0u], NULL);
    value[0u] = optimum.estimatedValue();
    leftIndex[0u] = optimum.leftIndex() + 1u;
    rightIndex[0u] = optimum.rightIndex() + 1u;
    // return cost as attribute
    ret.attr("cost") = optimum.costs();
  }
  
  return(ret);
}

/*************
* fitBandDynamicProgram
* fits the data with a dynamic program restricted to the confidence intervals, compute confidence bands
****************/
List fitBandDynamicProgram(Data * const data, IntervalSystem * const intervalSystem) {
  unsigned int n = data -> getN();
  std::vector<SingleBounds> intersectedBounds(n);
  
  // computation of left condidence intervals
  // intersectedBounds[i] will contain the intersected bounds from i to left.top() - 1
  std::stack<unsigned int> left;  // left confidence intervals starting with n, will be reordered
  left.push(n);
  
  unsigned int i = n;
  bool feasible = false;
  while (i > 0u) {
    --i;
    checkUserInterrupt();
    
    if (feasible) {
      intersectedBounds[i].intersect(intersectedBounds[i + 1u]);
    }
    
    data -> reset();
    for (unsigned int j = i; j < left.top(); ++j) {
      data -> addRight(j);
      
      if (intervalSystem -> isInIntervalSystem(i, j)) {
        intersectedBounds[i].intersect(data -> computeSingleBounds());
      }
    }    
    
    if (intersectedBounds[i].feasible()) {
      feasible = true;
    } else {
      feasible = false;
      left.push(i + 1u);
      intersectedBounds[i].reset();
      ++i;
    }
  }
  
  // reordering the left limits of the confidence intervals
  unsigned int K = left.size() - 1u;
  IntegerVector leftConfInt = IntegerVector(K);
  
  for (unsigned int k = 0u; k < K; ++k) {
    leftConfInt[k] = left.top();
    left.pop();
  }
  
  IntegerVector rightConfInt = IntegerVector(K);
  NumericVector lowerBand = NumericVector(n);
  NumericVector upperBand = NumericVector(n);
  std::vector<LocalOptimum> optimum(n);    // optimum[i]: current local solution to reach i
  // intersectedBounds[i] contains the intersected bounds from i to left.top() - 1
  
  if (K > 0u) {
    // segment [0, leftConfInt[0] - 1]
    data -> reset();
    for (int i = 0; i < leftConfInt[0u]; ++i) {
      data -> addRight(i);
      lowerBand[i] = intersectedBounds[0u].lower();
      upperBand[i] = intersectedBounds[0u].upper();
    }
    optimum[leftConfInt[0u] - 1] = data -> computeLocalOptimum(0, leftConfInt[0u] - 1,
                                                               intersectedBounds[0u], NULL);
    
    // first segment [0, j] until no constant solution is feasible
    for (unsigned int j = leftConfInt[0u]; j < n; ++j) {
      checkUserInterrupt();
      data -> reset();
      unsigned int i = j + 1u;
      while (i > 0u) {
        --i;
        data -> addLeft(i);
        if (intervalSystem -> isInIntervalSystem(i, j)) {
          intersectedBounds[0u].intersect(data -> computeSingleBounds());
        }
      }
      
      if (!intersectedBounds[0u].feasible()) {
        // j cannot reached without a change-point
        rightConfInt[0u] = j;
        break;
      }
      
      lowerBand[j] = std::min(intersectedBounds[0u].lower(), intersectedBounds[j].lower());
      upperBand[j] = std::max(intersectedBounds[0u].upper(), intersectedBounds[j].upper());
      
      optimum[j] = data -> computeLocalOptimum(0u, j, intersectedBounds[0u], NULL);
    }    
    
    // dynammic programming
    for (unsigned int k = 1u; k < K; ++k) {
      data -> reset();
      
      // step backward for right index leftConfInt[k] - 1
      for (unsigned int i = rightConfInt[k - 1u]; i < static_cast<unsigned int>(leftConfInt[k]); ++i) {
        data -> addRight(i);
        lowerBand[i] = intersectedBounds[rightConfInt[k - 1u]].lower();
        upperBand[i] = intersectedBounds[rightConfInt[k - 1u]].upper();
      }
      
      optimum[leftConfInt[k] - 1u] = data -> computeLocalOptimum(rightConfInt[k - 1u], leftConfInt[k] - 1u,
                                                                 intersectedBounds[rightConfInt[k - 1u]],
                                                                 &optimum[rightConfInt[k - 1u] - 1u]);
      
      for (unsigned int i = rightConfInt[k - 1u] - 1u; i >= static_cast<unsigned int>(leftConfInt[k - 1u]);
        --i) {
        data -> addLeft(i);
        optimum[leftConfInt[k] - 1u].update(data -> computeLocalOptimum(i, leftConfInt[k] - 1u,
                                                                        intersectedBounds[i], 
                                                                        &optimum[i - 1u]));
      }
      
      // step forward
      for (unsigned int j = leftConfInt[k]; j < n; ++j) {
        checkUserInterrupt();
        data -> reset();
        
        for (unsigned int i = j; i >= static_cast<unsigned int>(rightConfInt[k - 1u]); --i) {
          data -> addLeft(i);
          if (intervalSystem -> isInIntervalSystem(i, j)) {
            intersectedBounds[rightConfInt[k - 1u]].intersect(data -> computeSingleBounds());
          }
        }        
        
        if (!intersectedBounds[rightConfInt[k - 1u]].feasible()) {
          // no constant solution on [startPointRightIndex, ..., j] possible
          rightConfInt[k] = j;
          break;
        }
        
        lowerBand[j] = std::min(intersectedBounds[rightConfInt[k - 1u]].lower(),
                                intersectedBounds[j].lower());
        upperBand[j] = std::max(intersectedBounds[rightConfInt[k - 1u]].upper(),
                                intersectedBounds[j].upper());
        
        optimum[j] = data -> computeLocalOptimum(rightConfInt[k - 1u], j,
                                                 intersectedBounds[rightConfInt[k - 1u]],
                                                 &optimum[rightConfInt[k - 1u] - 1u]);
        
        // backward step
        for (unsigned int i = rightConfInt[k - 1u] - 1; i >= static_cast<unsigned int>(leftConfInt[k - 1u]);
          --i) {
          // find constant solution over [i, ..., j]
          intersectedBounds[i].intersect(intersectedBounds[i + 1u]); 
          data -> addLeft(i);
          if (intervalSystem -> isInIntervalSystem(i, j)) {
            intersectedBounds[i].intersect(data -> computeSingleBounds());
          }
          
          if (!intersectedBounds[i].feasible()) {
            // no constant solution on [i, ..., j] possible
            break;
          }
          
          optimum[j].update(data -> computeLocalOptimum(i, j, intersectedBounds[i], &optimum[i - 1u]));
        }
      }
    }
    
    // last segment
    data -> reset();
    for (unsigned int i = n - 1u; i >= static_cast<unsigned int>(rightConfInt[K - 1u]); --i) {
      data -> addLeft(i);
      lowerBand[i] = intersectedBounds[rightConfInt[K - 1u]].lower();
      upperBand[i] = intersectedBounds[rightConfInt[K - 1u]].upper();
    }    
    optimum[n - 1u] = data -> computeLocalOptimum(rightConfInt[K - 1u], n - 1u, 
                                                  intersectedBounds[rightConfInt[K - 1u]],
                                                  &optimum[rightConfInt[K - 1u] - 1u]);
    
    for (unsigned int i = rightConfInt[K - 1u] - 1u; i >= static_cast<unsigned int>(leftConfInt[K - 1u]);
      --i) {
      checkUserInterrupt();
      data -> addLeft(i);
      optimum[n - 1u].update(data -> computeLocalOptimum(i, n - 1u, intersectedBounds[i], &optimum[i - 1u]));
    }
  } else {
    for (unsigned int i = 0u; i < n; ++i) {
      lowerBand[i] = intersectedBounds[0u].lower();
      upperBand[i] = intersectedBounds[0u].upper();
    }
  }
  
  // return result
  IntegerVector leftIndex = IntegerVector(K + 1u);
  IntegerVector rightIndex = IntegerVector(K + 1u);
  NumericVector value = NumericVector(K + 1u);
  List ret = List::create(_["leftIndex"] = leftIndex, _["rightIndex"] = rightIndex, 
                          _["value"] = value,
                          _["leftConfInt"] = leftConfInt, _["rightConfInt"] = rightConfInt,
                          _["lowerBand"] = lowerBand, _["upperBand"] = upperBand);
  
  for (unsigned int i = 0; i < K; ++i) {
    ++leftConfInt[i];   // C to R-style
    ++rightConfInt[i];  // C to R-style
  }
  
  // backtracking
  if (K > 0u) {
    LocalOptimum localOptimum = optimum[n - 1u];  
    unsigned int k = K;
    do {
      value[k] = localOptimum.estimatedValue();
      leftIndex[k] = localOptimum.leftIndex() + 1;
      rightIndex[k] = localOptimum.rightIndex() + 1; // + 1 for R-style index
      
      if (k > 0u) {
        localOptimum = *localOptimum.lastSegment();
      } else {
        break;
      }
      --k;
    } while (true);
    // return cost as attribute
    ret.attr("cost") = optimum[n - 1u].costs();
  } else {
    LocalOptimum optimum = data -> computeLocalOptimum(0, n - 1u, intersectedBounds[0u], NULL);
    value[0] = optimum.estimatedValue();
    leftIndex[0] = optimum.leftIndex() + 1u;
    rightIndex[0] = optimum.rightIndex() + 1u;
    // return cost as attribute
    ret.attr("cost") = optimum.costs();
  }
  
  return(ret);
}
