.parametricFamily <- function(family, y = NULL, n, nq = 2L^(as.integer(log2(n) + 1e-12) + 1L) - 1L, ...) {
  if (is.null(family)) {
    family <- "gauss"
  }
  family <- match.arg(family, c("gauss", "mDependentPS", "hsmuce"))
  
  if (!is.numeric(n) || length(n) != 1 || !is.finite(n)) {
    stop("number of observations 'n' must be a single positive integer")
  }
  
  if (!is.integer(n)) {
    n <- as.integer(n + 1e-6)
  }
  
  if (n < 1L) {
    stop("number of observations 'n' must be a single positive integer")
  }
  
  if (!is.numeric(nq) || length(nq) != 1 || !is.finite(nq)) {
    stop("nq must be a single integer greather or equal than the number of observations 'n'")
  }
  
  if (!is.integer(nq)) {
    nq <- as.integer(nq + 1e-6)
  }
  
  if (nq < n) {
    stop("nq must be a single finite integer greather or equal than the number of observations 'n'")
  }
  
  if (is.null(y)) {
    data <- list(family = family, n = n, nq = nq)
  } else {
    data <- list(y = NULL, family = family, n = n, nq = nq)
  }
  
  switch(family,
         "gauss" = .familyGauss(data = data, y = y, ...),
         "mDependentPS" = .familyMDependentPS(data = data, y = y, ...),
         "hsmuce" = .familyHsmuce(data = data, y = y, ...),
         stop("unknown family"))
}

.familyGauss <- function(data, y, sd = NULL) {
  if (!is.null(y)) {
    if (!is.double(y) || any(!is.finite(y))) {
      stop("for family 'gauss' y must be a double vector containing only finite values")
    }
    
    if (length(y) != data$n) {
      stop("y must be of length n")
    }
    
    data$y <- y
  }
  
  if (is.null(sd)) {
    if (is.null(y)) {
      data$sd <- 1
    } else {
      data$sd <- sdrobnorm(y, supressWarningResultNA = TRUE)
      if (is.na(data$sd)) {
        stop("number of observations is too small for computing 'sd' by 'sdrobnorm'")
      }
    }
  } else {
    if (!is.numeric(sd) || length(sd) != 1 || !is.finite(sd) || sd <= 0) {
      stop("sd must be a single positive finite numeric")
    }
    
    data$sd <- as.numeric(sd)
  }
  data$type <- 0L
  data$argumentsList <- list(sd = data$sd)
  data$possibleLengths <- 1:data$n
  data$defaultLengths <- 1:data$n
  if (data$n <= 1e3L) {
    data$defaultIntervalSystem <- "all"
  } else {
    data$defaultIntervalSystem <- "dyaLen"
  }
  data$defaultPenalty <- "sqrt"
  
  data$key <- digest::sha1(list("gauss"), digits = 6)
  data$rand.gen <- function(data) rnorm(data$n)
  data$testGeneratedData <- function(data) {
    !is.numeric(data$y) || length(data$y) != data$n || any(!is.finite(data$y))
  }
  data$errorMessageGeneratedData <- paste("rand.gen must be a function with a single",
                                          "argument data returning for parametric family 'gauss'",
                                          "a finite numeric vector",
                                          "of length equal to the number of observations",
                                          sep = " ")
  data$MC <- function(data, n) {
    data$n <- n
    data$nq <- n
    data$possibleLengths <- 1:n
    data$defaultLengths <- 1:n
    if (n <= 1e3L) {
      data$defaultIntervalSystem <- "all"
    } else {
      data$defaultIntervalSystem <- "dyaLen"
    }
    data$sd <- 1
    data$argumentsList <- list(sd = data$sd)
    data
  }

  data
}

.familyMDependentPS <- function(data, y, covariances = NULL, correlations = NULL, filter = NULL, sd = NULL) {
  if (!is.null(y)) {
    if (!is.double(y) || !all(is.finite(y))) {
      stop("for family 'mDependentPS' y must be a double vector containing only finite values")
    }
    
    if (length(y) != data$n) {
      stop("y must be of length n")
    }
    
    data$y <- y
  }
  
  if (is.null(covariances)) {
    if (is.null(correlations)) {
      if (is.null(filter)) {
        stop("covariances, correlations or filter must be given for family 'mDependentPS'")
      } else {
        if (!is(filter, "lowpassFilter")) {
          stop("filter must be an object of class 'lowpassFilter'")
        }
        
        correlations <- filter$acf
      }
    } else {
      if (!is.numeric(correlations) || !all(is.finite(correlations))) {
        stop("correlations must be a finite numeric vector")
      }
      
      if (correlations[1] != 1 || any(abs(correlations[-1]) > 1) ||
          correlations[length(correlations)] == 0) {
        stop("correlations must be a correlation vector, ",
             "i.e. the first element must be 1, the absolute value of every other element must be ",
             "smaller or equal and the last element should not be zero")
      }
    }
    
    if (is.null(sd)) {
      if (is.null(y)) {
        sd <- 1
      } else {
        sd <- sdrobnorm(y, lag = length(correlations), supressWarningResultNA = TRUE)
        if (is.na(sd)) {
          stop("number of observations is too small for computing 'sd' by 'sdrobnorm'")
        }
      }
    } else {
      if (!is.numeric(sd) || length(sd) != 1 || !is.finite(sd) || sd <= 0) {
        stop("sd must be a single positive finite numeric")
      }
      
      sd <- as.numeric(sd)
    }
    
    covariances <- sd^2 * correlations
  } else {
    if (!is.numeric(covariances)  || !all(is.finite(covariances))) {
      stop("covariances has to be a finite numeric vector")
    }
    
    if (covariances[1] <= 0 || any(abs(covariances[-1]) > covariances[1]) ||
        covariances[length(covariances)] == 0) {
      stop("covariances has to be a covariance vector, ",
           "i.e. the first element has to be positive, the absolute value of every other element has to be ",
           "smaller or equal and the last element should not be zero")
    }
  }
  data$type <- 10L
  data$covariances <- covariances
  data$argumentsList <- list(covariances = data$covariances)
  data$possibleLengths <- 1:data$n
  data$defaultLengths <- 1:data$n
  data$defaultIntervalSystem <- "dyaLen"
  data$defaultPenalty <- "sqrt"
  data$ma <- .computeMA(data$covariances)
  data$rand.gen <- .rand.genMDependent
  data$testGeneratedData = function(data) {
    !is.numeric(data$y) || length(data$y) != data$n || !all(is.finite(data$y))
  }
  data$errorMessageGeneratedData = paste("rand.gen must be a function with a single", 
                                         "argument data returning for parametric family",
                                         "'mDependentPS' a finite numeric vector",
                                         "of length equal to the number of observations",
                                         sep = " ")
  data$key <- digest::sha1(list("mDependentPS", data$ma), digits = 6)
  
  data$MC <- function(data, n) {
    data$n <- n
    data$nq <- n
    data$possibleLengths <- 1:n
    data$defaultLengths <- 1:n
    data$covariances <- data$covariances / data$covariances[1]
    data$argumentsList <- list(covariances = data$covariances)
    data
  }

  data
}

.familyHsmuce <- function(data, y) {
  if (!is.null(y)) {
    if (!is.double(y) || !all(is.finite(y))) {
      stop("for family 'hsmuce' y must be a double vector containing only finite values")
    }
    
    if (length(y) != data$n) {
      stop("y must be of length n")
    }
    
    data$y <- y
  }
  
  data$type <- 20L
  data$argumentsList <- NULL
  data$possibleLengths <- 2:data$n
  data$defaultLengths <- 2:data$n
  data$defaultIntervalSystem <- "dyaPar"
  data$defaultPenalty <- "weights"
  data$rand.gen <- function(data) rnorm(data$n)
  data$testGeneratedData <- function(data) {
    !is.numeric(data$y) || length(data$y) != data$n || !all(is.finite(data$y))
  }
  data$errorMessageGeneratedData <- paste("rand.gen must be a function with a single",
                                          "argument data returning for parametric family",
                                          "'hsmuce' a finite numeric vector",
                                          "of length equal to the number of observations",
                                          sep = " ")
  data$key <- digest::sha1(list("hsmuce"), digits = 6)
  
  data$MC <- function(data, n) {
    data$n <- n
    data$nq <- n
    data$possibleLengths <- 2:n
    data$defaultLengths <- 2:n
    data
  }

  data
}

.rand.genMDependent <- function(data) {
  kern <- c(1, data$ma)
  z <- rnorm(data$n + length(kern) - 1, sd = 1)
  .convolve(z, kern) / sqrt(sum(kern^2))
}

.computeMA <- function(cov) {
  N <- length(cov) - 1
  alpha <- matrix(rep(0, N * N), nrow = N, ncol = N)
  vs <- rep(cov[1], N + 1)
  
  for(i in 1:N){
    alpha[i, 1:i] <- rep(0, i)
    alpha[i, i] <- cov[i + 1] / vs[1]
    if(i > 1) {
      for(k in 1:(i - 1)) {
        js <- 0:(k - 1)
        alpha[i, i - k] <- (cov[i - k + 1] - sum(alpha[i, i - js] * alpha[k, k - js] * vs[js + 1])) /
          vs[k + 1]
      }
    }
    js <- 0:(i - 1)
    vs[i + 1] <- vs[i + 1] - sum(alpha[i, i - js]^2 * vs[js + 1])
  }
  
  alpha[N, ]
}
