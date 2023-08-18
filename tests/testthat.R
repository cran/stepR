library(testthat)
library(stepR)

savePathRcache <- R.cache::getCacheRootPath()
R.cache::setCacheRootPath(path = file.path(R.cache::getCacheRootPath(), "test"))

tryCatch(
  test_check("stepR"),
  finally = {
    unlink(R.cache::getCacheRootPath(), force = TRUE, recursive = TRUE)
    R.cache::setCacheRootPath(savePathRcache)
  }
)


