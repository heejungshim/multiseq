test_that("multiseq gives sensible answer for baseline in simple settings",{ 
  set.seed(100)
  lambda = c(1,10,100)
  x = list()
  x.m=list()
  for(i in 1:length(lambda)){
    x[[i]] = rpois(1024,rep(lambda[i],1024))
    x.m[[i]] = multiseq(x[[i]])
    testthat::expect_less_than(mean((x.m[[i]]$baseline.mean-log(lambda[i]))^2), 0.01)
  }
})
