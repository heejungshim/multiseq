test_that("multiseq gives sensible answer for effect in simple settings",{ 
  set.seed(100)
  lambda = c(1,10,100)
  effect = c(2,4,8)
  for(i in 1:length(lambda)){
    xx = rbind(rpois(1024,rep(lambda[i],1024)),rpois(1024,rep(lambda[i],1024)),rpois(1024,rep(lambda[i],1024)))
    yy = rbind(rpois(1024,rep(effect[i]*lambda[i],1024)),rpois(1024,rep(effect[i]*lambda[i],1024)),rpois(1024,rep(effect[i]*lambda[i],1024)))
    yy.m = multiseq(rbind(xx,yy),g=c(0,0,0,1,1,1))
    testthat::expect_less_than(mean((yy.m$effect.mean-log(effect[i]))^2), 0.01)
  }
})



