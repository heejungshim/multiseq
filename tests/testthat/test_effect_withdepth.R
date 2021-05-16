test_that("multiseq gives sensible answer for effect in simple settings, with readdepth",{ 
  set.seed(100)
  lambda = 1
  depth=c(1,1,1,100,100,100)
  x = rbind(rpois(1024,lambda*depth[1]),rpois(1024,lambda*depth[2]),rpois(1024,lambda*depth[3]))
  y = rbind(rpois(1024,lambda*depth[4]),rpois(1024,lambda*depth[5]),rpois(1024,lambda*depth[6]))
  x.m = multiseq(rbind(x,y),read.depth=depth*100*1024,g=c(0,0,0,1,1,1))
  x.m2 = multiseq(rbind(x,y),read.depth=depth*10*1024,g=c(0,0,0,1,1,1))
  testthat::expect_less_than(mean((x.m$effect.mean^2)), 0.05)
  testthat::expect_equal(mean((x.m$effect.mean)), mean((x.m2$effect.mean)),tolerance=1e-2)
  }
)
