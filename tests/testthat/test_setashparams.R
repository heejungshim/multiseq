test_that("default ash parameters are correct", {
  temp = readRDS("default.ashparam.RDS")
  expect_equivalent(setAshParam(list()),temp)
})
