
if ( requireNamespace("tinytest", quietly=TRUE) ){
  
  expect_lt <- function(a,b) {
    
    tinytest::tinytest(
      result = a < b,
      call   = sys.call(sys.parent(1)),
      diff   = paste("b - a = ", b - a)
    )
    
  }
  
  expect_output <- function(x, pattern) {
    
    ans  <- paste(capture.output(eval(x)), collapse = "\n")
    test <- grepl(pattern, ans)
    
    tinytest::tinytest(
      result = test,
      call   = sys.call(sys.parent(1)),
      diff   = paste("Expected:", pattern, ". Returned:", ans)
    )
  }
  
  tinytest::test_package("fmcmc")
}

