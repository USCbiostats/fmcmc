f <- file()

p <- new_progress_bar(100, file = f, symbol = "/")
for (i in 1:100) {
  Sys.sleep(.5/100)
  p(i)
}

dat <- paste(readLines(f), collapse = "\n")
close(f)

expect_true(grepl("////////", dat, fixed = TRUE))
