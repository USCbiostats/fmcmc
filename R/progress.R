
#' Progress bar
#' 
#' A simple progress bar. This function is used in [MCMC] when the function has
#' been called with a single processor.
#' 
#' @param n Integer. Number of steps.
#' @param probs Double vector. Quantiles where to put marks
#' @param width Integer. Width of the bar in characters.
#' @param symbol Character. Single character symbol to print each bar.
#' @param ... Further arguments passed to [cat()] such as `file` and `append`.
#' @return A function that can be included at the interior of a loop to 
#' mark the progress of the loop. It receives a single argument, `i`,
#' which is the number of the current step.
#' @examples 
#' 
#' x <- new_progress_bar(20)
#' for (i in 1:20) {
#'   Sys.sleep(2/20)
#'   x(i)
#' }
#' 
#' @export
#' 
new_progress_bar <- function(
  n,
  probs  = c(0, .25, .5, .75, 1),
  width  = getOption("width", 80),
  symbol = "/",
  ...
  ) {
  
  
  # Figuring out the header
  where <- ceiling(stats::quantile(1:width, probs))
  
  line <- c(rep(" ", width), "\n")
  for (i in seq_along(where)) {
    
    # Compiling text
    
    if (i == length(where)) {
      line[where[i]] <- sprintf("%.0f%%|", probs[i]*100)
      line[(where[i] - 4):(where[i] - 1)] <- ""
    } else {
      line[where[i]] <- sprintf("|%.0f%%", probs[i]*100)
      line[(where[i] + 1):(where[i] + nchar(line[where[i]]) - 1 )] <- ""
      
    }
    
  }
  
  line <- paste(
    paste(line, collapse = ""),
    paste(rep("-", width), collapse = ""), sep=""
    )
  
  if (any(diff(where) < 4))
    stop("Some marks on the bar are overlapping.", call. = FALSE)

  # How many splits
  nsplits <- min(n, width)
  
  # What is the bar-width of each bar
  bwidth <- diff(floor(seq(from = 0, to = width, length.out = nsplits + 1)))
  
  # What is the number of iterations to count
  freq <- cumsum(diff(floor(seq(from = 0, to = n, length.out = nsplits + 1))))
  bar <- function(w) cat(rep(symbol, w), sep = "", ...)

  first_done <- FALSE
  head_position <- 1L
  function(i) {
    
    # If it is the first one, then start
    if (!first_done) {
      cat("\n", line, "\n", sep = "", ...)
      first_done <<- TRUE
    }
    
    # Intermediate steps are marked by this
    if (i == freq[head_position]) {
      bar(bwidth[head_position])
      head_position <<- head_position + 1
    }
    
    # if it is the last step, then write a new line
    if (i == freq[nsplits])
      cat("\n", ...)

  }
  
}
