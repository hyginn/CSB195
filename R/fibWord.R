# fibWord.R
#
# Construct a Fibonacci Word of a given length
#
# cf. https://en.wikipedia.org/wiki/Fibonacci_word

fibWord <- function(N = 8 * 13, mode = "str") {

  s <- integer(N+2)
  i <- 1
  iLast <- i+1

  while (iLast < N+1) {
    if (s[i] == 1) { s[iLast]             <- 0;   iLast <- iLast+1
    } else         { s[iLast:(iLast + 1)] <- 1:0; iLast <- iLast+2
    }
    i <- i+1
  }
  if (mode == "str") {                # string
    s <- paste(s[1:N], collapse = "")
  } else if (mode == "chr") {         # character
    s <- as.character(s[1:N])
  } else {                            # numeric
    s <- s[1:N]
  }

  return(s)
}

if (FALSE) {

  fibWord(1)
  fibWord(3)
  fibWord(5)
  fibWord(8)
  fibWord(13)
  fibWord(21)
  fibWord(1000)

  s <- fibWord()
  iStart <- ((0:7)*13)+1
  iStop  <- iStart+12
  cat(sprintf("\n%s", substring(s, iStart, iStop)))

}


# == fibNumbers =====
fibNumbers <- function(N) {
  if (N < 1)  { return(NaN) }
  if (N == 1) { return(1) }
  if (N == 2) { return(1) }
  v <- c(1,1, integer(N-2))
  for (i in seq_len(N)[-(1:2)]) {
    num <- v[i-2] + v[i-1]
    if (num > .Machine$integer.max) {
      stop("Result exceeds .Machine$integer.max for N > 46.")
    } else {
      v[i] <- num
    }
  }
  return(v)
}

if (FALSE) {
  i <- -1
  fibNumbers((i <- i+1))
}


# [END]
