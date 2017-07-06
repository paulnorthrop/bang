# Function to test whether the values in x are whole numbers

is_wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  return(abs(x - round(x)) < tol)
}
