#' clusters c(0, 5, 100, 95) ===>>> beta = c(0, 0, 0, 0, 0, 100, 100, ..., 100)
#' @export
#' @examples
#' get_data(200, 100, c(0, 5, 100, 95))
get_data <- function(n, p, clusters, sigma = 0, CovarMatrix = "I") {
  stopifnot(length(clusters) %% 2 == 0)
  if (is.character(CovarMatrix)) {
    if (CovarMatrix == "I") {
      CovarMatrix <- diag(p)
    } else if (CovarMatrix == "BlokDiag") {
      stopifnot(p %% 10 == 0)

      kor <- 0.8
      Covar1 <- matrix(kor, 10, 10)
      diag(Covar1) <- 1
      C1 <- diag(p / 10)
      CovarMatrix <- kronecker(C1, Covar1)
      diag(CovarMatrix) <- 1
    } else {
      stop('CovarMatrix jako napis musi byc "I" albo "BlokDiag".')
    }
  }
  CovarMatrixChol <- t(chol(CovarMatrix))

  X <- matrix(stats::rnorm(n * p), n, p) %*% t(CovarMatrixChol)
  eps <- sigma * matrix(stats::rnorm(n), n, 1)

  beta <- rep(0, p)
  j <- 1
  for (i in 1:(length(clusters) / 2)) {
    beta[j:(j + clusters[i * 2] - 1)] <- rep(clusters[i * 2 - 1], clusters[i * 2])
    j <- j + clusters[i * 2]
  }
  stopifnot(j == p + 1)

  Y <- X %*% beta + eps
  C <- t(X) %*% X / n

  return(list("X" = X, "Y" = Y, "C" = C, "beta" = beta))
}

get_lambda <- function(type = c("B", "K"), X, Y, p, q) {
  if (type == "K") {
    lambda <- sort(abs(t(X) %*% Y), decreasing = TRUE)
  } else {
    seq1 <- seq(1:p)
    lambda <- stats::qnorm(1 - q * seq1 / 2 / p)
  }

  lambda
}
