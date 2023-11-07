#' clusters = c(0, 5, 100, 95) ===>>> beta = c(0, 0, 0, 0, 0, 100, 100, ..., 100)
#' @export
#' @examples
#' get_data(200, 100, c(0, 5, 100, 95))
#' l <- get_data(200, 10, c(100, 5, 0, 5))
#' l <- get_data(200, 15, c(0, 5, -100, 5, 100, 5), sigma = 0, CovarMatrix = "I")
#'
#' X <- l$X
#' Y <- l$Y
#' n <- l$n
#' p <- l$p
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

  return(list("X" = X, "Y" = Y, "C" = C, "beta" = beta, "n" = n, "p" = p))
}

#' @export
get_lambda <- function(type = c("B", "K", "I", "InvI"), X, Y, p, q, n) {
  if (type == "K") {
    lambda <- sort(abs(t(X) %*% Y), decreasing = TRUE)
  } else if (type == "I"){ # TODO(Uprość)
    b_OLS <- solve(t(X) %*% X) %*% t(X) %*% Y

    C <- t(X) %*% X / n

    losowy_pomysl <- TRUE
    if (losowy_pomysl){
      b_OLS <- abs(b_OLS)
    }

    # Tu trzeba pomnożyć macierzą, gdzie każdy wiersz to b_OLS:
    B <- matrix(rep(b_OLS, p), nrow = p, byrow = TRUE)

    C_falka <- t(B) %*% C %*% B

    lambda <- sort(rowSums(C_falka), decreasing = TRUE)

    # lambda może być ujemna teraz. Trzba, żeby nie:
    lambda <- ifelse(lambda > 0, lambda, 0)
  } else if (type == "InvI"){
    lambdaI <- get_lambda(type = "I", X, Y, p, q, n)
    lambda <- sort(1/lambdaI, decreasing = TRUE)
  } else {
    seq1 <- seq(1:p)
    lambda <- stats::qnorm(1 - q * seq1 / 2 / p)
  }

  lambda
}
