.onUnload <- function(libpath) {
  library.dynam.unload("Rsubbotools", libpath)
}

newton <- function(x, f, dfdx = NULL, dfdx_h = 1e-4, tol = 1e-4, max_iter = 100) {
  # numerical derivative
  df <- function(x, h) {
    (f(x + h) - f(x)) / h
  }

  # if analitical derivative is not provided,
  # calculates one numerically with step
  # dfdx_h
  if (is.null(dfdx)) {
    dfdx <- function(x) {
      df(x, h = dfdx_h)
    }
  }

  print(paste0("iteration: 0; x value: ", x, "; f(x) value: ", round(f(x), 4)))

  # starts iteration
  iter <- 1
  x_vec <- numeric(max_iter)

  # initial value of the function
  f_x <- f(x)
  x_vec[iter] <- x - f_x / dfdx(x)
  f_x <- f(x_vec[iter])
  print(paste0("iteration: ", iter, "; x value: ", round(x_vec[iter], 4), "; f(x) value: ", round(f_x, 4)))
  iter <- iter + 1

  # checks convergence of the algorithm
  # (receives value after first iteration)
  check_convergence <- f_x
  count_warning <- 0


  while (abs(f_x) > tol && iter <= max_iter) {
    x_vec[iter] <- x_vec[iter - 1] - f_x / dfdx(x_vec[iter - 1])
    f_x <- f(x_vec[iter])
    print(paste0("iteration: ", iter, "; x value: ", round(x_vec[iter], 4), "; f(x) value: ", round(f_x, 4)))
    iter <- iter + 1
    if (is.na(f_x)) break

    # saves value from 5 interactions behind, to check convergence
    if (iter %% 5 == 0) {
      check_convergence <- f_x
    }

    if (abs(f_x) > abs(check_convergence)) {
      warning("Algorithm is not converging.")
      count_warning <- count_warning + 1
    }

    if (count_warning == 3) {
      stop("Algorithm is not converging.")
    }
  }

  return(list(x = x_vec[1:(iter - 1)], fx = f_x))
}



steffensen <- function(x, f, dfdx = NULL, dfdx_h = 1e-4, tol = 1e-4, max_iter = 100) {
  # numerical derivative
  df <- function(x, h) {
    (f(x + h) - f(x)) / h
  }

  # if analitical derivative is not provided,
  # calculates one numerically with step
  # dfdx_h
  if (is.null(dfdx)) {
    dfdx <- function(x) {
      df(x, h = dfdx_h)
    }
  }

  print(paste0("iteration: 0; x value: ", x, "; f(x) value: ", round(f(x))))

  # starts iteration
  iter <- 1
  x_vec <- numeric(max_iter)

  # initial value of the function
  f_x <- f(x)


  # first Newton iteration
  x_vec[iter] <- x - f_x / dfdx(x)
  f_x <- f(x_vec[iter])
  print(paste0("iteration: ", iter, "; x value: ", round(x_vec[iter], 4), "; f(x) value: ", round(f_x, 4)))
  iter <- iter + 1

  # checks convergence of the algorithm
  # (receives value after first iteration)
  check_convergence <- f_x
  count_warning <- 0

  while (abs(f_x) > tol && iter <= max_iter) {
    # calculates Newton value
    x <- x_vec[iter - 1] - f_x / dfdx(x_vec[iter - 1])

    if (iter < 3) {
      # starts with Newton's method
      x_vec[iter] <- x
    } else {
      # This is Steffensen's Method

      # some useful dereference
      xlag_1 <- x_vec[iter - 1]
      xlag_2 <- x_vec[iter - 2]
      denom <- (x - 2 * xlag_1 + xlag_2)

      # if denominator is null, uses Newton method
      if (denom == 0) {
        x_vec[iter] <- x
      } else {
        # uses accelerated method
        x_vec[iter] <- xlag_2 - (xlag_1 - xlag_2)^2 / denom
      }
    }

    f_x <- f(x_vec[iter])
    print(paste0("iteration: ", iter, "; x value: ", round(x_vec[iter], 4), "; f(x) value: ", round(f_x, 4)))
    iter <- iter + 1
    if (is.na(f_x)) break

    # saves value from 5 interactions behind, to check convergence
    if (iter %% 5 == 0) {
      check_convergence <- f_x
    }

    if (abs(f_x) > abs(check_convergence)) {
      warning("Algorithm is not converging.")
      count_warning <- count_warning + 1
    }

    if (count_warning == 3) {
      stop("Algorithm is not converging.")
    }
  }

  return(list(x = x_vec[1:(iter - 1)], fx = f_x))
}


# generates Subbotin Matrix of Covariances
subbo_covar_r <- function(subbo_test, N) {
  b <- subbo_test$dt[1, "coef"]
  a <- subbo_test$dt[2, "coef"]

  tmp <- log(b) + digamma(1 + 1. / b)
  dim <- dim(subbo_test$matrix)[[1]]
  I <- matrix(rep(0, 9), ncol = 3)

  #* b - b */
  I[1, 1] <- (tmp * tmp + trigamma(1 + 1 / b) * (1 + 1 / b) - 1) / b^3
  #* b - a */
  I[1, 2] <- -tmp / (a * b)
  I[2, 1] <- -tmp / (a * b)

  #* a - a */
  I[2, 2] <- b / (a * a)

  if (dim == 3) {
    #* b - m */
    I[1, 3] <- 0
    I[3, 1] <- 0

    #* a - m */
    I[2, 3] <- 0
    I[3, 2] <- 0

    #* m - m */
    I[3, 3] <- (b^(-2 / b + 1) * gamma(2 - 1 / b)) /
      (gamma(1 + 1 / b) * a * a)
  }

  # generates inverse
  I <- solve(I)
  # this is what Bottazzi does with GSL
  # library("Matrix")
  # I <- expand(lu(I))
  # I <- solve(I$U) %*% solve(I$L) %*% I$P

  #* set the var-covar matrix */
  for (i in 2:dim) {
    for (j in 1:(i - 1)) {
      I[i, j] <- I[i, j] / sqrt(I[i, i] * I[j, j])
    }
  }

  for (i in 1:dim) {
    for (j in i:dim) {
      I[i, j] <- I[i, j] / N
    }
  }
  return(I)
}

# generates Subbotin Matrix of Covariances
subboa_covar_r <- function(subbo_test, N) {
  B0 <- function(x) {
    return(x^(1. / x) * gamma(1. + 1. / x))
  }

  B1 <- function(x) {
    tmp <- 1. + 1. / x
    return(x^(1. / x - 1.) * gamma(tmp) * (log(x) + digamma(tmp)))
  }

  B2 <- function(x) {
    dt1 <- 1. + 1. / x
    dt2 <- log(x)
    dt3 <- digamma(dt1)
    dt4 <- trigamma(dt1)

    return(x^(1. / x - 2.) * gamma(dt1) * (dt2 * dt2 + 2 * dt2 * dt3 + dt3 * dt3 + dt4))
  }

  dB0dx <- function(x) {
    return(B0(x) / (x * x) - B1(x) / x)
  }

  dB0dx2 <- function(x) {
    dt1 <- x * x

    return(-B0(x) * (3. * x - 1.) / (dt1 * dt1) + 2. * B1(x) * (x - 1.) / (dt1 * x) + B2(x) / dt1)
  }

  M_EULER <- 0.57721566490153286060651209008
  bl <- subbo_test$dt[1, "coef"]
  br <- subbo_test$dt[2, "coef"]
  al <- subbo_test$dt[3, "coef"]
  ar <- subbo_test$dt[4, "coef"]

  A <- al * B0(bl) + ar * B0(br)
  B0l <- B0(bl)
  B0r <- B0(br)
  B1l <- B1(bl)
  B1r <- B1(br)
  B2l <- B2(bl)
  B2r <- B2(br)
  dB0ldx <- dB0dx(bl)
  dB0rdx <- dB0dx(br)
  dB0ldx2 <- dB0dx2(bl)
  dB0rdx2 <- dB0dx2(br)

  dim <- dim(subbo_test$matrix)[[1]]
  I <- matrix(rep(0, dim[1]^2), ncol = dim[1])

  # bl - bl
  I[1, 1] <- al * (dB0ldx2 - al * dB0ldx * dB0ldx / A +
    B2l / bl - 2 * B1l / (bl * bl) + 2 * B0l / bl^(3)) / A

  # bl - br
  I[1, 2] <- I[2, 1] <- -al * ar * dB0ldx * dB0rdx / (A * A)

  # bl - al
  I[1, 3] <- I[3, 1] <- dB0ldx / A - al * B0l * dB0ldx / (A * A) - B1l / A

  # bl - ar
  I[1, 4] <- I[4, 1] <- -al * B0r * dB0ldx / (A * A)

  # br - br
  I[2, 2] <- ar * (dB0rdx2 - ar * dB0rdx * dB0rdx / A +
    B2r / br - 2 * B1r / (br * br) + 2 * B0r / br^(3)) / A

  # br - al
  I[2, 3] <- I[3, 2] <- -ar * B0l * dB0rdx / (A * A)

  # br - ar
  I[2, 4] <- I[4, 2] <- dB0rdx / A - ar * B0r * dB0rdx / (A * A) - B1r / A

  # al - al
  I[3, 3] <- -B0l * B0l / (A * A) + (bl + 1) * B0l / (al * A)

  # al - ar
  I[3, 4] <- I[4, 3] <- -B0l * B0r / (A * A)

  # ar - ar
  I[4, 4] <- -B0r * B0r / (A * A) + (br + 1) * B0r / (ar * A)

  if (dim == 5) {
    dt1l <- gamma(2. - 1 / bl)
    dt1r <- gamma(2. - 1 / br)
    dt2l <- bl^(1. - 1. / bl)
    dt2r <- br^(1. - 1. / br)

    # bl - m
    I[1, 5] <- I[5, 1] <- (log(bl) - M_EULER) / (bl * A)

    # br - m
    I[2, 5] <- I[5, 2] <- -(log(br) - M_EULER) / (br * A)

    # al - m
    I[3, 5] <- I[5, 3] <- -bl / (al * A)

    # ar - m
    I[4, 5] <- I[5, 4] <- br / (ar * A)

    # m - m
    I[5, 5] <- (dt1l * dt2l / al + dt1r * dt2r / ar) / A
  }
  print(I)
  print(paste0(I[1, 1]))

  # generates inverse
  I <- solve(I)
  print(I)
  # this is what Bottazzi does with GSL
  # library("Matrix")
  # I <- expand(lu(I))
  # I <- solve(I$U) %*% solve(I$L) %*% I$P

  #* set the var-covar matrix */
  for (i in 2:dim) {
    for (j in 1:(i - 1)) {
      I[i, j] <- I[i, j] / sqrt(I[i, i] * I[j, j])
    }
  }
  print(I)

  for (i in 1:dim) {
    for (j in i:dim) {
      I[i, j] <- I[i, j] / N
    }
  }
  print(I)

  return(I)
}

