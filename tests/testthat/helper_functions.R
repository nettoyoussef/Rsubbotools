library(usethis)

test_function <- function(b, func1, func2, msg_param) {
  lapply(1:length(b), function(i) {
    # there is numerical difference, so we use equal instead of identical
    # to use the numerical threshold tolerance
    expect_equal(func1, func2, label = paste0(msg_param, b[i]))
  })
}

test_function_rng_seed <- function(func, param_list) {
  x <- compare_RNG(func, param_list)
  expect_equal(x[["x"]], x[["y"]])
}

test_function_rng <- function(size, rng_function_1, rng_function_2) {
  lapply(seq_len(size), function(i) {
    set.seed(i)
    x <- rng_function_1
    set.seed(i)
    y <- rng_function_2
    expect_equal(x, y)
  })
}

compare_RNG <- function(rng_function, params_list) {
  set.seed(1)
  x <- do.call(rng_function, params_list)
  set.seed(1)
  y <- do.call(rng_function, params_list)

  return(list(x = x, y = y))
}

generate_datasets <- function(b) {
  set.seed(1)
  sample <- data.frame(x = rpower(10^6, b = b))
  return(sample)
}

trunc_digit <- function(x, digits) {
  x <- x * 10^digits
  x <- trunc(x)
  return(x / 10^digits)
}

trunc_c_format <- function(vector, format) {
  sapply(seq_along(vector), function(x) {
    vector[x] <<- as.numeric(sprintf(format, vector[x]))
  })

  return(vector)
}


check_subbo <- function(subbo_test, orig_value) {
  items <- names(subbo_test) # c("dt", "log-likelihood", "matrix")

  # since the Subbotools routine generally outputs %10.4g significant
  # digits, we have two options:
  # either we trunc the numbers on the R routine or we
  # use the c routine to generate the same approximation
  # I opted for the later
  # To guarantee the output (since subbo not always do the 10.4g pattern)
  # I pass both to the same filter
  subbo_test$dt["coef"] <- trunc_c_format(subbo_test$dt[["coef"]], "%10.4g")
  orig_value$dt["coef"] <- trunc_c_format(orig_value$dt[["coef"]], "%10.4g")

  if ("std_error" %in% colnames(subbo_test$dt)) {
    subbo_test$dt["std_error"] <-
      trunc_c_format(subbo_test$dt[["std_error"]], "%10.4g")
    orig_value$dt["std_error"] <-
      trunc_c_format(orig_value$dt[["std_error"]], "%10.4g")
  }

  subbo_test$`log-likelihood` <-
    trunc_c_format(subbo_test$`log-likelihood`, "%10.4g")
  orig_value$`log-likelihood` <-
    trunc_c_format(orig_value$`log-likelihood`, "%10.4g")

  if (!is.null(subbo_test$matrix)) {
    colnames_matrix <- colnames(subbo_test$matrix)
    col_matrix <- ncol(subbo_test$matrix)
    subbo_test$matrix <-
      matrix(trunc_c_format(as.vector(subbo_test$matrix), "%.4f"),
        ncol = col_matrix
      )
    colnames(subbo_test$matrix) <- colnames_matrix

    orig_value$matrix <-
      matrix(trunc_c_format(as.vector(orig_value$matrix), "%.4f"),
        ncol = col_matrix
      )
    colnames(orig_value$matrix) <- colnames_matrix
  }

  lapply(items, function(z) {
    print(z)
    expect_identical(
      subbo_test[z], orig_value[z],
      label = paste0("parameter ", z, "=", subbo_test[z]),
      expected.label = paste0("parameter ", z, "=", orig_value[z])
    )
  })
}

check_fits <- function(orig_value, b_param, fit_function) {
  data_test <- generate_datasets(b_param)
  subbo_test <- fit_function(data_test$x)

  if (!is.null(subbo_test$matrix)) {
    # since the original Subbotools package doesnt' return the
    # variances, we have to transform the principal axis of the Covariance
    # Matrix to NA
    size <- dim(subbo_test$matrix)
    sapply(1:size, function(x) {
      subbo_test$matrix[x, x] <<- NA
    })
  }

  check_subbo(subbo_test, orig_value)
}

generate_orig_dt <- function(
    coef,
    log_likelihood,
    std_error = NULL,
    matrix = NULL,
    distribution = "subbofit") {
  # for Symmetric Laplace
  if (length(coef) == 2 && distribution == "laplafit") {
    param <- c("m", "a")
  }

  # for Subbotin
  if (length(coef) == 3 && distribution == "subbofit") {
    param <- c("b", "a", "m")
  }

  # for Asymmetric Laplace
  if (length(coef) == 3 && distribution == "alaplafit") {
    param <- c("m", "al", "ar")
  }

  # for less Asymmetric Subbotin
  if (length(coef) == 4 && distribution == "subbolafit") {
    param <- c("bl", "br", "a", "m")
  }

  # for Skewed Exponential
  if (length(coef) == 4 && distribution == "sepfit") {
    param <- c("mu", "si", "la", "al")
  }

  # for Asymmetric Subbotin
  if (length(coef) == 5 && distribution == "subboafit") {
    param <- c("bl", "br", "al", "ar", "m")
  }

  if (!is.null(matrix)) {
    matrix <- t(matrix(matrix, ncol = length(coef)))
  }

  dt <-
    data.frame(
      param = param,
      coef = coef
    )

  if (!is.null(std_error)) {
    dt$std_error <- std_error
  }


  func_list <-
    list(
      dt = dt,
      "log-likelihood" = log_likelihood
    )

  if (!is.null(matrix)) {
    func_list$matrix <- matrix
  }

  return(func_list)
}

sample_datasets <- function() {
  # This function create sample datasets to be used for testing the functions
  # exported by the package. For a longer description, see the file
  # test-fits.R

  dir.create(paste0(getwd(), "/data"))
  path <- paste0(getwd(), "/data")

  sublaplace <- generate_datasets(b = .5)
  laplace <- generate_datasets(b = 1)
  subnormal <- generate_datasets(b = 1.5)
  normal <- generate_datasets(b = 2)
  supernormal <- generate_datasets(b = 2.5)

  datasets <- c(
    "sublaplace",
    "laplace",
    "subnormal",
    "normal",
    "supernormal"
  )

  lapply(datasets, function(x) {
    write.table(get(x), paste0(path, "/", x, ".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  })


  # saves the data to be available on the package
  # to be
  usethis::use_data(
    sublaplace,
    laplace,
    subnormal,
    normal,
    supernormal,
    overwrite = TRUE
    # saves the data to be available internally on sysdata.rda
    # ,internal = TRUE
  )
}

r_inc_lower_gamma <- function(p, b) {
  return(gamma(b) * pgamma(p, shape = b, scale = 1))
}

r_inv_inc_lower_gamma <- function(p, b) {
  return(qgamma(p / gamma(b), shape = b, scale = 1, lower.tail = 1, log.p = 0))
}
