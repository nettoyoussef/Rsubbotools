
generate_datasets <- function(b){
    set.seed(1)
    sample  <- data.frame(x=rpower(10^6, 2, b))
    return(sample)
}

trunc_digit <- function(x, digits){
    x <- x*10^digits
    x <- trunc(x)
    return(x/10^digits)
}

trunc_c_format <- function(vector, format){

    sapply(seq_along(vector), function(x){
        vector[x] <<- as.numeric(sprintf(format,vector[x]))
    })

    return(vector)
}


check_subbo <- function(subbo_test, orig_value){

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

    if("std_error" %in% colnames(subbo_test$dt)){
        subbo_test$dt["std_error"] <-
            trunc_c_format(subbo_test$dt[["std_error"]], "%10.4g")
        orig_value$dt["std_error"] <-
            trunc_c_format(orig_value$dt[["std_error"]], "%10.4g")
    }

    subbo_test$`log-likelihood` <-
        trunc_c_format(subbo_test$`log-likelihood`, "%10.4g")
    orig_value$`log-likelihood` <-
        trunc_c_format(orig_value$`log-likelihood`, "%10.4g")

    if(!is.null(subbo_test$matrix)){
        colnames_matrix <- colnames(subbo_test$matrix)
        col_matrix <- ncol(subbo_test$matrix)
        subbo_test$matrix <-
            matrix(trunc_c_format(as.vector(subbo_test$matrix), "%.4f")
                  ,ncol = col_matrix)
        colnames(subbo_test$matrix) <- colnames_matrix

        orig_value$matrix <-
            matrix(trunc_c_format(as.vector(orig_value$matrix), "%.4f")
                  ,ncol = col_matrix)
        colnames(orig_value$matrix) <- colnames_matrix
    }

    lapply(items, function(z){
        print(z)
        expect_identical(
            subbo_test[z], orig_value[z]
           ,label = paste0("parameter ", z, "=", subbo_test[z])
           ,expected.label = paste0("parameter ", z, "=", orig_value[z])
        )
    })
}

check_fits <- function(orig_value, b_param, fit_function){

    data_test   <- generate_datasets(b_param)
    subbo_test <- fit_function(data_test$x)

    if(!is.null(subbo_test$matrix)){
        # since the original Subbotools package doesnt' return the
        # variances, we have to transform the principal axis of the Covariance
        # Matrix to NA
        size <- dim(subbo_test$matrix)
        sapply(1:size, function(x){
            subbo_test$matrix[x, x] <<- NA
        })
    }

    check_subbo(subbo_test, orig_value)

}

generate_orig_dt <- function(
                             coef
                           , log_likelihood
                           , std_error    = NULL
                           , matrix       = NULL
                           , distribution = "subbofit"
                             ){

    # for Symmetric Laplace
    if(length(coef) == 2 && distribution == "laplafit"){
        param <-  c("m", "a")

    }

    # for Subbotin
    if(length(coef) == 3 && distribution == "subbofit"){
        param <-  c("b", "a", "m")

    }

    # for Asymmetric Laplace
    if(length(coef) == 3 && distribution == "alaplafit"){
        param <-  c("m", "al", "ar")
    }

    # for less Asymmetric Subbotin
    if(length(coef) == 4 && distribution == "subbolafit"){
        param <-  c("bl", "br", "a", "m")
    }

    # for Skewed Exponential
    if(length(coef) == 4 && distribution == "sepfit"){
        param <-  c("mu", "si", "la", "al")
    }

    # for Asymmetric Subbotin
    if(length(coef) == 5 && distribution == "subboafit"){
        param <-  c("bl", "br", "al", "ar", "m")
    }

    if(!is.null(matrix)){
        matrix <- t(matrix(matrix, ncol = length(coef)))
    }

    dt <-
        data.frame(
            param     = param
           ,coef      = coef
        )

    if(!is.null(std_error)){
        dt$std_error <- std_error
    }


    func_list <-
        list(
            dt = dt
           ,'log-likelihood' = log_likelihood
        )

    if(!is.null(matrix)){
        func_list$matrix <- matrix
    }

    return(func_list)
}

