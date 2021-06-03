.onUnload <- function (libpath) {
  library.dynam.unload("Rsubbotools", libpath)
}

newton <- function(x, f, dfdx = NULL, dfdx_h = 1e-4, tol = 1e-4, max_iter = 100){

    # numerical derivative
    df <- function(x, h){
        (f(x+h)-f(x))/h
    }

    # if analitical derivative is not provided,
    # calculates one numerically with step
    # dfdx_h
    if(is.null(dfdx)){
        dfdx <- function(x){
            df(x, h = dfdx_h)
        }
    }

    print(paste0("iteration: 0; x value: ", x, "; f(x) value: ", round(f(x), 4)))

    # starts iteration
    iter = 1
    x_vec <- numeric(max_iter)

    # initial value of the function
    f_x <- f(x)
    x_vec[iter] <- x - f_x/dfdx(x)
    f_x <- f(x_vec[iter])
    print(paste0("iteration: ", iter, "; x value: ", round(x_vec[iter],4), "; f(x) value: ", round(f_x,4)))
    iter = iter +1

    # checks convergence of the algorithm
    # (receives value after first iteration)
    check_convergence <- f_x
    count_warning <- 0


    while(abs(f_x)>tol && iter <= max_iter){
        x_vec[iter] <- x_vec[iter-1] - f_x/dfdx(x_vec[iter-1])
        f_x <- f(x_vec[iter])
        print(paste0("iteration: ", iter, "; x value: ", round(x_vec[iter],4), "; f(x) value: ", round(f_x,4)))
        iter = iter +1
        if(is.na(f_x)) break

        # saves value from 5 interactions behind, to check convergence
        if(iter%%5 == 0){
            check_convergence <- f_x

        }

        if(abs(f_x) > abs(check_convergence)){
            warning("Algorithm is not converging.")
            count_warning = count_warning + 1
        }

        if(count_warning == 3){
            stop("Algorithm is not converging.")
        }


    }

    return(list(x=x_vec[1:(iter-1)], fx = f_x))
}



steffensen <- function(x, f, dfdx = NULL, dfdx_h = 1e-4, tol = 1e-4, max_iter = 100){

    # numerical derivative
    df <- function(x, h){
        (f(x+h)-f(x))/h
    }

    # if analitical derivative is not provided,
    # calculates one numerically with step
    # dfdx_h
    if(is.null(dfdx)){
        dfdx <- function(x){
            df(x, h = dfdx_h)
        }
    }

    print(paste0("iteration: 0; x value: ", x, "; f(x) value: ", round(f(x))))

    # starts iteration
    iter = 1
    x_vec <- numeric(max_iter)

    # initial value of the function
    f_x <- f(x)


    # first Newton iteration
    x_vec[iter] <- x - f_x/dfdx(x)
    f_x <- f(x_vec[iter])
    print(paste0("iteration: ", iter, "; x value: ", round(x_vec[iter],4), "; f(x) value: ", round(f_x,4)))
    iter = iter +1

    # checks convergence of the algorithm
    # (receives value after first iteration)
    check_convergence <- f_x
    count_warning <- 0

    while(abs(f_x)>tol && iter <= max_iter){

        # calculates Newton value
        x <- x_vec[iter-1] - f_x/dfdx(x_vec[iter-1])

        if(iter<3){
            # starts with Newton's method
            x_vec[iter] <- x
        }else{

            # This is Steffensen's Method

            # some useful dereference
            xlag_1 <- x_vec[iter-1]
            xlag_2 <- x_vec[iter-2]
            denom <- (x - 2*xlag_1 + xlag_2)

            # if denominator is null, uses Newton method
            if(denom == 0){
                x_vec[iter] <- x

            }else{
                # uses accelerated method
                x_vec[iter] <- xlag_2 - (xlag_1 - xlag_2)^2/denom

            }
        }

        f_x <- f(x_vec[iter])
        print(paste0("iteration: ", iter, "; x value: ", round(x_vec[iter],4), "; f(x) value: ", round(f_x,4)))
        iter = iter +1
        if(is.na(f_x)) break

        # saves value from 5 interactions behind, to check convergence
        if(iter%%5 == 0){
            check_convergence <- f_x

        }

        if(abs(f_x) > abs(check_convergence)){
            warning("Algorithm is not converging.")
            count_warning = count_warning + 1
        }

        if(count_warning == 3){
            stop("Algorithm is not converging.")
        }

    }

    return(list(x=x_vec[1:(iter-1)], fx = f_x))
}


# generates Subbotin Matrix of Covariances
subbo_covar_r <- function(subbo_test, N){

    b <- subbo_test$dt[1, "coef"]
    a <- subbo_test$dt[2, "coef"]

    tmp <- log(b)+digamma(1+1./b)
    dim <- dim(subbo_test$matrix)[[1]]
    I <- matrix(rep(0,9), ncol = 3)

    #* b - b */
    I[1,1] = (tmp*tmp+trigamma(1+1/b)*(1+1/b)-1)/b^3
    #* b - a */
    I[1,2] = -tmp/(a*b);
    I[2,1] = -tmp/(a*b);

    #* a - a */
    I[2,2] = b/(a*a);

    if(dim == 3){

        #* b - m */
        I[1,3] = 0;
        I[3,1] = 0;

        #* a - m */
        I[2,3] = 0;
        I[3,2] = 0;

        #* m - m */
        I[3,3] = (b^(-2/b+1)*gamma(2-1/b))/
            (gamma(1+1/b)*a*a);
    }

    # generates inverse
    I <- solve(I)
    # this is what Bottazzi does with GSL
    # library("Matrix")
    # I <- expand(lu(I))
    # I <- solve(I$U) %*% solve(I$L) %*% I$P

    #* set the var-covar matrix */
    for(i in 2:dim){
        for(j in 1:(i-1)){
            I[i,j] = I[i,j]/sqrt(I[i,i]* I[j,j]);
        }
    }

    for(i in 1:dim){
        for(j in i:dim){
            I[i,j] = I[i,j]/N;
        }
    }
    return(I)
}
