## functions used for supplement

##' Final size for model 1
##' @param n
##' @param q
##' @param Q
##' @param truncated
##' @param total.vector.length
##' @param log
##' @return
##' @author Andrew Azman
fs_probs_mod1 <- function(n,
                          q,
                          Q,
                          truncated=TRUE,
                          total.vector.length,
                          log){

    A = matrix(0,nrow=n+1,ncol=n+1)
    Y = matrix(0,nrow=n+1,ncol=1)
    for(k in 0:n){
        Y[k+1] = choose(n,k)
        for(m in 0:k){
            A[k+1,m+1] <- choose(n-m,k-m)/(Q^(n-k)*q^(m*(n-k)))
        }
    }
    probs = solve(A,Y)
    rownames(probs) = paste((0:n))

    if (truncated){
        probs <- probs[-1]/(1-probs[1])
    }

    if (log) probs <- log(probs)

    if(!missing(total.vector.length)) {
        add.more <- total.vector.length - length(probs)
        if (add.more > 0) probs <- c(probs,rep(0,add.more))
        if (add.more < 0) probs <- probs[-((length(probs) + add.more + 1):length(probs))]
    }
    return(probs)
}


##' .. content for \description{} (no empty lines) ..
##' gamma moment generating function (mean = beta/n^alpha, shape = shape)
##' @param x
##' @param shape
##' @param beta
##' @param alpha
##' @param n
##' @return
GammaMGF <- function(x,shape,beta,alpha,n){
    (shape/(shape + x*(beta/n^(alpha))))^shape
}

#includes random infectiousness
fs_probs_mod2 <- function(n,shape,beta,alpha,total.vector.length,log){

    A = matrix(0,nrow=n,ncol=n)
    Y = matrix(0,nrow=n,ncol=1)
    for(k in 0:(n-1)){
        Y[k+1] = choose(n-1,k)
        for(m in 0:k){
            A[k+1,m+1] <- choose(n-m-1,k-m)/(GammaMGF(n-1-k,shape,beta,alpha,n)^(m + 1))
        }
    }
    probs = solve(A,Y)
    rownames(probs) = paste((1:n))

    if (log) probs <- log(probs)

    if(!missing(total.vector.length)) {
        add.more <- total.vector.length - length(probs)
        if (add.more > 0) probs <- c(probs,rep(0,add.more))
        if (add.more < 0) probs <- probs[-((length(probs) + add.more + 1):length(probs))]
    }
    return(probs)
}


##' final size probabilities for model with household size dependent SITP and no random infectiousness
##' (Used in revised paper)
##' @param n household size
##' @param q
##' @param alpha household size depednence
##' @param total.vector.length how long of a vector do we want to output
##' @param log T/F
##' @param truncated do we want the truncated or full form?
##' @return
fs_probs_mod3 <- function(n,q,alpha,total.vector.length,log,truncated=TRUE){
    if (truncated){
        A = matrix(0,nrow=n,ncol=n)
        Y = matrix(0,nrow=n,ncol=1)
        for(k in 0:(n-1)){
            Y[k+1] = choose(n-1,k)
            for(m in 0:k){
                A[k+1,m+1] <- choose(n-m-1,k-m)/plogis(q/n^alpha)^((m+1)*(n-1-k))
            }
        }

        probs = solve(A,Y)
        rownames(probs) = paste((1:n))

        if (log) {
            probs <- log(probs)
            if (any(is.nan(probs))){
                warning("probs went below 0") #going below zero due to truncation error
                probs[is.nan(probs)] <- -50
            }
        }

    if(!missing(total.vector.length)) {
        add.more <- total.vector.length - length(probs)
        if (add.more > 0) probs <- c(probs,rep(0,add.more))
        if (add.more < 0) probs <- probs[-((length(probs) + add.more + 1):length(probs))]
    }
    } else {
        A = matrix(0,nrow=n+1,ncol=n+1)
        Y = matrix(0,nrow=n+1,ncol=1)
        for(k in 0:n){
            Y[k+1] = choose(n,k)

            for(m in 0:k){
                A[k+1,m+1] <- choose(n-m,k-m)/plogis(q/n^alpha)^((n-k)*(1+m))
            }
        }

        probs = solve(A,Y)
        rownames(probs) = paste((0:n))

        if (log) {
            probs <- log(probs)
            if (any(is.nan(probs))){
                warning("probs went below 0") #going below zero due to truncation error
                probs[is.nan(probs)] <- -50
            }
        }

        if(!missing(total.vector.length)) {
            add.more <- total.vector.length - length(probs)
            if (add.more > 0) probs <- c(probs,rep(0,add.more))
            if (add.more < 0) probs <- probs[-((length(probs) + add.more + 1):length(probs))]
        }

    }

    return(probs)
}

##'
##' @param params
##' @param table
##' @return
##' @author Andrew Azman
lik1 <- function(params,table){

    Q <- params["Q"]
    q <- params["q"]

    #initialize log-liklihood
    ll <- 0

    # deal with row and columns that are NA (ie.. missing household size)
    rem.cols <- is.na(colnames(table))
    rem.rows <- is.na(rownames(table))
    if (any(rem.cols)) table <- table[,-which(rem.cols)]
    if (any(rem.rows)) table <- table[-which(rem.rows),]

    #convert row and column names to integers
    sec.cases <- as.numeric(rownames(table))
    hhl.sizes <- as.numeric(colnames(table))

    # go through each column of the table
    for (i in 1:ncol(table)){
        ll <- ll + sum(table[,i] * fs_probs_mod1(hhl.sizes[i],q,Q,total.vector.length = nrow(table),log=TRUE))
    }
    return(-ll)
}

lik1.1 <- function(params,table1,table2){

    Q <- params["Q"]
    q1 <- params["q1"]
    q2 <- params["q2"]
    #initialize log-liklihood
    ll.1 <- ll.2 <- 0

    # deal with row and columns that are NA (ie.. missing household size)
    rem.cols1 <- is.na(colnames(table1))
    rem.rows1 <- is.na(rownames(table1))
    if (any(rem.cols1)) table1 <- table1[,-which(rem.cols1)]
    if (any(rem.rows1)) table1 <- table1[-which(rem.rows1),]

    rem.cols2 <- is.na(colnames(table2))
    rem.rows2 <- is.na(rownames(table2))
    if (any(rem.cols2)) table2 <- table2[,-which(rem.cols2)]
    if (any(rem.rows2)) table2 <- table2[-which(rem.rows2),]


    #convert row and column names to integers
    sec.cases1 <- as.numeric(rownames(table1))
    hhl.sizes1 <- as.numeric(colnames(table1))
    sec.cases2 <- as.numeric(rownames(table2))
    hhl.sizes2 <- as.numeric(colnames(table2))

    # go through each column of the table
    for (i in 1:ncol(table1)){
        ll.1 <- ll.1 + sum(table1[,i] * fs_probs_mod1(hhl.sizes1[i],q1,Q,total.vector.length = nrow(table1),log=TRUE))
    }

    for (i in 1:ncol(table2)){
        ll.2 <- ll.2 + sum(table2[,i] * fs_probs_mod1(hhl.sizes2[i],q2,Q,total.vector.length = nrow(table2),log=TRUE))
    }

    return(-(ll.1 + ll.2))
}

##'
##' @param params
##' @param table1
##' @param table2
##' @return
##' @author Andrew Azman
lik1.2 <- function(params,table1,table2){

    Q1 <- params["Q1"]
    Q2 <- params["Q2"]
    q1 <- params["q1"]
    q2 <- params["q2"]
    #initialize log-liklihood
    ll.1 <- ll.2 <- 0

    # deal with row and columns that are NA (ie.. missing household size)
    rem.cols1 <- is.na(colnames(table1))
    rem.rows1 <- is.na(rownames(table1))
    if (any(rem.cols1)) table1 <- table1[,-which(rem.cols1)]
    if (any(rem.rows1)) table1 <- table1[-which(rem.rows1),]

    rem.cols2 <- is.na(colnames(table2))
    rem.rows2 <- is.na(rownames(table2))
    if (any(rem.cols2)) table2 <- table2[,-which(rem.cols2)]
    if (any(rem.rows2)) table2 <- table2[-which(rem.rows2),]


    #convert row and column names to integers
    sec.cases1 <- as.numeric(rownames(table1))
    hhl.sizes1 <- as.numeric(colnames(table1))
    sec.cases2 <- as.numeric(rownames(table2))
    hhl.sizes2 <- as.numeric(colnames(table2))

    # go through each column of the table
    for (i in 1:ncol(table1)){
        ll.1 <- ll.1 + sum(table1[,i] * fs_probs_mod1(hhl.sizes1[i],q1,Q1,total.vector.length = nrow(table1),log=TRUE))
    }

    for (i in 1:ncol(table2)){
        ll.2 <- ll.2 + sum(table2[,i] * fs_probs_mod1(hhl.sizes2[i],q2,Q2,total.vector.length = nrow(table2),log=TRUE))
    }

    return(-(ll.1 + ll.2))
}

##' no community transmission prob and only 1 q
##' @param params
##' @param table1
##' @return
##' @author Andrew Azman
lik1.3 <- function(params,table1){

    Q1 <- 1
    q1 <- params["q1"]
    #initialize log-liklihood
    ll.1 <- 0

    # deal with row and columns that are NA (ie.. missing household size)
    rem.cols1 <- is.na(colnames(table1))
    rem.rows1 <- is.na(rownames(table1))
    if (any(rem.cols1)) table1 <- table1[,-which(rem.cols1)]
    if (any(rem.rows1)) table1 <- table1[-which(rem.rows1),]

    #convert row and column names to integers
    sec.cases1 <- as.numeric(rownames(table1))
    hhl.sizes1 <- as.numeric(colnames(table1))

    # go through each column of the table
    for (i in 1:ncol(table1)){
        ll.1 <- ll.1 + sum(table1[,i] * fs_probs_mod1(hhl.sizes1[i],q1,Q1,total.vector.length = nrow(table1),log=TRUE))
    }

    return(-ll.1)
}

##'
##' @param params
##' @param table1
##' @param table2
##' @return
##' @author Andrew Azman
lik1.4 <- function(params,table1,table2){

    Q1 <- 1
    Q2 <- 1
    q1 <- params["q1"]
    q2 <- params["q2"]
    #initialize log-liklihood
    ll.1 <- ll.2 <- 0

    # deal with row and columns that are NA (ie.. missing household size)
    rem.cols1 <- is.na(colnames(table1))
    rem.rows1 <- is.na(rownames(table1))
    if (any(rem.cols1)) table1 <- table1[,-which(rem.cols1)]
    if (any(rem.rows1)) table1 <- table1[-which(rem.rows1),]

    rem.cols2 <- is.na(colnames(table2))
    rem.rows2 <- is.na(rownames(table2))
    if (any(rem.cols2)) table2 <- table2[,-which(rem.cols2)]
    if (any(rem.rows2)) table2 <- table2[-which(rem.rows2),]


    #convert row and column names to integers
    sec.cases1 <- as.numeric(rownames(table1))
    hhl.sizes1 <- as.numeric(colnames(table1))
    sec.cases2 <- as.numeric(rownames(table2))
    hhl.sizes2 <- as.numeric(colnames(table2))

    # go through each column of the table
    for (i in 1:ncol(table1)){
        ll.1 <- ll.1 + sum(table1[,i] * fs_probs_mod1(hhl.sizes1[i],q1,Q1,total.vector.length = nrow(table1),log=TRUE))
    }

    for (i in 1:ncol(table2)){
        ll.2 <- ll.2 + sum(table2[,i] * fs_probs_mod1(hhl.sizes2[i],q2,Q2,total.vector.length = nrow(table2),log=TRUE))
    }

    return(-(ll.1 + ll.2))
}

##'
##' @param params
##' @param table
##' @return
##' @author Andrew Azman
lik2 <- function(params,table){

    shape <- exp(params["shape"])
    beta <- exp(params["beta"])
    alpha <- exp(params["alpha"])

    #initialize log-liklihood
    ll <- 0

    # deal with row and columns that are NA (ie.. missing household size)
    rem.cols <- is.na(colnames(table))
    rem.rows <- is.na(rownames(table))
    if (any(rem.cols)) table <- table[,-which(rem.cols)]
    if (any(rem.rows)) table <- table[-which(rem.rows),]

    #convert row and column names to integers
    sec.cases <- as.numeric(rownames(table))
    hhl.sizes <- as.numeric(colnames(table)) # since we are dealing with a truncted situation

    # go through each column of the table
    for (i in 1:ncol(table)){
        ll <- ll + sum(table[,i] * fs_probs_mod2(hhl.sizes[i],
                                                 shape = shape,
                                                 beta = beta,
                                                 alpha = alpha,
                                                 total.vector.length = nrow(table),
                                                 log=TRUE))
    }
    return(-ll)
}

##'
##' @param params
##' @param table1
##' @param table2
##' @return
##' @author Andrew Azman
lik2.1 <- function(params,table1,table2){

    shape <- exp(params["shape"])
    beta1 <- exp(params["beta1"])
    beta2 <- exp(params["beta2"])
    alpha <- exp(params["alpha"])

        #initialize log-liklihood
    ll.1 <- ll.2 <- 0

    # deal with row and columns that are NA (ie.. missing household size)
    rem.cols1 <- is.na(colnames(table1))
    rem.rows1 <- is.na(rownames(table1))
    if (any(rem.cols1)) table1 <- table1[,-which(rem.cols1)]
    if (any(rem.rows1)) table1 <- table1[-which(rem.rows1),]

    rem.cols2 <- is.na(colnames(table2))
    rem.rows2 <- is.na(rownames(table2))
    if (any(rem.cols2)) table2 <- table2[,-which(rem.cols2)]
    if (any(rem.rows2)) table2 <- table2[-which(rem.rows2),]

    #convert row and column names to integers
    sec.cases1 <- as.numeric(rownames(table1))
    hhl.sizes1 <- as.numeric(colnames(table1))
    sec.cases2 <- as.numeric(rownames(table2))
    hhl.sizes2 <- as.numeric(colnames(table2))

                                        # go through each column of the table
    for (i in 1:ncol(table1)){
        ll.1 <- ll.1 + sum(table1[,i] * fs_probs_mod2(hhl.sizes1[i],
                                                     shape = shape,
                                                     beta = beta1,
                                                     alpha = alpha,
                                                     total.vector.length = nrow(table1),
                                                     log=TRUE))
    }

    for (i in 1:ncol(table2)){
        ll.2 <- ll.2 + sum(table2[,i] * fs_probs_mod2(hhl.sizes2[i],
                                                     shape = shape,
                                                     beta = beta2,
                                                     alpha = alpha,
                                                     total.vector.length = nrow(table2),
                                                     log=TRUE))
    }
    return(-(ll.1 + ll.2))
}

##'
##' @title
##' @param params
##' @param table1
##' @param table2
##' @return
##' @author Andrew Azman
lik2.2 <- function(params,table1,table2){

    shape <- exp(params["shape"])
    beta1 <- exp(params["beta1"])
    beta2 <- exp(params["beta2"])
    alpha <- 0

        #initialize log-liklihood
    ll.1 <- ll.2 <- 0

    # deal with row and columns that are NA (ie.. missing household size)
    rem.cols1 <- is.na(colnames(table1))
    rem.rows1 <- is.na(rownames(table1))
    if (any(rem.cols1)) table1 <- table1[,-which(rem.cols1)]
    if (any(rem.rows1)) table1 <- table1[-which(rem.rows1),]

    rem.cols2 <- is.na(colnames(table2))
    rem.rows2 <- is.na(rownames(table2))
    if (any(rem.cols2)) table2 <- table2[,-which(rem.cols2)]
    if (any(rem.rows2)) table2 <- table2[-which(rem.rows2),]

    #convert row and column names to integers
    sec.cases1 <- as.numeric(rownames(table1))
    hhl.sizes1 <- as.numeric(colnames(table1))
    sec.cases2 <- as.numeric(rownames(table2))
    hhl.sizes2 <- as.numeric(colnames(table2))

                                        # go through each column of the table
    for (i in 1:ncol(table1)){
        ll.1 <- ll.1 + sum(table1[,i] * fs_probs_mod2(hhl.sizes1[i],
                                                     shape = shape,
                                                     beta = beta1,
                                                     alpha = alpha,
                                                     total.vector.length = nrow(table1),
                                                     log=TRUE))
    }

    for (i in 1:ncol(table2)){
        ll.2 <- ll.2 + sum(table2[,i] * fs_probs_mod2(hhl.sizes2[i],
                                                     shape = shape,
                                                     beta = beta2,
                                                     alpha = alpha,
                                                     total.vector.length = nrow(table2),
                                                     log=TRUE))
    }

    return(-(ll.1 + ll.2))
}


##'
##' @param params
##' @param table1
##' @param table2
##' @return
##' @author Andrew Azman
lik2.3 <- function(params,table1,table2){

    shape <- 0.94 #from Fraser
    beta1 <- exp(params["beta1"])
    beta2 <- exp(params["beta2"])
    alpha <- exp(params["alpha"])

        #initialize log-liklihood
    ll.1 <- ll.2 <- 0

    # deal with row and columns that are NA (ie.. missing household size)
    rem.cols1 <- is.na(colnames(table1))
    rem.rows1 <- is.na(rownames(table1))
    if (any(rem.cols1)) table1 <- table1[,-which(rem.cols1)]
    if (any(rem.rows1)) table1 <- table1[-which(rem.rows1),]

    rem.cols2 <- is.na(colnames(table2))
    rem.rows2 <- is.na(rownames(table2))
    if (any(rem.cols2)) table2 <- table2[,-which(rem.cols2)]
    if (any(rem.rows2)) table2 <- table2[-which(rem.rows2),]

    #convert row and column names to integers
    sec.cases1 <- as.numeric(rownames(table1))
    hhl.sizes1 <- as.numeric(colnames(table1))
    sec.cases2 <- as.numeric(rownames(table2))
    hhl.sizes2 <- as.numeric(colnames(table2))

                                        # go through each column of the table
    for (i in 1:ncol(table1)){
        ll.1 <- ll.1 + sum(table1[,i] * fs_probs_mod2(hhl.sizes1[i],
                                                      shape = shape,
                                                      beta = beta1,
                                                      alpha = alpha,
                                                      total.vector.length = nrow(table1),
                                                      log=TRUE))
    }

    for (i in 1:ncol(table2)){
        ll.2 <- ll.2 + sum(table2[,i] * fs_probs_mod2(hhl.sizes2[i],
                                                     shape = shape,
                                                     beta = beta2,
                                                     alpha = alpha,
                                                     total.vector.length = nrow(table2),
                                                     log=TRUE))
    }

    return(-(ll.1 + ll.2))
}


##'
##' @param params
##' @param table1
##' @param table2
##' @return
##' @author Andrew Azman
lik2.4 <- function(params,table1,table2){

    shape <- 0.94 #from Fraser
    beta1 <- exp(params["beta1"])
    beta2 <- exp(params["beta2"])
    alpha <- 0

        #initialize log-liklihood
    ll.1 <- ll.2 <- 0

    # deal with row and columns that are NA (ie.. missing household size)
    rem.cols1 <- is.na(colnames(table1))
    rem.rows1 <- is.na(rownames(table1))
    if (any(rem.cols1)) table1 <- table1[,-which(rem.cols1)]
    if (any(rem.rows1)) table1 <- table1[-which(rem.rows1),]

    rem.cols2 <- is.na(colnames(table2))
    rem.rows2 <- is.na(rownames(table2))
    if (any(rem.cols2)) table2 <- table2[,-which(rem.cols2)]
    if (any(rem.rows2)) table2 <- table2[-which(rem.rows2),]

    #convert row and column names to integers
    sec.cases1 <- as.numeric(rownames(table1))
    hhl.sizes1 <- as.numeric(colnames(table1))
    sec.cases2 <- as.numeric(rownames(table2))
    hhl.sizes2 <- as.numeric(colnames(table2))

                                        # go through each column of the table
    for (i in 1:ncol(table1)){
        ll.1 <- ll.1 + sum(table1[,i] * fs_probs_mod2(hhl.sizes1[i],
                                                      shape = shape,
                                                      beta = beta1,
                                                      alpha = alpha,
                                                      total.vector.length = nrow(table1),
                                                      log=TRUE))
    }

    for (i in 1:ncol(table2)){
        ll.2 <- ll.2 + sum(table2[,i] * fs_probs_mod2(hhl.sizes2[i],
                                                     shape = shape,
                                                     beta = beta2,
                                                     alpha = alpha,
                                                     total.vector.length = nrow(table2),
                                                     log=TRUE))
    }

    return(-(ll.1 + ll.2))
}


##'
##' @param params
##' @param table
##' @return
##' @author Andrew Azman
lik3 <- function(params,table){

    q <- params["q"]
    alpha <- 0

    #initialize log-liklihood
    ll <- 0

    # deal with row and columns that are NA (ie.. missing household size)
    rem.cols <- is.na(colnames(table))
    rem.rows <- is.na(rownames(table))
    if (any(rem.cols)) table <- table[,-which(rem.cols)]
    if (any(rem.rows)) table <- table[-which(rem.rows),]

    #convert row and column names to integers
    sec.cases <- as.numeric(rownames(table))
    hhl.sizes <- as.numeric(colnames(table)) # since we are dealing with a truncted situation

    # go through each column of the table
    for (i in 1:ncol(table)){
        ll <- ll + sum(table[,i] * fs_probs_mod3(hhl.sizes[i],
                                                 q=q,
                                                 alpha=alpha,
                                                 total.vector.length = nrow(table),
                                                 log=TRUE))
    }
    return(-ll)
}

lik3.1 <- function(params,table1,table2){

    q1 <- params["q1"]
    q2 <- params["q2"]
    alpha <- 0

    ## initialize log-liklihood
    ll.1 <- ll.2 <- 0

    # deal with row and columns that are NA (ie.. missing household size)
    rem.cols1 <- is.na(colnames(table1))
    rem.rows1 <- is.na(rownames(table1))
    if (any(rem.cols1)) table1 <- table1[,-which(rem.cols1)]
    if (any(rem.rows1)) table1 <- table1[-which(rem.rows1),]

    rem.cols2 <- is.na(colnames(table2))
    rem.rows2 <- is.na(rownames(table2))
    if (any(rem.cols2)) table2 <- table2[,-which(rem.cols2)]
    if (any(rem.rows2)) table2 <- table2[-which(rem.rows2),]

    #convert row and column names to integers
    sec.cases1 <- as.numeric(rownames(table1))
    hhl.sizes1 <- as.numeric(colnames(table1))
    sec.cases2 <- as.numeric(rownames(table2))
    hhl.sizes2 <- as.numeric(colnames(table2))

                                        # go through each column of the table
    for (i in 1:ncol(table1)){
        ll.1 <- ll.1 + sum(table1[,i] * fs_probs_mod3(hhl.sizes1[i],
                                                      q=q1,
                                                      alpha=alpha,
                                                      total.vector.length = nrow(table1),
                                                      log=TRUE))
    }

    for (i in 1:ncol(table2)){
        ll.2 <- ll.2 + sum(table2[,i] * fs_probs_mod3(hhl.sizes2[i],
                                                      q=q2,
                                                      alpha=alpha,
                                                      total.vector.length = nrow(table2),
                                                      log=TRUE))
    }

    return(-(ll.1 + ll.2))
}

##' Model used in revised paper
##' @param params
##' @param table1
##' @param table2
##' @return
##' @author Andrew Azman
lik3.2 <- function(params,table1,table2){

    q1 <- params["q1"]
    q2 <- params["q2"]
    alpha <- params["alpha"]

    ## initialize log-liklihood
    ll.1 <- ll.2 <- 0

    # deal with row and columns that are NA (ie.. missing household size)
    rem.cols1 <- is.na(colnames(table1))
    rem.rows1 <- is.na(rownames(table1))
    if (any(rem.cols1)) table1 <- table1[,-which(rem.cols1)]
    if (any(rem.rows1)) table1 <- table1[-which(rem.rows1),]

    rem.cols2 <- is.na(colnames(table2))
    rem.rows2 <- is.na(rownames(table2))
    if (any(rem.cols2)) table2 <- table2[,-which(rem.cols2)]
    if (any(rem.rows2)) table2 <- table2[-which(rem.rows2),]

                                        #convert row and column names to integers
    if (is.null(colnames(table1)) | is.null(colnames(table2))){
        hhl.sizes1 <- hhl.sizes2 <- 2:6
    } else {
        hhl.sizes1 <- as.numeric(colnames(table1))
        hhl.sizes2 <- as.numeric(colnames(table2))
    }

    if (is.null(rownames(table1)) | is.null(rownames(table2))){
        sec.cases1 <- sec.cases2 <- 0:5
    } else {
        sec.cases1 <- as.numeric(rownames(table1))
        sec.cases2 <- as.numeric(rownames(table2))
    }
                                        # go through each column of the table

    for (i in 1:ncol(table1)){
        ll.1 <- ll.1 + sum(table1[,i] * fs_probs_mod3(hhl.sizes1[i],
                                                      q=q1,
                                                      alpha=alpha,
                                                      total.vector.length = nrow(table1),
                                                      log=TRUE))
    }

    ## this is a little bit of a sloppy hack
    ## sometime we feed a table with all zeros to fit a one parameter model
    if (sum(table2) == 0){
        ll.2 <- 0
    } else {
        for (i in 1:ncol(table2)){
            ll.2 <- ll.2 + sum(table2[,i] * fs_probs_mod3(hhl.sizes2[i],
                                                          q=q2,
                                                          alpha=alpha,
                                                          total.vector.length = nrow(table2),
                                                          log=TRUE))
        }
    }

    #if (is.nan(-(ll.1 + ll.2))) recover()

    return(-(ll.1 + ll.2))
}

longs <- function(B,
                  Q,
                  trunc=FALSE,   #trunc is for the truncated
                  V1=FALSE){    #V1=TRUE indicates the NC version of the distribution NOTE THIS IS INCORRECT DO NOT USE!!!
                                        #zero-class case
  x <- matrix(rep(0,7*6),nrow=7,ncol=6)
  rownames(x) <- paste("m",0:6,"k",sep="")
  colnames(x) <- paste("k=",1:6,sep="")
  x[1,1] <- B       #start off with households with 1 person
  x[2,1] <- 1-x[1,1]
  for (k in 2:6){
    for (j in 0:k){
      if (j == 0) {  #for households with zero cases we know it is B^k
        x[1,k] <- B^k
      }
      else if (j == k){
        x[j+1,k] <- 1-sum(x[(1:j),k])
#        if (x[j+1,k] == 0) print(paste("Warning imputed zero probability for household size=",k,"-",x[j+1,k]))
      }
      else {
        if (V1){   #this version used in NC 2007
          x[j+1,k] <- choose(k,j)*x[j+1,j]*B^(k-j)*Q^(j*(k-1))
#          if (x[j+1,k] <= 0) print(x[j+1,k])

        }
        else {   ##this version used in Longini 1982
          x[(j+1),k] <- choose(k,j)*x[(j+1),j]*B^(k-j)*Q^(j*(k-j))
#          if (x[j+1,k] == 0) print(paste("Warning imputed zero probability for household size=",k))
                  }
      }
    }
  }
  if (trunc) {
    zero.class <- x[1,]
    new.x <- x[-1,]
    for (i in 1:length(zero.class)){
      new.x[,i] <- new.x[,i]/(1-zero.class[i])
    }
    return(new.x)
  }
  else {
    return(x)
  }
}


# profile likelihood for best model
prof.lik.3.2 <- function(params,table1,table2,prof,prof.param){
    if (prof == "q1") {
        q1 <- prof.param
        q2 <- params["q2"]
        alpha <- params["alpha"]
    } else if (prof == "q2") {
        q2 <- prof.param
        q1 <- params["q1"]
        alpha <- params["alpha"]
    } else if (prof == "alpha"){
        q2 <- params["q2"]
        q1 <- params["q1"]
        alpha <- prof.param
    } else {
        q1 <- params["q1"]
        q2 <- params["q2"]
        alpha <- params["alpha"]
    }

    ## initialize log-liklihood
    ll.1 <- ll.2 <- 0

    # deal with row and columns that are NA (ie.. missing household size)
    rem.cols1 <- is.na(colnames(table1))
    rem.rows1 <- is.na(rownames(table1))
    if (any(rem.cols1)) table1 <- table1[,-which(rem.cols1)]
    if (any(rem.rows1)) table1 <- table1[-which(rem.rows1),]

    rem.cols2 <- is.na(colnames(table2))
    rem.rows2 <- is.na(rownames(table2))
    if (any(rem.cols2)) table2 <- table2[,-which(rem.cols2)]
    if (any(rem.rows2)) table2 <- table2[-which(rem.rows2),]

    #convert row and column names to integers
    sec.cases1 <- as.numeric(rownames(table1))
    hhl.sizes1 <- as.numeric(colnames(table1))
    sec.cases2 <- as.numeric(rownames(table2))
    hhl.sizes2 <- as.numeric(colnames(table2))

                                        # go through each column of the table
    for (i in 1:ncol(table1)){
        ll.1 <- ll.1 + sum(table1[,i] * fs_probs_mod3(hhl.sizes1[i],
                                                      q=q1,
                                                      alpha=alpha,
                                                      total.vector.length = nrow(table1),
                                                      log=T))
    }

    for (i in 1:ncol(table2)){
        ll.2 <- ll.2 + sum(table2[,i] * fs_probs_mod3(hhl.sizes2[i],
                                                      q=q2,
                                                      alpha=alpha,
                                                      total.vector.length = nrow(table2),
                                                      log=TRUE))
    }
    if (is.nan(-(ll.1 + ll.2))) recover()
    return(-(ll.1 + ll.2))
}

GetProfLik <- function(params,prof,prof.start,table1,table2,prof.range=.1){
    proflik <- array(dim=c(100,3))
    #prof.seq <- seq(prof.start*(1-prof.range),prof.start*(1+prof.range),length=50)
    if (prof == "alpha"){
        prof.seq <- seq(-1.6,.1,length=100) #for alpha
    } else {
        prof.seq <- seq(0.08,4,length=100) #for q2
    }

    print(prof.seq)
    cur.param.start <- params
    for (i in seq_along(prof.seq)){
       # print(i)
       # print(cur.param.start)
        tmp <- optim(cur.param.start,
                     fn=prof.lik.3.2,
                     table1=table1,
                     table2=table2,
                     prof=prof,
                     prof.param=prof.seq[i],
                     )
        proflik[i,] <- c(tmp$par,tmp$value)
        #cur.param.start <- proflik[i,1:2]
        #names(cur.param.start) <- names(params)
    }

    print(proflik)
    min.ind <- which.min(proflik[,3])
    NLL <- proflik[min.ind,3]
    prof.lower <- proflik[1:min.ind,3]
    prof.Bvec.l <- prof.seq[1:min.ind]
    l <- approx(prof.lower,prof.Bvec.l,xout=NLL+qchisq(0.95,1)/2)
    l.1 <- approx(prof.lower,proflik[1:min.ind,1],xout=NLL+qchisq(0.95,1)/2)
    l.2 <- approx(prof.lower,proflik[1:min.ind,2],xout=NLL+qchisq(0.95,1)/2)

    prof.upper <- proflik[min.ind:nrow(proflik),3]
    prof.Bvec.u <- prof.seq[min.ind:nrow(proflik)]
    u <- approx(prof.upper,prof.Bvec.u,xout=NLL+qchisq(0.95,1)/2,yright=1)
    u.1 <- approx(prof.upper,proflik[min.ind:nrow(proflik),1],xout=NLL+qchisq(0.95,1)/2,yright=1)
    u.2 <- approx(prof.upper,proflik[min.ind:nrow(proflik),2],xout=NLL+qchisq(0.95,1)/2,yright=1)

    return(c("lower"=l$y,
             "upper"=u$y,
             "u1" = u.1$y,
             "u2" = u.2$y,
             "l1" = l.1$y,
             "l2" = l.2$y))
       }


GetEstimates <- function(){
    out <- vector("list",4)
    #first for AB NoTemp
    q1.1 <- GetProfLik(params=c("alpha"=-0.7070,"q2"=.5233),prof="q1",prof.start=0.4672,table1=hhl.table.A.NoTemp,table2=hhl.table.B.NoTemp)
    q2.1 <- GetProfLik(params=c("alpha"=-0.7070,"q1"=.4672),prof="q2",prof.start=0.5223,table1=hhl.table.A.NoTemp,table2=hhl.table.B.NoTemp)
    alpha.1 <- GetProfLik(params=c("q1"=.2,"q2"=.2),prof="alpha",prof.start=-0.7,table1=hhl.table.A.NoTemp[,-c(6:8)],table2=hhl.table.B.NoTemp[,-c(6:7)])
    out[[1]] <- rbind(q1.1,q2.1,alpha.1)

    ##then for AB Temp
    q1.2 <- GetProfLik(params=c("alpha"=-0.7070,"q2"=.5233),prof="q1",prof.start=0.4672,table1=hhl.table.Cont,table2=hhl.table.B)
    q2.2 <- GetProfLik(params=c("alpha"=-0.7070,"q1"=.4672),prof="q2",prof.start=0.5233,table1=hhl.table.A,table2=hhl.table.B)
    alpha.2 <- GetProfLik(params=c("q1"=.2,"q2"=.2),prof="alpha",prof.start=-0.7070,table1=hhl.table.A,table2=hhl.table.B)
    out[[2]] <- rbind(q1.2,q2.2,alpha.2)

    #then for ContInter NoTemp
    q1.3 <- GetProfLik(params=c("alpha"=-0.73,"q2"=.55),prof="q1",prof.start=0.45,table1=hhl.table.Cont.NoTemp,table2=hhl.table.Inter.NoTemp)
    q2.3 <- GetProfLik(params=c("alpha"=-0.7228,"q1"=.45),prof="q2",prof.start=0.5150,table1=hhl.table.Cont.NoTemp,table2=hhl.table.Inter.NoTemp)
    alpha.3 <- GetProfLik(params=c("q1"=.2,"q2"=.2),prof="alpha",prof.start=-0.7,table1=hhl.table.Cont.NoTemp,table2=hhl.table.Inter.NoTemp)
    out[[3]] <- rbind(q1.3,q2.3,alpha.3)

                                        #then for ContInter Temp
    q1.4 <- GetProfLik(params=c("alpha"=-0.73,"q2"=.55),prof="q1",prof.start=0.45,table1=hhl.table.Cont,table2=hhl.table.Inter)
    q2.4 <- GetProfLik(params=c("alpha"=-0.7228,"q1"=.45),prof="q2",prof.start=0.5150,table1=hhl.table.Cont,table2=hhl.table.Inter)
    alpha.4 <- GetProfLik(params=c("q1"=.2,"q2"=.2),prof="alpha",prof.start=-0.7,table1=hhl.table.Cont,table2=hhl.table.Inter)
    out[[4]] <- rbind(q1.4,q2.4,alpha.4)
    return(out)
}

# get expected final contigency tables

make.expected.tables <- function(){
    pars <- array(dim=c(2,3))
    pars[1,] <- optim(c("q1"=.85,"q2"=.85,"alpha"=0),fn=lik3.2,table1=hhl.table.A.NoTemp,table2=hhl.table.B.NoTemp)$par
    pars[2,] <- optim(c("q1"=.85,"q2"=.85,"alpha"=0),fn=lik3.2,table1=hhl.table.Cont.NoTemp,table2=hhl.table.Inter.NoTemp)$par

    tabs.out <- list()
    pvals <- c()
    ts <- c()

    tabs.out$A <- round(make.obs.table(hhl.table.A.NoTemp,q=pars[1,1],alpha=pars[1,3]),1)
    chisq.A <- c((hhl.table.A.NoTemp[1:6,][,1:5] - tabs.out$A)^2/tabs.out$A)
    chisq.A <- chisq.A[-which(is.nan(chisq.A))]
    ts[1] <- sum(chisq.A)
    pvals[1] <- 1 - pchisq(ts[1],df=length(chisq.A)  - 3 - 1)

    tabs.out$B <- round(make.obs.table(hhl.table.B.NoTemp,q=pars[1,2],alpha=pars[1,3]),1)
    chisq.B <- c((hhl.table.B.NoTemp[1:6,][,1:5] - tabs.out$B)^2/tabs.out$B)
    chisq.B <- chisq.B[-which(is.nan(chisq.B))]
    ts[2] <- sum(chisq.B)
    pvals[2] <- 1 - pchisq(ts[2],df=length(chisq.B)  - 3 - 1)

    tabs.out$Cont <- round(make.obs.table(hhl.table.Cont.NoTemp,q=pars[2,1],alpha=pars[2,3]),1)
    chisq.Cont <- c((hhl.table.Cont.NoTemp[1:6,][,1:5] - tabs.out$Cont)^2/tabs.out$Cont)
    chisq.Cont <- chisq.Cont[-which(is.nan(chisq.Cont))]
    ts[3] <- sum(chisq.Cont)
    pvals[3] <- 1 - pchisq(ts[3],df=length(chisq.Cont)  - 3 - 1)

    tabs.out$Inter <- round(make.obs.table(hhl.table.Inter.NoTemp,q=pars[2,1],alpha=pars[2,3]),1)
    chisq.Inter <- c((hhl.table.Inter.NoTemp[1:6,][,1:5] - tabs.out$Inter)^2/tabs.out$Inter)
    chisq.Inter <- chisq.Inter[-which(is.nan(chisq.Inter))]
    ts[4] <- sum(chisq.Inter)
    pvals[4] <- 1 - pchisq(ts[4],df=length(chisq.Inter)  - 3 - 1)

    return(list(tabs.out=tabs.out,pvals=pvals,ts=ts))
}

make.obs.table <- function(table,q=0.5150,alpha=-0.7228){
    ps <- NULL
    for (i in 2:6){
        p.temp <- fs_probs_mod3(n=i,
                                q=q,
                                alpha=alpha,
                                total.vector.length = 6,
                                log=FALSE,
                                truncated=TRUE)
        ps <- cbind(ps,p.temp)
    }
    ret <- ps
    for (i in 1:5) ret[,i] <- colSums(table)[i]*ps[,i]
    return(ret)
}
