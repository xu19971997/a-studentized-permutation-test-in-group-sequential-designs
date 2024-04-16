##################################################
###### Two stages balanced group sequential design, Type I error rate
##################################################
###### Authors: Long-Hao Xu, Tobias Mütze
###### Email: long-hao.xu@med.uni-goettingen.de    OR    longhao.xu@outlook.com
###################################

rm(list = ls())

###################################
## Define test statistics
getZ <- function(data_t,data_c){
    n_t <- length(data_t)
    n_c <- length(data_c)

    S <- ( mean(data_t) - mean(data_c) ) / sqrt( var(data_t)/n_t + var(data_c)/n_c )
    return(S)
}

## Define function that calculates test statistic for each row of two matrices
getZMat <- function(matT, matC) {
  nT <- dim(matT)[2]
  nC <- dim(matC)[2]
  rowmeanT <- rowMeans(matT)
  rowmeanC <- rowMeans(matC)
  varT <- (rowSums(matT^2) - nT*rowmeanT^2) / (nT-1)
  varC <- (rowSums(matC^2) - nC*rowmeanC^2) / (nC-1)
  
  teststat <- (rowmeanT - rowmeanC) / sqrt(varT/nT + varC/nC)
  return(teststat)
}

## Define the degree of freedom for Welch t-approximation
getdf <- function(data_t, data_c){
  n_t <- length(data_t)
  n_c <- length(data_c)
  n <- n_t + n_c

  sigma2_t <- var(data_t)
  sigma2_c <- var(data_c)

  v <- (sigma2_t/n_t + sigma2_c/n_c)^2/(sigma2_t^2/(n_t^2*(n_t-1)) + sigma2_c^2/(n_c^2*(n_c-1)))
  return(v)
}

## Solve the critical values based on the permutation distribution
g_t_1 <- function(t,Z_1){
    (1/length(Z_1))*( sum( Z_1>t ) )
}
g_t_2 <- function(t,Z_1,Z_2,u_1){
    (1/length(Z_1))*( sum( (Z_1<u_1) & (Z_2>t) ) )
}

## Define alpha-spending function
alpha_spending_P <- function(t) min(alpha*log(1+(exp(1)-1)*t),alpha)
alpha_spending_OF <- function(t) min(2*sided*(1-pnorm( qnorm( (1-alpha/(2*sided)) ) /sqrt(t) ) ), alpha)

## Prevent program interruptions, and record the number of interruptions
i_f_p_asP <- 0
i_f_p_asOF <- 0


###################################

N_a_p_asP <- 0      # Proposed method asP
N_a_norm_asP <- 0   # Traditional method (using normal distribution) asP
N_a_t_app_asP <- 0  # t-approximation asP
N_a_p_asOF <- 0     # Proposed method asOF
N_a_norm_asOF <- 0  # Traditional method (using normal distribution) asOF
N_a_t_app_asOF <- 0 # t-approximation asOF

n <- 5              # The sample size per treatment at any given stage
nperm <- 10000      # The number of permutations
nsim <- 10000       # The number of simulations
alpha <- 0.025      # Type I error rate
sided <- 1          # one-side hypothesis testing OR two-side hypothesis testing

###################################

for (ii in 1:nsim){
    ## Generate data
    data_t <- rexp(2*n,1) #rlnorm(2*n,meanlog=0,sdlog=1) #rchisq(2*n,df=3) #rnorm(2*n,mean=0,sd=1) 
    data_c <- rexp(2*n,1) #rlnorm(2*n,meanlog=0,sdlog=1) #rchisq(2*n,df=3) #rt(2*n,df=5)      
    
    ##
    vec <- c(getZ(data_t[1:n],data_c[1:n]),getZ(data_t,data_c))

    ######
    ###### Proposed method
    ######
    # create matrix that will be used to permutate data
    simMat <- t(sapply(1:nperm, FUN = function(x){sample(1:(2*n), replace = FALSE)}))
    
    # pooled data from stages 1 and 2
    dataPoolS1 <- c(data_t[1:n],data_c[1:n])
    dataPoolS2 <- c(data_t[(n+1):(2*n)], data_c[(n+1):(2*n)])
    
    # Permute the pooled data
    dataS1Perm <- matrix(dataPoolS1[simMat], nrow = nperm, byrow = FALSE)
    dataS2Perm <- matrix(dataPoolS2[simMat], nrow = nperm, byrow = FALSE)
    
    # Calculate the test statistics
    T1perm <- getZMat(matT = dataS1Perm[, 1:n], matC = dataS1Perm[, (n+1):(2*n)])
    T2perm <- getZMat(matT = cbind(dataS1Perm[, 1:n], dataS2Perm[, 1:n]), 
                      matC = cbind(dataS1Perm[, (n+1):(2*n)], dataS2Perm[, (n+1):(2*n)]))

    T_new_p <- cbind(T1perm, T2perm) 

    #### Solve the critical values based on the permutation distribution
    ## Proposed method, asP, directly using 1/2 as the information level
    r <- tryCatch({
        u_1_asP <- uniroot(function(x){ g_t_1(x,T_new_p[,1])-alpha_spending_P(1/2) }, c(1, 6))$root
        u_2_asP <- uniroot(function(x){ g_t_2(x,T_new_p[,1],T_new_p[,2],u_1_asP)-alpha+alpha_spending_P(1/2) }, c(1, 6))$root
        if ( vec[1]>u_1_asP ) {N_a_p_asP <- N_a_p_asP + 1}                    
        if ( vec[1]<u_1_asP & vec[2]>u_2_asP ) {N_a_p_asP <- N_a_p_asP + 1}  
        r <- 0
    },error=function(e){ 1 })
    if ( r > 0 ) { i_f_p_asP <- i_f_p_asP + 1 } ## record the number of interruptions

    ## Proposed method, asOF, directly using 1/2 as the information level
    r <- tryCatch({
        u_1_asOF <- uniroot(function(x){ g_t_1(x,T_new_p[,1])-alpha_spending_OF(1/2) }, c(1, 25))$root
        u_2_asOF <- uniroot(function(x){ g_t_2(x,T_new_p[,1],T_new_p[,2],u_1_asOF)-alpha+alpha_spending_OF(1/2) }, c(1, 25))$root
        if ( vec[1]>u_1_asOF ) {N_a_p_asOF <- N_a_p_asOF + 1}                    
        if ( vec[1]<u_1_asOF & vec[2]>u_2_asOF ) {N_a_p_asOF <- N_a_p_asOF + 1}  
        r <- 0
    },error=function(e){ 1 })
    if ( r > 0 ) {i_f_p_asOF <- i_f_p_asOF + 1 } ## record the number of interruptions
    ######
    ######
    ######

    ###### traditional method, asP
    u_1_norm_asP <- 2.156999
    u_2_norm_asP <- 2.200977
    ## These values are derived from the classical group sequential design using 1/2 as the information level. 
    ## (Pocock type alpha-spending function, two-stage design, one-side hypothesis testing, alpha=0.025)
    if ( vec[1]>u_1_norm_asP ) {N_a_norm_asP <- N_a_norm_asP + 1}  
    if ( vec[1]<u_1_norm_asP & vec[2]>u_2_norm_asP ) {N_a_norm_asP <- N_a_norm_asP + 1} 

    ###### t-approximation, asP
    df_1 <- getdf(data_t[1:n], data_c[1:n])
    df_2 <- getdf(data_t, data_c)
    u_1_t_app_asP <- qt(pnorm(2.156999), df_1)
    u_2_t_app_asP <- qt(pnorm(2.200977), df_2)
    if ( vec[1]>u_1_t_app_asP ) {N_a_t_app_asP <- N_a_t_app_asP + 1}
    if ( vec[1]<u_1_t_app_asP & vec[2]>u_2_t_app_asP ) {N_a_t_app_asP <- N_a_t_app_asP + 1}

    ###### traditional method, asOF
    u_1_norm_asOF <- 2.962588
    u_2_norm_asOF <- 1.968596
    ## These values are derived from the classical group sequential design using 1/2 as the information level. 
    ## (O’Brien-Fleming type alpha-spending function, two-stage design, one-side hypothesis testing, alpha=0.025)
    if ( vec[1]>u_1_norm_asOF ) {N_a_norm_asOF <- N_a_norm_asOF + 1} 
    if ( vec[1]<u_1_norm_asOF & vec[2]>u_2_norm_asOF ) {N_a_norm_asOF <- N_a_norm_asOF + 1} 

    ###### t-approximation, asOF
    df_1 <- getdf(data_t[1:n], data_c[1:n])
    df_2 <- getdf(data_t, data_c)
    u_1_t_app_asOF <- qt(pnorm(2.962588), df_1)
    u_2_t_app_asOF <- qt(pnorm(1.968596), df_2)
    if ( vec[1]>u_1_t_app_asOF ) {N_a_t_app_asOF <- N_a_t_app_asOF + 1} 
    if ( vec[1]<u_1_t_app_asOF & vec[2]>u_2_t_app_asOF ) {N_a_t_app_asOF <- N_a_t_app_asOF + 1}

}

#### Results -- empirical Type I error rate
N_a_norm_asP/nsim                    # Traditional method (using normal distribution) asP
N_a_t_app_asP/nsim                   # t-approximation asP
N_a_p_asP/(nsim-i_f_p_asP)           # Proposed method asP

N_a_norm_asOF/nsim                   # traditional method (using normal distribution) asOF
N_a_t_app_asOF/nsim                  # t-approximation asOF
N_a_p_asOF/(nsim-i_f_p_asOF)         # Proposed method asOF


#### Check the number of interruptions
i_f_p_asP                            # Proposed method asP
i_f_p_asOF                           # Proposed method asOF

####
####
# In cases where the sample size is extremely small (e.g., n=5), 
# an error may occur when computing the critical values using the permutation distribution. 
# To prevent any interruption of the program, we utilize the “tryCatch” function.
#
# This error, characterized by “f() values at end points not of opposite sign,” 
# stems from the use of the “uniroot” function. 
# The discreteness of the permutation distribution when dealing with very small sample sizes 
# may hinder the “uniroot” function’s ability to find a solution, resulting in an error.
#
# For n=5, there is approximately a 0.15% probability of encountering this error for the O’Brien-Fleming type alpha-spending function 
# and a 0.05% probability for the Pocock type alpha-spending function. 
# However, for n=10, the occurrence of this error is nearly negligible. 
# It is important to note that 
# this error does not impact the accuracy of the final results presented in the manuscript.
####
####


