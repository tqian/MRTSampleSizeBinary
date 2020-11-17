# Calculate sample size with given power
MRT_bin_samplesize.F <- function(avail_pattern,  ## E[I_t]  # TQ: will assume this is vector of length T
                                 f_t,             ## f(t)    # TQ: will assume this f_t object is matrix of dimension T * p
                                 g_t,             ## g(t)    # TQ: will assume this g_t object is matrix of dimension T * q
                                 beta,            ## beta_0
                                 alpha,           ## alpha_0
                                 p_t,             ## p_t (randomization probability)
                                 gamma,          ## gamma (type I error rate),
                                 b)               ## b (Type II error rate)
{
    
    p <- length(beta)
    q <- length(alpha)
    
    M_and_Sigma <- computeMSigma(avail_pattern, f_t, g_t, beta, alpha, p_t)
    M.Matrix <- M_and_Sigma$M
    Sigma.Matrix <- M_and_Sigma$Sigma
    
    ## Setting up the function that we will ultimately solve to get the sample size
    powerF <- function(n){
        ## Lambda as a function of x (i.e., the sample size)
        lambda <- function(x){
            as.numeric(x * t(beta) %*% solve(solve(M.Matrix) %*% Sigma.Matrix %*% t(solve(M.Matrix))) %*% beta)
        }
        
        righthand.side <- pf(q=qf(p=(1-gamma), df1=p, df2=n-q-p), df1=p, df2=n-q-p, ncp=lambda(n))
        lefthand.side = b
        return(righthand.side - lefthand.side)
    }
    
    sample_size <- uniroot(powerF, lower=p+q+1, upper=1000000)$root
    return(sample_size)
}


computeMSigma <- function(avail_pattern,  ## E[I_t]  # TQ: will assume this is vector of length T
                          f_t,             ## f(t)    # TQ: will assume this f_t object is matrix of dimension T * p
                          g_t,             ## g(t)    # TQ: will assume this g_t object is matrix of dimension T * q
                          beta,            ## beta_0
                          alpha,           ## alpha_0
                          p_t              ## p_t
) {
    p <- length(beta)
    q <- length(alpha)
    
    ## The M and Sigma matrices (needed to compute lambda)
    M.Matrix <- matrix(data=0, nrow=length(beta), ncol=length(beta))
    Sigma.Matrix <- matrix(data=0, nrow=length(beta), ncol=length(beta))
    for (i in 1:length(p_t)){## For each decision point (T = total # of decision points is taken as the length of p_t, for now)
        
        # breaking down steps to identify bug and improve robustness of code
        this_f_t <- as.numeric(f_t[i, ])
        this_g_t <- as.numeric(g_t[i, ])
        
        stopifnot(length(this_f_t) == p)
        stopifnot(length(this_g_t) == q)
        
        this_f_t_times_beta <- sum(this_f_t * beta) # both are vectors, so use this way to calculate inner product
        this_g_t_times_alpha <- sum(this_g_t * alpha)
        
        this_f_t_f_t <- outer(this_f_t, this_f_t) # this is f_t %*% f_t^T
        
        ThisM <- as.numeric(avail_pattern[i] * (exp(p_t[i] * this_f_t_times_beta)) * exp(this_g_t_times_alpha) * (1-p_t[i]) * p_t[i]) * this_f_t_f_t
        stopifnot("matrix" %in% class(ThisM) & all(dim(ThisM) == c(p, p))) # added a check of dimension
        
        M.Matrix <- M.Matrix + ThisM  ## A running sum so that we end up with each entry being the sum of that entry across all time points
        
        ThisSigma <- as.numeric(avail_pattern[i] * (exp(2 * p_t[i] * this_f_t_times_beta)) * exp(this_g_t_times_alpha) * (1-p_t[i]) * p_t[i] * ((1-p_t[i]) * exp(-1 * this_f_t_times_beta) + p_t[i] - exp(this_g_t_times_alpha))) * this_f_t_f_t
        Sigma.Matrix <- Sigma.Matrix + ThisSigma
        stopifnot("matrix" %in% class(ThisSigma) & all(dim(ThisSigma) == c(p, p)))
    }
    
    return(list(M = M.Matrix, Sigma = Sigma.Matrix))
    
}

MRT_bin_samplesize.F()

