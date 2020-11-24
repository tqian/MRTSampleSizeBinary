# Calculate sample size with given power

## E[I_t]  # TQ: will assume this is vector of length T
## f(t)    # TQ: will assume this f_t object is matrix of dimension T * p
## g(t)    # TQ: will assume this g_t object is matrix of dimension T * q
## beta_0
## alpha_0
## p_t (randomization probability)
## gamma (type I error rate),
## b (Type II error rate)
## exact # returns exact n if true, else returns ceiling of sample size
calculate_mrt_bin_samplesize_f <- function(avail_pattern,  
                                 f_t,             
                                 g_t,             
                                 beta,            
                                 alpha,           
                                 p_t,             
                                 gamma,          
                                 b,
                                 exact=FALSE)               
{
    
    p <- length(beta)
    q <- length(alpha)
    
    m_and_sigma <- compute_m_sigma(avail_pattern, f_t, g_t, beta, alpha, p_t)
    m_matrix <- m_and_sigma$m
    sigma_matrix <- m_and_sigma$sigma
    
    # Set up the function that we will ultimately solve to get the sample size
    power_f <- function(n){
        
        right_hand_side <- pf(q=qf(p=(1-gamma), df1=p, df2=n-q-p), 
                              df1=p, 
                              df2=n-q-p, 
                              ncp=compute_ncp(n, beta, m_matrix, sigma_matrix))
        left_hand_side <- b
        return(right_hand_side - left_hand_side)
    }
    
    # find min n to achieve desired power
    sample_size <- uniroot(power_f, lower=p+q+1, upper=1000000)$root
    
    # round up
    if(exact == FALSE){
        sample_size <- ceiling(sample_size)
    }
    
    return(sample_size)
}

# helper function for calculate_mrt_bin_samplesize_f
# computes the non-centrality parameter for the F distribution when solving for
# sample size
compute_ncp <- function(x, beta, m_matrix, sigma_matrix){
    as.numeric(x * t(beta) %*%
                   solve(solve(m_matrix) %*%
                             sigma_matrix %*%
                             t(solve(m_matrix))) %*% 
                   beta)
}


# helper function for calculate_mrt_bin_samplesize_f
## E[I_t]  # TQ: will assume this is vector of length T
## f(t)    # TQ: will assume this f_t object is matrix of dimension T * p
## g(t)    # TQ: will assume this g_t object is matrix of dimension T * q
## beta_0
## alpha_0
## p_t
compute_m_sigma <- function(avail_pattern,  
                          f_t,             
                          g_t,            
                          beta,            
                          alpha,           
                          p_t              
) {
    p <- length(beta)
    q <- length(alpha)
    
    ## The M and Sigma matrices (needed to compute lambda)
    m_matrix <- matrix(data=0, nrow=length(beta), ncol=length(beta))
    sigma_matrix <- matrix(data=0, nrow=length(beta), ncol=length(beta))
    
    # For each decision point 
    # (T = total # of decision points is taken as the length of p_t, for now)
    for (i in 1:length(p_t)){
        
        # breaking down steps to identify bug and improve robustness of code
        this_f_t <- as.numeric(f_t[i, ])
        this_g_t <- as.numeric(g_t[i, ])
        
        stopifnot(length(this_f_t) == p)
        stopifnot(length(this_g_t) == q)
        
        # both are vectors, so use this way to calculate inner product
        this_f_t_times_beta <- sum(this_f_t * beta) 
        this_g_t_times_alpha <- sum(this_g_t * alpha)
        
        this_f_t_f_t <- outer(this_f_t, this_f_t) # this is f_t %*% f_t^T
        
        this_m <- as.numeric(avail_pattern[i] * 
                                 (exp(p_t[i] * this_f_t_times_beta)) * 
                                 exp(this_g_t_times_alpha) *
                                 (1-p_t[i]) *
                                 p_t[i]) * this_f_t_f_t
        
        # added a check of dimension
        stopifnot("matrix" %in% class(this_m) & all(dim(this_m) == c(p, p))) 
        
        # A running sum so that we end up with each entry being the sum of that
        # entry across all time points
        m_matrix <- m_matrix + this_m 
        
        this_sigma <- as.numeric(avail_pattern[i] * 
                                     (exp(2 * p_t[i] * this_f_t_times_beta)) *
                                     exp(this_g_t_times_alpha) * 
                                     (1-p_t[i]) * 
                                     p_t[i] * 
                                     ((1-p_t[i]) * 
                                        exp(-1 * this_f_t_times_beta) +
                                        p_t[i] - exp(this_g_t_times_alpha))) * 
                    this_f_t_f_t
        
        
        sigma_matrix <- sigma_matrix + this_sigma
        stopifnot("matrix" %in% class(this_sigma) & 
                      all(dim(this_sigma) == c(p, p)))
    }
    
    return(list(m = m_matrix, sigma = sigma_matrix))
    
}

