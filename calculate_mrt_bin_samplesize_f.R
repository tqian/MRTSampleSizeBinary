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
    
    # round up if non-exact size is requested
    if(exact == FALSE){
        sample_size <- ceiling(sample_size)
    }
    
    return(sample_size)
}



