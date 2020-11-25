#' Returns sample size needed to achieve a specified power in the context of a
#' binary MRT.
#'
#' @param avail_pattern A vector of length T that is the average availability at
#'   each time point
#' @param f_t           Defines marginal excursion effect MEE(t) under
#'   alternative together with beta
#' @param g_t           Defines success probability null curve together with
#'   alpha
#' @param beta          Defines marginal excursion effect MEE(t) under
#'   alternative together with g_t
#' @param alpha         Defines success probability null curve together with f_t
#' @param p_t           Randomization probability at each time point
#' @param gamma         Desired Type I error
#' @param b             Desired Type II error
#' @param exact         Determines if exact n or ceiling will be returned
#'
#' @return              Power of the test.
#' @export
#'
#' @examples
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



