
if (0) {
    total_T <- 30
    p10 <- 0.5
    pT0 <- 0.5
    p11 <- 0.6
    pT1 <- 0.6
    alpha_shape <- "constant"
    beta_shape <- "constant"
    rand_prob <- rep(0.4, total_T)
    avail_pattern <- rep(0.7, total_T)
    typeIerror <- 0.05
    
    calculate_sample_size_binary_mrt_wrapper(p10, pT0, p11, pT1,  total_T,
                                         alpha_shape, 
                                         beta_shape, 
                                         rand_prob, 
                                         avail_pattern,
                                         type_I_error,
                                         0.8)
    
    calculate_power_binary_mrt_wrapper(p10, pT0, p11, pT1,  total_T,
                                         alpha_shape, 
                                         beta_shape, 
                                         rand_prob, 
                                         avail_pattern,
                                         type_I_error,
                                         35)
}

# Calculate sample size with given power: wrapper function
calculate_sample_size_binary_mrt_wrapper <- function(p10, pT0, p11, pT1,
                                                 total_T,
                                                 alpha_shape = c("constant", "loglinear"),
                                                 beta_shape = c("constant", "loglinear"),
                                                 rand_prob,  ## p_t
                                                 avail_pattern, ## E[I_t]  # TQ: will assume this is vector of length T
                                                 type_I_error,
                                                 power) {
    alpha_shape <- match.arg(alpha_shape)
    beta_shape <- match.arg(beta_shape)
    
    alpha_beta_ate <- compute_alpha_beta_from_prob(p10, pT0, p11, pT1, total_T, alpha_shape, beta_shape)
    
    if (alpha_shape == "constant") {
        
        g_t <- matrix(1, nrow = total_T, ncol = 1)
        
    } else if (alpha_shape == "loglinear") {
        
        g_t <- cbind(1, 1:total_T)
    }
    
    if (beta_shape == "constant") {
        
        f_t <- matrix(1, nrow = total_T, ncol = 1)
        
    } else if (beta_shape == "loglinear") {
        
        f_t <- cbind(1, 1:total_T)
    }
    
    sample_size <- mrt_bin_sample_size_f(avail_pattern = avail_pattern,
                                        f_t = f_t,
                                        g_t = g_t,
                                        beta = alpha_beta_ate$beta,
                                        alpha = alpha_beta_ate$alpha,
                                        p_t = rand_prob,
                                        gamma = typeIerror,
                                        b = 1 - power)
    return(ceiling(sample_size))
}


# Calculate sample size with given power
mrt_bin_sample_size_f <- function(avail_pattern,  ## E[I_t]  # TQ: will assume this is vector of length T
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
    
    m_and_sigma <- compute_m_sigma(avail_pattern, f_t, g_t, beta, alpha, p_t)
    m_matrix <- m_and_sigma$m      #TD : might need to change $M to $m but not sure -----------
    sigma_matrix <- m_and_sigma$sigma
    
    ## Setting up the function that we will ultimately solve to get the sample size
    power_f <- function(n){
        ## Lambda as a function of x (i.e., the sample size)
        lambda <- function(x){
            
            as.numeric(x * t(beta) %*% solve(solve(M.Matrix) %*% 
                                                 Sigma.Matrix %*%
                                                 t(solve(M.Matrix))) %*%
                                                  beta)
        }
        
        right_hand_side <- pf(q = qf(p = (1-gamma),
                                     df1 = p, 
                                     df2 = n-q-p), 
                                     df1 = p, 
                                     df2 = n-q-p,
                                     ncp = lambda(n))
        left_hand_side = b
        return(right_hand_side - left_hand_side)
    }
    
    sample_size <- uniroot(powerF, lower=p+q+1, upper=1000000)$root
    return(sample_size)
}

# Calculate sample size with given sample size: wrapper function
calculate_power_binary_mrt_wrapper <- function(p10, pT0, p11, pT1,
                                      total_T,
                                      alpha_shape = c("constant", "loglinear"),
                                      beta_shape = c("constant", "loglinear"),
                                      rand_prob,  ## p_t
                                      avail_pattern, ## E[I_t]  # TQ: will assume this is vector of length T
                                      type_I_error,
                                      sample_size) {
    
    alpha_shape <- match.arg(alpha_shape)
    beta_shape <- match.arg(beta_shape)
    
    alpha_beta_ate <- compute_alpha_beta_from_prob(p10, pT0, p11, pT1, total_T, alpha_shape, beta_shape)
    
    if (alpha_shape == "constant") {
        
        g_t <- matrix(1, nrow = total_T, ncol = 1)
        
    } else if (alpha_shape == "loglinear") {
        
        g_t <- cbind(1, 1:total_T)
    }
    
    if (beta_shape == "constant") {
        
        f_t <- matrix(1, nrow = total_T, ncol = 1)
        
    } else if (beta_shape == "loglinear") {
        
        f_t <- cbind(1, 1:total_T)
    }
    
    power <- mrt_bin_power_f(avail_pattern = avail_pattern,
                             f_t = f_t,
                             g_t = g_t,
                             beta = alpha_beta_ate$beta,
                             alpha = alpha_beta_ate$alpha,
                             p_t = rand_prob,
                             gamma = type_I_error,
                             sample_size = sample_size)
    return(power)
}

# Calculate the power with given sample size
mrt_bin_power_f <- function(avail_pattern,  ## E[I_t]  # TQ: will assume this is vector of length T
                            f_t,             ## f(t)    # TQ: will assume this f_t object is matrix of dimension T * p
                            g_t,             ## g(t)    # TQ: will assume this g_t object is matrix of dimension T * q
                            beta,            ## beta_0
                            alpha,           ## alpha_0
                            p_t,             ## p_t (randomization probability)
                            gamma,          ## gamma (type I error),
                            sample_size)     ## sample size
{
    
    p <- length(beta)
    q <- length(alpha)
    
    m_and_sigma <- compute_m_sigma(avail_pattern, f_t, g_t, beta, alpha, p_t)
    m_matrix <- m_and_sigma$m
    sigma_matrix <- m_and_sigma$Sigma
    
    ## Setting up the function that we will ultimately solve to get the sample size
    n <- sample_size
    lambda <- as.numeric(n * t(beta) %*%
                             solve(solve(m_matrix) %*% 
                                       sigma_matrix %*% 
                                       t(solve(M.Matrix))) %*% 
                                       beta)
         b <- pf(q = qf(p = (1-gamma),
                 df1 = p, 
                 df2 = n-q-p),
                 df1 = p,
                 df2 = n-q-p,
                 ncp = lambda)
    
    return(1 - b)
}


compute_m_sigma <- function(avail_pattern,  ## E[I_t]  # TQ: will assume this is vector of length T
                          f_t,             ## f(t)    # TQ: will assume this f_t object is matrix of dimension T * p
                          g_t,             ## g(t)    # TQ: will assume this g_t object is matrix of dimension T * q
                          beta,            ## beta_0
                          alpha,           ## alpha_0
                          p_t              ## p_t
) {
    p <- length(beta)
    q <- length(alpha)
    
    ## The M and Sigma matrices (needed to compute lambda)
    m_matrix <- matrix(data=0, 
                       nrow = length(beta), 
                       ncol = length(beta))
    
    sigma_matrix <- matrix(data = 0,
                           nrow = length(beta), 
                           ncol = length(beta))
    for (i in 1:length(p_t)){## For each decision point (T = total # of decision points is taken as the length of p_t, for now)
        
        # breaking down steps to identify bug and improve robustness of code
        this_f_t <- as.numeric(f_t[i, ])
        this_g_t <- as.numeric(g_t[i, ])
        
        stopifnot(length(this_f_t) == p)
        stopifnot(length(this_g_t) == q)
        
        this_f_t_times_beta <- sum(this_f_t * beta) # both are vectors, so use this way to calculate inner product
        this_g_t_times_alpha <- sum(this_g_t * alpha)
        
        this_f_t_f_t <- outer(this_f_t, this_f_t) # this is f_t %*% f_t^T
        
        this_m <- as.numeric(avail_pattern[i] * (exp(p_t[i] * this_f_t_times_beta)) * exp(this_g_t_times_alpha) * (1-p_t[i]) * p_t[i]) * this_f_t_f_t
        stopifnot("matrix" %in% class(ThisM) & all(dim(ThisM) == c(p, p))) # added a check of dimension
        
        m_matrix <- m_matrix + this_m  ## A running sum so that we end up with each entry being the sum of that entry across all time points
        
        this_sigma <- as.numeric(avail_pattern[i] *
                                     (exp(2 * p_t[i] * this_f_t_times_beta)) *
                                     exp(this_g_t_times_alpha) * (1-p_t[i]) * p_t[i] * ((1-p_t[i]) * exp(-1 * this_f_t_times_beta) + p_t[i] - exp(this_g_t_times_alpha))) * this_f_t_f_t
        sigma_matrix <- sigma_matrix + this_sigma
        stopifnot("matrix" %in% class(this_sigma) & all(dim(ThisSigma) == c(p, p)))
    }
    
    return(list(m = m_matrix, sigma = sigma_matrix))
    
}

## For computing the generative model from user inputs
compute_alpha_beta_from_prob <- function(p10, pT0, p11, pT1, total_T,
                                         alpha_shape = c("constant", "loglinear"),
                                         beta_shape = c("constant", "loglinear")) {
    
    #log(p11) = alpha_0 + alpha_1 + beta_0 + beta_1 
    #log(p10) = alpha_0 + alpha_1
    #log(pT1) = alpha_0 + total_T * alpha_1 + beta_0 + total_T * beta_1
    #log(pT0) = alpha_0 + total_T * alpha_1)
    
    #alpha_0 + alpha_1 = log(p10) => alpha_0 = log(p10) - alpha_1
    #log(pt0) = log(p10) - alpha_1 + total_T * alpha_1 = (total_T - 1)*alpha_1 + log(p10) => alpha_1 = log(pt0/p10) / (total_T - 1)
    
    alpha_shape <- match.arg(alpha_shape)
    beta_shape <- match.arg(beta_shape)
    
    alpha_1 <- log(pT0/p10) / (total_T - 1)
    alpha_0 <- log(p10) - alpha_1
    
    #beta_0 = log(p11) - (alpha_1 + alpha_0 + beta_1)
    #beta_0 = log(pT1) - alpha_0 - total_T * alpha_1 - total_T * beta_1
    
    eq <- function(beta_1) {
        log(pT1) - alpha_0 - total_T * alpha_1 - total_T * beta_1 - (log(p11) - (alpha_1 + alpha_0 + beta_1))
    }
    
    beta_1 <- uniroot(eq, interval = c(-4, 5), extendInt = "yes")$root
    beta_0 <- log(p11) - (alpha_1 + alpha_0 + beta_1)
    ate <-  sum(exp(alpha_0 + alpha_1 * (1:total_T) + beta_0 + beta_1 * (1:total_T))) /
        sum(exp(alpha_0 + alpha_1 * (1:total_T)))
    
    if (alpha_shape == "constant") {
        alpha <- alpha_0
    } else if (alpha_shape == "loglinear") {
        alpha <- c(alpha_0, alpha_1)
    }
    
    if (beta_shape == "constant") {
        beta <- beta_0
    } else if (beta_shape == "loglinear") {
        beta <- c(beta_0, beta_1)
    }
    
    return(list(alpha = alpha, beta = beta, ate = ate))
}