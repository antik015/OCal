posterior_projection <-
function(y, x, phi_sq, nmcmc, burnin, sig = sigT, m, prior = "BART", proj_type = "infinite", computer_model_functions, nugget = 0, alpha = 0.95, MH_sd = 0.1, prior_sd)
{
  t1 = Sys.time()
  n = length(x)
  
  #if(proj_type == "infinite")
  #{
    ## Quadrature points
    quad_points = gaussLegendre(m, 0, 1)
    x_t = quad_points$x
    w_t = quad_points$w
  #}
  
  q = length(computer_model_functions)
  p = q - 1
  f = computer_model_functions[[1]]
  for(j in 2:q)
  {
    assign(paste("g", j - 1, sep = ""), computer_model_functions[[j]])
  }
  
  l2_error = function(t)
  {
    f_t = f(t, x)
    e_sq = sum((y - f_t)^2)
    return(e_sq)
  }
  theta_tilde = optim(rep(0, p), l2_error)$par
  
  ## Initial values 
  theta = theta_tilde + rnorm(p, sd = MH_sd)
  f_theta = f(theta, x)
  if(prior == "GP")
  {
    calib_cov = function(x1, x2, phi_sq = phi_sq, psi = 1/2)
    {
      r = phi_sq*(1+ abs(x1 - x2)/psi)*(exp(-abs(x1 - x2)/psi))
      return(r)
    }
    cov_mat = function(x, phi_sq = phi_sq)
    {
      n = length(x)
      K = matrix(0, n, n)
      for(i in 1:(n-1))
      {
        for(j in (i+1):n)
        {
          K[i,j] = K[j,i] = calib_cov(x[i], x[j], phi_sq)
        }
      }
      diag(K) = phi_sq
      return(K)
    }
    K_all = cov_mat(c(x, x_t), phi_sq) + nugget*diag(n+m)
    K = K_all[1:n, 1:n]
    K_t = K_all[(n+1):(n+m), (n+1):(n+m)]
    K12 = K_all[1:n, (n+1):(n+m)]
    K_inv = solve(K)
    
    P_mat = solve(K + sig^2*diag(n))
    Q_mat = solve(K_inv + (1/sig)^2*diag(n))
    R_mat = K_t - t(K12)%*%P_mat%*%K12
    
  }
  ## Storage
  thetaout = matrix(0, p, nmcmc - burnin)
  MH_count = 0
  
  ## Start MCMC
  for(ii in 1:nmcmc)
  {
    
    ## Update bias function
    y_tilde = y - f_theta
    
    if(prior == "GP" && proj_type == "infinite")
    {
      post_cov = (Q_mat + t(Q_mat))/2
      post_mean = post_cov%*%(y_tilde/sig^2)
      b = rmvnorm(1, post_mean, post_cov)
      
      test_cov = (R_mat + t(R_mat))/2
      test_mean = t(K12)%*%P_mat%*%y_tilde
      b_test = rmvnorm(1, test_mean, test_cov)
      
      G = matrix(0, length(x_t), p)
      eta = rep(0, p)
      g_mat = matrix(0, n, p)
      for(j in 1:p)
      {
        eta[j] = sum(do.call(paste("g", j, sep = ""), list(theta_tilde, x_t))*b_test*w_t)
        G[,j] = do.call(paste("g", j, sep = ""), list(theta_tilde, x_t))
        g_mat[,j] = do.call(paste("g", j, sep = ""), list(theta_tilde, x))
      }
      G = diag(sqrt(w_t))%*%G
      lambda = solve(t(G)%*%G)%*%eta
      b_star = as.vector(b - t(lambda)%*%t(g_mat))
    }
    
    if(prior == "BART" && proj_type == "infinite")
    {
      invisible(capture.output(bfit <- wbart(x, y_tilde, x_t, sigest = sig, ndpost = 1, nskip = 0, nkeeptrain = 1, nkeeptest = 1, nkeeptestmean = 1, nkeeptreedraws = 1)))
      b = bfit$yhat.train
      b_test = bfit$yhat.test
      
      G = matrix(0, length(x_t), p)
      eta = rep(0, p)
      g_mat = matrix(0, n, p)
      for(j in 1:p)
      {
        eta[j] = sum(do.call(paste("g", j, sep = ""), list(theta_tilde, x_t))*b_test*w_t)
        G[,j] = do.call(paste("g", j, sep = ""), list(theta_tilde, x_t))
        g_mat[,j] = do.call(paste("g", j, sep = ""), list(theta_tilde, x))
      }
      G = diag(sqrt(w_t))%*%G
      lambda = solve(t(G)%*%G)%*%eta
      b_star = as.vector(b - t(lambda)%*%t(g_mat))
    }
    
    if(prior == "GP" && proj_type == "finite")
    {
      post_cov = (Q_mat + t(Q_mat))/2
      post_mean = post_cov%*%(y_tilde/sig^2)
      b = rmvnorm(1, post_mean, post_cov)
      
      g_mat = matrix(0, n, p)
      for(j in 1:p)
      {
        g_mat[,j] = do.call(paste("g", j, sep = ""), list(theta_tilde, x))
      }
      post_cov_half = sqrtm(post_cov)
      b_star = as.vector(post_cov_half$B%*%(diag(n) - g_mat%*%solve(t(g_mat)%*%g_mat)%*%t(g_mat))%*%post_cov_half$Binv%*%t(b))
    }

    if(prior == "BART" && proj_type == "finite")
    {
      invisible(capture.output(bfit <- wbart(x, y_tilde, x_t, sigest = sig, ndpost = 5*n, nskip = 5*n, nkeeptrain = 5*n, nkeeptest = 5*n, nkeeptestmean = 5*n, nkeeptreedraws = 5*n)))
      b = bfit$yhat.train
      b_test = bfit$yhat.test
      #post_cov = linearShrinkLWEst(b)
      post_cov = alpha*cov(b) + (1-alpha)*diag(n)
      
      g_mat = matrix(0, n, p)
      for(j in 1:p)
      {
        g_mat[,j] = do.call(paste("g", j, sep = ""), list(theta_tilde, x))
      }
      
      post_cov_half = sqrtm(post_cov)
      b_star = as.vector(post_cov_half$B%*%(diag(n) - g_mat%*%solve(t(g_mat)%*%g_mat)%*%t(g_mat))%*%post_cov_half$Binv%*%b[1,])
    }
    
    ## Update theta via Random-Walk Metropolis-Hastings ##
    theta_cand = theta + rnorm(p, sd = MH_sd)
    f_theta_cand = f(theta_cand, x)
    y_star = y - b_star
    num = sum(dnorm(y_star, mean = f_theta_cand, sd = sig, log = T)) + sum(dnorm(theta_cand , mean = 0, sd = prior_sd, log = T))
            + sum(dnorm(theta, mean = theta_cand, sd = MH_sd, log = T))
    
    den = sum(dnorm(y_star, mean = f_theta, sd = sig, log = T)) + sum(dnorm(theta, mean = 0, sd = prior_sd, log = T))
            + sum(dnorm(theta_cand, mean = theta, sd = MH_sd, log = T))
    if(log(runif(1))< min(0, num - den)) {
      theta = theta_cand
      MH_count = MH_count + 1
    }
    if(ii%%50 == 0 && MH_count/ii > 0.3)
    {
      MH_sd = MH_sd*exp(min(0.01, ii^(-0.5)))
    }else if(ii%%50 == 0 && MH_count/ii < 0.3){
      MH_sd = MH_sd*exp(-min(0.01, ii^(-0.5)))
    }
    
    f_theta = f(theta, x)
    
    if(ii%%100 == 0) print(ii)
    
    if(ii > burnin)
    {
      thetaout[,ii - burnin] = theta
    }
  }
  t2 = Sys.time() - t1 
  result = list("thetaout" = thetaout, "time" = t2, "MH_acceptance" = MH_count/nmcmc)
  return(result)
}
