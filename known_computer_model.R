
library(pracma)
library(mvtnorm)
library(BART)
library(cvCovEst)


posterior_projection = function(y, x, phi_sq, nmcmc, burnin, sig = sigT, m, prior = "BART", proj_type = "infinite", computer_model_functions, nugget = 0, alpha = 0.95, MH_sd = 0.1, prior_sd)
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

ogp_model1 = function(y, x, prior_sd, sigT, nmcmc, burnin)
{
  t1 = Sys.time()
  f = function(t, x)
  {
    return(t*x)
  }
  
  ## Computer model derivative
  g = function(t, x)
  {
    return(x)
  }
  
  calib_cov = function(x1, x2, phi_sq  = 1, psi = 1/2)
  {
    r = phi_sq*(1+ abs(x1 - x2)/psi)*(exp(-abs(x1 - x2)/psi))
    return(r)
  }
  
  l2_error = function(t)
  {
    f_t = f(t, x)
    e_sq = sum((y - f_t)^2)
    return(e_sq)
  }
  p = 1
  theta_tilde = optim(rep(0, p), l2_error)$par
  ## Plumee covariance ##
  
  f1 = function(x1, x2)
  {
    r = g(theta_tilde, x1)*g(theta_tilde, x2)*calib_cov(x1, x2)
  }
  H = integral2(f1, 0, 1, 0, 1)$Q
  
  ogp_calib_cov = function(x1, x2)
  {
    h = function(x)
    {
      f2 = function(u)
      {
        l = g(theta_tilde, u)*calib_cov(x, u)
        return(l)
      }
      r = integral(f2, 0, 1)
      return(r)
    }
    q = calib_cov(x1, x2) - h(x1)*h(x2)/H
    return(q)
  }
  
  cov_mat_ogp = function(x, phi_sq = 1)
  {
    n = length(x)
    K = matrix(0, n, n)
    for(i in 1:(n-1))
    {
      for(j in (i+1):n)
      {
        K[i,j] = K[j,i] = ogp_calib_cov(x[i], x[j])
      }
    }
    diag(K) = phi_sq
    return(K)
  }
  K_ogp = cov_mat_ogp(x)
  
  marg_cov_ogp = K_ogp + sigT^2*diag(n)
  
  theta_ogp_post_var = (1/prior_sd^2 + t(x)%*%solve(marg_cov_ogp)%*%x)^(-1)
  theta_ogp_post_mean = theta_ogp_post_var*(t(y)%*%solve(marg_cov_ogp)%*%x)
  
  #theta_ogp_post_mean = (t(x)%*%solve(marg_cov_ogp)%*%y)/t(x)%*%solve(marg_cov_ogp)%*%x 
  #theta_ogp_post_var = 1/(t(x)%*%solve(marg_cov_ogp)%*%x) 
  thetaout = rnorm(nmcmc - burnin, theta_ogp_post_mean, theta_ogp_post_var)
  t2 = Sys.time() -t1
  return(list("thetaout" = thetaout, "theta_mean" = theta_ogp_post_mean, "theta_var" = theta_ogp_post_var^2, "time" = t2))
}


ogp_model2 = function(y, x, prior_sd, MH_sd, sigT, nmcmc, burnin)
{
  n = length(x)
  t1 = Sys.time()
  ## Computer model ##
  f = function(t, x)
  {
    r = 7*sin(2*pi*t[1] - pi)^2 + 2*((2*pi*t[2] - pi)^2*sin(2*pi*x - pi)) 
    return(r)
  }
  
  ## Computer model derivatives ##
  
  g1 = function(t, x)
  {
    r = 28*pi*sin(2*pi*t[1] - pi)*cos(2*pi*t[1] - pi)
    return(r)
  }
  
  g2 = function(t, x)
  {
    r = 8*pi*(2*pi*t[2] - pi)*sin(2*pi*x - pi)
    return(r)
  }
  
  calib_cov = function(x1, x2, phi_sq  = 1, psi = 1/2)
  {
    r = phi_sq*(1+ abs(x1 - x2)/psi)*(exp(-abs(x1 - x2)/psi))
    return(r)
  }
  l2_error = function(t)
  {
    f_t = f(t, x)
    e_sq = sum((y - f_t)^2)
    return(e_sq)
  }
  p = 2
  theta_tilde = optim(rep(0, p), l2_error)$par
  H = matrix(0, 2, 2)
  H11 = function(x1, x2)
  {
    r = g1(theta_tilde, x1)*g1(theta_tilde, x2)*calib_cov(x1, x2)
    return(r)
  }
  H12 = function(x1, x2)
  {
    r = g1(theta_tilde, x1)*g2(theta_tilde, x2)*calib_cov(x1, x2)
    return(r)
  }
  H22 = function(x1, x2)
  {
    r = g2(theta_tilde, x1)*g2(theta_tilde, x2)*calib_cov(x1, x2)
    return(r)
  }
  H[1,1] = integral2(H11, 0, 1, 0, 1)$Q
  H[1,2] = H[2,1] = integral2(H12, 0, 1, 0, 1)$Q
  H[2,2] = integral2(H22, 0, 1, 0, 1)$Q
  
  ogp_calib_cov = function(x1, x2)
  {
    h1 = function(x)
    {
      f22 = function(u)
      {
        l = g1(theta_tilde, u)*calib_cov(x, u)
        return(l)
      }
      r = integral(f22, 0, 1)
      return(r)
    }
    
    h2 = function(x)
    {
      f23 = function(u)
      {
        l = g2(theta_tilde, u)*calib_cov(x, u)
        return(l)
      }
      r = integral(f23, 0, 1)
      return(r)
    }
    q = calib_cov(x1, x2) - c(h1(x1), h2(x1))%*%solve(H)%*%c(h1(x2), h2(x2))
    return(q)
  }
  
  cov_mat_ogp = function(x)
  {
    n = length(x)
    K = matrix(0, n, n)
    for(i in 1:(n-1))
    {
      for(j in (i+1):n)
      {
        K[i,j] = K[j,i] = ogp_calib_cov(x[i], x[j])
      }
    }
    diag(K) = phi_sq
    return(K)
  }
  K_ogp = cov_mat_ogp(x)
  marg_cov_ogp = K_ogp + sigT^2*diag(n)
  
  
  ## Start MCMC for theta ##
  theta = theta_tilde + rnorm(p, sd = MH_sd)
  f_theta = f(theta, x)
  thetaout = matrix(0, p, nmcmc - burnin)
  MH_count = 0
  for(ii in 1:nmcmc)
  {
    theta_cand = theta + rnorm(p, sd = MH_sd)
    f_theta_cand = f(theta_cand, x)
    num = dmvnorm(as.vector(y), mean = as.vector(f_theta_cand), sigma = marg_cov_ogp, log = T) + sum(dnorm(theta_cand , mean = 0, sd = prior_sd, log = T))
    + sum(dnorm(theta, mean = theta_cand, sd = MH_sd, log = T))
    
    den = dmvnorm(as.vector(y), mean = as.vector(f_theta), sigma = marg_cov_ogp, log = T) + sum(dnorm(theta, mean = 0, sd = prior_sd, log = T))
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
}

posterior_projection_q = function(Y, x, phi_sq, nmcmc, burnin, Sig = SigT, m, prior = "BART", computer_model_functions, nugget = 0, alpha = 0.95, MH_sd = 0.1, prior_sd)
{
  t1 = Sys.time()
  n = length(x)
  
  #if(proj_type == "infinite")
  #{
  ## Quadrature points
  quad_points = gaussLegendre(m, 0, 1)
  x_t = quad_points$x
  n_t = length(x_t)
  w_t = quad_points$w
  #}
  
  q = length(computer_model_functions)
  p = length(computer_model_functions[[1]]) - 1

  
  for(k in 1:q)
  {
    assign(paste("f", k, sep = ""), computer_model_functions[[k]][[1]])
    
    for(j in 1:p)
    {
      assign(paste("g", k, j, sep = ""), computer_model_functions[[k]][[j+1]])
    }
  }
  
  l2_error = function(t)
  {
    f_t = matrix(0, n, q)
    e_sq = 0
    for(k in 1:q)
    {
      f_t[,k] = eval(parse(text = paste("f", k,"(", "t", ",", "x", ")", sep="")))
      e_sq = e_sq + sum((Y[,k] - f_t[,k])^2)
    }
    return(e_sq)
  }
  theta_tilde = optim(rep(0, p), l2_error)$par
  
  ## Initial values 
  theta = theta_tilde + rnorm(p, sd = MH_sd)
  f_theta = matrix(0, n, q)
  B = matrix(0, n, q)
  for(k in 1:q)
  {
    f_theta[,k] = do.call(paste("f", k, sep = ""), list(theta, x))
  }
  SigT_inv = solve(SigT)
  sd_mat = sqrtm(SigT_inv)
  SigT_half_inv = sd_mat$B
  SigT_half = sd_mat$Binv
  sigT = 1
  Y1 = Y%*%SigT_half_inv
  f_theta1 = f_theta%*%SigT_half_inv
  
  log_post = function(y, f_theta, B)
  {
    M = y - f_theta - B
    r = sum(as.vector(M), log = T)
  }
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
    K_all = cov_mat(c(x, x_t), phi_sq = phi_sq) + nugget*diag(n+m)
    K = K_all[1:n, 1:n]
    K_t = K_all[(n+1):(n+m), (n+1):(n+m)]
    K12 = K_all[1:n, (n+1):(n+m)]
    K_inv = solve(K)
    
    P_mat = solve(K + sigT^2*diag(n))
    Q_mat = solve(K_inv + (1/sigT)^2*diag(n))
    R_mat = K_t - t(K12)%*%P_mat%*%K12
    
  }
  ## Storage
  thetaout = matrix(0, p, nmcmc - burnin)
  MH_count = 0
  
  ## Start MCMC
  for(ii in 1:nmcmc)
  {
    
    ## Update bias function
    Y_tilde = Y1 - f_theta1
    
    if(prior == "GP")
    {
      B_test = matrix(0, n_t, q)
      for(k in 1:q)
      {
        post_cov = (Q_mat + t(Q_mat))/2
        post_mean = post_cov%*%(Y_tilde[,k]/sigT^2)
        B[,k] = rmvnorm(1, post_mean, post_cov)
        test_cov = (R_mat + t(R_mat))/2
        test_mean = t(K12)%*%P_mat%*%Y_tilde[,k]
        B_test[,k] = rmvnorm(1, test_mean, test_cov)
      }
      
      G = array(0, c(length(x_t), p, q))
      eta = rep(0, p)
      g_mat = array(0, c(length(x), p, q))
      for(j in 1:p)
      {
        for(k in 1:q)
        {
          eta[j] = eta[j] + sum(do.call(paste("g", k, j, sep = ""), list(theta_tilde, x_t))*B_test[,k]*w_t)
          G[,j, k] =   do.call(paste("g", k, j, sep = ""), list(theta_tilde, x_t))
          g_mat[,j,k] = do.call(paste("g", k, j, sep = ""), list(theta_tilde, x))
        }
      }
      Q = matrix(0, p, p)
      for(k in 1:q)
      {
        G[,,k] = diag(sqrt(w_t))%*%G[,,k]
        Q = Q + t(G[,,k])%*%G[,,k]
      }
      lambda = solve(Q)%*%eta
      B_star = matrix(0, n, q)
      for(k in 1:q)
      {
         B_star[,k] = B[,k] - t(lambda)%*%t(g_mat[,,k])
      }
    }
    if(prior == "BART")
    {
      B_test = matrix(0, n_t, q)
      for(k in 1:q)
      {
        invisible(capture.output(bfit <- wbart(x, Y_tilde[,k], x_t, sigest = 1, ndpost = 1, nskip = 0, nkeeptrain = 1, nkeeptest = 1, nkeeptestmean = 1, nkeeptreedraws = 1)))
        B[,k] = bfit$yhat.train
        B_test[,k] = bfit$yhat.test
      }
      
      G = array(0, c(length(x_t), p, q))
      eta = rep(0, p)
      g_mat = array(0, c(length(x), p, q))
      for(j in 1:p)
      {
        for(k in 1:q)
        {
          eta[j] = eta[j] + sum(do.call(paste("g", k, j, sep = ""), list(theta_tilde, x_t))*B_test[,k]*w_t)
          G[,j, k] =   do.call(paste("g", k, j, sep = ""), list(theta_tilde, x_t))
          g_mat[,j,k] = do.call(paste("g", k, j, sep = ""), list(theta_tilde, x))
        }
      }
      Q = matrix(0, p, p)
      for(k in 1:q)
      {
        G[,,k] = diag(sqrt(w_t))%*%G[,,k]
        Q = Q + t(G[,,k])%*%G[,,k]
      }
      lambda = solve(Q)%*%eta
      B_star = matrix(0, n, q)
      for(k in 1:q)
      {
        B_star[,k] = B[,k] - t(lambda)%*%t(g_mat[,,k])
      }
    }
    
    
    ## Update theta via Random-Walk Metropolis-Hastings ##
    theta_cand = theta + rnorm(p, sd = MH_sd)
    f_theta_cand = matrix(0, n, q)
    for(k in 1:q)
    {
      f_theta_cand[,k] = do.call(paste("f", k, sep = ""), list(theta_cand, x))
    }
    f_theta_cand1 = f_theta_cand%*%SigT_half_inv
    Y_star = Y1 - B_star
    #Y_star = Y - B_star
    num = sum(dnorm(as.vector(Y_star), mean = as.vector(f_theta_cand1), sd = sigT, log = T)) + sum(dnorm(theta_cand , mean = 0, sd = prior_sd, log = T)) + sum(dnorm(theta, mean = theta_cand, sd = MH_sd, log = T))
    
    den = sum(dnorm(as.vector(Y_star), mean = as.vector(f_theta1), sd = sigT, log = T)) + sum(dnorm(theta, mean = 0, sd = prior_sd, log = T)) + sum(dnorm(theta_cand, mean = theta, sd = MH_sd, log = T))
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
    
    for(k in 1:q)
    {
      f_theta[,k] = do.call(paste("f", k, sep = ""), list(theta, x))
    }
    f_theta1 = f_theta%*%SigT_half_inv
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


# n = 100
# real_type = 2
# sigT = 0.02
# sim_data = gen_data(n, real_type, sigT)
# y = sim_data$y
# x = sim_data$x
# computer_model_functions = list()
# computer_model_functions[[1]] = sim_data$f
# computer_model_functions[[2]] = sim_data$g
# computer_model_functions[[3]] = sim_data$g2
# 
# 
# prior = "GP"
# proj_type = "infinite"
# nmcmc = 5000
# burnin = 1000
# phi_sq = 1
# m = 30
# nugget = 0.02
# alpha = 0.95
# MH_sd = 0.2 ## Reasonable for model 2 is 0.001, for model 1 0.05
# prior_sd = 10
# res_proj = posterior_projection(y, x, phi_sq, nmcmc, burnin, sigT, m, prior, proj_type, computer_model_functions, nugget, alpha, MH_sd, prior_sd)
# 
# res_ogp = ogp_model1(y,x, prior_sd)
