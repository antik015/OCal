posterior_projection_q <-
function(Y, x, phi_sq, nmcmc, burnin, Sig = SigT, m, prior = "BART", computer_model_functions, nugget = 0, alpha = 0.95, MH_sd = 0.1, prior_sd)
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
