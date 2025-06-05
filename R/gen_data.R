gen_data <-
function(n, real_type, sigT)
{
  if(real_type == 1)
  {
    ## Computer model
    f = function(t, x)
    {
      return(t*x)
    }
    
    # f = function(x, t)
    # {
    #   return(t*x)
    # }
    
    ## Computer model derivative
    g = function(t, x)
    {
      return(x)
      #return(x*sin(t*x))
    }
    
    ## Physical process
    y_R = function(x)
    {
      r = 4*x + x*sin(5*x)
      return(r)
    }
    ## Design points
    x = runif(n)
    y = y_R(x) + sigT*rnorm(n)
    result = list("y" = y, "x" = x, "f" = f, "g" = g)
  }
  
  if(real_type == 2)
  {
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
    
    ## Physical process
    y_R = function(x)
    {
      r = 7*sin(2*pi*0.2 - pi)^2 + 2*((2*pi*0.3 - pi)^2*sin(2*pi*x - pi)) 
      return(r)
    }
    
    ## Design points
    x = runif(n)
    y = y_R(x) + sigT*rnorm(n)
    result = list("y" = y, "x" = x, "f" = f, "g1" = g1, "g2" = g2)
  }
  
  if(real_type == 3)
  {
    ## Computer model ##
    f = function(t, x)
    {
      #r = 10*pi*sin(x[1]*x[2]) + 20*(x[3] - 0.5)^2 + t[1]*x[4] + t[2]*x[5] 
      r = t[1]*pi*sin(x[1]*x[2]) + t[2]*(x[3] - 0.5)^2 + 10*x[4] + 5*x[5] 
      return(r)
    }
    
    ## Computer model derivatives ##
    
    g1 = function(t, x)
    {
      #n = nrow(x)
      #r = x[,4]
      r = pi*sin(x[,1]*x[,2]) 
      return(r)
    }
    
    g2 = function(t, x)
    {
      #r = x[,5]
      r = (x[,3] - 0.5)^2
      return(r)
    }
    
    ## Physical process
    y_R = function(x)
    {
      r = 10*pi*sin(x[1]*x[2]) + 20*(x[3] - 0.5)^2 + 10*x[4] + 5*x[5] 
      return(r)
    }
    
    ## Design points
    if(real_type <=2){
      x = runif(n)
      y = y_R(x) + sigT*rnorm(n)
    }else if(real_type == 3){
      p = 5
      x = matrix(runif(n*p), n, p)
      y = rep(0, n)
     for(i in 1:n)
     {
       y[i] = y_R(x[i,]) + sigT*rnorm(1)
     }
    }
    result = list("y" = y, "x" = x, "f" = f, "g1" = g1, "g2" = g2)
  }
  
  return(result)
  
}
