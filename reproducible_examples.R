### Reproducible examples ###
install.packages("pracma")
install.packages("BART")
install.packages("BASS")
install.packages("mvtnorm")
install.packages("ggplot2")
install.packages("latex2exp")
install.packages("cvCovEst")
library(pracma)
library(BART)
library(BASS)
library(mvtnorm)
library(cvCovEst)
library(OCal)


### Code to reproduce figures ###
library(ggplot2)
library(latex2exp)
font_size = 22
My_Theme = theme(axis.title.x = element_text(size = font_size),axis.text.x = element_text(size = font_size),
                 axis.title.y = element_text(size = font_size), axis.text.y = element_text(size = font_size), 
                 legend.text = element_text(size = font_size), legend.title = element_text(size = font_size), plot.title = element_text(size = 22))

### Figure 1(a) ###
n = 100
sigT = 0.2
real_type = 1

sim_data = gen_data(n, real_type, sigT)
y = sim_data$y
x = sim_data$x
computer_model_functions = list()
computer_model_functions[[1]] = sim_data$f
computer_model_functions[[2]] = sim_data$g
nmcmc = 5000
burnin = 1000
phi_sq = 1
m = 15
nugget = 0
alpha = 0.95
MH_sd = 0.2 
prior_sd = 10
phi_sq = 1
# res1 = posterior_projection(y, x, phi_sq, nu = 1/2, psi = 1/2, nmcmc, burnin, sigT, m, "GP", "infinite", computer_model_functions, nugget, alpha, MH_sd, prior_sd)
# res2 = posterior_projection(y, x, phi_sq, nu = 3/2, psi = 1/2, nmcmc, burnin, sigT, m, "GP", "infinite", computer_model_functions, nugget, alpha, MH_sd, prior_sd)
# res3 = posterior_projection(y, x, phi_sq, nu = 5/2, psi = 1/2, nmcmc, burnin, sigT, m, "GP", "infinite", computer_model_functions, nugget, alpha, MH_sd, prior_sd)
res4 = posterior_projection(y, x, phi_sq, nu = 3/2, psi = 1/2, nmcmc, burnin, sigT, m, "BART", "infinite", computer_model_functions, nugget, alpha, MH_sd, prior_sd)

df = data.frame(theta = t(res4$thetaout))
pl1 = ggplot(df) + geom_line(aes(x = 1:4000, y = theta), linewidth = 0.5) + xlab("MCMC iterations") + ylab(TeX("$\\theta$"))
pl1 = pl1 + geom_hline(yintercept = 3.56, col = "red", linewidth = 0.5) 
pl1 = pl1 + ggtitle("Model 1") + My_Theme
pl1


### Figure 1(b) and 1(c) ###

real_type = 2
sim_data = gen_data(n, real_type, sigT)
y = sim_data$y
x = sim_data$x
computer_model_functions = list()
computer_model_functions[[1]] = sim_data$f
computer_model_functions[[2]] = sim_data$g1
computer_model_functions[[3]] = sim_data$g2
nmcmc = 5000
burnin = 1000
phi_sq = 1
m = 15
nugget = 0
alpha = 0.95
MH_sd = 0.001 
prior_sd = 10
phi_sq = 1
# res1 = posterior_projection(y, x, phi_sq, nu = 1/2, psi = 1/2, nmcmc, burnin, sigT, m, "GP", "infinite", computer_model_functions, nugget, alpha, MH_sd, prior_sd)
# res2 = posterior_projection(y, x, phi_sq, nu = 3/2, psi = 1/2, nmcmc, burnin, sigT, m, "GP", "infinite", computer_model_functions, nugget, alpha, MH_sd, prior_sd)
# res3 = posterior_projection(y, x, phi_sq, nu = 5/2, psi = 1/2, nmcmc, burnin, sigT, m, "GP", "infinite", computer_model_functions, nugget, alpha, MH_sd, prior_sd)
res4 = posterior_projection(y, x, phi_sq, nu = 3/2, psi = 1/2, nmcmc, burnin, sigT, m, "BART", "infinite", computer_model_functions, nugget, alpha, MH_sd, prior_sd)

df = data.frame(theta1 = res4$thetaout[1,], theta2 = res4$thetaout[2,])
pl2 = ggplot(df) + geom_line(aes(x = 1:4000, y = theta1), linewidth = 0.5) + xlab("MCMC iterations") + ylab(TeX("$\\theta_1$"))
pl2 = pl2 + geom_hline(yintercept = 0.2, col = "red", linewidth = 0.5) 
pl2 = pl2 + ggtitle("Model 2") + My_Theme
pl2

pl3 = ggplot(df) + geom_line(aes(x = 1:4000, y = theta2), linewidth = 0.5) + xlab("MCMC iterations") + ylab(TeX("$\\theta_2$"))
pl3 = pl3 + geom_hline(yintercept = 0.3, col = "red", linewidth = 0.5) 
pl3 = pl3 + ggtitle("Model 2") + My_Theme
pl3


### Figure 2(a) and 2(b) ###

real_type = 3
sim_data = gen_data(n, real_type, sigT)
y = sim_data$y
x = sim_data$x
computer_model_functions = list()
computer_model_functions[[1]] = sim_data$f
computer_model_functions[[2]] = sim_data$g1
computer_model_functions[[3]] = sim_data$g2
nmcmc = 5000
burnin = 1000
sig = sigT
prior_sd = 10
MH_sd = 0.05
prior = "BART"
proj_type = "infinite"
nugget = 0
m = 200
res1 = posterior_projection_hDim(y, x, phi_sq, psi = 1/2, nu = 1/2, m, nmcmc, burnin, sigT, prior, proj_type, computer_model_functions, nugget,  MH_sd, prior_sd)

df = data.frame(theta1 = res1$thetaout[1,], theta2 = res1$thetaout[2,])
pl4 = ggplot(df) + geom_line(aes(x = 1:4000, y = theta1), linewidth = 0.5) + xlab("MCMC iterations") + ylab(TeX("$\\theta_1$"))
pl4 = pl4 + geom_hline(yintercept = 10, col = "red", linewidth = 0.5) 
pl4 = pl4 + ggtitle("Model 3") + My_Theme
pl4

pl5 = ggplot(df) + geom_line(aes(x = 1:4000, y = theta2), linewidth = 0.5) + xlab("MCMC iterations") + ylab(TeX("$\\theta_2$"))
pl5 = pl5 + geom_hline(yintercept = 20, col = "red", linewidth = 0.5) 
pl5 = pl5 + ggtitle("Model 3") + My_Theme
pl5