#code to plot out posterior from ABC produced in matlab
setwd('/home/harrison/Documents/mRNA')

theta = read.csv('ABC_theta_read_from_R.dat', header = FALSE)

pairs(theta)

abc_theta <- matrix(unlist(theta), ncol = 6, byrow = FALSE)

contour(abc_theta[,1:2])
#filled.contour(abc_theta)

nu1 = abc_theta[,1]
nu2 = abc_theta[,2]
lambda = abc_theta[,3]
omega1 = abc_theta[,4]
omega2 = abc_theta[,5]
phi = abc_theta[,6]

par(mfrow=c(3,2))
hist(nu1)
hist(nu2)
hist(omega1)
hist(omega2)
hist(lambda)
hist(phi)

cov(theta)
