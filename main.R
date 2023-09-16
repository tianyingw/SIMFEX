#############################################################
####### Code for logit model
#############################################################
source('functions.R')
# 1.generate data------------
## parameters values
lam_real = 1 # Box(X,lam_real) is normally distributed
mean_box_x = 5 # true mean of Box(X,lam_real)
s2_box_x = 1 # true variance of Box(X,lam_real)
s2_u = 0.8 # true variance of U
beta_real = c(-0.42,log(1.5)) # c(beta_0,beta1)
## sample size
n = 1000 
r = 2
## generate dataset
set.seed(123)
X_int = Inv_Box(rnorm(n, mean = mean_box_x, sd = sqrt(s2_box_x)), lam_real) # X has no replicates
U_int = matrix(rnorm(n*r, 0, sqrt(s2_u)), n, r)
W_int = matrix(0, n, r)
for(j in 1:r){W_int[,j] = Inv_Box(Box(X_int,lam_real) + U_int[,j], lam_real)} 
w = W_int[,1] 
## generate response y
pr = H_fn(beta_real[1] + beta_real[2]*X_int)
y = c()
for(i in 1:n){y[i] = rbinom(1,size = 1,prob = pr[i])}

# 2.estimate nuisance parameters, misclassification matrix A_box and probability vector p------------
## define cut points
J = 5 # categorize W into five categories based on quintiles
C = rep(0, J-1)
for(j in 1:(J-1)){C[j] = quantile(w,j/J)} 
## estimate lambda in Box-Cox transformation (two methods are okay.)
lam_hat = optimize(fn_likeli,lower = 0,upper = 2)$minimum 
if(FALSE){
  library(MASS)
  b <- boxcox(lm(w ~ 1))
  lam_hat <- b$x[which.max(b$y)]  
}
## estimate mean_box_x, s2_u, s2_box_x
box_W_int = matrix(0,n,r)
for(j in 1:r){
  box_W_int[,j] = Box(W_int[,j],lam_hat)
}
box_W_int = na.omit(box_W_int)
row_mean_box_w = apply(box_W_int, 1, mean) 
mean_box_x_hat = mean(row_mean_box_w) 
s2_box_w = apply(box_W_int, 1, var) 
s2_u_hat = mean(s2_box_w) 
s2_box_x_hat = max(mean((row_mean_box_w - mean_box_x_hat) ^ 2) - s2_u_hat/r, 
                   0.2*(mean((row_mean_box_w - mean_box_x_hat)^2))) 
## estimate misclassification matrix A_box and probability vector p
a = rep(0,(J+1)); a[1] = min(w); a[2:J] = C; a[J+1] = max(w)
Ap = fn_box(mean_box_x_hat,s2_box_x_hat,s2_u_hat,lam_hat,a)
A_box = Ap$A_box 
p_hat = Ap$p_hat  

# 3.estimate theta_J-theta1 using our psimex method------------
## input
w = w # observations with measurement error, a n-dim vector
y = y # response, a n-dim vector
A = A_box # a misclassification matrix with J rows and J columns
p = p_hat # a J-dim probability vector
B.boot = 50 # Bootstrap number for variance estimation
## run our psimex method
fn_psimex(w,y,A,p,B.boot)
## output
## the point estimator, variance estimator and p-value for theta_J-theta_1

