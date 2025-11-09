# Over = 1-pnorm = at least
# Mean and stdev - All possible means of a sample
# No need to assume normal distribution of salaries.
# pnorm(mu+x,mu,stdev)-pnorm(mu-x,mu,stdev)
# Formulae (for reference)
# test statitics: s_p=sqrt(((n1-1)*s1^2+(n2-1)*s2^2)/(n1+n2-2)
# t = x1-x2/(sp*sqrt(1/n1+1/n2))

pool_test <- function(x1,x2,s1,s2,n1,n2,alpha) {
  sp = sqrt(((n1-1)*s1^2+(n2-1)*s2^2)/(n1+n2-2))
  t=(x1-x2)/(sp*sqrt(1/n1+1/n2))
  mu1 = (x1-x2)+qt(alpha/2,n1+n2-2)*sp*sqrt(1/n1+1/n2)
  mu2 = (x1-x2)-qt(alpha/2,n1+n2-2)*sp*sqrt(1/n1+1/n2)
  return(c(sp,t,mu1,mu2))
}
# sp,t,mu1,mu2,alpha = pool_test()
unpool_test <- function(x1,x2,s1,s2,n1,n2,alpha) {
  t=(x1-x2)/(sqrt(s1^2/n1+s2^2/n2))
  df = (s1^2/n1+s2^2/n2)^2/((s1^2/n1)^2/(n1-1)+(s2^2/n2)^2/(n2-1))
  mu1 = (x1-x2)+qt(alpha/2,df)*sqrt(s1^2/n1+s2^2/n2)
  mu2 = (x1-x2)-qt(alpha/2,df)*sqrt(s1^2/n1+s2^2/n2)
  return(c(t,df,mu1,mu2))
}
# sp,t,mu1,mu2,alpha = pool_test()

## CHPT 12
# E = z_(alpha/2)*sqrt(p*(1-p)/n)
# n = 0.25*(z_(a/2)/E)^2
# If z-test: p+-E
# p = x/n
# z test if: x > 5 and n-x > 5
p_crit <- function(p,alpha,n){
  ci1 = p+qnorm(alpha/2)*sqrt(p*(1-p)/n)
  ci2 = p-qnorm(alpha/2)*sqrt(p*(1-p)/n)
  E = qnorm(alpha/2)*sqrt(p*(1-p)/n)
  returnqb(c(ci1,ci2,E))
}
# n = (z/E)^2*p*(1-p)
# if Sample is 10 or more, can use plus 4
# Central Limit theorem
p4_crit <- function(x,alpha,n){
  p = (x+2)/(n+4)
  ci1 = p+qnorm(alpha/2)*sqrt(p*(1-p)/(n+4))
  ci2 = p-qnorm(alpha/2)*sqrt(p*(1-p)/(n+4))
  E = qnorm(alpha/2)*sqrt(p*(1-p)/(n+4))
  return(c(ci1,ci2,E))
}
# If n*p0 > 5, n*(1-p0) > 5, can use one-proportion
# z = (p-p0)/sqrt(p0*(1-p0)/n)
p0_test <- function(x, n, p0, alpha, tail) {
  p <- x / n
  tolerance <- 1e-10  # Small tolerance for floating-point comparisons
  if ((n * p0 >= 5 - tolerance) && (n * (1 - p0) >= 5 - tolerance)) {
    ret <- "good"
    z <- (p - p0) / sqrt(p0 * (1 - p0) / n)
    P <- pnorm(z)
    if (tail != 0) {
      za <- qnorm(alpha)
    } else {
      za <- qnorm(alpha / 2)
    }
  } else {
    ret <- "bad"
    return(c(p, ret, n * p0, n * (1 - p0)))
  }
  return(c(p, z, sprintf("Area = %f", P), sprintf("z_alpha = %f", za)))
}
px2_test <- function(x1,x2,n1,n2,alpha,alpha_conf,tail) {
  p1 = x1/n1;
  p2 = x2/n2;
  pp = (x1+x2)/(n1+n2);
  tolerance <- 1e-10  # Small tolerance for floating-point comparisons
  if ((x1 >= 5 - tolerance) && (x2 >= 5 - tolerance) && (n1-x1 >= 5 - tolerance) && (n2-x2 >= 5 - tolerance)) {
    ret <- "good"
    z <- (p1-p2)/(sqrt(pp*(1-pp))*sqrt(1/n1+1/n2))
    if (tail == -1) {
      za <- qnorm(alpha)
      P <- pnorm(z)
    } else if (tail == 0) {
      za <- qnorm(alpha / 2)
      P <- 2*pnorm(z)
    } else if (tail == 1) {
      za <- qnorm(alpha)
      P <- 1-pnorm(z)
    }
    ci1 <- (p1-p2)+qnorm(alpha_conf/2)*sqrt(p1*(1-p1)/n1+p2*(1-p2)/n2)
    ci2 <- (p1-p2)-qnorm(alpha_conf/2)*sqrt(p1*(1-p1)/n1+p2*(1-p2)/n2)
  } else {
    ret <- "bad"
    return(c(p1,p2,pp,ret))
  }
  return(c(p1,p2,pp, z, sprintf("Area = %f", P), sprintf("z_alpha = %f", za),sprintf("Conf_int = %f,%f",ci1,ci2)))
}
## CHPT 13
# Chi-sq = qchisq(1-a,df)
# Conditions: 
#1) 2 of them 1 or greater
#2) <=20% freq <5
# chisq.test(x=observed, p=expected)
# Associated values: cond distributions of one variable are not identical
# dof = (r-1)*(c-1)
# E = R*C/n
# chisq.test(x=c(aA,aB,bA,bB),p=c(aA/sum,aB/sum,bA/sum,bB/sum))
# Ptest= 1-pchisq(X-squared,dof)
chi_ass <- function(x,y) {
  el = array(0,length(x))
  p = array(0,2*length(x))
  for(i in 1:length(x)) {
    el[i] = x[i]+y[i]
  }
  frq = c(el*sum(x)/(sum(x)+sum(y)),el*sum(y)/(sum(x)+sum(y)))
  p = frq/(sum(x)+sum(y))
  print(frq)
  print(p)
  chisq.test(x=c(x,y),p=p)
}
chi_ass4 <- function(x,y,z,w) {
  el = array(0,length(x))
  p = array(0,4*length(x))
  for(i in 1:length(x)) {
    el[i] = x[i]+y[i]+z[i]+w[i]
  }
  summ = 1/(sum(x)+sum(y)+sum(z)+sum(w))
  frq = c(el*sum(x)*summ,el*sum(y)*summ,el*sum(z)*summ,el*sum(w)*summ)
  p = frq/(sum(x)+sum(y)+sum(z)+sum(w))
  print(frq)
  print(p)
  chisq.test(x=c(x,y,z,w),p=p)
}
chi_ass3 <- function(x,y,z) {
  el = array(0,length(x))
  p = array(0,3*length(x))
  for(i in 1:length(x)) {
    el[i] = x[i]+y[i]+z[i]
  }
  summ = 1/(sum(x)+sum(y)+sum(z))
  frq = c(el*sum(x)*summ,el*sum(y)*summ,el*sum(z)*summ)
  p = frq/(sum(x)+sum(y)+sum(z))
  print(frq)
  print(p)
  chisq.test(x=c(x,y,z),p=p)
}
## Chpt 14
regr <- function(x,y) {
 b1 <- (x %*% y - sum(x)*sum(y)/length(x))/(sum(x^2)-sum(x)^2/length(x))
 print(b1)
 b0 = mean(y)-b1*mean(x)
 print(b0)
}
# b1 - slope; b0 - intercept
# For big data: excel.
#slope(),intercept(),RSQ(),CORREL()
## Chpt 16
# df = (num,den)
# F-dist: 2 df.
#SSTR = sum(ni*(xi_av-x_av)^2)
#SSE =sum((ni-1)*si^2)
# MSTR =
#MSE =
custom_aov <- function(...) {
  i <- 1
  data <- c()
  groups <- c()
  for(el in list(...)) {
    data <- c(data,el)
    groups <- c(groups,rep(LETTERS[i],length(el)))
    i<- i+1
  }
  print(length(data))
  model <- aov(data~groups)
  options(digits = 8)
  summ <- summary(model)
  print(summ)
}
mean_aov <- function(n,xa,s){
  xav <- (n %*% xa)/sum(n)
  SSTR <- n%*%(xa-c(xav))^2
  SSE <- (n-1)%*%s^2
  MSTR <- SSTR/(length(n)-1)
  MSE <- SSE/(sum(n)-length(n))
  Fv <- MSTR/MSE
  print(Fv)
}
