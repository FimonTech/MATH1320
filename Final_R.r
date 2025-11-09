#### CASE I ####
# To compare two datasets, we don't know if the standard deviation is the same, 
# so let's use the non-pooled t-test:
unpool_test <- function(x1,x2,s1,s2,n1,n2,alpha) {
  t=(x1-x2)/(sqrt(s1^2/n1+s2^2/n2))
  df = (s1^2/n1+s2^2/n2)^2/((s1^2/n1)^2/(n1-1)+(s2^2/n2)^2/(n2-1))
  mu1 = (x1-x2)+qt(alpha/2,df)*sqrt(s1^2/n1+s2^2/n2)
  mu2 = (x1-x2)-qt(alpha/2,df)*sqrt(s1^2/n1+s2^2/n2)
  return(c(t,df,mu1,mu2))
}
survived <- c(687,703,709,715,721,723,723,726,728,728,728,729,
             730,730,733,733,735,736,739,741,741,741,741,743,
             749,751,752,752,755,756,766,767,769,770,780)
perished <- c(659,689,702,703,709,713,720,720,726,726,729,731,
              736,737,738,738,739,743,744,745,752,752,754,765)
# H1: mu1 != mu2, then 2*pt(t,df):
x1 <- mean(survived);
x2 <- mean(perished);
s1 <- sd(survived);
s2 <- sd(perished);
alpha <- 0.1 # Assumption.
result <- unpool_test(x1,x2,s1,s2,length(survived),length(perished),alpha)
names(result) <- c("t-distribution","Degrees of Freedom","x1-x2 - E","x1-x2 + E")
result
if (-abs(result[1]) < qt(alpha/2,result[2])) {
  print("The result indicates, with 90% confidence, there's no difference.")
} else {
  print("The result indicates, with 90% confidence, there is difference.")
}

alpha <- 0.05
result <- unpool_test(x1,x2,s1,s2,length(survived),length(perished),alpha)
names(result) <- c("t-distribution","Degrees of Freedom","x1-x2 - E","x1-x2 + E")
result
if (-abs(result[1]) < qt(alpha/2,result[2])) {
  print("The result indicates,with 95% confidence, there's no difference.")
} else {
  print("The result indicates, with 95% confidence, there is difference.")
}
x_diff <- x1-x2
error <- (result[4]-result[3])/2
sprintf("Based on data, with 95%% confidence, the difference is: %.2f+-%.2f",round(x_diff,2),round(error,2))
#### CASE II ####
regr <- function(x,y) {
  b1 <- (x %*% y - sum(x)*sum(y)/length(x))/(sum(x^2)-sum(x)^2/length(x))
  print(b1)
  b0 = mean(y)-b1*mean(x)
  print(b0)
  return(c(b0,b1))
}
velocity <- c(170,290,-130,-70,-185,-220,200,290,270,200,300,
              -30,650,150,500,920,450,500,500,960,500,850,800,1090)
distance <- c(0.032,0.034,0.214,0.263,0.275,0.275,0.45,0.5,0.5,
              0.63,0.8,0.9,0.9,0.9,0.9,1,1.1,1.1,1.4,1.7,2,2,2,2)
result <- regr(velocity,distance)
sprintf("The linear equation for the relationship is: %.3f+x*%.3f",result[1],result[2])
# To check for correlation, let's find the RSQ.
correl <- function(x,y) {
  r <- (x%*%y - sum(x)*sum(y)/length(x))/sqrt((sum(x^2)-sum(x)^2/length(x))*(sum(y^2)-sum(y)^2/length(y)))
  return(r)
}
result <- correl(velocity,distance)
sprintf("The coefficient of determination (r^2) is %.3f",result^2)
