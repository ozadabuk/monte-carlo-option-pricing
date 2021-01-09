library(pracma)
tic()
# X ~ N(0, 1)
# r(t) <- r(t-1) + 0.18 * (0.086 - r(t-1)) * delta + 0.02 * sqrt(delta) * X
r_path<-function(n){
  X <- matrix(rnorm(n*d, 0, 1),nrow=n) 
  r_t <- rep(0.07, n)
  r <- matrix(r_t)
  for(i in 2:d){
    r_t <- r_t + 0.18 * (0.086 - r_t) * delta + 0.02 * sqrt(delta) * X[,i]
    r <- cbind(r, r_t)
  }
  return(r)
}

S_path<-function(n, r){
  X<-matrix(rnorm(n*d, 0, 1),nrow=n)
  RM = t(apply(r, 1, cumsum)) # apply cumsum function to matrix R
  BM<-sqrt(delta)*t(apply(X,1,cumsum)) # apply cumsum function to matrix X and scale by sqrt(delta)
  S<-S_0*exp(sweep(sigma*BM + RM*delta,MARGIN=2,(k*sigma^2/2)*delta,'-'))

  return(S)
}

# Oguz price path
# S_path<-function(n, r){
#   Z<-matrix(rnorm(n*(d), 0, 1),nrow=n)
#   S_t <- S_0 * exp((r[,1] - (0.5 * sigma^2)) * delta + sigma * Z[,1]* sqrt(delta))
#   S<-matrix(S_t)
# 
#   for(i in 2:(d)){
#     S_t <- S_t * exp((r[,i] - (0.5 * sigma^2)) * delta + sigma * Z[,i]* sqrt(delta))
#     S<-cbind(S, S_t)
#   }
#   return(S)
# }

### 12 time steps

S_0 <- 50
K <- 50
r_0 <- 0.07
sigma <- 0.13
epsilon <- 0.05
MT <- 1 # 1 year
d <- 12 #monthly steps
delta <- MT / d
k <- seq(1,d)

# create first sample to find actual sample size
n.sigma <- 10000
r <- r_path(n.sigma) # get an interest rate path
S <- S_path(n.sigma, r) # get a Vasicek stock price path

euro_call_payoffs <- pmax(S[,d] - K, 0) * exp(-delta * rowSums(r)) # Vasicek call payoffs

N_Call <- ceiling((2.58 * 1.1 * sd(euro_call_payoffs) / epsilon)^2) # find sample size required for error tolerance (Vasicek)

r <- r_path(N_Call) # generate new interest rate path with new sample size
S <- S_path(N_Call, r) # generate new Vasicek price path with new sample size

euro_call_payoffs_monthly <- pmax(S[,d] - K, 0) * exp(-delta * rowSums(r)) # Vasicek call payoffs
euro_call_price_monthly <- mean(euro_call_payoffs_monthly) # Vasicek call price
error_call_price_monthly <- 2.58 * sd(euro_call_payoffs_monthly) / sqrt(N_Call) # error for Vasicek call price


### 52 time steps

d <- 52 # weekly time steps
delta <- MT / d
k<-seq(1,d)

n.sigma <- 1000
r <- r_path(n.sigma) # get an interest rate path
S <- S_path(n.sigma, r) # get a Vasicek stock price path

euro_call_payoffs <- pmax(S[,d] - K, 0) * exp(-delta * rowSums(r)) # Vasicek call payoffs

N_Call <- ceiling((2.58 * 1.1 * sd(euro_call_payoffs) / epsilon)^2) # find sample size required for error tolerance (Vasicek)

r <- r_path(N_Call) # generate new interest rate path with new sample size
S <- S_path(N_Call, r) # generate new Vasicek price path with new sample size

euro_call_payoffs_weekly <- pmax(S[,d] - K, 0) * exp(-delta * rowSums(r)) # Vasicek call payoffs
euro_call_price_weekly <- mean(euro_call_payoffs_weekly) # Vasicek call price
error_call_price_weekly <- 2.58 * sd(euro_call_payoffs_weekly) / sqrt(N_Call) # error for Vasicek call price

# true European call price
ExactEuroCall <- S_0*pnorm((log(S_0/K)+(r_0+sigma^2/2)*MT)/(sigma*sqrt(MT))) - K*exp(-r_0*MT)*pnorm((log(S_0/K)+(r_0-sigma^2/2)*MT)/(sigma*sqrt(MT)))

cat("Vasicek call price for monthly steps - price:", euro_call_price_monthly, " error:", error_call_price_monthly,"\n")
cat("Vasicek call price for weekly steps - price:", euro_call_price_weekly, " error:", error_call_price_weekly, "\n")
cat("True European call price from Black Scholes Model:", ExactEuroCall, "\n")
toc()
