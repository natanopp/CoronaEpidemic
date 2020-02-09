## CoronaVirus Analysis 
## Developer: Learning Machines
## Replicated by N. Pimonsathian
## Version: 0.1 
## Date: 09-02-2020

## Load and install packages ---------------------------------------------------
install.packages("yaml")
install.packages("deSolve")
library(deSolve)

## Data Input ------------------------------------------------------------------
Infected <- c(45, 62, 121, 198, 291, 
              440, 571, 830, 1287, 1975, 
              2744, 4515, 5974, 7711, 9692, 
              11791, 14380, 17205, 20440)
Day <- 1:(length(Infected))
N <- 1400000000 # population of mainland china

## Preliminary EDA -------------------------------------------------------------
old <- par(mfrow = c(1, 2))
plot(Day, Infected, type ="b")
plot(Day, Infected, log = "y")
abline(lm(log10(Infected) ~ Day))
title("Confirmed Cases 2019-nCoV China", outer = TRUE, line = -2)

## Define the SIR Function -----------------------------------------------------
SIR <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- -beta/N * I * S
    dI <- beta/N * I * S - gamma * I
    dR <- gamma * I
    list(c(dS, dI, dR))
  })
}

## Implement the SIR Function --------------------------------------------------
init <- c(S = N-Infected[1], I = Infected[1], R = 0)
RSS <- function(parameters) {
  names(parameters) <- c("beta", "gamma")
  out <- ode(y = init, times = Day, func = SIR, parms = parameters)
  fit <- out[ , 3]
  sum((Infected - fit)^2)
}

Opt <- optim(c(0.5, 0.5), RSS, method = "L-BFGS-B", 
             lower = c(0, 0), upper = c(1, 1)) # optimize with some sensible conditions
Opt$message

Opt_par <- setNames(Opt$par, c("beta", "gamma"))
Opt_par

t <- 1:70 # time in days
fit <- data.frame(ode(y = init, times = t, func = SIR, parms = Opt_par))
col <- 1:3 # colour

# Visualize the estimations 
matplot(fit$time, fit[ , 2:4], type = "l", 
        xlab = "Day", ylab = "Number of subjects", 
        lwd = 2, lty = 1, col = col)
matplot(fit$time, fit[ , 2:4], type = "l", 
        xlab = "Day", ylab = "Number of subjects", 
        lwd = 2, lty = 1, col = col, log = "y")

points(Day, Infected)
legend("bottomright", c("Susceptibles", "Infecteds", "Recovereds"), 
       lty = 1, lwd = 2, col = col, inset = 0.05)
title("SIR model 2019-nCoV China", outer = TRUE, line = -2)

## Calculate the basic reproduction number (R0) --------------------------------
par(old)
R0 <- setNames(Opt_par["beta"] / Opt_par["gamma"], "R0")

fit[fit$I == max(fit$I), "I", drop = FALSE] # height of pandemic
max(fit$I) * 0.02 # max deaths with supposed 2% mortality rate









