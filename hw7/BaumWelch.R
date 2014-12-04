ForwardBackwardRescaled <- function(x, E, P, nu) {
  numStates <- nrow(E)
  n <- length(x)
  backward <- matrix(0, nrow=n, ncol=numStates)
  backwardNorm <- rep(0, n)
  backward[n,] <- 1/numStates
  backwardNorm[n] <- numStates
  for (t in (n-1):1) {
    for (i in 1:numStates) {
      backward[t,i] <- sum(P[i,] * E[,x[t+1]] * backward[t+1,])
    }
    backwardNorm[t] <- sum(backward[t,])
    backward[t,] <- backward[t,] / backwardNorm[t]
  }
  forward <- matrix(0, nrow=n, ncol=numStates)
  forwardNorm <- rep(0, n)
  forward[1,] <- E[,x[1]] * nu
  forwardNorm[1] <- sum(forward[1,])
  forward[1,] <- forward[1,] / forwardNorm[1]
  for (t in 1:(n-1)) {
    for (i in 1:numStates) {
      forward[t+1,i] <- E[i,x[t+1]] * sum(forward[t,] * P[,i])
    }
    forwardNorm[t+1] <- sum(forward[t+1,])
    forward[t+1,] <- forward[t+1,] / forwardNorm[t+1]
  }
  list(backward=backward, backwardNorm=backwardNorm,
       forward=forward, forwardNorm=forwardNorm)
}
# probability that t-th state takes on value i
GetProbabilityRescaled <- function(forwardBackwardRescaled, t, i) {
  numerator <- forwardBackwardRescaled$forward[t,i] * forwardBackwardRescaled$backward[t,i]
  return(numerator/sum(forwardBackwardRescaled$forward[t,] * forwardBackwardRescaled$backward[t,]))
}

BaumWelch <- function(x, E, P, nu, updateNu=F) {
  n <- length(x)
  ll = NULL
  previous_ll = 0
  iter = 1
  while (TRUE) {
    P.hat <- matrix(0, nrow=2, ncol=2)
    E.hat <- matrix(0, nrow=2, ncol=16)
    fb <- ForwardBackwardRescaled(x, E, P, nu)
    log.lik <- sum(log(fb$forwardNorm))
    ll = cbind(ll, log.lik)
    if (abs(log.lik - previous_ll) < 0.01) {
      break;
    }
    previous_ll = log.lik
    cat(paste("Iteration ", iter, ": ", log.lik, "\n", sep=""))
    # don't update nu
    if (updateNu) {
      nu <- c(GetProbabilityRescaled(fb, 1, 1), GetProbabilityRescaled(fb, 1, 2))
      nu <- nu / sum(nu)
    }
    iter = iter + 1

    for (i in 1:2) {
      E.hat[i,x[1]] <- E.hat[i,x[1]] + GetProbabilityRescaled(fb, 1, i)
    }
    for (t in 2:n) {
      g <- matrix(0, nrow=2, ncol=2)
      for (i in 1:2) {
        E.hat[i,x[t]] <- E.hat[i,x[t]] + GetProbabilityRescaled(fb, t, i)
        g[i,] <- fb$backward[t,] * E[,x[t]] * P[i,] * fb$forward[t-1,i]
        #for (j in 1:2) {
        #  g[i,j] <- fb$backward[t,j] * E[j,x[t]] * P[i,j] * fb$forward[t-1,i]
        #}
      }
      g <- g / sum(g)
      P.hat <- P.hat + g
    }
    for (i in 1:2) {
      P[i,] <- P.hat[i,] / sum(P.hat[i,])
      E[i,] <- E.hat[i,] / sum(E.hat[i,])
    }
  }
  list(E=E, P=P, nu=nu, ll=ll)
}

Viterbi <- function(x, E, P, nu) {
  numStates <- nrow(E)
  n <- length(x)
  d <- matrix(0, nrow=n, ncol=numStates)
  f <- matrix(0, nrow=n, ncol=numStates)
  for (i in 1:numStates) {
    d[1,] <- log(nu) + log(E[,x[1]])
  }
  for (t in 2:n) {
    for (i in 1:numStates) {
      inner <- d[t-1,] + log(P[,i])
      maximizer <- which.max(inner)
      f[t,i] <- maximizer
      d[t,i] <- inner[maximizer] + log(E[i,x[t]])
    }
  }
  viterbi <- rep(0,n)
  viterbi[n] <- which.max(d[n,])
  for (t in (n-1):1) {
    viterbi[t] <- f[t+1,viterbi[t+1]]
  }
  viterbi
}
GetMovingAverage <- function(sequence, window, weighted=F) {
  n <- length(sequence) - window

  if (weighted) {
    sd <- (window/2)/2
    w <- rep(0, window)
    if ((window %% 2) == 0) {
      half <- window/2
      w[1:half] <- dnorm((half-1):0, 0, sd)
      w[(half+1):window] <- w[half:1]
    } else {
      half <- floor(window/2)
      w[1:(half+1)] <- dnorm(half:0, 0, sd)
      w[(half+2):window] <- w[half:1]
    }
  }

  avg <- rep(0, n)
  for (i in 1:n) {
    if (weighted) {
      avg[i] <- weighted.mean(sequence[i:(i+window-1)], w)
    } else {
      avg[i] <- mean(sequence[i:(i+window-1)])
    }
  }
  avg
}
GetMovingAverageDates <- function(dates, window) {
  dates[floor(window/2):(length(dates) - floor(window/2) - 1)]
}

set.seed(200)
wind <- scan("wind.txt", quiet=T)

# MLE estimates of transition and emission probabilities
best_ll = -10000000000000000000000000000000
best_mles = 0
for (i in 1:5) {
  P <- matrix(runif(2 * 2), nrow=2, ncol=2)
  E <- matrix(runif(2 * 16), nrow=2, ncol=16)
  for (i in 1:2) {
    P[i,] <- P[i,] / sum(P[i,])
    E[i,] <- E[i,] / sum(E[i,])
  }
nu <- c(1/2, 1/2)
  mles <- BaumWelch(wind, E, P, nu, updateNu=T)
  if (mles$ll[length(mles$ll)] > best_ll) {
    best_ll = mles$ll[length(mles$ll)]
    best_mles = mles
  }
}

# Plot wind directions
library(ggplot2)
PlotDirections <- function(emissions, title) {
  directions <- c("N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE", "S",
                  "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW")
  dir.factor <- factor(directions, levels=directions)
  probs.frame <- data.frame(direction=dir.factor, probability=emissions)
  ggplot(probs.frame, aes(x=direction, y=probability)) +
    geom_bar(stat="identity", width=1) + coord_polar(start=-pi/16) +
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank()) + ggtitle(title)
}
plot1 <- PlotDirections(best_mles$E[1,], "State 1 Wind Directions")
plot2 <- PlotDirections(best_mles$E[2,], "State 2 Wind Directions")
library(gridExtra)
grid.arrange(plot1, plot2, ncol=2)

# latent sequence
viterbi <- Viterbi(wind, best_mles$E, best_mles$P, best_mles$nu)
week <- 24 * 7
month <- week * 4
d <- as.Date("1985-5-1")
dates <- d + 0:(length(viterbi)/24 - 1)
dates <- rep(dates, rep(24, length(dates)))
move.avg <- GetMovingAverage(viterbi, month, weighted=T)
move.avg.dates <- GetMovingAverageDates(dates, month)
plot(move.avg.dates, move.avg, type='l', main="Monthly Moving Average Latent State", xlab="", ylab="Monthly Moving Average", xaxt="n")
labDates <- seq(head(move.avg.dates, 1), tail(move.avg.dates, 1), by = "2 months")
axis.Date(side=1, move.avg.dates, at=labDates, format="%b %y", las=2)

plot(1:length(best_mles$ll), best_mles$ll, xlab="Iteration", ylab="Log likelihood",main="Log likelihood per iteration")

# a fun "video"
# for (x in c(10, 50, 100, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000)) {
#   move.avg <- GetMovingAverage(viterbi, x, weighted=T)
#   move.avg.dates <- GetMovingAverageDates(dates, x)
#   plot(move.avg.dates, move.avg, type='l',
#        main=paste("W = ", x, sep=""), xlab="", ylab="Moving Average", xaxt='n')
#   labDates <- seq(head(move.avg.dates, 1), tail(move.avg.dates, 1), by = "2 months")
#   axis.Date(side=1, move.avg.dates, at=labDates, format="%b %y", las=2)
# }
