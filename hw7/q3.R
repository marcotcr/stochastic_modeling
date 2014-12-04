ForwardBackwardRescaled <- function(Y, P, E, nu) {
  n = length(Y)
  num_states <- nrow(P)
  b = matrix(1, n, num_states)
  a = matrix(0, n, num_states)
  c = rep(0, n)
  for (i in 1:num_states) {
    a[1,i] = nu[i] * E[i, Y[1]]
  }
  c[1] = sum(a[1,])
  a[1,] = a[1,] / c[1]
  for (t in 1:(n-1)) {
    for (i in 1:num_states) {
      a[t+1,i] = E[i, Y[t+1]] * sum(a[t,] * P[,i])
    }
    c[t+1] = sum(a[t+1,])
    a[t+1,] = a[t+1,] / c[t+1]
  }
  for (t in (n-1):1 ) {
    for (i in 1:num_states ) {
      b[t,i] = sum(P[i,] * E[, Y[t+1]] * b[t+1,]) / c[t+1]
    }
  }
  return (list(a=a,b=b, c=c))
}
# Probability of state t (vector). fb is the return of the function above
Gamma <- function(fb, t) {
  return(fb$a[t,] * fb$b[t,])
}
# Probability of state t-1 = i, t=j. fb is the return of the function above
Csi <- function(fb, t, i, j, E, P) {
  1 / fb$c[t] * fb$a[t - 1, i] * E[j, Y[t]] * P[i,j] * fb$b[t,j]
}
BaumWelch <- function(Y, theta0) {
  n = length(Y)
  num_states <- nrow(P)
  nu = theta0$nu
  P = theta0$P
  E = theta0$E
  previous_logpy = 0
  while (TRUE) {
    z = ForwardBackwardRescaled(Y, P, E, nu)
    logpy = sum(log(z$c))
    nu = Gamma(z, 1)
    new_P = P
    # Computing P
    for (i in 1:num_states) {
      for (j in 1:num_states) {
        new_P[i,j] = sum(Csi(z, 2:n, i, j, E, P))
      }
      new_P[i,] = new_P[i,] / sum(new_P[i,])
    }
    P = new_P
    # Computing E
    norm = apply(Gamma(z, 1:n),2, sum)
    for (i in 1:16) {
      E[,i] = apply((Y[1:n] == i) * Gamma(z, 1:n), 2, sum) / norm
    }
    if (abs(logpy - previous_logpy) < .01) {
      break;
    }
    previous_logpy = logpy
    print(logpy)
  }
  return (list(nu=nu,P=P,E=E))
}

Viterbi <- function(Y, P, E, nu) {
  n = length(Y)
  num_states <- nrow(P)
  d = matrix(0, n, num_states)
  f = matrix(0, n, num_states)
  d[1,] = log(nu) + log(E[,Y[1]])
  for (t in 2:n) {
    d[t,] = apply(d[t-1,] + log(P), 2, max) + log(E[, Y[t]])
    f[t,] = apply(d[t-1,] + log(P), 2, which.max)
  }
  y = rep(0, n)
  y[n] = which.max(d[n,])
  for (t in (n-1):1) {
    y[t] = f[t+1,y[t+1]]
  }
  return (y)
}
P = matrix(c(.5,.5,.2,.8), nrow=2, byrow=TRUE)
E = matrix(c(rep(1/32, 32)), nrow=2)
#
P <- matrix(c(1:2, 2:1), nrow=2, byrow=T)
E <- matrix(c(1:16, 16:1), nrow=2, byrow=T)
#P <- matrix(runif(2 * 2), nrow=2, ncol=2)
#E <- matrix(runif(2 * 16), nrow=2, ncol=16)
for (i in 1:2) {
  P[i,] <- P[i,] / sum(P[i,])
  E[i,] <- E[i,] / sum(E[i,])
}
nu <- c(1/2, 1/2)
Y = scan('set3')
theta0 = list(nu=nu,P=P,E=E)
z = BaumWelch(Y, theta0)
viterbi <- Viterbi(Y, z$P, z$E, z$nu)
hours_in_month = 24 * 28
month = rep(0, length(Y) / hours_in_month)
j = 1;
for (i in seq(1, length(Y), hours_in_month)) {
  month[j] = mean(viterbi[i:(i+hours_in_month)])
  j = j + 1
}
fb = ForwardBackwardRescaled(Y, z$P, z$E, z$nu)


#z = ForwardBackwardRescaled(Y,P,E, nu)

