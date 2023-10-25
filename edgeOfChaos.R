# tocID <- "edgeOfChaos.R"
#
#
# Purpose: Chaotic Sytems explorations for CSB195
#         (Computational Biology Foundations).
#
#
# Version: 1.1
# Date:    2020-08 - 2023-10
# Author:  boris.steipe@utoronto.ca
#
# Versions:
#   1.1  Generalize to other chaotic systems. Add Henon Map and Double
#          Pendulum. Start the Neural Network Exploration.
#   1.0  Logistic map code
#
# ToDo:
# Notes:
#
# Tools to plot the logistic map
#
# ==============================================================================


#TOC> ==========================================================================
#TOC>
#TOC>   Section  Title                                                   Line
#TOC> -----------------------------------------------------------------------
#TOC>   1        Introduction                                              39
#TOC>   2        The Logistic Map                                          42
#TOC>   2.1.1          Logistic Map - Supporting functions                 95
#TOC>   2.2        Visualizing the Logistic Map                           192
#TOC>   3        The Henon map                                            200
#TOC>   4        Three-species Lotka-Volterra (in a chaotic regime)       244
#TOC>   5        The Double Pendulum                                      250
#TOC>   6        A Neural Network at the Edge of Chaos                    332
#TOC>
#TOC> ==========================================================================


# =    1  Introduction  ========================================================


# =    2  The Logistic Map  ====================================================

# The logistic map is a simple recurrence equation:
#
#   x_n+1 <- r * x_n * (1 − x_n)
#
#   Where:
#     x_n:    the current state (between 0 and 1).
#     r:      an arbitrarily chosen parameter (typically between 0 and 4). It
#               turns out that r controls the behavior of the map.
#     x_n+1:  is the next state.
#
#   Let's examine the behavior of this relation:

runLogistic <- function(x1, r, N) {
  v <- numeric(N)
  v[1] <- x1
  for (i in 2:N) {
    x <- v[i-1]
    v[i] <- r * x * (1.0 - x)
  }
  return(v)
}

myCols <- colorRampPalette(c("#e1e6ff", "#D2BED8", "#9f3352"))(10)
N <- 50
x1 <- 0.25
plot(1:N, runLogistic(x1 = x1, r = 0.1, N = N),
     type = "l", xlim = c(0, N), ylim <- c(0, 1),
     main = "Behaviour of the logistic map recurrence",
     xlab = "number of iterations",
     ylab = x,
     col = myCols[1])
points((1:N) + 0.0, runLogistic(x1 = x1, r = 0.999, N = N),
       type = "l", xlim = c(0, N), ylim <- c(0, 1), col = myCols[2])
points((1:N) + 0.0, runLogistic(x1 = x1, r = 1.1, N = N),
       type = "l", xlim = c(0, N), ylim <- c(0, 1), col = myCols[3])
points((1:N) + 0.0, runLogistic(x1 = x1, r = 2.0, N = N),
       type = "l", xlim = c(0, N), ylim <- c(0, 1), col = myCols[4])
points((1:N) + 0.0, runLogistic(x1 = x1, r = 2.9, N = N),
       type = "l", xlim = c(0, N), ylim <- c(0, 1), col = myCols[5])
points((1:N) + 0.0, runLogistic(x1 = x1, r = 3.2, N = N),
       type = "l", xlim = c(0, N), ylim <- c(0, 1), col = myCols[6])
points((1:N) + 0.2, runLogistic(x1 = x1, r = 3.4, N = N),
       type = "l", xlim = c(0, N), ylim <- c(0, 1), col = myCols[7])
points((1:N) + 0.4, runLogistic(x1 = x1, r = 3.5, N = N),
       type = "l", xlim = c(0, N), ylim <- c(0, 1), col = myCols[8])
points((1:N) + 0.6, runLogistic(x1 = x1, r = 3.572, N = N),
       type = "l", xlim = c(0, N), ylim <- c(0, 1), col = myCols[9])
points((1:N) + 0.8, runLogistic(x1 = x1, r = 3.572, N = N),
       type = "l", xlim = c(0, N), ylim <- c(0, 1), col = myCols[10])


# ===   2.1.1  Logistic Map - Supporting functions

logisticIteration <- function(r, xMin, xMax, x, nInit, nRun, nPoints) {
  # Compute one evolution of logistic map equation for a given r at some
  # x, by default 0.25.
  # break() when a value is repeated and multiply resulting histogram
  # to give correct overall number of trials
  IMAXMATCH <- 200  # tradeoff between early breaks (good) and overhead of
                    # finding matches (bad)

  br <- seq(xMin, xMax, length.out = (nPoints + 1))

  for (i in 1:nInit) {  # burn in from starting value
    x <- r * x * (1.0 - x)   # Logistic Map equation
  }
  v <- numeric(nRun)
  v[1] <- x
  for (i in 2:nRun) {
    x <- v[i-1]
    v[i] <- r * x * (1.0 - x)
    # if i > IMAXMATCH, no point in matching, there are likely going to be
    # very many unique values, so its faster to just compute them out.
    # Otherwise: match(), and break() if the values cycle.
    if (i < IMAXMATCH && !is.na(match(v[i], v[1:(i-1)]))) {
      break()
    }
  }
  v <- v[1:i]                       # use only computed values
  v <- v[(v >= xMin) & (v <= xMax)] # discard values outside domain, otherwise
                                    # hist() will complain
  h <- hist(v[1:i], breaks = br, plot = FALSE, warn.unused = FALSE)$counts
  h <- h * (nRun / i)
  return(h)
}


logisticImage <- function(rMin, rMax,
                          xMin, xMax,
                          nR, nX, x0,
                          nInit, nRun) {
  # helper function
  img <- matrix(numeric(nR * nX), nrow = nX)
  vR <- seq(rMin, rMax, length.out = nR)

  for (i in 1:nR) {
    pBar(i, nR)
    img[ ,i] <- logisticIteration(r = vR[i],
                                  xMin = xMin, xMax = xMax,
                                  x = x0,
                                  nInit = nInit, nRun = nRun,
                                  nPoints = nX)
  }
  return(img)
}


getGrid <- function(i, a, n) {
  # Helper function
  x <- seq(i, a, length.out = (n + 1))
  x <- x + ((x[2] - x[1]) / 2)
  return( x[-length(x)] )
}


plotLogisticMap <- function(rMin = 2.7,
                            rMax = 4.0,
                            xMin = 0.0,
                            xMax = 1.0,
                            nR = 300,
                            nX = 200,
                            x0 = 0.25,
                            nInit = 1000,
                            nRun = 100000,
                            pal = c(hcl.colors(63), "#FFFFFF"),
                            zMin = 0.97) {
  # Plot the logistic map given the parameters
  img <- logisticImage(rMin = rMin, rMax = rMax,
                       xMin = xMin, xMax = xMax,
                       nR = nR, nX = nX, x0 = x0,
                       nInit = nInit, nRun = nRun)

  img <- t(1 - (img / max(img)))
  img <- ifelse(img < zMin, zMin, img)

  image(x = getGrid(rMin, rMax, nR),
        y = getGrid(xMin, xMax, nX),
        z = img,
        zlim = c(zMin, 1.0),
        col = pal,
        xlab = expression(italic("r")),
        ylab = expression(italic("x")),
        cex.axis = 0.8,
        main = sprintf("Logistic map between %3.2f and %3.2f", rMin, rMax),
        useRaster = TRUE)

}

# ==   2.2  Visualizing the Logistic Map  ======================================
if (FALSE) {
  plotLogisticMap(nR = 100, nX = 100, zMin = 0.97)
  plotLogisticMap(rMin = 3.562, rMax = 3.575, xMin = 0.888, xMax = 0.894,
                  nR = 1000, nX = 1000, zMin = 0.995)

}

# =    3  The Henon map  =======================================================

# The Henon map is defined by the following two coupled
# recurrence relations:
#
#   x_n+1 <- y_n +1 - a* x_n^2
#   y_n+1 = b* x_n
#

henonMap <- function(a, b, N, x0, y0) {
  x <- numeric(N)
  y <- numeric(N)
  x[1] <- x0
  y[1] <- y0
  for (i in 2:N) {
    x[i] <- y[i - 1] + 1 - a * x[i - 1]^2
    y[i] <- b * x[i - 1]
  }

  return(list(x = x, y = y))
}

# Parameters:
a <- 1.4
b <- 0.3
N <- 5000  # number of iterations
x0 <- 0    # initial condition for x
y0 <- 0    # initial condition for y

v <- henonMap(a, b, N, x0, y0)

# Plot the results
plot(v$x, v$y,
     pch = 20, cex = 0.7,
     col = colorRampPalette(c("#e1e6ff",
                              "#D2BED844",
                              "#9f335244"),
                            alpha = TRUE)(N),
     main = "Henon Map (phase space)",
     xlab = "x_n", ylab="y_n")




# =    4  Three-species Lotka-Volterra (in a chaotic regime)  ==================

# TBD



# =    5  The Double Pendulum  =================================================

# One simple system with a direct physical interpretation is the double
# pendulum. While the equations governing its motion are a bit more involved due
# to their nonlinear nature, the system itself is easy to understand: it's just
# one pendulum attached to the end of another. The double pendulum is a classic
# example of a chaotic system in physics. When you set it into motion, it's
# nearly impossible to predict its behavior after a short amount of time, even
# though it's completely deterministic.

# Physical Interpretation:
#   Imagine a pendulum (the "upper pendulum") hanging from a fixed point. Now,
#   imagine attaching a second pendulum (the "lower pendulum") to the end of the
#   first one. Even this simple mechanical system will exhibit chaotic behavior
#   for certain initial conditions.

# Biological Interpretation:
#   Chaotic dynamics can be found in various biological systems. One classic
#   example is the beating of the human heart. Under certain conditions, the
#   heart can enter chaotic rhythms (like fibrillation). While the heart's
#   electrical system is vastly more complicated than a double pendulum, the
#   idea is similar: a deterministic system that can exhibit unpredictable
#   behavior.

# (Example code originally contributed by ChatGPT-4, refactored for clarity)

if (! requireNamespace("deSolve")) {
  install.packages("deSolve")
}


library(deSolve)

doublePendulum <- function(t, state, parameters) {
  g <- 9.81

  theta1 <- state[1]
  omega1 <- state[2]
  theta2 <- state[3]
  omega2 <- state[4]

  deltaTheta <- theta2 - theta1
  denom1 <- (parameters[1]-parameters[2]*cos(delta_theta)^2)
  denom2 <- (2-parameters[2]*cos(delta_theta)^2)

  dy1 <- omega1
  dy2 <- ((g*(sin(theta2)*cos(deltaTheta)-parameters[2]*sin(theta1))-
             (sin(deltaTheta)*(parameters[2]*omega2^2+omega1^2*cos(deltaTheta))))/denom1)
  dy3 <- omega2
  dy4 <- ((g*(sin(theta1)*cos(deltaTheta)-sin(theta2))+
             (sin(deltaTheta)*(omega1^2+parameters[2]*omega2^2*cos(deltaTheta))))/denom2)

  return(list(c(dy1, dy2, dy3, dy4)))
}

parameters <- c(2, 0.5)  # Adjust these for different behaviors
state <- c(pi/2, 0, pi, 0)  # initial conditions: (theta1, omega1, theta2, omega2)

N <- 50
times <- seq(0, N, length.out = 5000)

result <- deSolve::ode(y = state, times = times,func = double_pendulum,
                       parms = parameters, method = "rk4")

# Trajectories: Theta over time
myCols = colorRampPalette(c("#91a6ff", "#D2BED8", "#9f3352"))(N)
plot(result[, "time"], result[, "1"],
     type = "l",
     xlab = "Time", ylab = "Theta1 / 2",
     main = "Phase Plot of Double Pendulum",
     col = myCols[1])
points(result[, "time"], result[, "3"], type = "l", col = myCols[N])

# Position in phase space
plot(result[, "1"], result[, "3"],
     xlab = "Theta1", ylab = "Theta2",
     main = "Phase Space Plot of Double Pendulum",
     col=myCols)
rectXY <- c(-3.5, -1.5, 2.5, 4.5)  # Zoom in here
rect(rectXY[1], rectXY[2], rectXY[3], rectXY[4], border = "#ffdd77")

# Zoomed in:
plot(result[, "1"], result[, "3"],
     xlim = rectXY[c(1,3)],
     ylim = rectXY[c(2,4)],
     xlab = "Theta1", ylab = "Theta2",
     main = "Phase Plot of Double Pendulum",
     col=myCols)





# =    6  A Neural Network at the Edge of Chaos  ===============================

# Let's design a simple neural network "at the edge of order and chaos".
# Assume the input is an oscillator (and maybe some directed perturbation).
#
# Let's build a network that has a complex output that can have chaotic
# behaviour but a well defined attractor domain. Let's take frequency and
# amplitude as the relevant phase space.
#
# How ? ...
#







# [END]
