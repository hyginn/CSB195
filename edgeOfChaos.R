# tocID <- "edgeOfChaos.R"
#
#
# Purpose: Chaotic Sytems explorations for CSB195
#         (Computational Biology Foundations).
#
#
# Version: 1.2
# Date:    2020-08 - 2023-10
# Author:  boris.steipe@utoronto.ca
#
# Versions:
#   1.2  Annotate the logistic map operations more. Allow different lengths
#          and masses for the coupled pendulums.
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
#TOC>   Section  Title                                                       Line
#TOC> ---------------------------------------------------------------------------
#TOC>   1        Introduction                                                  43
#TOC>   2        The Logistic Map                                              46
#TOC>   2.1.1          Logistic Map - Supporting functions                    182
#TOC>   2.2        Visualizing the Logistic Map                               279
#TOC>   3        The Henon map                                                308
#TOC>   4        Three-species Lotka-Volterra (in a chaotic regime)           355
#TOC>   5        The Double Pendulum                                          361
#TOC>   6        A Neural Network at the Edge of Chaos                        486
#TOC> 
#TOC> ==========================================================================


# =    1  Introduction  ========================================================


# =    2  The Logistic Map  ====================================================

# The logistic map is a simple recurrence relation:
#
#   x_n+1 <- r * x_n * (1 − x_n)
#
#   Where:
#     x_n:    the current state (between 0 and 1).
#     r:      an arbitrarily chosen parameter (typically between 0 and 4). It
#               turns out that r controls the behavior of the map.
#     x_n+1:  is the next state.
#
# (N.b. There are a number of related terms: "relation", "map", "equation",
#  "transformation", "operator", ... all of which describe ways to associate
#  one "thing" with another. The term "relation" in "recurrence relation" can
#  be a bit misleading: here it  refers to an _equation_ that defines a
#  sequence in a recursive manner.  The term "map" refers to the function's
#   ability to "map" or "send" one value (from the "domain") to another
#   value (in the "codomain"). This is common terminology in the study of
#   dynamical systems. For the logistic map, the iteration "maps" each value
#   in the interval [0,1] to some other value in [0,1].
#
#   Let's examine the behavior of this relation:

runLogistic <- function(x1, r, N) {
  # This function iterates the Logistic Map from the starting value
  # x1 for N iterations using the parameter r. It returns every value
  # x that is produced, in a vector.

  x <- numeric(N)                       # Initialize the vector of results
  x[1] <- x1                            # x1 is the first value
  for (n in 1:(N - 1)) {                # for N steps ...
    x[n+1] <- r * x[n] * (1.0 - x[n])   #  ... apply the relation.
  }
  return(x)
}

# Different starting values x1 converge to the same value quite quickly ...
N <- 50
plot(1:N, runLogistic(x1 = 0.1, r = 0.99, N = N),
     type = "l", xlim = c(0, N), ylim <- c(0, 0.6),
     main = "Logistic map: convergence",
     xlab = "number of iterations",
     ylab = "x",
     col = "#9f3352")
for (x1 in seq(0.1, 0.6, by = 0.1)) {
  points((1:N), runLogistic(x1, 0.99, N), type = "l", col = "#9f3352")
}
# ... but that's not what makes the Logistic Map interesting. The relation
# behaves quite remarkably for different choices of r !

# Define a set of colors from a gradient
myCols <- colorRampPalette(c("#e1e6ff", "#D2BED8", "#9f3352"))(10)

N <- 50
x1 <- 0.25   # This is not a particularly special value

plot(1:N, runLogistic(x1 = x1, r = 0.1, N = N),
     pch = 20, xlim = c(0, N), ylim = c(0, 1),
     main = "Logistic map: Dependence on r",
     xlab = "number of iterations",
     ylab = "x",
     col = myCols[1])
#                         Watch what happens for different values of r:
points((1:N), runLogistic(x1 = x1, r = 0.999, N = N),  # r just below 1
       pch = 20, col = myCols[2])
points((1:N), runLogistic(x1 = x1, r = 1.1, N = N),    # r = 1.1
       pch = 20, col = myCols[3])
points((1:N), runLogistic(x1 = x1, r = 2.0, N = N),    # 2.0
       pch = 20, col = myCols[4])
points((1:N), runLogistic(x1 = x1, r = 2.9, N = N),    # 2.9
       pch = 20, col = myCols[5])
points((1:N), runLogistic(x1 = x1, r = 3.2, N = N),
       pch = 20, col = myCols[6])
points((1:N), runLogistic(x1 = x1, r = 3.4, N = N),
       pch = 20, col = myCols[7])
points((1:N), runLogistic(x1 = x1, r = 3.5, N = N),
       pch = 20, col = myCols[8])
points((1:N), runLogistic(x1 = x1, r = 3.572, N = N),
       pch = 20, col = myCols[9])
points((1:N), runLogistic(x1 = x1, r = 3.67, N = N),
       pch = 20, col = myCols[10])

# Let's look again at the behaviour of the last six values of r - I'll draw
# them apart vertically.

plot(1:N, runLogistic(x1 = x1, r = 2.9, N = N),
     pch = 20, xlim = c(0, N), ylim = c(0, 6),
     main = "Logistic map: Dependence on r",
     xlab = "number of iterations",
     ylab = "x",
     yaxt = "n",
     col = myCols[5])
abline(h = 1:5)
points((1:N), runLogistic(x1 = x1, r = 3.2, N = N) + 1,
       pch = 20, col = myCols[6])
points((1:N), runLogistic(x1 = x1, r = 3.4, N = N) + 2,
       pch = 20, col = myCols[7])
points((1:N), runLogistic(x1 = x1, r = 3.5, N = N) + 3,
       pch = 20, col = myCols[8])
points((1:N), runLogistic(x1 = x1, r = 3.572, N = N) + 4,
       pch = 20, col = myCols[9])
points((1:N), runLogistic(x1 = x1, r = 3.67, N = N) + 5,
       pch = 20, col = myCols[10])
text(rep(5, 6), (0:4) + 0.1,
     labels = sprintf("r = %.3f",c(2.9, 3.2, 3.4, 3.5, 3.572, 3.67)))

# At r = 2.9, instead of quickly moving to a converging value, we see some
# oscillation in x values - but they do converge. However, at r = 3.2, we get
# divergence. At r = 3.4, there is a noticeable oscillation superimposed
# on the diverging values ... which gets stronger at r = 3.5, and
# at r = 3.572, we get oscillations on top of the oscillations.

# Let's examine the distribution of the x values at r = 3.67, by creating
# a million of them:

x <- runLogistic(x1 = x1, r = 3.67, N = 1000000) # a million
h <- hist(x, breaks = 150, plot = FALSE)         # get the histogram values

# define 50 colors from the hcl palette from light to dark
myPal <- rev(c(hcl.colors(49), "#ffffff"))
# assign a color from the palette to each histogram bar
myCols <- myPal[round(1 + ((h$counts / max(h$counts)) * 49))]
# plot the histogram
hist(x, breaks = 150, col = myCols, border = myCols)

# To explore this behavior of r more clearly, we can produce very fine-
# grained histograms of x values, stack them vertically and plot them
# side by side. A single stack looks like this:

plot(rep(1, length(h$mids)), h$mids, pch = 15, col = myCols)

# This requires a few supporting functions, you can ignore them and skip to the
# next section where we plot the map.


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
  # Helper function
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


plotLogisticMap <- function(rMin = 2.5,
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

  image(x = getGrid(rMin, rMax, nR),   # This determines the number of pixels
        y = getGrid(xMin, xMax, nX),   # that are being computed
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

  # Low resolution plot of a large domain of r.
  plotLogisticMap(nR = 200, nX = 200, zMin = 0.97)

  # High resolution plot of the region in the tiny interval:
  rI <- 3.562
  rA <- 3.575
  xI <- 0.888
  xA <- 0.894
  abline(h = c(xI, xA), col = "#AA0000")
  abline(v = c(rI, rA), col = "#AA0000")

  #  Compute four million pixels (2000 * 2000).
  plotLogisticMap(rMin = rI, rMax = rA, xMin = xI, xMax = xA,
                  nR = 2000, nX = 2000, zMin = 0.995)

  #  To actually see the plot at full resolution, and appreciate its
  #  intricate structure, click on Zoom when it is done and look at it on a
  #  large monitor, or export it and zoom in.

  # If we would make the interval even smaller, we would get similar structures.
  # The x values can lie arbitrarily close to each other. The relation in this
  # domain is exquisitely sensitive to the values of r. "Sensitive to initial
  # conditions" - that is the definition of chaotic behaviour.

}

# =    3  The Henon map  =======================================================

# The Henon map is another "classical" chaotic system defined by the
# following _two_ coupled recurrence relations:
#
#   x_n+1 <- y_n +1 - a* x_n^2
#   y_n+1 = b* x_n

# We can look at its behaviour, to introduce the idea of a "phase space".
# A phase space plot plots one

henonMap <- function(a, b, N, x0, y0) {
  x <- numeric(N)
  y <- numeric(N)
  x[1] <- x0
  y[1] <- y0
  for (n in 1:N-1) {
    x[n+1] <- y[n] + 1 - a * x[n]^2
    y[n+1] <- b * x[n]
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
     col = colorRampPalette(c("#81a6ff44",
                              "#a2BED844",
                              "#cf335244"),
                            alpha = TRUE)(N),
     main = "Henon Map (phase space)",
     xlab = "x_n", ylab="y_n")

# Note that these points do not fall on this curve consecutively, but at
# widely different points in our sequence of 5,000 iterations.


# =    4  Three-species Lotka-Volterra (in a chaotic regime)  ==================

# TBD



# =    5  The Double Pendulum  =================================================

# One simple system with a direct physical interpretation is the double
# pendulum. The equations governing its motion are a bit more involved, but the
# system itself is easy to understand: it's just one pendulum attached to the
# end of another. The double pendulum is a classic example of a chaotic system
# in physics. When you set it into motion, it is nearly impossible to predict
# its behavior after a short amount of time, even though it is completely
# deterministic.

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
#   behavior. The double pendulum is just a metaphor in this case, but the
#   dynamics we get from the interacting oscillators in the heart can look
#   very similar to such a simple system. The we can ask: what aspects
#   of the actual heartbeat can be mapped to the simple physical (or
#   computational) model.

# (Example code originally contributed by ChatGPT-4, refactored for style
# and clarity, and commented)

if (! requireNamespace("deSolve")) {  # differential equation solver
  install.packages("deSolve")
}

doublePendulum <- function(t,        # timepoint
                           state,    # four numbers: initial conditions
                           parm) {   # four numbers: lengths and masses

  # This is the function to compute the state of a pendulum at timepoint t
  # as a coupled differential equation. deSolve::ode() solves this equation
  # at the requested timepoints to compuite the trajectori for the two
  # pendulums over time.

  g <- 9.81  # (earth standard surface gravitational constant: 9.1 m/s^2)

  # Initial conditions
  th1 <- state[1]  # angle with the vertical of first pendulum
  om1 <- state[2]  # angular velocity of first pendulum (rad/sec)
  th2 <- state[3]
  om2 <- state[4]

  # Parameters
  l1 <- parm[1]    # length of first pendulum
  m1 <- parm[2]    # mass of first pendulum
  l2 <- parm[3]
  m2 <- parm[4]

  deltaTheta <- th2 - th1
  denom1 <- ((m1 + m2) * l1) - (m2 * l1 * cos(deltaTheta)^2)
  denom2 <- (l2 / l1) * denom1

  dy1 <- om1
  dy2 <- ((m2 * l2 * w2^2 * sin(deltaTheta) * cos(deltaTheta)
           + m2 * g * sin(y3) * cos(deltaTheta)
           + m2 * l2 * w2^2 * sin(deltaTheta)
           - (m1 + m2) * g * sin(y1))
          / denom1)

  dy3 <- om2
  dy4 <- ((-l1 / l2 * w1^2 * sin(deltaTheta) * cos(deltaTheta)
           + (m1 + m2) * g * sin(y1) * cos(deltaTheta)
           - (m1 + m2) * l1 * w1^2 * sin(deltaTheta)
           - (m1 + m2) * g * sin(y3))
          / denom2)

  return(list(c(dy1, dy2, dy3, dy4)))
}

# Initial state: (th1, om1, th2, om2)
state <- c(pi/2,   # angle from the vertical: 90°
           0,      # angular velocity zero: pendulum initially at rest
           pi,     # 180°
           0)      # ... at rest

# Parameters
parm <- c(1.8,     # first pendulum: 1.8 meters long ...
          0.5,     #   ... with a mass of 0.5 kg
          0.8,     # 0.8m
          0.4)     # 0.4 kg

N <- 100  # (seconds)
times <- seq(0, N, length.out = 5000) # 5,000 timepoints from 0 to 100 seconds

# Solve this
result <- deSolve::ode(y = state, times = times,func = double_pendulum,
                       parms = parm, method = "rk4")

# Trajectories: Theta over time
myCols = colorRampPalette(c("#91a6ff", "#D2BED8", "#9f3352"))(N)
plot(result[, "time"], result[, "1"],
     type = "l",
     ylim = c(-90, 90),
     xlab = "Time (sec)", ylab = "Theta1 , 2 (rad)",
     main = "Phase Plot of Double Pendulum",
     col = myCols[1])
points(result[, "time"], result[, "3"], type = "l", col = myCols[N])

# Position in phase space
plot(result[, "1"], result[, "3"],
     xlab = "Theta1", ylab = "Theta2",
     main = "Phase Space Plot of Double Pendulum",
     col=myCols)
rectXY <- c(82, 37, 95, 48)  # Zoom in here
rect(rectXY[1], rectXY[2], rectXY[3], rectXY[4], lwd = 1.5, border = "#00ddff")

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
