# tocID <- "edgeOfChaos.R"
#
#
# Purpose: Chaotic Sytems explorations for CSB195
#         (Computational Biology Foundations).
#
#
# Version: 1.3
# Date:    2020-08 - 2023-10
# Author:  boris.steipe@utoronto.ca
#
# Versions:
#   1.3  Fix faulty double pendulum code. Add better plotting: three panels
#          in one window, power spectrum analysis, plot legends, explanations.
#   1.2  Annotate the logistic map operations more. Allow different lengths
#          and masses for the coupled pendulums. Comment the code.
#   1.1  Generalize to other chaotic systems. Add Henon Map and Double
#          Pendulum. Start the Neural Network Exploration.
#   1.0  Logistic map code
#
# ToDo:
# Notes:
#
# ==============================================================================


#TOC> ==========================================================================
#TOC>
#TOC>   Section  Title                                                       Line
#TOC> ---------------------------------------------------------------------------
#TOC>   1        Required packages                                             44
#TOC>   2        Introduction                                                  79
#TOC>   3        The Logistic Map                                              82
#TOC>   3.1.1          Logistic Map - Supporting functions                    218
#TOC>   3.2        Visualizing the Logistic Map                               316
#TOC>   4        The Henon map                                                346
#TOC>   5        Three-species Lotka-Volterra (in a chaotic regime)           395
#TOC>   6        The Double Pendulum                                          401
#TOC>   7        A Neural Network at the Edge of Chaos                        658
#TOC>
#TOC> ==========================================================================


# =    1  Required packages  ===================================================

# Check that required packages have been installed. Install if needed. Get
# information about the package contents.

# A package to compute ODEs ...
if (! requireNamespace("deSolve", quietly=TRUE)) {
  install.packages("deSolve")
}
# Package information:
#  library(help   = deSolve)     # basic information
#  browseVignettes("deSolve")    # available vignettes
#  data(package  = "deSolve")    # available datasets


# A package to support plotting. We need this  to plot lines with a colour
# gradient which we can't do in base R.
if (! requireNamespace("plotrix", quietly=TRUE)) {
  install.packages("plotrix")
}
# Package information:
#  library(help   = plotrix)     # basic information
#  data(package  = "plotrix")    # available datasets


# A package to support plotting 3D plots ...
if (! requireNamespace("plot3D", quietly = TRUE)) {
  install.packages("plot3D")
}
#  library(help   = plot3D)     # basic information
#  data(package  = "plot3D")    # available datasets




# =    2  Introduction  ========================================================


# =    3  The Logistic Map  ====================================================

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


# ===   3.1.1  Logistic Map - Supporting functions

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


plotLogisticMap <- function(rMin  = 2.5,
                            rMax  = 4.0,
                            xMin  = 0.0,
                            xMax  = 1.0,
                            nR    = 300,
                            nX    = 200,
                            x0    = 0.25,
                            nInit = 1000,
                            nRun  = 100000,
                            pal   = c(hcl.colors(63), "#FFFFFF"),
                            zMin  = 0.97) {
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


# ==   3.2  Visualizing the Logistic Map  ======================================

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

# =    4  The Henon map  =======================================================

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




# =    5  Three-species Lotka-Volterra (chaotic regime)  ==================

# TBD



# =    6  The Double Pendulum  =================================================

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

# Note: This example code was originally contributed by ChatGPT-4. The general
# layout of the input parameters and how to set up the function for the coupled
# differential equations and solve them with desolve::ode() was entirely correct
# - and that is very helpful. However, the details of the (complicated)
# first-order differential equations were wrong, even though they used the right
# elements of lengths, masses, and angular velocities, and desolve::ode()
# converged on solutions. The correct equations were taken from
# https://www.myphysicslab.com/pendulum/double-pendulum-en.html (actual code
# translated from the Javascript source at
# https://github.com/myphysicslab/myphysicslab/
# blob/master/src/sims/pendulum/DoublePendulumSim).

if (! requireNamespace("deSolve")) {  # differential equation solver
  install.packages("deSolve")
}

doublePendulum <- function(t,        # timepoint
                           state,    # current values of theta and omega
                           parm) {   # lenghts and masses
  # This function to computes the state of a double pendulum at timepoint t
  # given the two angles from the vertical and angular velocities, and the
  # parameters for length and mass of both pendulums. deSolve::ode() solves this
  # equation at the requested timepoints to compute the trajectories for the two
  # pendulums over time.

  g <- 9.81  # (earth standard surface gravitational constant: 9.81 m/s^2)

  # Initial conditions
  th1 <- state["th1"]  # angle with the vertical of first pendulum
  om1 <- state["om1"]  # angular velocity of first pendulum (rad/sec)
  th2 <- state["th2"]
  om2 <- state["om2"]

  # Parameters
  L1 <-   parm[["L1"]]    # length of first pendulum
  m1 <-   parm[["m1"]]    # mass of first pendulum
  L2 <-   parm[["L2"]]
  m2 <-   parm[["m2"]]

  dy1 <- om1

  dy2 <- -g*(2*m1+m2)*sin(th1)
  dy2 <- dy2 - g*m2*sin(th1-2*th2)
  dy2 <- dy2 - 2*m2*om2*om2*L2*sin(th1-th2)
  dy2 <- dy2 - m2*om1*om1*L1*sin(2*(th1-th2))
  dy2 <- dy2/(L1*(2*m1+m2-m2*cos(2*(th1-th2))))

  dy3 <- om2

  dy4 <- (m1+m2)*om1*om1*L1
  dy4 <- dy4 + g*(m1+m2)*cos(th1)
  dy4 <- dy4 + m2*om2*om2*L2*cos(th1-th2)
  dy4 <- dy4*2*sin(th1-th2)
  dy4 <- dy4/(L2*(2*m1+m2-m2*cos(2*(th1-th2))))

  return(list(c(dy1, dy2, dy3, dy4)))
}

deg2rad <- function(deg) {
  # convert degrees to radians
  return((deg / 180) * pi)
}

N <- 40  # (seconds)
times <- seq(0, N, length.out = 5000) # 5,000 timepoints from 0 to N seconds

# Initial state: (th1, om1, th2, om2)
state <- c(th1 = deg2rad( 90),   # angle from the vertical (degrees)
           om1 = 0,              # angular velocity - zero: pendulum at rest
           th2 = deg2rad(180),   # 180°
           om2 = 0)              # ... at rest

# Parameters
parm <- list(L1   = 1.8,     # first pendulum: 1.8 meters long ...
             m1   = 0.5,     #   ... with a mass of 0.5 kg
             L2   = 0.8,     # 0.8m
             m2   = 0.4)     # 0.4 kg

# =====================================================
state <- c(th1 = deg2rad( 5),
           om1 = 0,
           th2 = deg2rad(1),
           om2 = 0)

# Parameters
parm <- list(L1   = 1.8,
             m1   = 1.0,
             L2   = 0.1,
             m2   = 0.01)

# Solve this
result <- deSolve::ode(y = state,
                       times = times,
                       func = doublePendulum,
                       parms = parm,
                       method = "rk4")

# Trajectories: Theta over time

plotPendulum <- function(res) {
  # plot of the double pendulum results
  myCols <- colorRampPalette(
    c("#6d7c8a","#8d9199","#9c99a1","#a87d86","#c96772"))(nrow(res))

  oPar <- par(mfrow = c(1, 2))

  # Define a layout for the plots: the next three plots go into the
  # sections defined in the matrix.
  myMat <- matrix(c(1, 2, 1, 2, 1, 3), ncol=2, byrow = TRUE)

  # Set the layout
  layout(myMat, widths=c(1,1), heights=c(1,1))

  # First plot: y over time
  t <- result[, "time"]
  th1 <- result[, "th1"]  # theta of first pendulum
  th2 <- result[, "th2"]  # theta of second pendulum
  plot(t, th1,
       type = "l",
       ylim = c(min(c(th1, th2)), max(c(th1, th2))),
       xlab = "Time (sec)", ylab = "Theta 1 , Theta 2 (rad)",
       main = "Double Pendulum",
       col = myCols[1])
  points(t, th2, type = "l", col = myCols[length(myCols)])
  legend("bottomright",
         legend = c("Theta 1", "Theta 2"),
         lty = 1,
         bty = "n",
         col = c(myCols[1], myCols[length(myCols)]))

  # Second plot: phase space
  th1 <- result[, "th1"]
  th2 <- result[, "th2"]
  plot(th1, th2,
       xlab = "Theta1", ylab = "Theta2",
       main = "Phase Space Plot of Double Pendulum",
       type = "n",
       col = myCols)
  plotrix::color.scale.lines(th1, th2, col = myCols, lwd = 2)
  legend("bottomright",
         legend = c("t = 0",
                    sprintf("t = %d",
                            round(result[nrow(result), "time"]))),
         fill =  c(myCols[1], myCols[length(myCols)]),
         bty = "n")

  # Third plot: Power Spectra of theta 1 and theta 2
  # x-axis points are frequencies (discretized to the stated bandwidth)
  # y-axis values are the "power" of a particular frequency in the time-series,
  #     i.e. how much a harmonic function of the frequency at the x-axis
  #     contributes to the time series data.
  sp1 <- spectrum(th1, log = "no", plot = FALSE)
  sp2 <- spectrum(th2, log = "no", plot = FALSE)
  plot(sp1$freq, sp1$spec,
       xlim=c(0, 50/length(th1)),
       ylim = c(0, max(c(sp1$spec, sp2$spec))),
       xlab = sprintf("frequency (bandwith = %.5f)", sp1$bandwidth),
       ylab = "spectrum",
       main = "Power Spectrum of Double Pendulum",
       type = "n")
  points(sp1$freq, sp1$spec, type ="l", lwd=3, col=myCols[1])
  points(sp2$freq, sp2$spec, type="l", lwd=1.5, col=myCols[length(myCols)])
  legend("topright",
         legend = c("Theta 1", "Theta 2"),
         lty = 1,
         bty = "n",
         col = c(myCols[1], myCols[length(myCols)]))

  par(oPar)

}

N <- 20  # (seconds)
times <- seq(0, N, length.out = 5000) # 5,000 timepoints from 0 to N seconds

# Initial state: theta (angle) and omega (angular velocity)
state <- c(th1 = deg2rad( 15), om1 = 0, th2 = deg2rad(2), om2 = 0)

# Pendulum parameters: lengths and masses
parm <- list(L1 = 1.8, m1 = 1.0, L2 = 0.7, m2 = 0.1)

# Solve this
result <- deSolve::ode(y = state,
                       times = times,
                       func = doublePendulum,
                       parms = parm,
                       method = "rk4")

# Plot it (Open a new graphics window with dev.new() and make it large.)
plotPendulum(result)


# Next

state <- c(th1 = deg2rad( 45), om1 = 0, th2 = deg2rad(30), om2 = 0)
parm <- list(L1 = 1.6, m1 = 0.6, L2 = 1.2, m2 = 0.4)
result <- deSolve::ode(y = state,  times = times, func = doublePendulum,
                       parms = parm, method = "rk4")
plotPendulum(result)

# Next

state <- c(th1 = deg2rad( 70), om1 = 0, th2 = deg2rad(80), om2 = 0)
parm <- list(L1 = 1.4, m1 = 0.8, L2 = 1.4, m2 = 0.35)
result <- deSolve::ode(y = state,  times = times, func = doublePendulum,
                       parms = parm, method = "rk4")
plotPendulum(result)


# Next

N <- 20  # (seconds)
times <- seq(0, N, length.out = 5000)

state <- c(th1 = deg2rad( 171), om1 = 0, th2 = deg2rad(98), om2 = 0)
parm <- list(L1 = 0.87, m1 = 1.15, L2 = 0.92, m2 = 1.33)
result <- deSolve::ode(y = state,  times = times, func = doublePendulum,
                       parms = parm, method = "rk4")
plotPendulum(result)

N <- 200  # (seconds)
times <- seq(0, N, length.out = 10000)
result <- deSolve::ode(y = state,  times = times, func = doublePendulum,
                       parms = parm, method = "rk4")
plotPendulum(result)



# =    7  A Neural Network at the Edge of Chaos  ===============================

# Let's design a simple neural network "at the edge of order and chaos".
# Assume the input is an oscillator (and maybe some directed perturbation).
#
# Let's build a network that has a complex output that can have chaotic
# behaviour but a well defined attractor domain. Let's take frequency and
# amplitude as the relevant phase space.
#
# How ? ...
#

# = 1. Code glossary






# [END]
