#' Fit ellipse
#'
#' Given n points in 2D cartesian coords compute the fittest ellipse to that points.
#'
#' @param x vector of values of the points to fit the ellipse on the x axis
#' @param y vector of values of the points to fit the ellipse on the y axis
#'
#' @return list of parameters needed to represent and a ellipse like centroid, length of axis, angle of the ellipse and regression coefficients
#' @noRd
fit.ellipse <- function (x, y = NULL) {

  EPS <- 1.0e-8
  dat <- xy.coords(x, y)

  D1 <- cbind(dat$x * dat$x, dat$x * dat$y, dat$y * dat$y)
  D2 <- cbind(dat$x, dat$y, 1)
  S1 <- t(D1) %*% D1
  S2 <- t(D1) %*% D2
  S3 <- t(D2) %*% D2
  T <- -solve(S3) %*% t(S2)
  M <- S1 + S2 %*% T
  M <- rbind(M[3,] / 2, -M[2,], M[1,] / 2)
  evec <- eigen(M)$vec
  cond <- 4 * evec[1,] * evec[3,] - evec[2,]^2
  a1 <- evec[, which(cond > 0)]
  f <- c(a1, T %*% a1)
  names(f) <- letters[1:6]


  A <- matrix(c(2*f[1], f[2], f[2], 2*f[3]), nrow=2, ncol=2, byrow=T )
  b <- matrix(c(-f[4], -f[5]), nrow=2, ncol=1, byrow=T)
  soln <- solve(A) %*% b

  b2 <- f[2]^2 / 4

  center <- c(soln[1], soln[2])
  names(center) <- c("x", "y")

  num  <- 2 * (f[1] * f[5]^2 / 4 + f[3] * f[4]^2 / 4 + f[6] * b2 - f[2]*f[4]*f[5]/4 - f[1]*f[3]*f[6])
  den1 <- (b2 - f[1]*f[3])
  den2 <- sqrt((f[1] - f[3])^2 + 4*b2)
  den3 <- f[1] + f[3]

  semi.axes <- sqrt(c( num / (den1 * (den2 - den3)),  num / (den1 * (-den2 - den3)) ))

  # calculate the angle of rotation
  term <- (f[1] - f[3]) / f[2]
  angle <- atan(1 / term) / 2

  list(coef = f, center = center, major = max(semi.axes), minor = min(semi.axes), angle = unname(angle))
}

#' Sample points for a ellipse
#'
#' Given the set of parameters centroid, length of axis, angle of the ellipse and regression coefficients sample n points from that ellipse
#'
#' @param fit list of parameters centroid, length of axis, angle of the ellipse and regression coefficients
#' @param n is an integer that denotes the number of points to sample
#'
#' @return matrix nx3 where each column denotes the cartesian coordinates of each one of the n points that determines the ellipse
#' @noRd
get.ellipse <- function(fit, n = 360 )
{

  tt <- seq(0, 2*pi, length=n)
  sa <- sin(fit$angle)
  ca <- cos(fit$angle)
  ct <- cos(tt)
  st <- sin(tt)

  x <- fit$center[1] + fit$maj * ct * ca - fit$min * st * sa
  y <- fit$center[2] + fit$maj * ct * sa + fit$min * st * ca

  cbind(x = x, y = y)
}

#' Compute rotation matrix to rotate a given angle around a given vector
#'
#' Given a vector in cartesian coordinates and an angle compute rotation matrix to rotate that angle around the given vector
#'
#' @param u is a perpendicular vector to the vectors that are going to be aligned
#' @param theta is the rotation angle in radians
#'
#' @return a 3x3 rotation matrix
#' @noRd
rotateMatrix<-function(u, theta)
{
  rotationMatrix <- t(matrix(c(cos(theta)+u[1]^2*(1-cos(theta)), u[1]*u[2]*(1-cos(theta))-u[3]*sin(theta), u[1]*u[3]*(1-cos(theta))+u[2]*sin(theta),
                               u[2]*u[1]*(1-cos(theta))+u[3]*sin(theta), cos(theta)+u[2]^2*(1-cos(theta)), u[2]*u[3]*(1-cos(theta))-u[1]*sin(theta),
                               u[3]*u[1]*(1-cos(theta))-u[2]*sin(theta), u[3]*u[2]*(1-cos(theta))+u[1]*sin(theta), cos(theta)+u[3]^2*(1-cos(theta))),nrow=3,ncol=3))

  return(rotationMatrix)
}

#'  Aling vector with the Z axis
#'
#' Given a vector in cartesian coordinates compute the translation and rotatio matrices needed to place that vector on the origin and aligned with Z axis
#'
#' @param vector
#'
#' @return list of matrices. First matrix is the translation, second is rotation along X axis and the third one refers to the rotation along Y axis
#' @noRd
align3dZ <-function(vector)
{
  #Translacion al origen
  trans<-trans3d(-vector)
  itrans <- trans3d(vector);

  u <- vector/sqrt(sum(vector^2));
  uz <- u;  uz[1] = 0; uz[2] = 0;
  uzy <- u; uzy[1] = 0;


  modUz <- uz[3];
  modUzy <- sqrt(sum(uzy^2))
  radsX <- acos(modUz/modUzy);
  if (u[2] < 0) {
    radsX <- -radsX;
  }
  rotX <- rotation3dX(radsX);
  irotX <- rotation3dX(-radsX);

  modUzx <- 1;
  modUzBis <- modUzy;
  radsY <- acos(modUzBis / modUzx);
  if (u[1] < 0) {
    radsY = -radsY;
  }

  rotY <- rotation3dY(-radsY);
  irotY <- rotation3dY(radsY);
  return(list(trans, rotX, rotY))
}

#' Cartesian coords to spherical
#'
#' Given a vector with origin 0 cast its cartesian coordinates to spherical.
#'
#' @param x is a double value which denotes the value of the vector tip along the x axis
#' @param y is a double value which denotes the value of the vector tip along the y axis
#' @param z is a double value which denotes the value of the vector tip along the z axis
#'
#' @return a matrix Nx3 where first column is the azimuth angle, the second one is the elevation angle and the third one is the length of the vector
#' @noRd
cartesian2Spherical <- function(x,y,z)
{
  azimuth <- atan2(y,x)
  elevation <- atan2(z,sqrt(x^2 + y^2))
  r <- sqrt(x^2 + y^2 + z^2)
  return(cbind(azimuth,elevation,r))
}

#' Spherical coords to cartesian.
#'
#' Given a vector azimuth angle, elevation angle and vector length cast them to cartesian coordinates
#'
#' @param azimuth is a double value which denotes the azimuthal angle
#' @param elevation is a double value which denotes the elevation angle
#' @param r is a double value that denotes the length of the vector
#'
#' @return a matrix Nx3 where first column is the value along X axis, the second one is the value along Y axis and the third one is the value along Z axis
#' @noRd
spherical2Cartesian<-function(azimuth,elevation,r)
{
  x <- r * cos(elevation) * cos(azimuth)
  y <- r * cos(elevation) * sin(azimuth)
  z <- r * sin(elevation)
  return(cbind(x, y, z))
}

#' Compute cross product
#'
#' Compute cross product between two vectors
#'
#' @param a first vector
#' @param b second vector
#'
#' @return result which is a vector perpendicular to a and b
#' @noRd
vcrossp <- function( a, b ) {
  result <- matrix( NA, 1, 3 )
  result[1] <- a[2] * b[3] - a[3] * b[2]
  result[2] <- a[3] * b[1] - a[1] * b[3]
  result[3] <- a[1] * b[2] - a[2] * b[1]
  return(result)
}

#' Compute norm of a vector
#'
#' Compute the length of a given vector
#'
#' @param a it is a vector
#'
#' @return the length of the vector
#' @noRd
normv<- function(a)
{
  a<-as.numeric(a)
  return(sqrt(a%*%a))
}
