require(readxl)
require(functional)
require(ggplot2)

data <- read_xlsx("matrixAndVector.xlsx")
# Preparing data
vector <- as.vector(data$Vector)
matrix <- as.matrix(data[,3:ncol(data)])
matrix <- matrix[apply(matrix, 1, Compose(is.finite, all)),]
#write.csv(matrix(runif(60 * 1, min = 0, max = 1), ncol = 1, nrow = 60), file = "test1.csv")

#h <- as.numeric(readline(prompt = "Enter h: "))
#delta <- as.numeric(readline(prompt = "Enter delta: "))
h <- 0.0001
delta <- 0.0001
eps <- 0.00001

transponentMatrix <- t(matrix)
# Multiply transponented matrix with general matrix
MTGmatrix <- transponentMatrix %*% matrix
# Multiply transponented matrix with vector
MTVcoefficients <- transponentMatrix %*% vector

alpha1 <- 0.01
alpha2 <- 0.1

if (alpha1 > alpha2) {
  stop("alpha1 > alpha2")
}

E <- diag(1,ncol(MTGmatrix))
# Approximation to solve an ill-posed problem
psi <- matrix(c(0,0),nrow = ncol(MTGmatrix), ncol = 1)

x1 <- solve(MTGmatrix + alpha1 * E, MTVcoefficients + alpha1 * psi)
x2 <- solve(MTGmatrix + alpha2 * E, MTVcoefficients + alpha2 * psi)

if (det(MTGmatrix) == 0) {
  stop("det of matrix has zero value")
} else {
  x <- solve(MTGmatrix, MTVcoefficients)
  myu <- norm((matrix %*% x) - vector)
  delta1 <- 2 * norm(x) * h + 2 * delta + myu
}

phia1 <- norm(matrix %*% x1 - vector) - delta1
phia2 <- norm(matrix %*% x2 - vector) - delta1

delta11 <- h*norm(x1) + delta1
delta12 <- h*norm(x2) + delta1
flag <- FALSE

# Tikhonov method
while (TRUE)
{
  if (phia1 * phia2 > 0 && phia1 * phia2 < eps / 100) {
    break
  }
  else if (phia1 * phia2 > 0) {
    alpha1 <- alpha1 / 2
    alpha2 <- alpha2 * 2
  }
  else if (phia1 * phia2 < 0) {
    flag <- TRUE
    break
  }
  
  x1 <- solve(MTGmatrix + alpha1 * E, MTVcoefficients + alpha1 * psi)
  delta11 <- h * norm(x1) + delta11
  x2 <- solve(MTGmatrix + alpha2 * E, MTVcoefficients + alpha2 * psi)
  delta12 <- h * norm(x2) + delta12
  
  phia1 <- norm(matrix %*% x1 - vector) - delta11
  phia2 <- norm(matrix %*% x2 - vector) - delta12
}

# Half-division method
while (flag)
{
  if (abs(alpha1 - alpha2) / 2 <= eps) {
    break
  }
  
  m <- (alpha1 + alpha2) / 2
  xm <- solve(MTGmatrix + m * E, MTVcoefficients + m * psi)
  deltam <- h * norm(xm) + delta
  phiam <- norm(matrix %*% xm - vector) - deltam
  if (phia1 * phiam < 0) {
    alpha2 = m
  } else {
    alpha1 = m
  }
}

x1 <- solve(MTGmatrix + alpha1 * E, MTVcoefficients + alpha1 * psi)

graphicsData <- data.frame(matrix %*% x1, vector)
colnames(graphicsData) <- c("x", "y")

# Scatter diagram
ggplot(data = graphicsData, mapping = aes(x = graphicsData$x, y = graphicsData$y)) +
  theme_light() +
  geom_point(fill = "black", size = 1) +
  geom_smooth()
