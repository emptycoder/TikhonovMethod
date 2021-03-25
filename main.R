require(readxl)
require(functional)
require(ggplot2)

source("methods.R")

data <- read_xlsx("matrixAndVector.xlsx")
# Preparing data
vector <- as.vector(data$Vector)
matrix <- as.matrix(data[,3:ncol(data)])
matrix <- matrix[apply(matrix, 1, Compose(is.finite, all)),]
# Matrix generation
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

x1 <- TihonovMethod(MTGmatrix, MTVcoefficients, alpha1, alpha2, delta, h, eps)

graphicsData <- data.frame(matrix %*% x1, vector)
colnames(graphicsData) <- c("x", "y")

# Scatter diagram
ggplot(data = graphicsData, mapping = aes(x = graphicsData$x, y = graphicsData$y)) +
  theme_light() +
  geom_point(fill = "black", size = 1) +
  geom_smooth()
