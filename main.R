require(readxl)
require(functional)
require(ggplot2)

source("methods.R")

data <- read_xlsx("matrixAndVector.xlsx")
# Preparing data
vector <- as.vector(data$Vector)
vector2 <- as.vector(data$Vector2)

matrix <- as.matrix(data[,3:6])
matrix <- matrix[apply(matrix, 1, Compose(is.finite, all)),]
matrix2 <- as.matrix(data[,10:13])
matrix2 <- matrix2[apply(matrix2, 1, Compose(is.finite, all)),]
# Matrix generation
#write.csv(matrix(runif(60 * 1, min = 0, max = 1), ncol = 1, nrow = 60), file = "test1.csv")

#h <- as.numeric(readline(prompt = "Enter h: "))
#delta <- as.numeric(readline(prompt = "Enter delta: "))

transponentMatrix <- t(matrix)
# Multiply transponented matrix with general matrix
MTGmatrix <- transponentMatrix %*% matrix
# Multiply transponented matrix with vector
MTVcoefficients <- transponentMatrix %*% vector

transponentMatrix2 <- t(matrix2)
# Multiply transponented matrix with general matrix
MTGmatrix2 <- transponentMatrix2 %*% matrix2
# Multiply transponented matrix with vector
MTVcoefficients2 <- transponentMatrix2 %*% vector2

if (alpha1 > alpha2) {
  stop("alpha1 > alpha2")
}

x1 <- TihonovMethod(MTGmatrix, MTVcoefficients, alpha1 = 0.01, alpha2 = 0.1, delta = 0.0001, h = 0.0001, eps = 0.00001)
x2 <- TihonovMethod(MTGmatrix2, MTVcoefficients2, alpha1 = 0.01, alpha2 = 0.1, delta = 0.0001, h = 0.0001, eps = 0.00001)

graphicsData <- data.frame(vector, matrix %*% x1)
colnames(graphicsData) <- c("x", "y")

linearGraphicsData <- data.frame(vector2, matrix2 %*% x2)
colnames(linearGraphicsData) <- c("x", "y")

# Scatter diagram
ggplot(data = graphicsData, mapping = aes(x = graphicsData$x, y = graphicsData$y)) +
  theme_light() +
  geom_point(colour = "black", size = 1) +
  geom_abline()

ggplot(data = linearGraphicsData, mapping = aes(x = linearGraphicsData$x, y = linearGraphicsData$y)) +
  theme_light() +
  geom_point(colour = "red", size = 1) + 
  geom_abline()
