library(ggplot2)
data("iris")

ggplot(data = iris, aes(x=iris$Sepal.Length, y=iris$Petal.Length)) + geom_point()

#----------

# Library
library(ggplot2)
install.packages("hrbrthemes")
library(hrbrthemes)


# Dummy data
x <- LETTERS[1:20]
y <- paste0("var", seq(1,20))
data <- expand.grid(X=x, Y=y)
data$Z <- runif(400, 0, 5)

# Give extreme colors:
ggplot(data, aes(X, Y, fill= Z)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="blue") +
  theme_classic()

# Color Brewer palette
ggplot(data, aes(X, Y, fill= Z)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdPu") +
  theme_get()

# Color Brewer palette
ggplot(data, aes(X, Y, fill= Z)) + 
  geom_tile() +
  theme_bw()

