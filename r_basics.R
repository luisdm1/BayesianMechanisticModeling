1+1
3/4
6*4

.libPaths()

print('Hello World')


?mean 

mean

library(dplyr)

install.packages('lmer')

library(utils)
sessionInfo()

getwd()
setwd("C:/Users/luis_/GitHub/BayesianMechanisticModeling")


data("iris")

head(iris)
tail(iris)

summary(iris)

names(iris)

head(iris$Petal.Length)

summary(iris$Petal.Length)

typeof(iris)
str(iris)


plot(iris)

plot(iris$Petal.Length, iris$Petal.Width)


iris[4]

iris[4, ]

iris[4, 'Petal.Width']


x = c()
x = c(1,2,3,4)
length(x)

dim(iris)

result = 0
for (idx in 1:10){
  result = result + idx
}
print(result)
result

if (4<6){
  print('4<6')
} else{
  print('4>6')
}

quadratic = function(x){
  return(x^2)
}
quadratic(9)
quadratic(c(1,2,3,4,5,6))
quadratic(seq(1, 100))

?rnorm


x = rnorm(n=1000000, mean = 6, sd=10)

summary(x)
sd(x)

hist(x)

