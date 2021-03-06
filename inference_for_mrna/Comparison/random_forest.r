library(randomForest)
library(miscTools)
library(ggplot2)

#fit a random forest model now

# read parameter inputs from file
pp <- read.csv('params.csv',header=TRUE, stringsAsFactors=FALSE)

# import data from csv file
x <- read.csv(pp$filename,header=FALSE)
n <- nrow(x)
prop <- pp$prop
n.tree <- pp$ntree
train <- x[1:(prop*n),]
test <- x[(prop*n+1):n,]

sz <- pp$size_ss # size of summary statistics
m <- pp$size_theta

# start by fitting a simple linear model

models <- list()
responsenames <- paste("V", 1:m, sep='')
varnames <- paste("V", (m+1):(sz+m), sep='') 

for (y in responsenames){
	print(y)
  form <- formula(paste(y, "~ ."))  #, varnames))  # varnames))
  models[[y]] <- randomForest(form, data=train, ntree=n.tree)
  }

z <- test
for (j in 1:m){
  lm.predict <- predict(models[[responsenames[j]]], z)

  #loss based on mse
  loss <- mean((lm.predict/z[,j] - 1)^2) 
  line=paste("RF: ", loss)
  print(line)
  write(line,file="loss.csv",append=TRUE)

  #Rsquared
  # r2 <- rSquared(z[,j], z[,j] - predict(models[[responsenames[j]]], z))
  # print(r2)
 
plot(z[,j], lm.predict,
    main="RF regression predictions vs actual",
    xlab="Actual")
#NB currently the 7th piece of info is spurious, so should not get a decent (predictive) model for that

}
