library(randomForest)
library(miscTools)
library(ggplot2)

#fit a random forest model now



# import data from csv file
x <- read.csv('train2.csv',header=FALSE)
n <- nrow(x)/10
train <- x[1:(0.8*n),]
print(nrow(train))
test <- x[(0.8*n+1):n,]
sz1 <- 53
sz2 <- 21
sz <- sz1*sz2 # 1113 size of summary statistics

# start by fitting a simple linear model
# before moving on to the nnet model

models <- list()
responsenames <- paste("V", 1:7, sep='')
varnames <- paste("V", 8:(sz+7), sep='') 

for (y in responsenames){
  print(y)
  form <- formula(paste(y, "~ ."))  #, varnames))  # varnames))
  print(form)
  models[[y]] <- randomForest(form, data=train, ntree=20)
  }

z <- test
for (j in 1:7){
  lm.predict <- predict(models[[responsenames[j]]], z)

  #mse
  print(mean((lm.predict - z[,j])^2)) 
  #Rsquared
  r2 <- rSquared(z[,j], z[,j] - predict(models[[responsenames[j]]], z))
  print(r2)
 
plot(z[,j], lm.predict,
    main="RF regression predictions vs actual",
    xlab="Actual")
#NB currently the 7th piece of info is spurious, so should not get a decent (predictive) model for that 

#  p <- ggplot(aes(x=actual, y=pred),
#    data=data.frame(actual=z[,j], pred=lm.predict))
#  p + geom_point() +
#	geom_abline(color="red") +
#	ggtitle(paste("RandomForest Regression in R r^2=", r2, sep=""))
}
