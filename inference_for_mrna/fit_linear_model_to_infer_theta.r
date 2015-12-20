
# import data from csv file
x <- read.csv('train2.csv',header=FALSE)
n <- nrow(x)
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
  form <- formula(paste(y, "~ .")) #, varnames))  # varnames))
  models[[y]] <- lm(form, data=train) 
  }

z <- test
for (j in 1:7){
lm.predict <- predict(models[[responsenames[j]]],z)
 
# mean squared error: 21.89483
print(mean((lm.predict - z[,j])^2)) 
 
plot(z[,j], lm.predict,
    main="Linear regression predictions vs actual",
    xlab="Actual")
}
#NB currently the 7th piece of info is spurious, so should not get a decent (predictive) model for that 
