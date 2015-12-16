
# import data from csv file
x <- read.csv('train.csv',header=FALSE)
train <- x[1:800,]
test <- x[801:1000,]
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
  form <- formula(paste(y, "~ ."))  # varnames))
  models[[y]] <- lm(form, data=train) 
  }

#lm.fit <- lm(V1 ~ ., data=x)
for (j in 1:7){
lm.predict <- predict(models[[responsenames[j]]],test)
 
# mean squared error: 21.89483
print(mean((lm.predict - test[,j])^2)) 
 
plot(test[,j], lm.predict,
    main="Linear regression predictions vs actual",
    xlab="Actual")
}
#NB currently the 7th piece of info is spurious, so should not get a decent (predictive) model for that 
