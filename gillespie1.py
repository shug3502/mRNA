import numpy as np
import random

def matrixtolinearindexing(i,j,n,m):
	k = (i-1)*m+j
	return k


# initialise
h=0.05
D=1.0/400
n=20
m=20
X = [0]*n*m
X[0] = 10

final_time = 10

time=0

while time <final_time:
	#calculate propensities
	propensities = np.divide(np.multiply(X,4*D),(h**2))  #since von neumann neighbourhood
	print propensities
	
	temp_row = [2]
	temp_row.append(3 for i in range(1,(m-1)))
	temp_row.append(2)
	print temp_row
	
	coeffs = [2,[3]*(m-2),2]
	print coeffs
	for i in range(n-2):
		coeffs = coeffs.append([4]*m)
	coeffs = coeffs.append([2,[3]*(m-2),2])
	print coeffs	
	a0 = np.sum(propensities)
	#find time to next reaction
	rr = random()
	tau = -1/a0*np.log(rr)
	
	rr2 = random()
	#which reaction?
	ss=0
	while ss*D/(h*h) < rr2:
		ss+=1 
	
	
	
	
	time +=tau


