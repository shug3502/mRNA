import numpy as np
import random as rnd
from operator import add
from matplotlib.pylab import *

def matrixtolinearindexing(i,j,m,n):
	#check this
	k = j*m+i
	return k

def lineartomatrixindexing(k,m,n):
	i = k%m
	j = (k-i)/m
	return [i,j]

def is_boundary(m,n):
	is_boundary_cell = [1]*n
	temp_col = [0]*m
	temp_col[0] = 1
	temp_col[m-1] = 1

	for r in range(n-2):
		for q in range(m):
			is_boundary_cell.append(temp_col[q])
	for qq in range(m):
		is_boundary_cell.append(1)
	return is_boundary_cell

def calculate_valency(m,n):
	# for rectangular domain
	lmbda = [3]*m
	lmbda[0] = 2
	lmbda[m-1] = 2
	valency = lmbda #first column of matrix
	j = 1
	new_col = [x+1 for x in lmbda]
	while j<(n-1):
		for w in range(len(new_col)):
			valency.append(new_col[w])
		j+=1
	for v in range(m):
		valency.append(lmbda[v])
	return valency

def calculate_propensities(diff_matrix,valency,X):
	#alph = np.multiply(diff_matrix, np.multiply(valency, X))
    # Bug noticed involving getting stuck in corners where the propensity is lower.
    #attempted fix
    alph = np.multiply(diff_matrix, np.multiply(4, X))

    return alph


def restricted_neighbourhood_indices(cell,neighbour_indices):
	#restrict set of neighbours to take account of the boundary
	[i,j] = lineartomatrixindexing(cell,m,n)
	#print 'i,j', i, j
	if is_boundary_cell[cell]:
		#assumes here a square domain
		if [i,j]==[0,0]: #cell == 0:
			new_indices = [[1,0],[0,1]]
		elif j==0 and i>0 and i<m-1: #cell>0 and cell<n-1:
			new_indices = [[1,0],[-1,0],[0,1]]
		elif [i,j]==[m-1,0]: #cell == n-1:
			new_indices = [[-1,0],[0,1]]
		elif i==0 and j>0 and j<n-1: #cell>n-1 and cell<n**2-n and cell%n==0:
			new_indices = [[1,0],[0,1],[0,-1]]
		elif i==(m-1) and j>0 and j<n-1: #cell>n-1 and cell<n**2-n and cell%n!=0:
			new_indices = [[-1,0],[0,-1],[0,1]]
		elif [i,j]==[0,n-1]: #cell==n**2-n:
			new_indices = [[1,0],[0,-1]]
		elif j==(n-1) and i>0 and i<m-1: #cell>n**2-n and cell<n**2-1:
			new_indices = [[-1,0],[0,-1],[1,0]]
		elif [i,j]==[m-1,n-1]: #cell==n**2-1:
			new_indices = [[-1,0],[0,-1]]
		else:
			print 'woops'
	else:
		new_indices = neighbour_indices
	return new_indices



# initialise
h=0.05
D=1.0/400
d=D/(h*h)
n=11
m=11
X = [0]*n*m

#Initial condition
#randomaly pick a position within the nucleus?
X[50] = 8

Total_occupancy = np.array(X)
f = open('output.txt', 'w')
#X_plot = []
#X_plot.append(X)

diff_matrix = [d]*n*m #homogeneous diffusion
#alternatively include inhomogenous diffusion here
is_boundary_cell = is_boundary(m,n)


neighbour_indices = [[x,y] for x in [-1,0,1] for y in [-1,0,1] if abs(x)+abs(y)==1]
valency = calculate_valency(m,n)

final_time = 1000
time=0
num_jumps=0

while time <final_time:
	#calculate propensities
	alpha = calculate_propensities(diff_matrix, valency, X)
	#print alpha
	a0 = np.sum(alpha)

	#update time
	rr = rnd.random()
	tau = -1/a0*np.log(rr)

	#find position of reaction
	rr2 = rnd.random()
	ss=0
	k=-1
	while rr2>ss:
		k+=1
		ss += alpha[k]


	#Now we've found voxel to jump out of. To find target of jump:
	rr3 = rnd.random()
	sss=0
	a_index=-1
	#print "valency", valency[k]
	while rr3>sss:
		sss+=1.0/valency[k]
		a_index+=1
		#print "sss", sss, "a", a_index


	#print k
	#print lineartomatrixindexing(k,m,n)
	new_indices = restricted_neighbourhood_indices(k,neighbour_indices)
	print a_index
	print new_indices[a_index]
	temp_index = map(add, lineartomatrixindexing(k,m,n), new_indices[a_index])
	#print temp_index
	target_index = matrixtolinearindexing(temp_index[0],temp_index[1],m,n)
	#print "target", target_index

	X[k] -=1
	X[target_index] +=1
	num_jumps +=1


	time += tau
	#print X
	#X_temp = X
	#X_plot.append(X_temp)
	for h in range(len(X)):
		f.write(str(X[h]))
	f.write('\n')
	#print time

#print "X_plot", X_plot

#print X
f.close()
g = open('output.txt', 'r')
shape = (m,n)


for index in range(num_jumps):
	#print "reading...", '\n'
	X_temp = g.readline().strip()
	numbers =[int(e.strip()) for e in X_temp]
	X_temp2 = np.array(numbers)
	Total_occupancy = [x+y for x,y in zip(numbers, Total_occupancy)]
	X_temp2 = X_temp2.reshape(shape)
	#matshow(X_temp2, fignum=100+index, cmap=cm.gray)
X_end_temp = np.array(X)
X_end = X_end_temp.reshape(shape)
print X_end
matshow(X_end, fignum=99, cmap=cm.gray)
print 'occupancy', Total_occupancy
Total_occupancy = np.array(Total_occupancy)
Total_jumps_in = Total_occupancy.reshape(shape)
matshow(Total_jumps_in, fignum =150)
show()


print num_jumps
