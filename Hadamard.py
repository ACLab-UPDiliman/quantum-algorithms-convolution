from math import sqrt
from numpy import array

def H2x2():
	a = array([[1,1],[1,-1]])
	h = (1/sqrt(2)) * a
	return h

def X():
	a = array([[0,1],[1,0]])
	return a
