from numpy import array, dot, kron, zeros, ones, identity, outer
from scipy.fftpack import fft, ifft
from math import sqrt
import matplotlib.pyplot as plotter

T = ['a','a','b','b','b','a']
P = ['b','b','a']

P.reverse()
#a
encoding_vector_T = array([
							1,  
							1,
							0, 
							0,  
							0,
							1,  
							0,
							0
						])

encoding_vector_P = array([
							1, 
							0, 
							0, 
							0, 
							0, 
							0,
							0, 
							0])
#b
# encoding_vector_T = array([
# 							0, 
# 							0, 
# 							1,  
# 							1, 
# 							1, 
# 							0, 
# 							0,
# 							0
# 						])

# encoding_vector_P = array([
# 							0, 
# 							1, 
# 							1, 
# 							0, 
# 							0, 
# 							0,
# 							0, 
# 							0])


#a
# encoding_vector_T = array([
# 							sqrt(1.0/3),  
# 							sqrt(1.0/3),
# 							0, 
# 							0,  
# 							0,
# 							sqrt(1.0/3),  
# 							0,
# 							0
# 						])

# encoding_vector_P = array([
# 							1.0, 
# 							0, 
# 							0, 
# 							0, 
# 							0, 
# 							0,
# 							0, 
# 							0])

#b
# encoding_vector_T = array([
# 							0, 
# 							0, 
# 							sqrt(1.0/3),  
# 							sqrt(1.0/3), 
# 							sqrt(1.0/3), 
# 							0, 
# 							0,
# 							0
# 						])

# encoding_vector_P = array([
# 							0, 
# 							sqrt(1.0/2), 
# 							sqrt(1.0/2), 
# 							0, 
# 							0, 
# 							0,
# 							0, 
# 							0])
N = 6
M = 3
dim = N + M - 1

def DFT(a,t):
	m = sqrt(1.0/t) * fft(a)
	
	# m = fft(a)
	return m

def IDFT(a,t):
	m =  sqrt(1.0/t) * ifft(a)
	
	# m =  sqrt(1.0/t) * t * ifft(a)
	
	# m = ifft(a)
	return m

def normalize(a):
	dim = a.shape[0]
	b = get_probabilities(a)
	total = getTotalProbabilities(b)
	c = zeros(dim)
	for i in xrange(0,dim):
		c[i] = a[i]/sqrt(total)

	return c

def get_probabilities(a):
	dim = a.shape[0]
	b = zeros(dim)
	for i in xrange(0,dim):
		b[i] = a[i].real**2 + a[i].imag**2

	return b
 
def getTotalProbabilities(a):
	dim = a.shape[0]
	total = 0
	for i in xrange(0,dim):
		total = total + a[i]

	return total

def finalize(a):
	dim = a.shape[0]
	b = zeros(dim)
	print 'dim = ' + str(dim)
	for i in xrange(0,dim):
		b[i] = (a[i].real + M)/2

	return b

# circular convolute
def cconvolute(encoding_vector_T, encoding_vector_P):
	dim = encoding_vector_T.shape[0]
	DFT_T = DFT(encoding_vector_T, dim)
	DFT_P = DFT(encoding_vector_P, dim)
	DFT_T_x_DFT_P = DFT_T * DFT_P
	IDFT_DFT_T_x_DFT_P = IDFT(DFT_T_x_DFT_P, dim)
	probabilities = get_probabilities(IDFT_DFT_T_x_DFT_P)
	total_prob = getTotalProbabilities(probabilities)
	normalized = normalize(IDFT_DFT_T_x_DFT_P)
	normalized_probabilities = get_probabilities(normalized)
	total_normalized_prob = getTotalProbabilities(normalized_probabilities)
	finalized = finalize(IDFT_DFT_T_x_DFT_P)

	print 'DFT(text) = ' + str(DFT_T)
	print 'DFT(pattern) = ' + str(DFT_P)
	print 'DFT(text) * DFT(pattern) = ' + str(DFT_T_x_DFT_P)
	print 'IDFT(DFT(text) * DFT(pattern)) = ' + str(IDFT_DFT_T_x_DFT_P)
	# print 'Probabilities = ' + str(probabilities)
	# print 'Total Probabilities = ' + str(total_prob)
	# print 'Normalized amplitudes = ' + str(normalized)
	# print 'Normalized probabilities = ' + str(normalized_probabilities)
	# print 'Total normalized probabilities = ' + str(total_normalized_prob)
	# print 'Finalized vector = ' + str(finalized)

	plotter.title('Initial Probabilities of Basis States')
	plotter.ylabel('Probability')
	plotter.xlabel('Basis State')
	x = range(probabilities.shape[0])
	width = 0.10
	plotter.bar(x, probabilities, width, color="black", alpha=0.75)
	plotter.grid(True)
	plotter.show()


cconvolute(encoding_vector_T, encoding_vector_P)