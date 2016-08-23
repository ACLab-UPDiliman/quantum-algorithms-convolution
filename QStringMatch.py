from numpy import zeros, ones, identity, empty, array, dot, sort
from math import sqrt
import matplotlib.pyplot as plotter

from CEncoder import cencode_BIV_a, cencode_BIV_b
from QDeleter import qdelete_matrix, qdelete_arbitrary
from QEncoderFinal import qencode
from QConvolute import q_convolute
from QOperators import qft, iqft, euclid_norm, Euclid_Norm_Squared

def construct_U_M(convolute):
	dim = convolute.shape[0]
	match_less_mismatch = empty(dim,dtype=complex)
	for i in xrange(0,dim):
		match_less_mismatch[i] = complex(convolute[i].real * sqrt(N*M), convolute[i].imag * sqrt(N*M) )

	num = empty(dim)
	for i in xrange(0,dim):
		num[i] = (match_less_mismatch[i].real + M) / sqrt(N*M)

	denum = empty(dim)
	for i in xrange(0,dim):
		denum[i] = euclid_norm(num[i])

	m = zeros(shape=(dim,dim))
	for i in xrange(0,dim):
		m[i][i] = num[i] / denum[i]

	print '\nm: '
	print m

	return m

def constructU_z(convolute, encoding_state_T):
	dim = convolute.shape[0]
	z = identity(dim)
	qft_T = qft(encoding_state_T)
	sum = 0.0
	for i in xrange(0,dim):
		sum = sum + 1.0 / euclid_norm(qft_T[i])
	for j in xrange(M-1,N):
		numerator = (1.0/2.0) + (((M/dim)*sum)/(2*convolute[j]))
		denumerator = euclid_norm(numerator)
		z[j][j] = numerator/denumerator

	return z

def adjust(convolute, N, M):
	# Delete first and last M-1 basis states
	# c = qdelete_matrix(N,M,convolute,array([[M-1,N]]))
	c_array = []
	for i in xrange(M-1, N):
		c_array.append(i)

	tau_array = []
	for i in xrange(0, M-1):
		tau_array.append(i)
	for i in xrange(N,N+M-1):
		tau_array.append(i)

	# c_array = array([[M-1, N]])
	# tau_array = array([[0, M-1],[N, N+M-1]])
	c = qdelete_arbitrary(N+M-1, convolute, c_array, tau_array)
	return c

def getProbabilities(a):
	dim = a.shape[0]
	b = empty(dim, dtype=float)
	total_prob = 0.0
	for i in xrange(0,dim):
		b[i] = a[i].real**2 + a[i].imag**2
		total_prob = total_prob + b[i]
	return array([b,total_prob])

def prettyprint_probabilities(a):
	# map each item in array to an index
	mapping = [('index',int),('amp',float)]
	# copy array a to list b
	b = list(a)
	b_dim = len(b)
	print 'b_dim = ' + str(b_dim)

	# for each item in list b map an index to create a tuple, start from index 0
	for i in xrange(0,b_dim):
		b[i] = (i,b[i])

	# associate the metadata in mapping to each tuple in list b
	b = array(b, dtype=mapping)
	print 'b_dim = ' + str(b_dim)
	print 'b :'
	print b

	# copy the middle N-M+1 tuples from b to b
	# b = list(b[M-1:N])
	# print 'b_dim = ' + str(b_dim)
	# print 'b :'
	# print b

	# sort the tuples in list b according to each tuple's second element, the probability
	b = sort(b, kind='mergesort', order='amp')

	# reverse the sorting of tuples in list b into order from highest probability to lowest 
	b = b[::-1]
	b_dim = len(b)
	print 'b_dim = ' + str(b_dim)
	string = ''
	
	# construct a string for printing the tuple items in list b
	pad = ''
	for i in xrange(0,b_dim):
		if (b[i][0] + (2-M)) < 100:
			pad = pad + '0'
			if (b[i][0] + (2-M)) < 10:
				pad = pad + '0'
		print 'printing for text[ ' + str(b[i][0] + (2-M) - 1) + ' ]'
		
		if (b[i][0] + (2-M)) > 0 and ((b[i][0] + (2-M)) <= (N-M+1)):
			string = string + '\n[' + pad + str(b[i][0] + (2-M)) +'] : ' + str(T[b[i][0] + (2-M) - 1 : b[i][0] + 1]) + ' : ' + str(b[i][1])
		
		elif ((b[i][0] + (2-M)) > N):
			string = string + '\n[+' + str(b[i][0] + (2-M) - N) +'] : N/A : ' + str(b[i][1])	
		
		else:
			string = string + '\n[' + str(b[i][0] + (2-M)) +'] : N/A : ' + str(b[i][1])

		pad = ''
	
	string = string + '\n'

	return string


def prettyprint_probabilities_simple(a):
	# map each item in array to an index
	mapping = [('index',int),('amp',float)]
	# copy array a to list b
	b = list(a)
	b_dim = len(b)
	print 'b_dim = ' + str(b_dim)

	# for each item in list b map an index to create a tuple, start from index 0
	for i in xrange(0,b_dim):
		b[i] = (i,b[i])

	# associate the metadata in mapping to each tuple in list b
	b = array(b, dtype=mapping)
	print 'b_dim = ' + str(b_dim)
	print 'b :'
	print b

	# copy the middle N-M+1 tuples from b to b
	# b = list(b[M-1:N])
	# print 'b_dim = ' + str(b_dim)
	# print 'b :'
	# print b

	# sort the tuples in list b according to each tuple's second element, the probability
	b = sort(b, kind='mergesort', order='amp')

	# reverse the sorting of tuples in list b into order from highest probability to lowest 
	b = b[::-1]
	b_dim = len(b)
	print 'b_dim = ' + str(b_dim)
	string = ''
	
	# construct a string for printing the tuple items in list b
	pad = ''
	for i in xrange(0,b_dim):
		if (b[i][0] + (2-M)) < 100:
			pad = pad + '0'
			if (b[i][0] + (2-M)) < 10:
				pad = pad + '0'

		string = string + '\n[' + pad + str(b[i][0] + (2-M)) +'] : ' + str(b[i][1])
		pad = ''
	
	string = string + '\n'

	return string

def slice_probabilities(a):
	# copy array a to list b
	b = a[M-1:N]

	return b

#=== INPUT
# sample input text and pattern
T = [ 'a', 'a', 'b', 'b', 'b', 'a' ]
P = [ 'b', 'b', 'a' ]

# text = [
# 		'b','a','a','a','b','b','b','b','b','b',
# 		'b','b','b','b','b','b','b','b','b','b',
# 		'b','b','b','b','b','b','b','b','b','b',
# 		'b','b','b','b','b','b','b','b','b','b',
# 		'b','b','b','b','b','b','b','b','b','b',
# 		'b','b','b','b','b','b','b','b','b','b',
# 		'b','a','a','a','b','b','b','b','b','b',
# 		'b','b','b','b','b','b','b','b','b','b',
# 		'b','b','b','b','b','b','b','b','b','b',
# 		'b','b','b','b','b','b','b','b','b','b',
# 		'b','b','b','b','b','b','b','b','b','b',
# 		'b','b','b','b','b','b','b','b','b','b',
# 		'b','b','a','a','a'
# 	]
# pattern = ['b','a','a','a']

N = len(T)
M = len(P)
dim = N+M-1
mid_dim = N-M+1
print 'text: ' + str(T) + ', N = ' + str(N)
print 'pattern: ' + str(P) + ', M = ' + str(M)
print 'dim: ' + str(dim) + ', middle dim: ' + str(mid_dim)

#=== ENCODE
# encode text and pattern into quantum states
print '============ ENCODING PHASE START ============\n'
# BIV_encoded_strings = cencode_BIV_b(text,pattern)
# T_cencoded = BIV_encoded_strings[0]
# P_cencoded = BIV_encoded_strings[1]
# encoding_states = qencode(T_cencoded,P_cencoded)
# a
# encoding_states = array([[sqrt(1.0/3.0), sqrt(1.0/3.0),0.0,0.0,0.0,sqrt(1.0/3.0),0.0,0.0],[1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]])
# b
encoding_states = array([[0.0,0.0,sqrt(1.0/3.0),sqrt(1.0/3.0),sqrt(1.0/3.0),0.0,0.0,0.0],[0.0,sqrt(1.0/2.0),sqrt(1.0/2.0),0.0,0.0,0.0,0.0,0.0]])
# encoding_states = array([[0.0,0.0,1.0,1.0,1.0,0.0,0.0,0.0],[0.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0]])
prob_T = getProbabilities(encoding_states[0])
prob_P = getProbabilities(encoding_states[1])
print 'Encoding state for text: ' + prettyprint_probabilities_simple(encoding_states[0])
print 'Encoding state for pattern: ' + prettyprint_probabilities_simple(encoding_states[1]) + '\n'
print 'Probabilities for text: ' + prettyprint_probabilities_simple(prob_T[0])
print 'Total probability for text: ' + str(prob_T[1]) + '\n'
print 'Probabilities for pattern:' + prettyprint_probabilities_simple(prob_P[0])
print 'Total probability for pattern: ' + str(prob_P[1])
print '============ ENCODING PHASE END ==============\n'

#=== CONVOLUTE
print '============ QUANTUM CONVOLUTION PHASE START =\n'
convolute = q_convolute(encoding_states[0], encoding_states[1])
prob_convolute = getProbabilities(convolute)
print 'Convolute: ' + str(convolute) + '\n'
print 'Probabilities for convolute: ' + prettyprint_probabilities_simple(prob_convolute[0])
print 'Total probability for convolute: ' + str(prob_convolute[1])
print '============ QUANTUM CONVOLUTION PHASE END ===\n'

# # # corresponds to adding M to each element of vector z and dividing the 
# # # resulting value by 2 to get the match count at each index of text
print '============ ADJUSTMENT PHASE START ==========\n'
# # U_z = constructU_z(convolute, encoding_states[0])
# # vector_c = dot(U_z,convolute)
# # prob_c = getProbabilities(vector_c)
# # print 'Vector C: ' + str(vector_c) + '\n'
# # print 'Probabilities for Vector C: ' + prettyprint_probabilities(prob_c[0][0])
# # print 'Total probability for Vector C: ' + str(prob_c[0][1])
# # deletes the first and last M-1 basis states of resulting vector
# # and distributes the amplitude evenly into the retained basis states

# U_m = construct_U_M(convolute)
# matches = dot(U_m,convolute)
# prob_matches = getProbabilities(matches)
# print 'matches vector: '
# print matches
# print 'Probabilities for matches vector: ' + prettyprint_probabilities(prob_matches[0])
# print 'Total probability for matches vector: ' + str(prob_matches[1]) + '\n'

# trimmed_vector_c = adjust(convolute, N, M) 
# prob_trimmed_c = getProbabilities(trimmed_vector_c)
# print 'trimmed vector c: '
# print trimmed_vector_c
# print 'Probabilities for trimmed Vector C: ' + prettyprint_probabilities(prob_trimmed_c[0])
# print 'Total probability for trimmed Vector C: ' + str(prob_trimmed_c[1])
print '============ ADJUSTMENT PHASE END ============\n'

# b = slice_probabilities(prob_convolute[0][0])
# plotter.title('Trimmed Probability of Occurence for Each Basis State')
# plotter.ylabel('Probability')
# plotter.xlabel('Basis state index')
# x = range(prob_trimmed_c[0].shape[0])
# width = 0.10
# plotter.bar(x, prob_trimmed_c[0], width, color="black", alpha=0.75)
# plotter.grid(True)
# plotter.show()

#=== OUTPUT
# probabilities = zeros(dim)
# total_probability = 0.0
# for i in xrange(0,dim):
# 	probabilities[i] = vector_c[i].real**2 + vector_c[i].imag**2  #(S_ket_T[i].conjugate())*S_ket_T[i]
# 	total_probability = total_probability + probabilities[i]
# print 'Final probabilities of basis states'
# for i in xrange(0,dim):
# 	print '|' + str(i) + '> : ' + str(probabilities[i])
# print 'Total probability: ' + str(total_probability)

# probabilities = zeros(dim)
# total_probability = 0.0
# for i in xrange(0,dim):
# 	probabilities[i] = trimmed_vector_c[i].real**2 + trimmed_vector_c[i].imag**2  #(S_ket_T[i].conjugate())*S_ket_T[i]
# 	total_probability = total_probability + probabilities[i]
# print 'Final probabilities of basis states'
# for i in xrange(0,dim):
# 	print '|' + str(i) + '> : ' + str(probabilities[i])
# print 'Total probability: ' + str(total_probability)


