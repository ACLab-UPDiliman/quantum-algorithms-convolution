from numpy import array, dot, kron, zeros, ones, identity, outer, array, empty
from scipy.fftpack import fft, ifft
from math import sqrt, log, ceil
from Util import BinaryNumber, Bit, Symbol


def H2x2():
	a = array([[1,1],[1,-1]])
	h = (1/sqrt(2)) * a
	return h

def H2x2_kron_log_dim(dim):
	log_dim = int(ceil(log(dim,2)))
	# print 'log(' + str(dim) + ') = ' + str(log_dim)
	h = H2x2()
	for i in xrange(0,log_dim-1):
		h = kron(h,H2x2())
	return h

def QFT(a):
	dim = a.shape[0]
	factor = sqrt(1.0/dim)
	m = factor * fft(a)
	return m

def IQFT(a):
	dim = a.shape[0]
	factor = sqrt(1.0/dim)
	# print 'IQFT mult factor = ' + str(factor)
	# print 'a = ' + str(a)
	#NOTE: In IFFT, the elements of the matrix are multiplied by factor 1/dim
	#		so that DFT*IFFT = I. So, in the quantum setting, we need to compensate
	#		for the inherently included multiplicative factor 1/dim in each element of
	#		IFFT by factoring in dim. This will make our quantum version of IFFT,
	#		IQFT, be equal to (1/dim)*(IFFT) without the previous factor dim inside
	#		IFFT.
	m = factor * dim * ifft(a)
	return m

####### Quantum half adder components ########
def NOT(x):
	if x == 0:
		not_x = 1
	else:
		not_x = 0
	return not_x

def CNOT(control,target):
	if control == 1:
		if target == 0:
			target = 1
		else:
			target = 0
	return (control,target)

def CCNOT(control1,control2,target):
	if control1 == 1 and control2 == 1:
		if target == 0:
			target = 1
		else:
			target = 0
	return (control1,control2,target)

def bitToBitSum(x_i,y_i):
	return CNOT(x_i,y_i) # XOR

def bitToBitCarry(x_i,y_i,aux):
	return CCNOT(x_i,y_i,aux) # AND

# compute for total carry (previous carry + current carry) for a bit-to-bit full-adder
# input:
# 	- adder1_carry = current carry; carry from first half-adder (carry(x_i + y_i))
#	- adder2_carry = partial carry; carry from second half-adder (carry((x_i + y_i) + previous carry))
# 	- aux = auxilliary bit set to 0
# output:
#	- 
def bitToBitTotalCarry(adder1_carry,adder2_carry,aux):
	# carry from first half-adder OR carry from second half-adder
	ccnot_result = CCNOT(adder1_carry,adder2_carry,aux)
	cnot_result1 = CNOT(ccnot_result[0],ccnot_result[1])
	cnot_result2 = CNOT(cnot_result1[1],ccnot_result[2])
	return (cnot_result1[0],cnot_result2[0],cnot_result2[1])

def bitToBitHalfAdder(x_i,y_i,aux):
	carry_result = bitToBitCarry(x_i,y_i,aux)[2] # AND
	sum_result = bitToBitSum(x_i,y_i)[1] # XOR
	return (x_i,sum_result,carry_result)

def bitToBitFullAdder(x_i,y_i,previous_carry):
	# first half-adder
	# compute for sum of current bits x_i and y_i, and carry (current_carry) in first half-adder
	# input:
	# - x_i = i-th bit of x
	# - y_i = i-th bit of y
	# - aux1 = auxilliary bit set to 0
	# output:
	# - add_result1[0] = x_i
	# - add_result1[1] = x_i + y_i (XOR) (current sum)
	# - add_result1[2] = carry(x_i + y_i) (AND) (current carry)
	aux1 = 0
	add_result1 = bitToBitHalfAdder(x_i,y_i,aux1)

	# second half-adder
	# compute for sum and carry (current_carry) in second half-adder
	# input:
	# - previous_carry = carry from previous bit-to-bit adder
	# - current_sum = add_result1[1] = x_i + y_i (XOR) (current sum)
	# - aux2 = auxilliary bit set to 0
	# output:
	# - previous carry = add_result2[0]
	# - total sum = add_result2[1] (x_i + y_i) + previous_carry (XOR) (total sum of adding current sum and previous carry)
	# - partial carry = add_result2[2] = carry((x_i + y_i) + previous_carry) (AND) (carry from adding current sum and previous carry)
	current_sum = add_result1[1]
	aux2 = 0
	add_result2 = bitToBitHalfAdder(previous_carry,current_sum,aux2)
	total_sum = add_result2[1]
	
	# total carry
	# compute for total carry of adding x_i, y_i and previous carry
	# input:
	# - partial_carry = add_result2[2] = carry((x_i + y_i) + previous_carry) (AND) (carry from adding current sum and previous carry)
	# - current_carry = add_result1[2] = carry(x_i + y_i) (AND) (current carry)
	# - aux3 = auxilliary bit set to 0
	# output:
	# - carry_result[0] = partial_carry
	# - carry_result[1] = partial_carry XOR current_carry
	# - carry_result[2] = (partial_carry XOR current_carry) XOR (partial_carry AND current_carry) (total carry of the full-adder)
	partial_carry = add_result2[2]
	current_carry = add_result1[2]
	aux3 = 0
	carry_result = bitToBitTotalCarry(partial_carry,current_carry,aux3)
	total_carry = carry_result[2]

	return (total_sum, total_carry)	

# input:
# - x = binary number represented as array of 0s and 1s
# - y = binary number represented as array of 0s and 1s
# output:
# - z = binary number representing x+y expressed as an array (not list) of 0s and 1s
# note:
# - x and y must have same length
def binaryAdd(binary_number_1,binary_number_2):
	# Just for printing purposes.
	binary_number_1_reversed = list(binary_number_1)
	binary_number_1_reversed.reverse()
	binary_number_2_reversed = list(binary_number_2)
	binary_number_2_reversed.reverse()	
	# print 'binaryAdd >>',binary_number_1_reversed,'+',binary_number_2_reversed,'=',

	dim = len(binary_number_1)
	previous_carry = 0
	z = empty(dim+1,dtype=int)
	for i in xrange(0,dim):
		result = bitToBitFullAdder(binary_number_1[i],binary_number_2[i],previous_carry)
		
		# the resulting sum of ith bit of x with ith bit of y
		z[i] = result[0]

		# the resulting carry from bit-to-bit addition
		previous_carry = result[1]

	# the overflow carry of addition of binary numbers x and y
	z[dim] = previous_carry
	z_reversed = list(z)
	# print z_reversed
	return z.tolist()

# compute for the 2's complement of a binary number represented as an array of 0s and 1s
# input:
# 	- x = binary number to convert
# output:
#	- 2's complement of x represented as an array of 0s and 1s
def twosComplement(binary_number):
	binary_number_reversed = list(binary_number)
	binary_number_reversed.reverse()
	# print 'twosComplement >> 2\'s complement of binary number',binary_number_reversed,':',
	dim = len(binary_number)
	# 1's complement of each bit of x
	for i in xrange(0,dim):
		binary_number[i] = NOT(binary_number[i])

	# 2's complement of each bit of x
	y = zeros(dim,dtype=int)
	y[0] = 1
	binary_sum = binaryAdd(binary_number,y)
	binary_sum.pop()
	binary_number_reversed = list(binary_sum)
	binary_number_reversed.reverse()
	# print binary_number_reversed

	return binary_sum

##############################################

def Euclid_Norm(cplex):
	return sqrt((cplex.conjugate()*cplex).real)

def Euclid_Norm_Squared(cplex):
	return (cplex.conjugate()*cplex).real

############################################## CLASSES ##############################################		

class Hadamard:
	def __init__(self, dim):
		H = array([[1,1],[1,-1]])
		self.matrix = H
		for i in xrange(0,dim-1):
			self.matrix = kron(H,self.matrix)
		self.norm_factor = pow(1/sqrt(2),dim)

