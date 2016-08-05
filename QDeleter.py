from numpy import array, dot, kron, zeros, ones, identity, outer, ceil, floor
from scipy.fftpack import fft, ifft
from math import sqrt, pi, cos, sin, expm1, asin, acos
import matplotlib.pyplot as plotter

from QOperators import H2x2_kron_log_dim, Euclid_Norm_Squared

#Operator S from Yang and Xiaoping's quantum basis deletion algorithm for uniform superposition state input
def construct_S(N):
	S = array(
		[
			[ (N-2.0)/(2.0*N) - 1j*(sqrt(3.0)/2.0), (sqrt(N-1.0)/(2.0*N)) - 1j*(sqrt(3.0*(N-1.0))/(2.0*N)) ],
			[(sqrt(N-1.0))/N, (((-2.0)*N + 1.0)/(2.0*N)) - 1j*((sqrt(3.0))/(2.0*N)) ]
		])
	return S

#Operator S from Yang and Xiaoping's quantum basis deletion algorithm for non-uniform superposition state input
def construct_S_arbitrary(beta,phi):
	cos_beta = cos(beta)
	cos_beta_squared = cos_beta**2
	sin_beta = sin(beta)
	phase_factor = cos(phi) + 1j*sin(phi)
	a_00 = -1.0*phase_factor*(1.0 + (phase_factor-1.0)*cos_beta_squared)
	a_01 = -1.0*(phase_factor-1.0)*sin_beta*cos_beta
	a_10 = -1.0*phase_factor*(phase_factor-1.0)*sin_beta*cos_beta
	a_11 = -1.0*phase_factor + (phase_factor-1.0)*cos_beta_squared
	S = array([
			[a_00, a_01],
			[a_10, a_11]
		])
	return S

# x - number of basis states to retain
# y - number basis states to delete
# gamma - a list of ranges of basis state indices to retain in the superposition, i.e. ([[0,2],[5,8]]) means 
# 			basis states from |0> to |2> and |5> to |8>,
def construct_S_Matrix(dim,beta,phi,gamma_array):
	print '=============== construct_S_matrix start ================================='
	phase_factor = cos(phi) + 1j*sin(phi) - 1
	# print 'phase factor: ' + str(phase_factor)
	I = identity(dim)

	# I + phase factor * |0><0|
	ket_0_outer_product = zeros((dim,dim))
	ket_0_outer_product[0][0] = 1.0
	I_0 = I + (phase_factor*ket_0_outer_product)
	# print 'I_0 ='
	# print I_0


	# I + phase factor * sum(|i><i|) where i neq to tau
	ket_gamma_outer_product = zeros((dim,dim))
	gamma_dim = gamma_array.shape[0] # number of ranges in the list gamma
	for i in xrange(0,gamma_dim): # for each range indexed i
		for j in xrange(gamma_array[i][0],gamma_array[i][1]): # from starting index gamma[i][0] to ending index gamma[i][1]-1 of range i
			ket_gamma_outer_product[j][j] = 1.0
	I_gamma = I + (phase_factor*ket_gamma_outer_product)
	# print '|c><c| ='
	# print ket_gamma_outer_product
	# print 'I_gamma ='
	# print I_gamma 

	U = H2x2_kron_log_dim(dim)
	U_inv = U
	# print 'U = '
	# print U
	
	# S = -U_inv dot I_0 dot U dot I_gamma
	S = dot(U,I_gamma)
	S = dot(I_0,S)
	S = dot(U_inv,S)
	S = (-1)*S
	# print 'S = '
	# print S
	print '=============== construct_S_matrix end ================================='

	return S

# J is the number of iterations required in Yang and Xiaoping's quantum deletion algorithm
def getJ(beta):
	x = float(pi / (2.0*pi - 4.0*beta))
	y = float(1.0/2.0)
	raw_J = x-y
	ceil_J = ceil(raw_J)
	floor_J = floor(raw_J)
	print 'x: ' + str(x) + ', y: ' + str(y)
	print 'raw J: ' + str(raw_J) + ', floor J: ' + str(floor_J) + ', ceil J: ' + str(ceil_J)
	print str(int(raw_J)) + ' != ' + str(int(floor_J)) + ': ' + str(int(raw_J) != int(floor_J))
	if raw_J > floor_J or int(floor_J) == 0:
		print 'Returning ceil J.'
		return ceil_J
	else:
		print 'Returning floor J.'
		return floor_J

# cos(beta) is the amplitude of |c> in Yang and Xiaoping's quantum deletion algorithm
def getBeta(superposition_state, c_array):
	sum_gamma = 0.0
	string = ''
	for i in xrange(c_array[0][0],c_array[0][1]):
		norm_squared = Euclid_Norm_Squared(superposition_state[i])
		sum_gamma = sum_gamma + norm_squared
		string = string + '|a_' + str(i) + '|^2 = |' + str(superposition_state[i]) + '|^2' + ' = ' + str(norm_squared) + '\n'

	sqrt_of_sum_of_norms = sqrt(sum_gamma)
	print 'Sum of Norms = ' + str(sum_gamma)
	print 'Square root of sum of norms = sqrt(' + str(sum_gamma) + ') = ' + str(sqrt_of_sum_of_norms)
	beta = acos(sqrt_of_sum_of_norms)
	print 'Acos of Square root of sum of norms = ' + str(beta)
	print 'beta = acos(sqrt(\n' + string + '))'
	sin_beta_squared = (sin(beta))**2
	cos_beta_squared = (cos(beta))**2
	print 'sin^2(beta) = ' + str(sin_beta_squared)
	print 'cos^2(beta) = ' + str(cos_beta_squared)
	print 'sin^2(beta) + cos^2(beta) = ' + str(sin_beta_squared) + ' + ' + str(cos_beta_squared) + ' = ' + str(sin_beta_squared + cos_beta_squared)
	return beta

# cos(beta) is the amplitude of |c> in Yang and Xiaoping's quantum deletion algorithm
def getBeta_scattered_indices(superposition_state, c_array):
	sum_gamma = 0.0
	string = ''
	c_array_len = len(c_array)
	for i in xrange(0, c_array_len):
		norm_squared = (abs(superposition_state[c_array[i]]))**2
		sum_gamma = sum_gamma + norm_squared
		string = string + '|a_' + str(i) + '|^2 = |' + str(superposition_state[c_array[i]]) + '|^2' + ' = ' + str(norm_squared) + '\n'

	sqrt_of_sum_of_norms = sqrt(sum_gamma)
	print 'Sum of Norms = ' + str(sum_gamma)
	print 'Square root of sum of norms = sqrt(' + str(sum_gamma) + ') = ' + str(sqrt_of_sum_of_norms)
	beta = acos(sqrt_of_sum_of_norms)
	print 'Acos of Square root of sum of norms = ' + str(beta)
	print 'beta = acos(sqrt(\n' + string + '))'
	sin_beta_squared = (sin(beta))**2
	cos_beta_squared = (cos(beta))**2
	print 'sin^2(beta) = ' + str(sin_beta_squared)
	print 'cos^2(beta) = ' + str(cos_beta_squared)
	print 'sin^2(beta) + cos^2(beta) = ' + str(sin_beta_squared) + ' + ' + str(cos_beta_squared) + ' = ' + str(sin_beta_squared + cos_beta_squared)
	return beta

def getPhi(beta, J):
	x = sin(pi/((4.0*J)+2.0))
	y = cos(beta)
	print 'getPhi: x = ' + str(x)
	print 'getPhi: y = ' + str(y)
	print 'getPhi: x/y = ' + str(x/y)
	phi = 2.0*asin(x/y)

	return phi

def prettyprint(a):
	dim = a.shape[0]
	string = ''
	for i in xrange(0,dim):
		string = string + '|' + str(i) + '>: ' + str(a[i]) + '\n'

	return string

#Yang and Xiaoping's quantum basis deletion algorithm for non-uniform superposition state input
# def qdelete_arbitrary(dim, superposition_state, c_array, tau_array):
# 	print 'dim = ' + str(dim)
# 	print 'superposition_state size = ' + str(superposition_state.shape[0])
# 	probabilities = zeros(dim)
# 	total_probability = 0.0	
# 	for i in xrange(0,dim):
# 		probabilities[i] = superposition_state[i].real**2 + superposition_state[i].imag**2
# 		total_probability = total_probability + probabilities[i]

# 	print '=================================='
# 	print 'Probabilities prior to deletion'
# 	print prettyprint(probabilities)
# 	print 'Total Probability prior to deletion = ' + str(total_probability)
# 	print '=================================='

# 	#beta = getBeta(superposition_state, c_array)
# 	beta = getBeta_scattered_indices(superposition_state, c_array)

# 	# prepare amplitude of |c> which is cos(beta)
# 	cos_beta = cos(beta)

# 	# compute for normalizing factor for |c>
# 	c_norm_sum = 0.0
# 	for i in xrange(c_array[0][0],c_array[0][1]):
# 		c_norm_sum = c_norm_sum + Euclid_Norm_Squared(superposition_state[i])
# 	c_norm = sqrt(1.0/c_norm_sum)

# 	# construct |c> from superposition_state
# 	c_dim = c_array[0][1]-c_array[0][0]
# 	ket_c = zeros(c_dim, dtype=complex)
# 	c_count = 0
# 	for i in xrange(0, c_dim):
# 		ket_c[i] = complex( superposition_state[c_array[0][0] + c_count].real * c_norm, superposition_state[c_array[0][0] + c_count].imag * c_norm)
# 		c_count = c_count + 1

# 	# prepare amplitude of |tau> which is sin(beta)
# 	sin_beta = sin(beta)

# 	# compute for normalizing factor for |tau>
# 	tau_norm_sum = 0.0
# 	for i in xrange(tau_array[0][0], tau_array[0][1]):
# 		tau_norm_sum = tau_norm_sum + Euclid_Norm_Squared(superposition_state[i])
# 	tau_norm = sqrt(1.0/tau_norm_sum)

# 	# count total number of basis in |tau>
# 	tau_array_dim = tau_array.shape[0]
# 	tau_total_dim = 0
# 	for k in xrange(0, tau_array_dim):
# 		tau_total_dim = tau_total_dim + (tau_array[k][1] - tau_array[k][0])

# 	# construct |tau> from superposition_state
# 	ket_tau = zeros(tau_total_dim, dtype=complex)
# 	ket_count = 0
# 	for i in xrange(0, tau_array_dim):
# 		curr_tau_dim = tau_array[i][1] - tau_array[i][0] # size of current range
# 		curr_tau_count = 0
# 		for j in xrange(0, curr_tau_dim):
# 			ket_tau[ket_count] =  complex(superposition_state[tau_array[i][0] + curr_tau_count].real * tau_norm, superposition_state[tau_array[i][0] + curr_tau_count].imag * tau_norm)
# 			curr_tau_count = curr_tau_count + 1
# 			ket_count = ket_count + 1
	
# 	# prepare input superposition state in terms of basis states |c> and |tau>
# 	ket_gamma = array([cos_beta, sin_beta])
# 	ket_gamma_probabilities = array([ Euclid_Norm_Squared(ket_gamma[0]), Euclid_Norm_Squared(ket_gamma[1])])
# 	# plotter.title('Probability of Occurence for Each Basis State in |gamma>')
# 	# plotter.ylabel('Probability')
# 	# plotter.xlabel('Basis state index')
# 	# x = range(2)
# 	# width = 1/1.5
# 	# plotter.bar(x, ket_gamma_probabilities, width, color="blue", alpha=0.35)
# 	# plotter.grid(True)
# 	# plotter.show()

# 	# determine number of required iterations J
# 	J = getJ(beta)

# 	# determine phase factor phi
# 	phi = getPhi(beta,J)

# 	# construct unitary deletion operator S which operates on superposition 
# 	# |gamma> = cos(beta)|c> + sin(beta)|tau>
# 	S = construct_S_arbitrary(beta, phi)

# 	# Perfrom deletion of basis states with indices in tau_array form superposition_state.
# 	# After application of S|gamma> = |gamma^\prime>, the amplitudes of |c> and |tau> are
# 	# changed into new values.
# 	for i in xrange(0, int(J)):
# 		ket_gamma = dot(S,ket_gamma)
# 	print '|gamma^prime> = ' + str(ket_gamma[0]) + '|c> + ' + str(ket_gamma[1]) + '|tau>'
# 	print 'Probabilities of basis states of |gamma^prime>: |c> = ' + str(Euclid_Norm_Squared(ket_gamma[0])) + ', |tau> = ' + str(Euclid_Norm_Squared(ket_gamma[1]))
	

# 	ket_gamma_prime_probabilities = array([ Euclid_Norm_Squared(ket_gamma[0]), Euclid_Norm_Squared(ket_gamma[1])])
# 	# plotter.title('Probability of Occurence for Each Basis State in |gamma^prime>')
# 	# plotter.ylabel('Probability')
# 	# plotter.xlabel('Basis state index')
# 	# x = range(2)
# 	# width = 1/1.5
# 	# plotter.bar(x, ket_gamma_prime_probabilities, width, color="blue", alpha=0.35)
# 	# plotter.grid(True)
# 	# plotter.show()

# 	# construct new superposition_state from |gamma^prime>
# 	# copy amplitudes of |c> into appropriate basis states of superposition_state_prime
# 	superposition_state_prime = zeros(dim, dtype=complex)
# 	count = 0
# 	for i in xrange(c_array[0][0],c_array[0][1]):
# 		print '|gamma^prime>[0] x |c>[' + str(count) +'] = ' + str(Euclid_Norm_Squared(ket_gamma[0])) + ' x ' + str(ket_c[count]) + ' = ' + str(Euclid_Norm_Squared(ket_gamma[0])*ket_c[count])
# 		print 'superposition_state_prime['+ str(i) +'] = ' + str(Euclid_Norm_Squared(ket_gamma[0]) * ket_c[count])
# 		superposition_state_prime[i] =  complex(ket_c[count].real * abs(ket_gamma[0]), ket_c[count].imag * abs(ket_gamma[0]),)
# 		count = count + 1

# 	# copy amplitudes of |tau> into appropriate basis states of superposition_state_prime
# 	count = 0
# 	for i in xrange(0, tau_array_dim):
# 		for j in xrange(tau_array[i][0], tau_array[i][1]):
# 			print '|gamma^prime>[1] x |tau>[' + str(count) +'] = ' + str(Euclid_Norm_Squared(ket_gamma[1])) + ' x ' + str(ket_tau[count]) + ' = ' + str(Euclid_Norm_Squared(ket_gamma[1]) * ket_tau[count])
# 			print 'superposition_state_prime['+ str(j) +'] = ' + str(Euclid_Norm_Squared(ket_gamma[1]) * ket_tau[count])
# 			superposition_state_prime[j] = Euclid_Norm_Squared(ket_gamma[1]) * ket_tau[count]
# 			count = count + 1


# 	probabilities = zeros(dim)
# 	total_probability = 0.0
# 	for i in xrange(0,dim):
# 		probabilities[i] = superposition_state_prime[i].real**2 + superposition_state_prime[i].imag**2
# 		total_probability = total_probability + probabilities[i]

# 	print '=================================='
# 	print 'Probabilities after deletion'
# 	print prettyprint(probabilities)
# 	print 'Total Probability after deletion = ' + str(total_probability)
# 	print '=================================='

# 	# plotter.title('Probability of Occurence for Each Basis State')
# 	# plotter.ylabel('Probability')
# 	# plotter.xlabel('Basis state index')
# 	# x = range(dim)
# 	# width = 1/1.5
# 	# plotter.bar(x, probabilities, width, color="blue", alpha=0.35)
# 	# plotter.grid(True)
# 	# plotter.show()

# 	return superposition_state_prime


# c_array is a list of indices of states to amplify
# tau_array is a list of indices of states to delete
def qdelete_arbitrary(dim, superposition_state, c_array, tau_array):
	print 'dim = ' + str(dim)
	print 'superposition_state size = ' + str(superposition_state.shape[0])
	probabilities = zeros(dim)
	total_probability = 0.0	
	for i in xrange(0,dim):
		probabilities[i] = (abs(superposition_state[i]))**2
		total_probability = total_probability + probabilities[i]

	print '=================================='
	print 'Probabilities prior to deletion'
	print prettyprint(probabilities)
	print 'Total Probability prior to deletion = ' + str(total_probability)
	print '=================================='

	#beta = getBeta(superposition_state, c_array)
	beta = getBeta_scattered_indices(superposition_state, c_array)

	# prepare amplitude of |c> which is cos(beta)
	cos_beta = cos(beta)

	# compute for normalizing factor for |c>
	c_norm_sum = 0.0
	c_array_len = len(c_array)
	for i in xrange(0, c_array_len):
		c_norm_sum = c_norm_sum + (abs(superposition_state[c_array[i]]))**2
	c_norm = sqrt(1.0/c_norm_sum)

	# construct |c> from superposition_state
	ket_c = zeros(c_array_len, dtype=complex)
	for i in xrange(0, c_array_len):
		ket_c[i] = complex( superposition_state[c_array[i]].real * c_norm, superposition_state[c_array[i]].imag * c_norm)

	# prepare amplitude of |tau> which is sin(beta)
	sin_beta = sin(beta)

	# compute for normalizing factor for |tau>
	tau_norm_sum = 0.0
	tau_array_len = len(tau_array)
	for i in xrange(0, tau_array_len):
		tau_norm_sum = tau_norm_sum + (abs(superposition_state[tau_array[i]]))**2
	tau_norm = sqrt(1.0/tau_norm_sum)

	# construct |tau> from superposition_state
	ket_tau = zeros(tau_array_len, dtype=complex)
	for i in xrange(0, tau_array_len):
		ket_tau[i] =  complex(superposition_state[tau_array[i]].real * tau_norm, superposition_state[tau_array[i]].imag * tau_norm)
			
	# prepare input superposition state in terms of basis states |c> and |tau>
	ket_gamma = array([cos_beta, sin_beta])
	ket_gamma_probabilities = array([ (abs(ket_gamma[0]))**2, (abs(ket_gamma[1]))**2 ])

	# determine number of required iterations J
	J = getJ(beta)

	# determine phase factor phi
	phi = getPhi(beta,J)

	# construct unitary deletion operator S which operates on superposition 
	# |gamma> = cos(beta)|c> + sin(beta)|tau>
	S = construct_S_arbitrary(beta, phi)

	# Perfrom deletion of basis states with indices in tau_array form superposition_state.
	# After application of S|gamma> = |gamma^\prime>, the amplitudes of |c> and |tau> are
	# changed into new values.
	for i in xrange(0, int(J)):
		ket_gamma = dot(S,ket_gamma)
	print '|gamma^prime> = ' + str(ket_gamma[0]) + '|c> + ' + str(ket_gamma[1]) + '|tau>'
	print 'Probabilities of basis states of |gamma^prime>: |c> = ' + str((abs(ket_gamma[0]))**2) + ', |tau> = ' + str((abs(ket_gamma[1]))**2)
	
	ket_gamma_prime_probabilities = array([ (abs(ket_gamma[0]))**2, (abs(ket_gamma[1]))**2])
	
	# construct new superposition_state from |gamma^prime>
	# copy amplitudes of |c> into appropriate basis states of superposition_state_prime
	superposition_state_prime = zeros(dim, dtype=complex)
	for i in xrange(0, c_array_len):
		superposition_state_prime[c_array[i]] =  complex(ket_c[i].real * ket_gamma[0], ket_c[i].imag * ket_gamma[0])
		print 'superposition_state_prime['+ str(c_array[i]) +'] = ' + str(superposition_state_prime[c_array[i]])

	# copy amplitudes of |tau> into appropriate basis states of superposition_state_prime
	for j in xrange(0, tau_array_len):
		superposition_state_prime[tau_array[j]] = complex(ket_tau[j].real * ket_gamma[1], ket_tau[j].imag * ket_gamma[1]) 
		print 'superposition_state_prime['+ str(tau_array[j]) +'] = ' + str(superposition_state_prime[tau_array[j]])


	probabilities = zeros(dim)
	total_probability = 0.0
	for i in xrange(0,dim):
		probabilities[i] = abs(superposition_state_prime[i])**2
		total_probability = total_probability + probabilities[i]

	print '=================================='
	print 'Probabilities after deletion'
	print prettyprint(probabilities)
	print 'Total Probability after deletion = ' + str(total_probability)
	print '=================================='

	# plotter.title('Probability of Occurence for Each Basis State')
	# plotter.ylabel('Probability')
	# plotter.xlabel('Basis state index')
	# x = range(dim)
	# width = 1/1.5
	# plotter.bar(x, probabilities, width, color="blue", alpha=0.35)
	# plotter.grid(True)
	# plotter.show()

	return superposition_state_prime


# def qdelete_matrix(N,M):
def qdelete_matrix(N,M,superposition_state,gamma_array):
	print '=============== qdelete_matrix start ================================='
	dim = N+M-1
	range_size = gamma_array[0][1] - gamma_array[0][0]
	print 'Basis states to retain: |' + str(gamma_array[0][0]) + '> to |' + str(gamma_array[0][1]-1) + '>'

	beta = getBeta(superposition_state, gamma_array)
	J = getJ(beta)
	phi = getPhi(beta, J)
	print 'beta : ' + str(beta) + ', J : ' + str(J) + ', phi : ' + str(phi)
	
	S = construct_S_Matrix(dim,beta,phi,gamma_array)
	initial_probabilities = zeros(dim)
	for i in xrange(0,dim):
		initial_probabilities[i] = superposition_state[i].real**2 + superposition_state[i].imag**2

	print '=================================='
	print 'Probabilities prior to deletion'
	print prettyprint(initial_probabilities)
	print '=================================='

	S_dot_state = dot(S,superposition_state)

	probabilities = zeros(dim)
	total_probability = 0.0
	for i in xrange(0,dim):
		probabilities[i] = S_dot_state[i].real**2 + S_dot_state[i].imag**2  #(S_dot_state[i].conjugate())*S_dot_state[i]
		total_probability = total_probability + probabilities[i]

	print '=================================='
	print 'Probabilities after deletion'
	print prettyprint(probabilities)
	print 'Total Probabilities after deletion = ' + str(total_probability)
	print '=================================='

	# plotter.plot(probabilities)
	# plotter.title('Initial Probability of Occurence of Basis States')
	# plotter.ylabel('Probability')
	# plotter.xlabel('Basis state index')
	# x = range(dim)
	# width = 1/1.5
	# plotter.bar(x, initial_probabilities, width, color="blue", alpha=0.35)
	# plotter.axis([-0.5, 16, 0.0, 1.0])
	# plotter.grid(True)
	# plotter.show()

	plotter.title('Probability of Occurence for Each Basis State')
	plotter.ylabel('Probability')
	plotter.xlabel('Basis state index')
	x = range(dim)
	width = 1/1.5
	plotter.bar(x, probabilities, width, color="blue", alpha=0.35)
	# # plotter.axis([-0.5, 16, 0.0, 1.0])
	plotter.grid(True)
	plotter.show()
	print '=============== qdelete_matrix end ================================='

	return S_dot_state

# param1 is total number of basis states
# param2 is number of basis states to retain
# param3 is number of basis states to delete
def qdelete(N, c, t):
	phi = pi/3 
	general_phase = cos((pi - phi)/2.0) + 1j*sin((pi - phi)/2.0)
	U = construct_S(N)
	ket_psi_0 = array([ sqrt(c/N), sqrt(t/N) ])
	ket_psi_1 = dot( U, ket_psi_0 )
	phase_corrected = array([ ket_psi_1[0]*general_phase, ket_psi_1[1]*general_phase ])
	amplitude_basis_states_c = phase_corrected[0]/sqrt(c)
	total_probability_basis_states_c = (amplitude_basis_states_c.conjugate()*amplitude_basis_states_c)*c
	amplitude_basis_states_t = phase_corrected[1]/sqrt(t)
	total_probability_basis_states_t = (amplitude_basis_states_t.conjugate()*amplitude_basis_states_t)*t
	probabilities = array([ (phase_corrected[0].conjugate())*phase_corrected[0], (phase_corrected[1].conjugate())*phase_corrected[1] ])
	
	# print 'general_phase = ' + str(general_phase)
	# print 'ket_psi_0 = ' + str(ket_psi_0)
	# print 'ket_psi_1 = ' + str(ket_psi_1)
	# print 'phase_corrected = ' + str(phase_corrected)
	# print 'amplitude of each base in c: ' + str(amplitude_basis_states_c)
	# print 'amplitude of each base in t: ' + str(amplitude_basis_states_t)
	# print 'total probability of basis states in c: ' + str(total_probability_basis_states_c)
	# print 'total probability of basis states in t: ' + str(total_probability_basis_states_t)
	# print 'probabilities = ' + str(probabilities)
	# print 'sum of probabilities = ' + str(probabilities[0] + probabilities[1])

# superposition_state = array([ 0.10871478283-0.171893164536j,
# 								-0.0121481439569+0.0192079021238j,
# 								-0.00401718124347+0.00635172125156j,
# 								0.0173950351342-0.0275039655014j,
# 								-0.0441507429156+0.0698084540011j,
# 								0.116614954989-0.184384433501j,
# 								-0.438511331679+0.693347293952j,
# 								0.256102626843-0.404933807789j
# 							])
# c_array = array([[1,7]])
# qdelete_matrix(7,2,superposition_state,c_array)

# superposition_state = array([ 0.10871478283-0.171893164536j,
# 								-0.0121481439569+0.0192079021238j,
# 								-0.00401718124347+0.00635172125156j,
# 								0.0173950351342-0.0275039655014j,
# 								-0.0441507429156+0.0698084540011j,
# 								0.116614954989-0.184384433501j,
# 								-0.438511331679+0.693347293952j,
# 								0.256102626843-0.404933807789j
# 							])
# c_array = array([[1,7]])
# tau_array = array([[0,1],[7,8]])
# qdelete_arbitrary(7+2-1, superposition_state, c_array, tau_array)

# superposition_state = array([ sqrt(1.0/8.0),
# 								sqrt(1.0/8.0),
# 								sqrt(1.0/8.0),
# 								sqrt(1.0/8.0),
# 								sqrt(1.0/8.0),
# 								sqrt(1.0/8.0),
# 								sqrt(1.0/8.0),
# 								sqrt(1.0/8.0)
# 							])

# c_array = array([[1,7]])
# tau_array = array([[0,1],[7,8]])

# qdelete_arbitrary(7+2-1, superposition_state, c_array, tau_array)





