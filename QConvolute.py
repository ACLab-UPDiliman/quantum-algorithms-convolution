
from numpy import zeros, ones, empty, array, dot
from math import sqrt
import matplotlib.pyplot as plotter

from QOperators import QFT, IQFT, Euclid_Norm

# T = ['b','b','a']
# P = ['b','a']
# N = len(T)
# M = len(P)
# dim = N + M - 1
# P.reverse()
# encoding_state_T = (1.0/sqrt(N)) * array([1.0, 1.0, -1.0, 0.0])
# encoding_state_P = (1.0/sqrt(M)) * array([-1.0, 1.0, 0.0, 0.0])



def constructV(a):
	dim = a.shape[0]
	b = zeros(shape=(dim,dim),dtype=complex)
	for i in xrange(0,dim):
		if a[i] != 0.0+0.0j:
			norm = Euclid_Norm(a[i])
			b[i][i] = complex(a[i].real/norm, a[i].imag/norm)
			# b[i][i] = complex(a[i].real, a[i].imag)
		else:
		 	# b[i][i] = 1.0
		 	b[i][i] = complex(0.0, 1.0)
		 	# b[i][i] = 0.0
	return b

def constructV_Inv(a):
	dim = a.shape[0]
	b = zeros(shape=(dim,dim),dtype=complex)
	for i in xrange(0,dim):
		if a[i][i] != 0.0+0.0j:
			b[i][i] = a[i][i].conjugate()
	return b

def getProbabilities(a):
	dim = a.shape[0]
	b = empty(dim, dtype=float)
	for i in xrange(0,dim):
		b[i] = a[i].real**2 + a[i].imag**2

	return b

def getTotalProbabilities(a):
	dim = a.shape[0]
	total = 0.0
	for i in xrange(0,dim):
		total = total + a[i]

	return total

def getDFTFromQFT(a,strln):
	dim = a.shape[0]
	factor = 1.0/(sqrt(1.0/dim)*sqrt(1.0/strln))
	print 'DFT factor = ' + str(factor)
	b = empty(dim, dtype=complex)
	for i in xrange(0,dim):
		b[i] = complex(factor * a[i].real, factor * a[i].imag)

	return b

def getDFTT_dot_DFTP_from_VQFTP(v,vqft_p,dft_t,Pstrln):
	dim = vqft_p.shape[0]
	qft_P_factor = (1.0/(sqrt(1.0/dim)*sqrt(1.0/Pstrln)))
	b = empty(dim, dtype=complex)
	for i in xrange(0,dim):
		factor = ((1.0/v[i][i])*dft_t[i])*qft_P_factor
		b[i] = complex(factor*vqft_p[i].real, factor*vqft_p[i].imag)

	return b

def getIDFT_from_IQFT(v,iqft_v,dft_t,dft_p,Pstrln):
	dim = iqft_v.shape[0]
	qft_P_factor = 1.0/(sqrt(1.0/dim)*sqrt(1.0/Pstrln))
	iqft_factor = dim/(sqrt(1.0/dim))
	print 'qft_P_factor = ' + str(qft_P_factor)
	print 'iqft_factor = ' + str(iqft_factor)
	b = empty(dim, dtype=complex)
	for j in xrange(0,dim):
		factor = 0.0
		for i in xrange(0,dim):
			factor = factor + iqft_factor*(1.0/v[i][i])*dft_t[i]*qft_P_factor
		b[j] = complex(factor*iqft_v[j].real, factor*iqft_v[j].imag)

	return b

def prettyprint(a):
	dim = a.shape[0]
	string = ''
	for i in xrange(0,dim):
		string = string + '|' + str(i) + '>: ' + str(a[i]) + '\n'

	return string

def qconvolute(encoding_state_T, encoding_state_P):
	dim = encoding_state_T.shape[0]
	# print '================================================='
	# print 'encoding_state_T = ' + str(encoding_state_T)
	# print 'encoding_state_P = ' + str(encoding_state_P)
	# print 'encoding state\'s dimension = ' + str(dim)
	# print 'N = ' + str(N)
	# print 'M = ' + str(M)
	# print '=================================================\n'

	#====== Text =======
	# print '================================================='
	qft_T = QFT(encoding_state_T)
	print 'Amplitudes of QFT for T'
	print prettyprint(qft_T)
	# probabilities = getProbabilities(qft_T)
	# totalProbability = getTotalProbabilities(probabilities)
	# print 'Probabilities after QFT T = ' + str(probabilities)
	# print 'Total probabilities after QFT T = ' + str(totalProbability) + '\n'

	## dft_T = getDFTFromQFT(qft_T,N)
	## print 'DFT(T) derived from QFT_T = ' + str(dft_T)
	# print '=================================================\n'

	#====== Pattern ========
	# print '================================================='
	qft_P = QFT(encoding_state_P)
	print 'Amplitudes of QFT for P'
	print prettyprint(qft_P)
	# probabilities = getProbabilities(qft_P)
	# print 'Probabilities after QFT P = ' + str(probabilities)
	# print 'Total probabilities after QFT P = ' + str(totalProbability) + '\n'

	## dft_P = getDFTFromQFT(qft_P,M)
	## print 'DFT(P) derived from QFT_P = ' + str(dft_P)
	# print '=================================================\n'

	#====== V(QFT(P)) ========
	# V = constructV(qft_T)
	V = constructV(qft_P)
	print 'V: '
	print V

	# V_QFT_P = dot(V, qft_P)
	V_QFT_T = dot(V, qft_T)
	# print 'Amplitudes of V on QFT of P'
	print 'Amplitudes of V on QFT of T'
	# print prettyprint(V_QFT_P)
	print prettyprint(V_QFT_T)
	# probabilities = getProbabilities(V_QFT_P)
	# totalProbability = getTotalProbabilities(probabilities)
	# print 'Probabilities after V = ' + str(probabilities)
	# print 'Total probabilities after V = ' + str(totalProbability) + '\n'

	## dft_T_times_dft_P = getDFTT_dot_DFTP_from_VQFTP(V, V_QFT_P, dft_T, M)
	## print 'DFT(T)*DFT(P) = ' + str(dft_T_times_dft_P)

	#====== IQFT =======
	# IQFT_V_QFT_P = IQFT(V_QFT_P)
	IQFT_V_QFT_T = IQFT(V_QFT_T)
	print 'Amplitudes of QFT-Inv (convolution)'
	# print prettyprint(IQFT_V_QFT_P) + '\n'
	print prettyprint(IQFT_V_QFT_T) + '\n'

	# probabilities = getProbabilities(IQFT_V_QFT_P)
	probabilities = getProbabilities(IQFT_V_QFT_T)
	totalProbability = getTotalProbabilities(probabilities)
	print 'Convolution Probabilities Per Basis State'
	print prettyprint(probabilities)
	print 'Total probabilities after IQFT = ' + str(totalProbability) + '\n'
	print '================================================='

	plotter.title('Probability of Occurence for Convolute')
	plotter.ylabel('Probability')
	plotter.xlabel('Basis state index')
	x = range(dim)
	width = 0.15
	plotter.bar(x, probabilities, width, color="black", alpha=0.75)
	# # plotter.axis([-0.5, 16, 0.0, 1.0])
	plotter.grid(True)
	plotter.show()

	# return IQFT_V_QFT_P
	return IQFT_V_QFT_T

# encoding_state_T = array([ 0.20203050891-0.3194382825j,
# 							0.20203050891-0.3194382825j,
# 							0.20203050891-0.3194382825j,
# 							0.20203050891-0.3194382825j,
# 							0.20203050891-0.3194382825j,
# 							0.20203050891-0.3194382825j,
# 							-0.20203050891+0.3194382825j,
# 							-8.32667268469e-17+2.42861286637e-17j])
# encoding_state_P = array([ 0.707106781187-4.32978028118e-17j,
# 							-0.707106781187+4.32978028118e-17j,
# 							-5.55111512313e-17+0j,
# 							-5.55111512313e-17+0j,
# 							-5.55111512313e-17+0j,
# 							-5.55111512313e-17+0j,
# 							-5.55111512313e-17+0j,
# 							-5.55111512313e-17+0j ])

# qconvolute(encoding_state_T, encoding_state_P)

