# QEncoder sub-routine
from numpy import kron, dot, identity, array
from math import ceil, log,sqrt
from QOperators import ket0, ket1, H2x2, iqft, qft
from QDeleter import qdelete_matrix, qdelete_arbitrary

#============ methods ================
def F(T,i):
	if T[i] == 'a':
		return -1
	elif T[i] == 'b':
		return 1

def G(P,i):
	if P[i] == 'a':
		return -1
	elif P[i] == 'b':
		return 1

def U_f(N,M,T):
	dim = N+M-1
	u = identity(dim)
	for i in xrange(0, N):
		u[i,i] = F(T,i)
	return u

def U_g(N,M,P):
	dim = N+M-1
	u = identity(dim)
	for i in xrange(0, M):
		u[i,i] = G(P,i)
	return u

def construct_c_tau(a, dim_x):
	a_len = len(a)
	c_array = []
	tau_array = []
	for i in xrange(0,a_len):
		if a[i] == 1:
			c_array.append(i)
		else:
			tau_array.append(i)
	for i in xrange(a_len, a_len + dim_x - 1):
		tau_array.append(i)

	print 'construct_c_tau: c_array = ' + str(c_array)
	print 'construct_c_tau: tau_array = ' + str(tau_array)

	return array([array(c_array), array(tau_array)])

def qencode(T,P):
	P.reverse()
	dim_T = len(T)
	dim_P = len(P)
	dim = dim_T + dim_P - 1
	log_dim = int(ceil(log(dim,2)))
	print 'log(dim) = ' + str(log(dim))
	print 'int(ceil(log(dim))): log('+ str(dim) +') = ' + str(log_dim)

	# initial state of QRegT and QRegP, log(dim)-dimensional |0>
	qRegT_state = ket0()
	qRegP_state = ket0()
	for i in xrange(0,log_dim - 1):
		qRegT_state = kron(qRegT_state,ket0())
		qRegP_state = kron(qRegP_state,ket0())

	# print '===================================='
	# print 'QRegT\'s initial state: ' + str(qRegT_state)
	# print 'QRegP\'s initial state: ' + str(qRegP_state)
	# print '===================================='

	# preparation of log(dim)-dimensional superposition state for
	# QRegT and QRegP
	# H = H2x2()
	# for i in xrange(0,log_dim - 1):
	# 	H = kron(H, H2x2())

	# print 'H: '
	# print H

	# qRegT_state = dot(H,qRegT_state)
	# qRegP_state = dot(H,qRegP_state)

	# create an (N+M-1)-dimension uniform quantum superposition state
	qRegT_state = qft(qRegT_state)
	qRegP_state = qft(qRegP_state)
	# print '===================================='
	# print 'QRegT in uniform superposition state: ' + str(qRegT_state)
	# print 'QRegP in uniform superposition state: ' + str(qRegP_state)
	# print '===================================='

	# deletion of last dim_P-1 basis state of QRegT
	# qRegT_state = qdelete_matrix(dim_T,dim_P,qRegT_state,array([[0,dim_T]]))
	
	# c_array = array([[0,dim_T]])
	# tau_array = array([[dim_T,dim]])
	# qRegT_state = qdelete_arbitrary(dim, qRegT_state, c_array, tau_array)

	c_tau = construct_c_tau(T, dim_P)
	c_array = c_tau[0]
	tau_array = c_tau[1]
	qRegT_state = qdelete_arbitrary(dim, qRegT_state, c_array, tau_array)	

	# print '===================================='
	# print 'QRegT with deleted padding basis states: ' + str(qRegT_state)
	# deletion of last dim_T-1 basis state of QRegP
	# qRegP_state = qdelete_matrix(dim_P,dim_T,qRegP_state,array([[0,dim_P]]))
	
	# c_array = array([[0,dim_P]])
	# tau_array = array([[dim_P,dim]])
	# qRegP_state = qdelete_arbitrary(dim, qRegP_state, c_array, tau_array)

	c_tau = construct_c_tau(P, dim_T)
	c_array = c_tau[0]
	tau_array = c_tau[1]
	qRegP_state = qdelete_arbitrary(dim, qRegP_state, c_array, tau_array)		

	# print 'QRegP with deleted padding basis states: ' + str(qRegP_state)
	# print '===================================='

	# encoding of T into state of QRegT
	# U = U_f(dim_T,dim_P, T)
	# qRegT_state = dot(U,qRegT_state)
	# print '===================================='
	# print 'Encoding of T into state of QRegT: ' + str(qRegT_state)
	# encoding of P into state of QRegP

	# U = U_g(dim_T,dim_P, P)
	# qRegP_state = dot(U,qRegP_state)
	# print 'Encoding of P into state of QRegP: ' + str(qRegP_state)
	# print '===================================='

	# qRegT_state = array([0,0,sqrt(1.0/3),sqrt(1.0/3),sqrt(1.0/3),0,0,0])
	# qRegP_state = array([0,sqrt(1.0/2),sqrt(1.0/2),0,0,0,0,0])
	return array([qRegT_state,qRegP_state])

# ===================================
# T = ['b','b','a']
# P = ['b','a']

# qencode(T,P)



