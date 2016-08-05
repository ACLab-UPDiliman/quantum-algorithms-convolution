from Hadamard import H2x2, X
from numpy import array, dot, kron, zeros, ones, identity, outer
from scipy.fftpack import fft, ifft
from math import sqrt

#=== functions

def S_T(i, N, M):
	if i >= 0 and i <= (N-1):
		return -1
	elif i >= N and i <= (N + M - 2):
		return 1

def S_P(i, N, M):
	if i >= 0 and i <= (M-1):
		return -1
	elif i >= M and i <= (N + M - 2):
		return 1

def S_T_Matrix(N,M):
	m = identity(N + M - 1)
	for i in xrange(N, N + M -1):
		m[i,i] = 0
	return m

def S_P_Matrix(N,M):
	m = identity(N + M - 1)
	for i in xrange(M, N + M - 1):
		m[i,i] = 0
	return m
		
def A(N,M):
	I = identity(N + M - 1)
	ket_psi = (1/sqrt(N+M-1))*ones(N+M-1)
	bra_psi = (1/sqrt(N+M-1))*ones(N+M-1)
	outerproduct = outer(ket_psi, bra_psi)
	projector = 2 * outerproduct
	return (projector - I)

def normalize_T(a,N,M):
	norm_factor = sqrt(float(N+M-1)/float(N))
	m = norm_factor * a
	return m

def normalize_P(a,N,M):
	norm_factor = sqrt(float(N+M-1)/float(M))
	m = norm_factor * a
	return m

def F(i):
	if T[i] == 'a':
		return -1
	elif T[i] == 'b':
		return 1

def G(i):
	if P[i] == 'a':
		return -1
	elif P[i] == 'b':
		return 1

#== non-unitary
def U_f():
	u = identity(N + M - 1)
	for i in xrange(0, N):
		u[i,i] = F(i)
	return u

#== non-unitary
def U_g():
	u = identity(N + M - 1)
	for i in xrange(0, M):
		u[i,i] = G(i)
	return u

def normalize_qft(a):
	for i in xrange(0, N + M - 2):
		if a[i] != 0:
			a[i] = a[i]

#--- Hadamard
ket0 = array([1,0])
ket1 = array([0,1])
# H_ket0 = dot(H2x2(),ket0)
# H_ket1 = dot(H2x2(),ket1)
# print 'H ket0 = ' + str(H_ket0)
# print 'H ket1 = ' + str(H_ket1)

# #--- Kronecker Product
# ket0_kron_ket_0 = kron(ket0, ket0)
# print 'Ket0 kron Ket0 = ' + str(ket0_kron_ket_0)

# ket0_kron_ket_1 = kron(ket0, ket1)
# print 'Ket0 kron Ket1 = ' + str(ket0_kron_ket_1)

# H2x2_kron_H2x2 = kron(H2x2(), H2x2())
# print 'H2x2 kron H2x2 = ' + str(H2x2_kron_H2x2)

# H2x2_kron_H2x2_dot_ket0_kron_ket0 = dot(H2x2_kron_H2x2, ket0_kron_ket_0)
# print 'H2x2 kron H2x2 dot (ket0 kron ket0) = ' + str(H2x2_kron_H2x2_dot_ket0_kron_ket0)

T = ['b','b','a']
N = len(T)
P = ['b','a']
P.reverse()
print 'P reversed = ' + str(P)
M = len(P)
print 'Length of T = ' + str(N)
print 'Length of P = ' + str(M)

S_T_Matrix = S_T_Matrix(N,M)
print '4x4 S_T Matrix = '
print str(S_T_Matrix)

S_P_Matrix = S_P_Matrix(N,M)
print '4x4 S_P Matrix = '
print str(S_P_Matrix)

# A_Matrix = A(N,M)
# print 'A = '
# print str(A_Matrix)

U_f_Matrix = U_f()
print 'U_f Matrix = '
print str(U_f_Matrix)

U_g_Matrix = U_g()
print 'U_g Matrix = '
print str(U_g_Matrix)

#===== Text ======
T_qubit1 = ket0
T_qubit2 = ket0
T_ancilla = ket0

H_dot_T_qubit1 = dot(H2x2(),T_qubit1)
H_dot_T_qubit2 = dot(H2x2(),T_qubit2)
H_dot_T_qubit1_kron_H_dot_T_qubit2 = kron(H_dot_T_qubit1, H_dot_T_qubit2)
X_dot_T_ancilla = dot(X(),T_ancilla)
H_dot_X_dot_T_ancilla = dot(H2x2(),X_dot_T_ancilla)
S_T_Matrix_dot_H_dot_T_qubit1_kron_H_dot_T_qubit2 = dot(S_T_Matrix, H_dot_T_qubit1_kron_H_dot_T_qubit2)
normalized_state_T = normalize_T(S_T_Matrix_dot_H_dot_T_qubit1_kron_H_dot_T_qubit2, N, M)
encoding_state_T = dot(U_f_Matrix, normalized_state_T)
# A_dot_S_T_Matrix_dot_H_dot_T_qubit1_kron_H_dot_T_qubit2 = dot(A_Matrix, S_T_Matrix_dot_H_dot_T_qubit1_kron_H_dot_T_qubit2)
print '======== Encoding T into quantum superposition state ======='
print 'H_dot_T_qubit1 = ' + str(H_dot_T_qubit1)
print 'H_dot_T_qubit2 = ' + str(H_dot_T_qubit2)
print 'X_dot_T_ancilla = ' + str(X_dot_T_ancilla)
print 'H_dot_X_dot_T_ancilla = ' + str(H_dot_X_dot_T_ancilla)
print 'H_dot_T_qubit1_kron_H_dot_T_qubit2 = ' + str(H_dot_T_qubit1_kron_H_dot_T_qubit2)
print 'S_T_Matrix_dot_H_dot_T_qubit1_kron_H_dot_T_qubit2 = ' + str(S_T_Matrix_dot_H_dot_T_qubit1_kron_H_dot_T_qubit2)
print 'Normalized State for T = ' + str(normalized_state_T)
print 'Encoding State for T = ' + str(encoding_state_T)


#===== Pattern ======
P_qubit1 = ket0
P_qubit2 = ket0
P_ancilla = ket0

H_dot_P_qubit1 = dot(H2x2(),P_qubit1)
H_dot_P_qubit2 = dot(H2x2(),P_qubit2)
H_dot_P_qubit1_kron_H_dot_P_qubit2 = kron(H_dot_P_qubit1, H_dot_P_qubit2)
X_dot_P_ancilla = dot(X(),P_ancilla)
H_dot_X_dot_P_ancilla = dot(H2x2(),X_dot_P_ancilla)
S_P_Matrix_dot_H_dot_P_qubit1_kron_H_dot_P_qubit2 = dot(S_P_Matrix, H_dot_P_qubit1_kron_H_dot_P_qubit2)
normalized_state_P = normalize_P(S_P_Matrix_dot_H_dot_P_qubit1_kron_H_dot_P_qubit2, N, M)
encoding_state_P = dot(U_g_Matrix, normalized_state_P)

print '======== Encoding P into quantum superposition state ======='
print 'H_dot_P_qubit1 = ' + str(H_dot_P_qubit1)
print 'H_dot_P_qubit2 = ' + str(H_dot_P_qubit2)
print 'X_dot_P_ancilla = ' + str(X_dot_P_ancilla)
print 'H_dot_X_dot_P_ancilla = ' + str(H_dot_X_dot_P_ancilla)
print 'H_dot_P_qubit1_kron_H_dot_P_qubit2 = ' + str(H_dot_P_qubit1_kron_H_dot_P_qubit2)
print 'S_P_Matrix_dot_H_dot_P_qubit1_kron_H_dot_P_qubit2 = ' + str(S_P_Matrix_dot_H_dot_P_qubit1_kron_H_dot_P_qubit2)
print 'Normalized State for P = ' + str(normalized_state_P)
print 'Encoding State for P = ' + str(encoding_state_P)

#================================= Convolution =================================

