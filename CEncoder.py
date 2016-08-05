# -------------
#  CEncoder
# -------------

from numpy import zeros, array

def F(x):
	if x == 'a':
		return 1
	else:
		return 0

def G(x):
	if x == 'b':
		return 1
	else:
		return 0		

def cencode_BIV_a(a,b):
	
	dim_T = len(a)
	dim_P = len(b)
	BIV_encoded_T_a = zeros(dim_T)
	BIV_encoded_P_a = zeros(dim_P)

	for i in xrange(0,dim_T):
		BIV_encoded_T_a[i] = F(a[i])

	for i in xrange(0,dim_P):
		BIV_encoded_P_a[i] = F(b[i])

	print 'BIV_encoded_T_a: ' + str(BIV_encoded_T_a)
	print 'BIV_encoded_P_a: ' + str(BIV_encoded_P_a)
	return array([BIV_encoded_T_a.tolist(), BIV_encoded_P_a.tolist()])

def cencode_BIV_b(a,b):
	
	dim_T = len(a)
	dim_P = len(b)
	BIV_encoded_T_b = zeros(dim_T)
	BIV_encoded_P_b = zeros(dim_P)

	for i in xrange(0,dim_T):
		BIV_encoded_T_b[i] = G(a[i])

	for i in xrange(0,dim_P):
		BIV_encoded_P_b[i] = G(b[i])

	print 'BIV_encoded_T_b: ' + str(BIV_encoded_T_b)
	print 'BIV_encoded_P_b: ' + str(BIV_encoded_P_b)
	return array([BIV_encoded_T_b.tolist(), BIV_encoded_P_b.tolist()])
