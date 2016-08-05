
from numpy import identity, ones, zeros, dot, kron, array

# CONVENTION:
# 1. Maximum allowable length of text is 32.
# 	- implies that the fixed number of bits in any register is 5.
# 	- binary values with less than 5 bits will be padded by 0s as necessary.
registerBitCount = 1
fixedBinaryNumberLength = 5

def ket0():
	return array([1,0])

def ket1():
	return array([0,1])

class QuantumRegister(object):
	"""docstring for QuantumRegister"""
	def __init__(self,bitCount, state):
		# super(QuantumRegister, self).__init__()
		self.bitCount = bitCount
		self.state = state # list of QuantumState
		
	def toString(self):
		return 'Bit count: ' + str(self.bitCount) + ', State: {' + str(self.state.toString()) + '}'

class QuantumState(object):
	def __init__(self,amplitude,binaryNumberFormat,symbolFormat,basisVectorFormat):
		super(QuantumState, self).__init__()
		self.amplitude = amplitude # float
		self.binaryNumberFormat = binaryNumberFormat # list of 0s and 1s
		self.symbolFormat = symbolFormat # int or char
		self.basisVectorFormat = basisVectorFormat # list of 0s and single 1

	def toString(self):
		return 'Amplitude: '+ str(self.amplitude) + ', Binary number format: ' + str(self.binaryNumberFormat) + ', Symbol format: ' + str(self.symbolFormat) + ', Basis vector format: ' + str(self.basisVectorFormat)

	def toStringSimple(self):
		return '['+ str(self.amplitude) + ', ' + str(self.binaryNumberFormat) + ', ' + str(self.symbolFormat) + ', ' + str(self.basisVectorFormat) + ']'

class OneBitAdder(object):
	"""Unitary operator for adding a bit to another bit. It is implemented using a unitary matrix representation of a reversible one-bit full adder circuit."""
	def __init__(self):
		super(OneBitAdder, self).__init__()

		class bcolors:
			ONE = '\033[92m'
			ENDC = '\033[0m'

		# Operator 1 of one-bit full adder circuit, CCNOT gate
		operator_1 = identity(32,dtype=int)
		for i in xrange(0,2):
			for j in xrange(0,2):
				decimal_value_0 = pow(2,4)*i + pow(2,3)*1 + pow(2,2)*1 + pow(2,1)*0 + pow(2,0)*j
				decimal_value_1 = pow(2,4)*i + pow(2,3)*1 + pow(2,2)*1 + pow(2,1)*1 + pow(2,0)*j
				operator_1[decimal_value_0][decimal_value_1] = 1
				operator_1[decimal_value_1][decimal_value_0] = 1
				operator_1[decimal_value_0][decimal_value_0] = 0
				operator_1[decimal_value_1][decimal_value_1] = 0
		
		# Operator 2 of one-bit full adder circuit, CNOT gate
		operator_2 = identity(32,dtype=int)
		for i in xrange(0,2):
			for j in xrange(0,2):
				for k in xrange(0,2):
					decimal_value_0 = pow(2,4)*i + pow(2,3)*1 + pow(2,2)*0 + pow(2,1)*j + pow(2,0)*k
					decimal_value_1 = pow(2,4)*i + pow(2,3)*1 + pow(2,2)*1 + pow(2,1)*j + pow(2,0)*k
					operator_2[decimal_value_0][decimal_value_1] = 1
					operator_2[decimal_value_1][decimal_value_0] = 1
					operator_2[decimal_value_0][decimal_value_0] = 0
					operator_2[decimal_value_1][decimal_value_1] = 0

		# Operator 3 of one-bit full adder circuit, CCNOT gate
		operator_3 = identity(32,dtype=int)
		for i in xrange(0,2):
			for j in xrange(0,2):
					decimal_value_0 = pow(2,4)*1 + pow(2,3)*i + pow(2,2)*1 + pow(2,1)*j + pow(2,0)*0
					decimal_value_1 = pow(2,4)*1 + pow(2,3)*i + pow(2,2)*1 + pow(2,1)*j + pow(2,0)*1
					operator_3[decimal_value_0][decimal_value_1] = 1
					operator_3[decimal_value_1][decimal_value_0] = 1
					operator_3[decimal_value_0][decimal_value_0] = 0
					operator_3[decimal_value_1][decimal_value_1] = 0

		# Operator 4 of one-bit full adder circuit, CNOT gate
		operator_4 = identity(32,dtype=int)
		for i in xrange(0,2):
			for j in xrange(0,2):
				for k in xrange(0,2):
					decimal_value_0 = pow(2,4)*1 + pow(2,3)*i + pow(2,2)*0 + pow(2,1)*j + pow(2,0)*k
					decimal_value_1 = pow(2,4)*1 + pow(2,3)*i + pow(2,2)*1 + pow(2,1)*j + pow(2,0)*k
					operator_4[decimal_value_0][decimal_value_1] = 1
					operator_4[decimal_value_1][decimal_value_0] = 1
					operator_4[decimal_value_0][decimal_value_0] = 0
					operator_4[decimal_value_1][decimal_value_1] = 0

		# Operator 5 of one-bit full adder circuit, CNOT gate
		operator_5 = identity(32,dtype=int)
		for i in xrange(0,2):
			for j in xrange(0,2):
				for k in xrange(0,2):
					decimal_value_0 = pow(2,4)*i + pow(2,3)*j + pow(2,2)*k + pow(2,1)*1 + pow(2,0)*0
					decimal_value_1 = pow(2,4)*i + pow(2,3)*j + pow(2,2)*k + pow(2,1)*1 + pow(2,0)*1
					operator_5[decimal_value_0][decimal_value_1] = 1
					operator_5[decimal_value_1][decimal_value_0] = 1
					operator_5[decimal_value_0][decimal_value_0] = 0
					operator_5[decimal_value_1][decimal_value_1] = 0

		U_adder = dot(operator_2, operator_1)
		U_adder = dot(operator_3, U_adder)
		U_adder = dot(operator_4, U_adder)
		U_adder = dot(operator_5, U_adder)

		self.matrix = U_adder
		# for x in xrange(0,32):
		# 	for i in xrange(0,32):
		# 		if U_adder[x][i] == 1:
		# 			print bcolors.ONE + str(U_adder[x][i]) + bcolors.ENDC,
		# 		else:
		# 			print U_adder[x][i],
		# 	print

class FullAdder(object):
	"""docstring for FullAdder"""
	def __init__(self):
		super(FullAdder, self).__init__()
		self.operator = OneBitAdder().matrix

	def add(self, carryBitValue, bitValue1, bitValue2):
		# convert to basis vectors
		bitValue1BasisVector = ket0()
		if bitValue1 == 1:
			bitValue1BasisVector = ket1()

		bitValue2BasisVector = ket0()
		if bitValue2 == 1:
			bitValue2BasisVector = ket1()

		carryBitValueBasisVector = ket0()
		if carryBitValue == 1:
			carryBitValueBasisVector = ket1()	

		# prepare register state expressed as basis vector
		registerBasisVector = kron(carryBitValueBasisVector,
									kron(bitValue1BasisVector, 
											kron(bitValue2BasisVector,
													kron(ket0(),ket0()))))

		# perform binary addition 
		registerBasisVector = dot(self.operator, registerBasisVector).tolist()
		indexOf1 = registerBasisVector.index(1)
		
		# resulting sum; store the value of the sum into the second input bit parameter, bitValue2
		sumBitValue = 0
		for startIndex in [0,8,16,24]:
			if indexOf1 in range(startIndex, startIndex + 8):
				if indexOf1 in range(startIndex, startIndex + 4):
					sumBitValue = 0
				else:
					sumBitValue = 1

		# resulting carry
		carryBitValue = 0
		if indexOf1 % 2 == 1 :
			carryBitValue = 1
		# print 'newCarryBitValue:', newCarryBitValue
		
		additionResult = AdditionResult(sumBitValue, carryBitValue)

		return additionResult

class SymbolComparator(object):
	"""Unitary operator for comparing two symbols. Takes as input bit binary code for symbols x, y and result bit c.
	Operation puts bit c into 1 if x and y do not match and 0 if they match."""
	def __init__(self):
		class bcolors:
				ONE = '\033[92m'
				ENDC = '\033[0m'
		# super(self).__init__()
		alphabet = Alphabet()
		operatorBitCount = (2 * alphabet.codeBitCount) + 1
		matrixDimension = pow(2,alphabet.codeBitCount) * pow(2,alphabet.codeBitCount) * pow(2,1)
		self.matrix = identity(matrixDimension,dtype=int)
		for i in xrange(0,alphabet.codeBitCount):
			# CNOT gate
			operator = identity(matrixDimension, dtype=int)
			controlIndex = (2 * alphabet.codeBitCount)-i # 2*2-0=4
			targetIndex = alphabet.codeBitCount-i #2-0=2
			
			decimalIndex0 = pow(2,controlIndex)*1 + pow(2,targetIndex)*0
			decimalIndex1 = pow(2,controlIndex)*1 + pow(2,targetIndex)*1

			indexList = (ones(operatorBitCount, dtype=int)).tolist()
			indexList[controlIndex] = 0
			indexList[targetIndex] = 0
			
			# A list of decimal values corresponding to all permutation of bits
			# which are not currently control and target bits
			valuesList = [0]
			for j in xrange(0,operatorBitCount):
				if indexList[j] == 1:
					for k in xrange(0,len(valuesList)):
						valuesList.append(valuesList[k] + pow(2,j)*1)
			
			for j in xrange(0,len(valuesList)):
				operator[decimalIndex0 + valuesList[j]][decimalIndex1 + valuesList[j]] = 1
				operator[decimalIndex1 + valuesList[j]][decimalIndex0 + valuesList[j]] = 1
				operator[decimalIndex0 + valuesList[j]][decimalIndex0 + valuesList[j]] = 0
				operator[decimalIndex1 + valuesList[j]][decimalIndex1 + valuesList[j]] = 0
			self.matrix = dot(operator, self.matrix)

			# NOT gate
			decimalIndex0 = pow(2,targetIndex)*0
			decimalIndex1 = pow(2,targetIndex)*1

			indexList = (ones(operatorBitCount, dtype=int)).tolist()
			indexList[targetIndex] = 0

			# A list of decimal values corresponding to all permutation of bits
			# which are not currently control and target bits
			valuesList = [0]
			for j in xrange(0,operatorBitCount):
				if indexList[j] == 1:
					for k in xrange(0,len(valuesList)):
						valuesList.append(valuesList[k] + pow(2,j)*1)
			
			operator = identity(matrixDimension,dtype=int)
			for j in xrange(0,len(valuesList)):
				operator[decimalIndex0 + valuesList[j]][decimalIndex1 + valuesList[j]] = 1
				operator[decimalIndex1 + valuesList[j]][decimalIndex0 + valuesList[j]] = 1
				operator[decimalIndex0 + valuesList[j]][decimalIndex0 + valuesList[j]] = 0
				operator[decimalIndex1 + valuesList[j]][decimalIndex1 + valuesList[j]] = 0
			self.matrix = dot(operator, self.matrix)


		# MultipleControlNOT operator
		# The control bits will be the bits of the second symbol in which the result of
		# the previous operation was stored.
		indexList = (ones(operatorBitCount, dtype=int)).tolist()
		# target bits
		targetIndex = 0
		decimalIndex0 = pow(2,targetIndex)*0
		decimalIndex1 = pow(2,targetIndex)*1
		indexList[targetIndex] = 0

		for controlIndex in xrange(1,alphabet.codeBitCount+1):
			# control bits
			decimalIndex0 = decimalIndex0 + pow(2,controlIndex)*1
			decimalIndex1 = decimalIndex1 + pow(2,controlIndex)*1
			indexList[controlIndex] = 0
		# A list of decimal values corresponding to all permutation of bits
		# which are not currently control and target bits
		valuesList = [0]
		for j in xrange(0,operatorBitCount):
			if indexList[j] == 1:
				for k in xrange(0,len(valuesList)):
					valuesList.append(valuesList[k] + pow(2,j)*1)
		operator = identity(matrixDimension,dtype=int)
		for j in xrange(0,len(valuesList)):
			operator[decimalIndex0 + valuesList[j]][decimalIndex1 + valuesList[j]] = 1
			operator[decimalIndex1 + valuesList[j]][decimalIndex0 + valuesList[j]] = 1
			operator[decimalIndex0 + valuesList[j]][decimalIndex0 + valuesList[j]] = 0
			operator[decimalIndex1 + valuesList[j]][decimalIndex1 + valuesList[j]] = 0
		self.matrix = dot(operator, self.matrix)			


		# NOT operator
		targetIndex = 0
		decimalIndex0 = pow(2,targetIndex)*0
		decimalIndex1 = pow(2,targetIndex)*1
		indexList = (ones(operatorBitCount, dtype=int)).tolist()
		indexList[targetIndex] = 0
		# A list of decimal values corresponding to all permutation of bits
		# which are not currently control and target bits
		valuesList = [0]
		for j in xrange(0,operatorBitCount):
			if indexList[j] == 1:
				for k in xrange(0,len(valuesList)):
					valuesList.append(valuesList[k] + pow(2,j)*1)
			
		operator = identity(matrixDimension,dtype=int)
		for j in xrange(0,len(valuesList)):
			operator[decimalIndex0 + valuesList[j]][decimalIndex1 + valuesList[j]] = 1
			operator[decimalIndex1 + valuesList[j]][decimalIndex0 + valuesList[j]] = 1
			operator[decimalIndex0 + valuesList[j]][decimalIndex0 + valuesList[j]] = 0
			operator[decimalIndex1 + valuesList[j]][decimalIndex1 + valuesList[j]] = 0
		self.matrix = dot(operator, self.matrix)	

class OneBitComparator(object):
	"""docstring for OneBitComparator"""
	def __init__(self):
		super(OneBitComparator, self).__init__()
		matrixDimension = pow(2,1) * pow(2,1) * pow(2,1) * pow(2,1) # the first 2 multiplicands correspond to bits to compare; last two correspond to auxiliary bits
		self.matrix = identity(matrixDimension, dtype=int)
		
		# Operator 1: CNOT(3[1],1)
		dummyOperator = identity(matrixDimension, dtype=int)
		for bit2 in xrange(0,1):
			for bit0 in xrange(0,1):
				decimalValue0 = pow(2,3)*1 + pow(2,2)*bit2 + pow(2,1)*0 + pow(2,0)*bit0
				decimalValue1 = pow(2,3)*1 + pow(2,2)*bit2 + pow(2,1)*1 + pow(2,0)*bit0
				dummyOperator[decimalValue0][decimalValue1] = 1
				dummyOperator[decimalValue1][decimalValue0] = 1
				dummyOperator[decimalValue0][decimalValue0] = 0
				dummyOperator[decimalValue1][decimalValue1] = 0
		self.matrix = dot(self.matrix, dummyOperator)

		# Operator 2: CNOT(2[1],1)
		dummyOperator = identity(matrixDimension, dtype=int)
		for bit3 in xrange(0,1):
			for bit0 in xrange(0,1):
				decimalValue0 = pow(2,3)*bit3 + pow(2,2)*1 + pow(2,1)*0 + pow(2,0)*bit0
				decimalValue1 = pow(2,3)*bit3 + pow(2,2)*1 + pow(2,1)*1 + pow(2,0)*bit0
				dummyOperator[decimalValue0][decimalValue1] = 1
				dummyOperator[decimalValue1][decimalValue0] = 1
				dummyOperator[decimalValue0][decimalValue0] = 0
				dummyOperator[decimalValue1][decimalValue1] = 0
		self.matrix = dot(dummyOperator, self.matrix)

		# Operator 3: CCNOT(3[1], 1[1], 0)
		dummyOperator = identity(matrixDimension, dtype=int)
		for bit2 in xrange(0,1):
			decimalValue0 = pow(2,3)*1 + pow(2,2)*bit2 + pow(2,1)*1 + pow(2,0)*0
			decimalValue1 = pow(2,3)*1 + pow(2,2)*bit2 + pow(2,1)*1 + pow(2,0)*1
			dummyOperator[decimalValue0][decimalValue1] = 1
			dummyOperator[decimalValue1][decimalValue0] = 1
			dummyOperator[decimalValue0][decimalValue0] = 0
			dummyOperator[decimalValue1][decimalValue1] = 0
		self.matrix = dot(dummyOperator, self.matrix)

class FiveBitComparator(object):
	"""docstring for FiveBitComparator"""
	def __init__(self):
		super(FiveBitComparator, self).__init__()
		print 'FiveBitComparator <<'
		matrixDimension = pow(2,fixedBinaryNumberLength) * pow(2,fixedBinaryNumberLength) * pow(2,1) * pow(2,1) # first two sets of 5 bits correspond to binary numbers to compare; the other 2 bits are auxiliary bits for storing result
		self.matrix = identity(matrixDimension, dtype=int)

		# Operator 1: CNOT(11[1], 6); CNOT(controlBitIndex[controlBitValue], targetBitIndex)
		dummyOperator = identity(matrixDimension,dtype=int)
		for bit10 in xrange(0,2):
			for bit9 in xrange(0,2):
				for bit8 in xrange(0,2):
					for bit7 in xrange(0,2):
						for bit5 in xrange(0,2):
							for bit4 in xrange(0,2):
								for bit3 in xrange(0,2):
									for bit2 in xrange(0,2):
										for bit1 in xrange(0,2):
											for bit0 in xrange(0,2):
												decimalValue0 = pow(2,11)*1 + pow(2,10)*bit10 + pow(2,9)*bit9 + pow(2,8)*bit8 + pow(2,7)*bit7 + pow(2,6)*0 + pow(2,5)*bit5 + pow(2,4)*bit4 + pow(2,3)*bit3 + pow(2,2)*bit2 + pow(2,1)*bit1 + pow(2,0)*bit0
												decimalValue1 = pow(2,11)*1 + pow(2,10)*bit10 + pow(2,9)*bit9 + pow(2,8)*bit8 + pow(2,7)*bit7 + pow(2,6)*1 + pow(2,5)*bit5 + pow(2,4)*bit4 + pow(2,3)*bit3 + pow(2,2)*bit2 + pow(2,1)*bit1 + pow(2,0)*bit0
												dummyOperator[decimalValue0][decimalValue1] = 1
												dummyOperator[decimalValue1][decimalValue0] = 1
												dummyOperator[decimalValue0][decimalValue0] = 0
												dummyOperator[decimalValue1][decimalValue1] = 0
		print 'Operator 1: CNOT'
		self.matrix = dot(self.matrix, dummyOperator)

		# Operator 2: CNOT(6[1], 1)
		dummyOperator = identity(matrixDimension,dtype=int)
		for bit11 in xrange(0,2):
			for bit10 in xrange(0,2):
				for bit9 in xrange(0,2):
					for bit8 in xrange(0,2):
						for bit7 in xrange(0,2):
							for bit5 in xrange(0,2):
								for bit4 in xrange(0,2):
									for bit3 in xrange(0,2):
										for bit2 in xrange(0,2):
											for bit0 in xrange(0,2):
												decimalValue0 = pow(2,11)*bit11 + pow(2,10)*bit10 + pow(2,9)*bit9 + pow(2,8)*bit8 + pow(2,7)*bit7 + pow(2,6)*1 + pow(2,5)*bit5 + pow(2,4)*bit4 + pow(2,3)*bit3 + pow(2,2)*bit2 + pow(2,1)*0 + pow(2,0)*bit0
												decimalValue1 = pow(2,11)*bit11 + pow(2,10)*bit10 + pow(2,9)*bit9 + pow(2,8)*bit8 + pow(2,7)*bit7 + pow(2,6)*1 + pow(2,5)*bit5 + pow(2,4)*bit4 + pow(2,3)*bit3 + pow(2,2)*bit2 + pow(2,1)*1 + pow(2,0)*bit0
												dummyOperator[decimalValue0][decimalValue1] = 1
												dummyOperator[decimalValue1][decimalValue0] = 1
												dummyOperator[decimalValue0][decimalValue0] = 0
												dummyOperator[decimalValue1][decimalValue1] = 0		
		print 'Operator 2: CNOT'
		self.matrix = dot(dummyOperator, self.matrix)		

		# Operator 3: CCNOT(11[1], 1[1], 0)
		dummyOperator = identity(matrixDimension,dtype=int)
		for bit10 in xrange(0,2):
			for bit9 in xrange(0,2):
				for bit8 in xrange(0,2):
					for bit7 in xrange(0,2):
						for bit6 in xrange(0,2):
							for bit5 in xrange(0,2):
								for bit4 in xrange(0,2):
									for bit3 in xrange(0,2):
										for bit2 in xrange(0,2):
											decimalValue0 = pow(2,11)*1 + pow(2,10)*bit10 + pow(2,9)*bit9 + pow(2,8)*bit8 + pow(2,7)*bit7 + pow(2,6)*bit6 + pow(2,5)*bit5 + pow(2,4)*bit4 + pow(2,3)*bit3 + pow(2,2)*bit2 + pow(2,1)*1 + pow(2,0)*0
											decimalValue1 = pow(2,11)*1 + pow(2,10)*bit10 + pow(2,9)*bit9 + pow(2,8)*bit8 + pow(2,7)*bit7 + pow(2,6)*bit6 + pow(2,5)*bit5 + pow(2,4)*bit4 + pow(2,3)*bit3 + pow(2,2)*bit2 + pow(2,1)*1 + pow(2,0)*1
											dummyOperator[decimalValue0][decimalValue1] = 1
											dummyOperator[decimalValue1][decimalValue0] = 1
											dummyOperator[decimalValue0][decimalValue0] = 0
											dummyOperator[decimalValue1][decimalValue1] = 0		
		print 'Operator 3: CCNOT'
		self.matrix = dot(dummyOperator, self.matrix)								

		# Operator 4: CCNOT(10[1], 6[0], 5)
		dummyOperator = identity(matrixDimension,dtype=int)
		for bit11 in xrange(0,2):
			for bit9 in xrange(0,2):
				for bit8 in xrange(0,2):
					for bit7 in xrange(0,2):
						for bit4 in xrange(0,2):
							for bit3 in xrange(0,2):
								for bit2 in xrange(0,2):
									for bit1 in xrange(0,2):
										for bit0 in xrange(0,2):
											decimalValue0 = pow(2,11)*bit11 + pow(2,10)*1 + pow(2,9)*bit9 + pow(2,8)*bit8 + pow(2,7)*bit7 + pow(2,6)*0 + pow(2,5)*0 + pow(2,4)*bit4 + pow(2,3)*bit3 + pow(2,2)*bit2 + pow(2,1)*bit1 + pow(2,0)*bit0
											decimalValue1 = pow(2,11)*bit11 + pow(2,10)*1 + pow(2,9)*bit9 + pow(2,8)*bit8 + pow(2,7)*bit7 + pow(2,6)*0 + pow(2,5)*1 + pow(2,4)*bit4 + pow(2,3)*bit3 + pow(2,2)*bit2 + pow(2,1)*bit1 + pow(2,0)*bit0
											dummyOperator[decimalValue0][decimalValue1] = 1
											dummyOperator[decimalValue1][decimalValue0] = 1
											dummyOperator[decimalValue0][decimalValue0] = 0
											dummyOperator[decimalValue1][decimalValue1] = 0		
		print 'Operator 4: CCNOT'
		self.matrix = dot(dummyOperator, self.matrix)				

		# Operator 5: CCNOT(6[0], 5[1], 1)
		dummyOperator = identity(matrixDimension,dtype=int)
		for bit11 in xrange(0,2):
			for bit10 in xrange(0,2):
				for bit9 in xrange(0,2):
					for bit8 in xrange(0,2):
						for bit7 in xrange(0,2):
							for bit4 in xrange(0,2):
								for bit3 in xrange(0,2):
									for bit2 in xrange(0,2):
										for bit0 in xrange(0,2):
											decimalValue0 = pow(2,11)*bit11 + pow(2,10)*bit10 + pow(2,9)*bit9 + pow(2,8)*bit8 + pow(2,7)*bit7 + pow(2,6)*0 + pow(2,5)*1 + pow(2,4)*bit4 + pow(2,3)*bit3 + pow(2,2)*bit2 + pow(2,1)*0 + pow(2,0)*bit0
											decimalValue1 = pow(2,11)*bit11 + pow(2,10)*bit10 + pow(2,9)*bit9 + pow(2,8)*bit8 + pow(2,7)*bit7 + pow(2,6)*0 + pow(2,5)*1 + pow(2,4)*bit4 + pow(2,3)*bit3 + pow(2,2)*bit2 + pow(2,1)*1 + pow(2,0)*bit0
											dummyOperator[decimalValue0][decimalValue1] = 1
											dummyOperator[decimalValue1][decimalValue0] = 1
											dummyOperator[decimalValue0][decimalValue0] = 0
											dummyOperator[decimalValue1][decimalValue1] = 0		
		print 'Operator 5: CCNOT'
		self.matrix = dot(dummyOperator, self.matrix)						

		# Operator 6: CCCNOT(10[1], 6[0], 1[1], 0)
		dummyOperator = identity(matrixDimension,dtype=int)
		for bit11 in xrange(0,2):
			for bit9 in xrange(0,2):
				for bit8 in xrange(0,2):
					for bit7 in xrange(0,2):
						for bit5 in xrange(0,2):
							for bit4 in xrange(0,2):
								for bit3 in xrange(0,2):
									for bit2 in xrange(0,2):
										decimalValue0 = pow(2,11)*bit11 + pow(2,10)*1 + pow(2,9)*bit9 + pow(2,8)*bit8 + pow(2,7)*bit7 + pow(2,6)*0 + pow(2,5)*bit5 + pow(2,4)*bit4 + pow(2,3)*bit3 + pow(2,2)*bit2 + pow(2,1)*1 + pow(2,0)*0
										decimalValue1 = pow(2,11)*bit11 + pow(2,10)*1 + pow(2,9)*bit9 + pow(2,8)*bit8 + pow(2,7)*bit7 + pow(2,6)*0 + pow(2,5)*bit5 + pow(2,4)*bit4 + pow(2,3)*bit3 + pow(2,2)*bit2 + pow(2,1)*1 + pow(2,0)*1
										dummyOperator[decimalValue0][decimalValue1] = 1
										dummyOperator[decimalValue1][decimalValue0] = 1
										dummyOperator[decimalValue0][decimalValue0] = 0
										dummyOperator[decimalValue1][decimalValue1] = 0		
		print 'Operator 6: CCCNOT'
		self.matrix = dot(dummyOperator, self.matrix)								

		# Operator 7: CCNOT(9[1], 5[0], 4)
		dummyOperator = identity(matrixDimension,dtype=int)
		for bit11 in xrange(0,2):
			for bit10 in xrange(0,2):
				for bit8 in xrange(0,2):
					for bit7 in xrange(0,2):
						for bit6 in xrange(0,2):
							for bit3 in xrange(0,2):
								for bit2 in xrange(0,2):
									for bit1 in xrange(0,2):
										for bit0 in xrange(0,2):
											decimalValue0 = pow(2,11)*bit11 + pow(2,10)*bit10 + pow(2,9)*1 + pow(2,8)*bit8 + pow(2,7)*bit7 + pow(2,6)*bit6 + pow(2,5)*0 + pow(2,4)*0 + pow(2,3)*bit3 + pow(2,2)*bit2 + pow(2,1)*bit1 + pow(2,0)*bit0
											decimalValue1 = pow(2,11)*bit11 + pow(2,10)*bit10 + pow(2,9)*1 + pow(2,8)*bit8 + pow(2,7)*bit7 + pow(2,6)*bit6 + pow(2,5)*0 + pow(2,4)*1 + pow(2,3)*bit3 + pow(2,2)*bit2 + pow(2,1)*bit1 + pow(2,0)*bit0
											dummyOperator[decimalValue0][decimalValue1] = 1
											dummyOperator[decimalValue1][decimalValue0] = 1
											dummyOperator[decimalValue0][decimalValue0] = 0
											dummyOperator[decimalValue1][decimalValue1] = 0		
		print 'Operator 7: CCNOT'
		self.matrix = dot(dummyOperator, self.matrix)				

		# Operator 8: CCNOT(5[0], 4[1], 1)
		dummyOperator = identity(matrixDimension,dtype=int)
		for bit11 in xrange(0,2):
			for bit10 in xrange(0,2):
				for bit9 in xrange(0,2):
					for bit8 in xrange(0,2):
						for bit7 in xrange(0,2):
							for bit6 in xrange(0,2):
								for bit3 in xrange(0,2):
									for bit2 in xrange(0,2):
										for bit0 in xrange(0,2):
											decimalValue0 = pow(2,11)*bit11 + pow(2,10)*bit10 + pow(2,9)*bit9 + pow(2,8)*bit8 + pow(2,7)*bit7 + pow(2,6)*bit6 + pow(2,5)*0 + pow(2,4)*1 + pow(2,3)*bit3 + pow(2,2)*bit2 + pow(2,1)*0 + pow(2,0)*bit0
											decimalValue1 = pow(2,11)*bit11 + pow(2,10)*bit10 + pow(2,9)*bit9 + pow(2,8)*bit8 + pow(2,7)*bit7 + pow(2,6)*bit6 + pow(2,5)*0 + pow(2,4)*1 + pow(2,3)*bit3 + pow(2,2)*bit2 + pow(2,1)*1 + pow(2,0)*bit0
											dummyOperator[decimalValue0][decimalValue1] = 1
											dummyOperator[decimalValue1][decimalValue0] = 1
											dummyOperator[decimalValue0][decimalValue0] = 0
											dummyOperator[decimalValue1][decimalValue1] = 0		
		print 'Operator 8: CCNOT'
		self.matrix = dot(dummyOperator, self.matrix)						

		# Operator 9: CCCNOT(9[1], 5[0], 1[1], 0)
		dummyOperator = identity(matrixDimension,dtype=int)
		for bit11 in xrange(0,2):
			for bit10 in xrange(0,2):
				for bit8 in xrange(0,2):
					for bit7 in xrange(0,2):
						for bit6 in xrange(0,2):
							for bit4 in xrange(0,2):
								for bit3 in xrange(0,2):
									for bit2 in xrange(0,2):
										decimalValue0 = pow(2,11)*bit11 + pow(2,10)*bit10 + pow(2,9)*1 + pow(2,8)*bit8 + pow(2,7)*bit7 + pow(2,6)*bit6 + pow(2,5)*0 + pow(2,4)*bit4 + pow(2,3)*bit3 + pow(2,2)*bit2 + pow(2,1)*1 + pow(2,0)*0
										decimalValue1 = pow(2,11)*bit11 + pow(2,10)*bit10 + pow(2,9)*1 + pow(2,8)*bit8 + pow(2,7)*bit7 + pow(2,6)*bit6 + pow(2,5)*0 + pow(2,4)*bit4 + pow(2,3)*bit3 + pow(2,2)*bit2 + pow(2,1)*1 + pow(2,0)*1
										dummyOperator[decimalValue0][decimalValue1] = 1
										dummyOperator[decimalValue1][decimalValue0] = 1
										dummyOperator[decimalValue0][decimalValue0] = 0
										dummyOperator[decimalValue1][decimalValue1] = 0		
		print 'Operator 9: CCCNOT'
		self.matrix = dot(dummyOperator, self.matrix)										

		# Operator 10: CCNOT(8[1], 4[0], 3)
		dummyOperator = identity(matrixDimension,dtype=int)
		for bit11 in xrange(0,2):
			for bit10 in xrange(0,2):
				for bit9 in xrange(0,2):
					for bit7 in xrange(0,2):
						for bit6 in xrange(0,2):
							for bit5 in xrange(0,2):
								for bit2 in xrange(0,2):
									for bit1 in xrange(0,2):
										for bit0 in xrange(0,2):
											decimalValue0 = pow(2,11)*bit11 + pow(2,10)*bit10 + pow(2,9)*bit9 + pow(2,8)*1 + pow(2,7)*bit7 + pow(2,6)*bit6 + pow(2,5)*bit5 + pow(2,4)*0 + pow(2,3)*0 + pow(2,2)*bit2 + pow(2,1)*bit1 + pow(2,0)*bit0
											decimalValue1 = pow(2,11)*bit11 + pow(2,10)*bit10 + pow(2,9)*bit9 + pow(2,8)*1 + pow(2,7)*bit7 + pow(2,6)*bit6 + pow(2,5)*bit5 + pow(2,4)*0 + pow(2,3)*1 + pow(2,2)*bit2 + pow(2,1)*bit1 + pow(2,0)*bit0
											dummyOperator[decimalValue0][decimalValue1] = 1
											dummyOperator[decimalValue1][decimalValue0] = 1
											dummyOperator[decimalValue0][decimalValue0] = 0
											dummyOperator[decimalValue1][decimalValue1] = 0		
		print 'Operator 10: CCNOT'
		self.matrix = dot(dummyOperator, self.matrix)				

		# Operator 11: CCNOT(4[0], 3[1], 1)
		dummyOperator = identity(matrixDimension,dtype=int)
		for bit11 in xrange(0,2):
			for bit10 in xrange(0,2):
				for bit9 in xrange(0,2):
					for bit8 in xrange(0,2):
						for bit7 in xrange(0,2):
							for bit6 in xrange(0,2):
								for bit5 in xrange(0,2):
									for bit2 in xrange(0,2):
										for bit0 in xrange(0,2):
											decimalValue0 = pow(2,11)*bit11 + pow(2,10)*bit10 + pow(2,9)*bit9 + pow(2,8)*bit8 + pow(2,7)*bit7 + pow(2,6)*bit6 + pow(2,5)*bit5 + pow(2,4)*0 + pow(2,3)*1 + pow(2,2)*bit2 + pow(2,1)*0 + pow(2,0)*bit0
											decimalValue1 = pow(2,11)*bit11 + pow(2,10)*bit10 + pow(2,9)*bit9 + pow(2,8)*bit8 + pow(2,7)*bit7 + pow(2,6)*bit6 + pow(2,5)*bit5 + pow(2,4)*0 + pow(2,3)*1 + pow(2,2)*bit2 + pow(2,1)*1 + pow(2,0)*bit0
											dummyOperator[decimalValue0][decimalValue1] = 1
											dummyOperator[decimalValue1][decimalValue0] = 1
											dummyOperator[decimalValue0][decimalValue0] = 0
											dummyOperator[decimalValue1][decimalValue1] = 0		
		print 'Operator 11: CCNOT'
		self.matrix = dot(dummyOperator, self.matrix)						

		# Operator 12: CCCNOT(8[1], 4[0], 1[1], 0)
		dummyOperator = identity(matrixDimension,dtype=int)
		for bit11 in xrange(0,2):
			for bit10 in xrange(0,2):
				for bit9 in xrange(0,2):
					for bit7 in xrange(0,2):
						for bit6 in xrange(0,2):
							for bit5 in xrange(0,2):
								for bit3 in xrange(0,2):
									for bit2 in xrange(0,2):
										decimalValue0 = pow(2,11)*bit11 + pow(2,10)*bit10 + pow(2,9)*bit9 + pow(2,8)*1 + pow(2,7)*bit7 + pow(2,6)*bit6 + pow(2,5)*bit5 + pow(2,4)*0 + pow(2,3)*bit3 + pow(2,2)*bit2 + pow(2,1)*1 + pow(2,0)*0
										decimalValue1 = pow(2,11)*bit11 + pow(2,10)*bit10 + pow(2,9)*bit9 + pow(2,8)*1 + pow(2,7)*bit7 + pow(2,6)*bit6 + pow(2,5)*bit5 + pow(2,4)*0 + pow(2,3)*bit3 + pow(2,2)*bit2 + pow(2,1)*1 + pow(2,0)*1
										dummyOperator[decimalValue0][decimalValue1] = 1
										dummyOperator[decimalValue1][decimalValue0] = 1
										dummyOperator[decimalValue0][decimalValue0] = 0
										dummyOperator[decimalValue1][decimalValue1] = 0		
		print 'Operator 12: CCCNOT'
		self.matrix = dot(dummyOperator, self.matrix)										

		# Operator 13: CCNOT(7[1], 3[0], 2)
		dummyOperator = identity(matrixDimension,dtype=int)
		for bit11 in xrange(0,2):
			for bit10 in xrange(0,2):
				for bit9 in xrange(0,2):
					for bit8 in xrange(0,2):
						for bit6 in xrange(0,2):
							for bit5 in xrange(0,2):
								for bit4 in xrange(0,2):
									for bit1 in xrange(0,2):
										for bit0 in xrange(0,2):
											decimalValue0 = pow(2,11)*bit11 + pow(2,10)*bit10 + pow(2,9)*bit9 + pow(2,8)*bit8 + pow(2,7)*1 + pow(2,6)*bit6 + pow(2,5)*bit5 + pow(2,4)*bit4 + pow(2,3)*0 + pow(2,2)*0 + pow(2,1)*bit1 + pow(2,0)*bit0
											decimalValue1 = pow(2,11)*bit11 + pow(2,10)*bit10 + pow(2,9)*bit9 + pow(2,8)*bit8 + pow(2,7)*1 + pow(2,6)*bit6 + pow(2,5)*bit5 + pow(2,4)*bit4 + pow(2,3)*0 + pow(2,2)*1 + pow(2,1)*bit1 + pow(2,0)*bit0
											dummyOperator[decimalValue0][decimalValue1] = 1
											dummyOperator[decimalValue1][decimalValue0] = 1
											dummyOperator[decimalValue0][decimalValue0] = 0
											dummyOperator[decimalValue1][decimalValue1] = 0		
		print 'Operator 13: CCNOT'
		self.matrix = dot(dummyOperator, self.matrix)				

		# Operator 14: CCNOT(3[0], 2[1], 1)
		dummyOperator = identity(matrixDimension,dtype=int)
		for bit11 in xrange(0,2):
			for bit10 in xrange(0,2):
				for bit9 in xrange(0,2):
					for bit8 in xrange(0,2):
						for bit7 in xrange(0,2):
							for bit6 in xrange(0,2):
								for bit5 in xrange(0,2):
									for bit4 in xrange(0,2):
										for bit0 in xrange(0,2):
											decimalValue0 = pow(2,11)*bit11 + pow(2,10)*bit10 + pow(2,9)*bit9 + pow(2,8)*bit8 + pow(2,7)*bit7 + pow(2,6)*bit6 + pow(2,5)*bit5 + pow(2,4)*bit4 + pow(2,3)*0 + pow(2,2)*1 + pow(2,1)*0 + pow(2,0)*bit0
											decimalValue1 = pow(2,11)*bit11 + pow(2,10)*bit10 + pow(2,9)*bit9 + pow(2,8)*bit8 + pow(2,7)*bit7 + pow(2,6)*bit6 + pow(2,5)*bit5 + pow(2,4)*bit4 + pow(2,3)*0 + pow(2,2)*1 + pow(2,1)*1 + pow(2,0)*bit0
											dummyOperator[decimalValue0][decimalValue1] = 1
											dummyOperator[decimalValue1][decimalValue0] = 1
											dummyOperator[decimalValue0][decimalValue0] = 0
											dummyOperator[decimalValue1][decimalValue1] = 0		
		print 'Operator 14: CCNOT'
		self.matrix = dot(dummyOperator, self.matrix)						

		# Operator 15: CCCNOT(7[1], 3[0], 1[1], 0)
		dummyOperator = identity(matrixDimension,dtype=int)
		for bit11 in xrange(0,2):
			for bit10 in xrange(0,2):
				for bit9 in xrange(0,2):
					for bit8 in xrange(0,2):
						for bit6 in xrange(0,2):
							for bit5 in xrange(0,2):
								for bit4 in xrange(0,2):
									for bit2 in xrange(0,2):
										decimalValue0 = pow(2,11)*bit11 + pow(2,10)*bit10 + pow(2,9)*bit9 + pow(2,8)*bit8 + pow(2,7)*1 + pow(2,6)*bit6 + pow(2,5)*bit5 + pow(2,4)*bit4 + pow(2,3)*0 + pow(2,2)*bit2 + pow(2,1)*1 + pow(2,0)*0
										decimalValue1 = pow(2,11)*bit11 + pow(2,10)*bit10 + pow(2,9)*bit9 + pow(2,8)*bit8 + pow(2,7)*1 + pow(2,6)*bit6 + pow(2,5)*bit5 + pow(2,4)*bit4 + pow(2,3)*0 + pow(2,2)*bit2 + pow(2,1)*1 + pow(2,0)*1
										dummyOperator[decimalValue0][decimalValue1] = 1
										dummyOperator[decimalValue1][decimalValue0] = 1
										dummyOperator[decimalValue0][decimalValue0] = 0
										dummyOperator[decimalValue1][decimalValue1] = 0		
		print 'Operator 15: CCCNOT'
		self.matrix = dot(dummyOperator, self.matrix)										

class BinaryNumberComparator(object):
	"""docstring for BinaryNumberComparator"""
	def __init__(self):
		super(BinaryNumberComparator, self).__init__()
		self.operator = OneBitComparator()
		
	def compare(self, binaryNumber1, binaryNumber2):
		basisVector1 = convertToBasisVectorFromBinaryNumber( binaryNumber1 )
		basisVector2 = convertToBasisVectorFromBinaryNumber( binaryNumber2 )
		registerBasisVector = kron(basisVector1, kron(basisVector2, kron(ket0(), ket0())))
		registerBasisVector = dot(self.operator.matrix, registerBasisVector).tolist()
		indexOf1 = registerBasisVector.index(1)

		# identify value for first auxiliary bit
		equivalenceIndicator = '0'
		startingIndices = [(index * 4) for index in range(4)]
		if indexOf1 in startingIndices or (indexOf1-1) in startingIndices:
			equivalenceIndicator = '0'
		elif (indexOf1-2) in startingIndices or (indexOf1-3) in startingIndices:
			equivalenceIndicator = '1'

		# identify value for second auxiliary bit
		relationIndicator = '0'
		if indexOf1 % 2 == 1:
			relationIndicator = '1'

		return equivalenceIndicator + relationIndicator


class Alphabet(object):
	"""docstring for Alphabet"""
	def __init__(self):
		super(Alphabet, self).__init__()
		self.symbols = [Symbol('a'), Symbol('c'), Symbol('t'), Symbol('g')]
		self.codes = [[0,0],[1,0],[0,1],[1,1]]
		self.codeBitCount = 2
		self.size = len(self.symbols)
	
	def getSymbols(self):
		return self.symbols

	def getCodeBitCount(self):
		return self.codeBitCount

	def getBinaryCode(self,character):
		for i in xrange(0,self.size):
			if self.symbols[i].character == character:
				return self.codes[i]

	def toString(self):
		text = 'Sigma = {'
		for symbol in self.symbols:
			text = text + symbol + ', '
		text = text + '}'
		return text

class Symbol(object):
	"""A class for a symbol in an alphabet. A symbol is represented as a character."""
	def __init__(self, character):
		super(Symbol, self).__init__()
		self.character = character
	
	def toBinaryCodeFromCharacter(self):
		print 'Symbol character:',self.character
		alphabet = Alphabet()
		binaryCode = alphabet.getBinaryCode(self.character)
		print 'Symbol character:',self.character, ', Binary code:',binaryCode
		return binaryCode

	def toBasisVectorFromBinaryCode(self):
		basisVector = convertToBasisVectorFromSymbolBinaryCode(self.toBinaryCodeFromCharacter())
		return basisVector

	# Compare this symbol's character with another symbol's character.
	# INPUT:
	# 	symbol - the character to be compared with this symbol's character
	# OUTPUT:
	# 	The output is either a 0 for a match and 1 for a non-match of the two symbols.
	def compare(self,symbol):
		comparisonResult = compareSymbolBinaryCodes(self.toBinaryCodeFromCharacter(), symbol.toBinaryCodeFromCharacter())
		return comparisonResult

	def toString(self):
		return self.character

# A class for defining a bit.
# INPUT:
# 	value - the bit value, a 0 or a 1
class Bit(object):
	"""docstring for Bit"""
	def __init__(self, value):
		super(Bit, self).__init__()
		self.value = value
		
	def add(self, carryBit, bit):
		print 'BEFORE: Bit.add: carryBit.value =',carryBit.value,', self.value =',self.value,', bit.value =',bit.value
		additionResult = FullAdder().add(carryBit.value, self.value, bit.value)	
		print 'AFTER: Bit.add: carryBit.value =',additionResult.carryBitValue,', self.value =',additionResult.sumBitValue,', bit.value =',bit.value
		self.value = additionResult.sumBitValue
		return Bit(additionResult.carryBitValue)
		
	def toString(self):
		return str(self.value)


# A class for defining a binary number.
# INPUT:
#	binaryString - a list of 0s and 1s of the binary number to create. The purpose of the constructor is to pad the list with 0s to
#					fit the register which has, by convention, maximum of 5 bits.
class BinaryNumber(object):
	"""docstring for BinaryNumber"""
	def __init__(self, binaryString):
		super(BinaryNumber, self).__init__()
		self.value = zeroBits(fixedBinaryNumberLength)
		for (index,bitValue) in enumerate(binaryString): # create a list of Bit()s
			self.value[index].value = bitValue
	
	def singleBit(self):
		self.value = [self.value[0]]

	def add(self, otherBinaryNumber):
		print 'Add:',self.toString(),' +',otherBinaryNumber.toString()
		carryBit = Bit(0)
		for selfBit, otherBit in zip(self.value, otherBinaryNumber.value):
			carryBit = selfBit.add(carryBit, otherBit)

	# Check if binary number self.value is greater than input binary number.
	# INPUT:
	# 	otherBinaryNumber - the binaryNumber to compare to
	# OUTPUT:
	# 	True - binary number self.value is greater than input binary number
	# 	False - binary number self.value is less or equal to input binary number
	def isGreaterThan(self, otherBinaryNumber):
		binaryNumberComparator = BinaryNumberComparator()
		result = binaryNumberComparator.compare(self, otherBinaryNumber)
		print 'result:',result
		if result == '11':
			return True
		else:
			return False

	# Return the list of bit values for this binary number.
	def toBitValuesList(self):
		bitValuesList = []
		for bit in self.value:
			bitValuesList.append(bit.value)
		return bitValuesList

	def toString(self):
		text = ''
		for bit in self.value:
			text = text + bit.toString()
		return text

class AdditionResult(object):
	"""docstring for AdditionResult"""
	def __init__(self, sumBitValue, carryBitValue):
		super(AdditionResult, self).__init__()
		self.sumBitValue = sumBitValue
		self.carryBitValue = carryBitValue		

# Input:
#	amplitude - amplitude of state
# 	substring - an M-length list of log_2(Sigma)-length binary list representing symbols
# 	pattern - an M-length list of log_2(Sigma)-length binary list representing symbols
class QuantumStateForVerification(object):
	"""A quantum state representation for the verification phase of algorithm 2."""
	def __init__(self, amplitude,substring,pattern,auxiliary,thresholdDistance):
		self.amplitude = amplitude # float
		self.substringState = substring # M-length list of binary numbers (in log_2(M)-length list format) representing symbols in the substring
		self.patternState = pattern # M-length list of binary numbers (in log_2(M)-length list format) representing symbols in the pattern
		self.auxiliary = auxiliary # M-length list of 0s and 1s each corresponding to symbol in substring
		# self.mismatchCount = BinaryNumber([Bit(0)]) # log_2(M)-bit-length binary number representing count of mismatches in the substring
		self.mismatchCount = 0
		self.thresholdDistance = thresholdDistance # maximum number mismatches allowed for the substring
		self.isValid = 0

	def markIfValid(self):
		# mark this state with a negative amplitude if self.mismatchCount <= self.thresholdDistance
		# result is a binary list; 00 => mismatchCount == thresholdDistance, 10 => mismatchCount > thresholdDistance, 01 => mismatchCount < thresholdDistance
		# binaryComparisonResult = Util.compareBinaryNumbers(self.mismatchCount, self.thresholdDistance) 

		mismatchCountBinaryNumber = convertToBinaryNumberFromDecimalNumber(self.mismatchCount) # convert decimal value for mismatch count for this substring into binary number
		thresholdDistanceBinaryNumber = convertToBinaryNumberFromDecimalNumber(self.thresholdDistance) # convert decimal value for threshold distance into binary number
		print 'self.mismatchCount:',self.mismatchCount,', mismatchCountBinaryNumber:',mismatchCountBinaryNumber.toString()
		print 'thresholdDistanceBinaryNumber:',thresholdDistanceBinaryNumber.toString()
		mismatchCountIsGreaterThanThresholdDistance = False
		for (mismatchCountBit, thresholdDistanceBit) in zip(mismatchCountBinaryNumber.value, thresholdDistanceBinaryNumber.value):
			bit1 = BinaryNumber([mismatchCountBit.value])
			bit2 = BinaryNumber([thresholdDistanceBit.value])
			bit1.singleBit()
			bit2.singleBit()
			if bit1.isGreaterThan(bit2): # compare in binary the value for mismatch count for this substring itno threshold distance
				mismatchCountIsGreaterThanThresholdDistance = True
		if not mismatchCountIsGreaterThanThresholdDistance:
			self.amplitude = -1 * self.amplitude # mark this state with negative sign of amplitude if the substring it represents has mismatch count less or equal to threshold		


		

		
	def toString(self):
		text = 'amplitude: ' + str(self.amplitude) + ', substring: '
		for symbol in self.substringState:
			text = text + symbol.character
		text = text + ', pattern: '
		for symbol in self.patternState:
			text = text + symbol.character
		text = text + ', auxiliary: '
		for auxiliaryBit in self.auxiliary:
		 	text = text + str(auxiliaryBit)
		text = text + ', mismatchCount: ' + str(self.mismatchCount)
		text = text + ', thresholdDistance:' + str(self.thresholdDistance)
		return text


	def toStringSimple(self):
		return '['+ str(self.amplitude) + ', [' + str(self.substringState) + ', ' + str(self.patternState) + ', ' + str(self.auxiliary) + ', ' + str(self.mismatchCount) + ',' + str(thresholdDistance) +']]'

def zeroBits(bitCount):
	bits = []
	for i in xrange(0,bitCount):
		bits.append(Bit(0))
	return bits


# Compare two symbols expressed in their binary codes. Computation is performed using matrix operator.
# INPUT:
#	binaryCode1, binaryCode2 - binary codes of symbols to compare, each of same length
# OUTUT:
#   1 - the two symbols do not match
# 	0 - the two symbols match
def compareSymbolBinaryCodes(string1, string2):
	symbolComparator = SymbolComparator().matrix
	string1BasisVector = convertToBasisVectorFromSymbolBinaryCode(string1)
	string2BasisVector = convertToBasisVectorFromSymbolBinaryCode(string2)
	registerBasisVector = kron(string1BasisVector, kron(string2BasisVector, ket0())) # ket0() corresponds to the single auxiliary bit
	registerBasisVector = dot(symbolComparator, registerBasisVector).tolist()
	indexOf1 = registerBasisVector.index(1)
	comparisonResult = 0
	if indexOf1 % 2 == 1 :
		comparisonResult = 1
	print 'Comparison result:',comparisonResult
	return comparisonResult # a single bit value; 1 for non-match and 0 for match


# Get the state vector representation of a binary code for an alphabet symbol.
# INPUT:
#	binaryCode - the binary code of the alphabet symbol to convert to basis vector representation
#	dimension - the dimension of the basis vector to create
# OUTPUT:
#	A basis vector representation of the input binary code with dimension as specified
#	by input dimension.
def convertToBasisVectorFromSymbolBinaryCode(binaryCode):
	bitCount = len(binaryCode)
	decimalNumber = 0
	for i in xrange(0,bitCount):
		decimalNumber = decimalNumber + pow(2,i)*binaryCode[i]
	return convertToBasisVectorFromDecimalNumber(decimalNumber, Alphabet().size)

# Get the basis vector representation of a decimal number.
# INPUT:
#	binaryCode - the binary code to convert to state vector representation
#	dimension - the dimension of the state vector to create
# OUTPUT:
#	A state vector representation of the input binary code with dimension as specified
#	by input dimension.
def convertToBasisVectorFromDecimalNumber(decimalNumber, dimension):
	basisVector = zeros(dimension, dtype=int)
	basisVector[decimalNumber] = 1
	return basisVector.tolist()

def convertToBasisVectorFromBinaryNumber(binaryNumber):
	bitCount = len(binaryNumber.value)
	decimalNumber = 0
	print 'binaryNumber:',binaryNumber.toString()
	for bitIndex in xrange(0,bitCount):
		decimalNumber = decimalNumber + pow(2,bitIndex)*binaryNumber.value[bitIndex].value
	return convertToBasisVectorFromDecimalNumber(decimalNumber, pow(2,bitCount))

def convertToBinaryNumberFromDecimalNumber(decimalNumber):
	binaryString = []
	for x in xrange(0,fixedBinaryNumberLength):
		i = fixedBinaryNumberLength - x - 1
		if decimalNumber >= pow(2,i):
			binaryString.append(1)
			decimalNumber = decimalNumber - pow(2,i)
		else:
			binaryString.append(0)
	binaryString.reverse()
	return BinaryNumber(binaryString)

def convertToDecimalNumberFromBinaryNumber(binaryNumber):
	decimalNumber = 0
	for bitIndex in xrange(0,len(binaryNumber.value)):
		decimalNumber = decimalNumber + pow(2,bitIndex)*binaryNumber.value[bitIndex].value
	print 'binaryNumber:',binaryNumber.toString()
	print 'decimalNumber:',decimalNumber
	return decimalNumber


def prepareQuantumRegistersForVerificationPhase(substrings, pattern, threshold):
	
	patternQuantumRegisters = []
	for patternCharacter in pattern:
		patternCharacterSymbol = Symbol(patternCharacter)
		patternCharacterBinaryCode = patternCharacterSymbol.toBinaryCodeFromCharacter()
		patternCharacterBasisVector = patternCharacterSymbol.toBasisVectorFromBinaryCode()
		patternCharacterQuantumState = QuantumState(1.0, patternCharacterBinaryCode, patternCharacter, patternCharacterBasisVector)
		patternCharacterQuantumRegister = QuantumRegister(2, patternCharacterQuantumState)
		patternQuantumRegisters.append(patternCharacterQuantumRegister)


	mismatchesIndicatorQuantumRegister = []
	for i in xrange(0,len(pattern)):
		mismatchesIndicatorQuantumRegister.append(QuantumState(1.0, [0], 0, ket0())) 


	mismatchesCountQuantumRegister = []
	mismatchesCountBinaryNumber = convertToBinaryNumberFromDecimalNumber(0)
	for bit in mismatchesCountBinaryNumber.value:
		basisVector = ket0()
		if bit.value == 1:
			basisVector = ket1()
		mismatchesCountQuantumRegister.append(QuantumState(1.0, [bit.value], bit.value, basisVector))

	
	thresholdQuantumRegister = []
	thresholdBinaryNumber = convertToBinaryNumberFromDecimalNumber(threshold)
	for bit in thresholdBinaryNumber.value:
		basisVector = ket0()
		if bit.value == 1:
			basisVector = ket1()
		thresholdQuantumRegister.append(QuantumState(1.0, [bit.value], bit.value, basisVector))

	verificationPhaseQuantumRegisters = []
	for substring in substrings:
		print 'substring:',substring
		substringQuantumRegisters = []
		for substringCharacter in substring:
			print 'substring character:',substringCharacter
			substringCharacterSymbol = Symbol(substringCharacter)
			substringCharacterBinaryCode = substringCharacterSymbol.toBinaryCodeFromCharacter()
			substringCharacterBasisVector = substringCharacterSymbol.toBasisVectorFromBinaryCode()
			substringCharacterQuantumState = QuantumState(1.0, substringCharacterBinaryCode, substringCharacter, substringCharacterBasisVector)
			substringCharacterQuantumRegister = QuantumRegister(2, substringCharacterQuantumState)
			substringQuantumRegisters.append(substringCharacterQuantumRegister)
		# verificationPhaseQuantumRegister = []
		# verificationPhaseQuantumRegister.append(substringQuantumRegisters)
		# verificationPhaseQuantumRegister.append(patternQuantumRegisters)
		# verificationPhaseQuantumRegister.append(mismatchesIndicatorQuantumRegister)
		# verificationPhaseQuantumRegister.append(mismatchesCountQuantumRegister)
		# verificationPhaseQuantumRegister.append(thresholdQuantumRegister)
		# verificationPhaseQuantumRegisters.append(verificationPhaseQuantumRegister)

		verificationPhaseQuantumRegister = VerificationPhaseQuantumRegister(substringQuantumRegisters, patternQuantumRegisters, mismatchesIndicatorQuantumRegister, mismatchesCountQuantumRegister, thresholdQuantumRegister)
		verificationPhaseQuantumRegisters.append(verificationPhaseQuantumRegister)

	return verificationPhaseQuantumRegisters

class VerificationPhaseQuantumRegister(object):
	"""docstring for VerificationPhaseQuantumRegister"""
	def __init__(self, substringQuantumRegisters, patternQuantumRegisters, mismatchesIndicatorQuantumRegister, mismatchesCountQuantumRegister, thresholdQuantumRegister):
		super(VerificationPhaseQuantumRegister, self).__init__()
		self.substringQuantumRegister = substringQuantumRegisters # list of qubits
		self.patternQuantumRegisters = patternQuantumRegisters
		self.mismatchesIndicatorQuantumRegister = mismatchesIndicatorQuantumRegister
		self.mismatchesCountQuantumRegister = mismatchesCountQuantumRegister
		self.thresholdQuantumRegister = thresholdQuantumRegister
		


# substrings = ['ac','at']
# pattern = 'ac'
# verificationPhaseQuantumRegisters = prepareQuantumRegistersForVerificationPhase(substrings, pattern, 3)

# # for verificationPhaseQuantumRegister in verificationPhaseQuantumRegisters:
# # 	substringCharacterSymbolList = []
# # 	for substringCharacterQuantumRegister in verificationPhaseQuantumRegister[0]:
# # 		substringCharacterSymbolList.append(Symbol(substringCharacterQuantumRegister.state.symbolFormat))
# # 	patternCharacterSymbolList = []
# # 	for patternCharacterQuantumRegister in verificationPhaseQuantumRegister[1]:
# # 		patternCharacterSymbolList.append(Symbol(patternCharacterQuantumRegister.state.symbolFormat))
# # 	mismatchesBitList = []
# # 	for mismatchesBit in verificationPhaseQuantumRegister[2]
# # 		mismatchesBitList.append(mismatchesBit.state.binaryNumberFormat[0])
# # 	thresholdBitList = []
# # 	for thresholdBit in verificationPhaseQuantumRegister[4]
# # 		thresholdBitList.append(thresholdBit.state.decimalNumberFormat)

# # QuantumStateForVerification(1.0, substringCharacterSymbolList, patternCharacterSymbolList, mismatchesBitList, threshold)

# for (index,verificationPhaseQuantumRegister) in enumerate(verificationPhaseQuantumRegisters):
# 	substringQuantumRegisters = verificationPhaseQuantumRegister[0]
# 	print 'Superposition state', index
# 	print 'Substring register <<'
# 	for substringCharacterQuantumRegister in substringQuantumRegisters:
# 		print substringCharacterQuantumRegister.toString()
# 	print '>>'

# 	patternQuantumRegisters = verificationPhaseQuantumRegister[1]
# 	print 'Pattern register <<'
# 	for patternCharacterQuantumRegister in patternQuantumRegisters:
# 		print patternCharacterQuantumRegister.toString()
# 	print '>>'

# 	mismatchesQuantumRegister = verificationPhaseQuantumRegister[2]
# 	print 'Mismatches register <<'
# 	for mismatchesBitQuantumRegister in mismatchesQuantumRegister:
# 		print mismatchesBitQuantumRegister.toString()
# 	print '>>' 

# 	mismatchesCountQuantumRegister = verificationPhaseQuantumRegister[3]
# 	print 'Mismatches count register <<'
# 	for mismatchesCountBitQuantumRegister in mismatchesCountQuantumRegister:
# 		print mismatchesCountBitQuantumRegister.toString()
# 	print '>>' 	

# 	thresholdQuantumRegister = verificationPhaseQuantumRegister[4]
# 	print 'Threshold register <<'
# 	for thresholdBitQuantumRegister in thresholdQuantumRegister:
# 		print thresholdBitQuantumRegister.toString()

# 	print 
