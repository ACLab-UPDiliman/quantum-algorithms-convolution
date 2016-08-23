from numpy import kron, array, identity, ones, zeros, dot, matrix
from math import sqrt, ceil, log, pi, floor
import matplotlib.pyplot as plotter
import logging

maxPatternLength = 2
oneBitAdderBitCount = int(ceil(log(maxPatternLength, 2))) + 1

registerBitCount = 2 # translates to maximum pattern length of 2^registerBitCount characters

# Define logging tool.
logging.basicConfig(format='%(levelname)s:%(message)s',level=logging.INFO)
# ==================================================================== CLASSES ====================================================================

# Already added to Xcode version.
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

# Already added to Xcode version.
class Symbol(object):
	"""docstring for Symbol"""
	def __init__(self, character):
		super(Symbol, self).__init__()
		self.character = character

	def getCharacter(self):
		return self.character

	def getBinaryCode(self):
		if self.character == 'a':
			return [0,0]
		elif self.character == 'c':
			return [1,0]
		elif self.character == 't':
			return [0,1]
		elif self.character == 'g':
			return [1,1]

# Already added to Xcode version.
class Bit(object):
	"""docstring for Bit"""
	def __init__(self, state):
		super(Bit, self).__init__()
		self.state = state
		if state == 0:
			self.basisVector = ket0()
		else:
			self.basisVector = ket1()

	def getState(self):
		return self.state

	def getBasisVector(self):
		if self.state == 0:
			self.basisVector = ket0()
		else:
			self.basisVector = ket1()
		return self.basisVector

	def add(self, carryBit, bit):
		# print 'BEFORE: Bit.add: carryBit.value =',carryBit.state,', self.value =',self.state,', bit.value =',bit.state
		additionResult = FullAdder().add(carryBit.state, self.state, bit.state)	
		# print 'AFTER: Bit.add: carryBit.value =',additionResult.carryBitState,', self.value =',additionResult.sumBitState,', bit.value =',bit.state
		self.state = additionResult.sumBitState
		return Bit(additionResult.carryBitState)

	def isGreaterThan(self, otherBit):
		# print 'Bit.isGreaterThan: self.state:',self.state, ' otherBit.state',otherBit.state
		result = BitComparator().compare(self, otherBit)
		# print 'Bit.isGreaterThan: result:',result
		return result

	def toString(self):
		return str(self.state)

# Already added to Xcode version.
class BinaryNumber(object):
	"""docstring for BinaryNumber"""
	def __init__(self, binaryString):
		super(BinaryNumber, self).__init__()
		self.value = zeroBits(registerBitCount) # list of Bit
		for (index,bitState) in enumerate(binaryString): # create a list of Bit()s
			self.value[index].state = bitState
	
	def singleBit(self):
		self.value = [self.value[0]]

	def add(self, otherBinaryNumber):
		# print 'Add:',self.toString(),' +',otherBinaryNumber.toString()
		carryBit = Bit(0)
		for selfBit, otherBit in zip(self.value, otherBinaryNumber.value):
			carryBit = selfBit.add(carryBit, otherBit)

	def subtract(self, otherBinaryNumber):
		# Convert self.value to 2's complement representation
		otherTwosComplement = otherBinaryNumber.getTwosComplement()
		add(self, otherTwosComplement)



	def isGreaterThan(self, otherBinaryNumber):
		reversedSelfBitList = list(self.value)
		reversedSelfBitList.reverse()
		reversedOtherBitList = list(otherBinaryNumber.value)
		reversedOtherBitList.reverse()
		for (bit, otherBit) in zip(reversedSelfBitList, reversedOtherBitList):
			result = bit.isGreaterThan(otherBit)
			if result == '11':
				return True
		return False
		# binaryNumberComparator = BinaryNumberComparator()
		# result = binaryNumberComparator.compare(self, otherBinaryNumber)
		# print 'result:',result
		# if result == '11':
		# 	return True
		# else:
		# 	return False

	def toString(self):
		text = '['
		for bit in self.value:
			text = text + bit.toString()
		text = text + ']'
		return text

class Register(object):
	"""docstring for Register"""
	def __init__(self, bits):
		super(Register, self).__init__()
		self.bits = bits

	def getState(self):
		bitList = []
		for bit in self.bits:
			bitList.append(bit.getState())
		return bitList

	def getBitsBasisVector(self):
		basisVectorList = []
		for bit in self.bits:
			basisVectorList.append(bit.basisVector)
		return basisVectorList

	def getRegisterBasisVector(self):
		basisVector = self.bits[0].basisVector
		for index in xrange(1, len(self.bits)):
			basisVector = kron(self.bits[index].basisVector, basisVector)
		return basisVector

class SymbolRegister(object):
	"""docstring for Symbol Register"""
	def __init__(self, symbol, bits):
		super(SymbolRegister, self).__init__()
		self.symbol = symbol
		self.bits = bits

	def getSymbol(self):
		return self.symbol

	def getState(self):
		bitList = []
		for bit in self.bits:
			bitList.append(bit.state)
		return bitList

	def getBitsBasisVector(self):
		basisVectorList = []
		for bit in self.bits:
			basisVectorList.append(bit.basisVector)
		return basisVectorList

	def getRegisterBasisVector(self):
		basisVector = self.bits[0].basisVector
		for index in xrange(1, len(self.bits)):
			basisVector = kron(self.bits[index].basisVector, basisVector)
		return basisVector		

	def compare(self, otherRegister):
		symbolComparator = SymbolComparator()
		comparisonResult = symbolComparator.compare(self.getRegisterBasisVector(), otherRegister.getRegisterBasisVector())
		return comparisonResult

class FilteringPhaseRegister(object):
	"""docstring for FilteringPhaseRegister"""
	def __init__(self, amplitude, indexRegister, startingLocationRegister):
		super(FilteringPhaseRegister, self).__init__()
		self.amplitude = amplitude 
		
		# -index in text text
		# -valid values range from 0 to N-1
		self.indexRegister = indexRegister
		
		# -starting index of a possble solution substring in text text
		# -expressed as a list of bits representing the starting index with the rightmost bit as the most significant bit
		self.startingLocationRegister = startingLocationRegister 
		
	def getAmplitude(self):
		return self.amplitude

	def getIndexRegister(self):
		return self.indexRegister

	def getStartingLocationRegister(self):
		return self.startingLocationRegister

	# get the basis vector representation of the binary state of the indexRegister
	def getIndexRegisterBasisVector(self):
		return convertToBasisVectorFromBinaryNumber(BinaryNumber(self.indexRegister))

	# get the basis vector representation of the binary state of the locationRegister
	def getStartingLocationRegisterBasisVector(self):
		return convertToBasisVectorFromBinaryNumber(BinaryNumber(self.startingLocationRegister))

	# -Fetch location of first occurrence in pattern pattern of the symbol located in the index represented by the binary state of the indexRegister
	# -Return the location of first occurrence in pattern pattern as a list of bits
	def identifyLocationOfFirstOccurrence(self):
		# print 'FilteringPhaseRegister: identifyLocationOfFirstOccurrence: self.getIndexRegisterBasisVector():', self.getIndexRegisterBasisVector()
		# print 'FilteringPhaseRegister: identifyLocationOfFirstOccurrence: self.getStartingLocationRegisterBasisVector():', self.getStartingLocationRegisterBasisVector()
		basisVector = kron(self.getIndexRegisterBasisVector(), self.getStartingLocationRegisterBasisVector()) 
		# print 'FilteringPhaseRegister: identifyLocationOfFirstOccurrence: basisVector:', basisVector
		self.startingLocationRegister = convertToBinaryStringFromBasisVector(locationOperator.identifyLocationOfFirstOccurrence(basisVector).tolist(), registerBitCount) # list of 0s and 1s representing binary state of startingLocationRegister
		# print 'FilteringPhaseRegister: identifyLocationOfFirstOccurrence: new basisVector:', self.startingLocationRegister

	# -Compute for possible starting location of pattern pattern in text text given the index represented by the state of the indexRegister and the location of first occurrence
	#	in pattern pattern
	def computeStartingLocation(self):
		# get the negative representation of the location represented by the state of locationRegister
		startingLocationRegisterStateTwosComplement = getTwosComplement(self.startingLocationRegister)

		# startingLocationRegisterStateTwosComplementBinaryNumber = BinaryNumber(startingLocationRegisterStateTwosComplement)
		# self.startingLocationRegister = convertToBinaryStringFromBasisVector()
		#===== LAST EDIT AS OF JANUARY 22, 2016 ===== 

	def toString(self):
		text = 'Amplitude:', self.amplitude, ' Index register:', self.indexRegister, ' Starting location register:', self.startingLocationRegister
		return text


class VerificationPhaseRegister(object):
	"""docstring for VerificationPhaseRegister"""
	def __init__(self, amplitude, substringRegister, patternRegister, mismatchesIndicatorRegister, mismatchesCountRegister, thresholdRegister):
		super(VerificationPhaseRegister, self).__init__()
		self.amplitude = amplitude
		self.substringRegister = substringRegister # list of qubits
		self.patternRegister = patternRegister
		self.mismatchesIndicatorRegister = mismatchesIndicatorRegister
		self.mismatchesCountRegister = mismatchesCountRegister
		self.thresholdRegister = thresholdRegister

	def getAmplitude(self):
		return self.amplitude

	def setAmplitude(self,amplitude):
		self.amplitude = amplitude

	def getSubstringRegister(self):
		return self.substringRegister

	def getPatternRegister(self):
		return self.patternRegister

	def getMismatchesIndicatorRegister(self):
		return self.mismatchesIndicatorRegister

	def getMismatchesCountRegister(self):
		return self.mismatchesCountRegister

	def getThresholdRegister(self):
		return self.thresholdRegister

	def countMismatches(self):
		mismatchesCountBinaryNumber = BinaryNumber(self.mismatchesCountRegister)
		for indicatorBit in self.mismatchesIndicatorRegister:
			if indicatorBit == 1:
				indicatorBitBinaryNumber = BinaryNumber([indicatorBit])
				mismatchesCountBinaryNumber.add(indicatorBitBinaryNumber)
		for (index,bit) in enumerate(mismatchesCountBinaryNumber.value):
			self.mismatchesCountRegister[index] = bit.state

	def markIfMismatchesCountIsValid(self):
		# mark this state with a negative amplitude if self.mismatchCount <= self.thresholdDistance
		# result is a binary list; 00 => mismatchCount == thresholdDistance, 10 => mismatchCount > thresholdDistance, 01 => mismatchCount < thresholdDistance
		# binaryComparisonResult = Util.compareBinaryNumbers(self.mismatchCount, self.thresholdDistance) 
		isGreaterThan = BinaryNumber(self.mismatchesCountRegister).isGreaterThan(BinaryNumber(self.thresholdRegister))
		if not isGreaterThan:
			self.amplitude = -1 * self.amplitude

		# mismatchesCountBinaryNumber = BinaryNumber(self.mismatchesCountRegister)
		# thresholdDistanceBinaryNumber = BinaryNumber(self.thresholdRegister)
		# mismatchesCountIsGreaterThanThresholdDistance = False
		# for (mismatchCountBit, thresholdDistanceBit) in zip(mismatchCountBinaryNumber.value, thresholdDistanceBinaryNumber.value):
		# 	singeBitBinaryNumber1 = BinaryNumber([mismatchCountBit.state])
		# 	singeBitBinaryNumber2 = BinaryNumber([thresholdDistanceBit.state])
		# 	singeBitBinaryNumber1.singleBit()  
		# 	singeBitBinaryNumber2.singleBit()
		# 	if singeBitBinaryNumber1.isGreaterThan(singeBitBinaryNumber2): # compare in binary the value for mismatch count for this substring itno threshold distance
		# 		mismatchCountIsGreaterThanThresholdDistance = True
		# if not mismatchCountIsGreaterThanThresholdDistance:
		# 	self.amplitude = -1 * self.amplitude # mark this state with negative sign of amplitude if the substring it represents has mismatch count less or equal to threshold		

	def compareMismatchCountWithThresholdDistance(self):
		comparisonResult = BinaryNumber(self.mismatchesCountRegister).isGreaterThan(BinaryNumber(self.thresholdRegister))
		return comparisonResult

	def toString(self):
		text = 'Amplitude: ' + str(self.amplitude)
		text = text + ', Substring register: '
		for symbolRegister in self.substringRegister:
			text = text + str(symbolRegister.getState())
		text = text + ', Pattern register: '
		for symbolRegister in self.patternRegister:
			text = text + str(symbolRegister.getState())
		text = text + ', Mismatches indicator register: ' + str(self.mismatchesIndicatorRegister) + ', Mismatches count register: ' + str(self.mismatchesCountRegister) + ', Threshold register: ' + str(self.thresholdRegister)
		return text

# Already encorporated to method FullAdder.init() the construction of a 1-bit full adder. This class is already deprecated.
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

# Already added to Xcode version.
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

# Already added to Xcode version.
class AdditionResult(object):
	"""docstring for AdditionResult"""
	def __init__(self, sumBitState, carryBitState):
		super(AdditionResult, self).__init__()
		self.sumBitState = sumBitState
		self.carryBitState = carryBitState

class SymbolComparator(object):
	"""docstring for SymbolComparator"""
	def __init__(self):
		super(SymbolComparator, self).__init__()
		class bcolors:
				ONE = '\033[92m'
				ENDC = '\033[0m'
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

	def compare(self, symbol1BasisVector, symbol2BasisVector):
		basisVector = kron(symbol1BasisVector, kron(symbol2BasisVector, ket0()))
		basisVector = dot(self.matrix, basisVector).tolist()
		indexOf1 = basisVector.index(1)
		comparisonResult = 0
		if indexOf1 % 2 == 1 :
			comparisonResult = 1
		# print 'Comparison result:',comparisonResult
		return comparisonResult

# Already encorporated to method BitComparator.init() the construction of this operator. This class is then already deprecated.
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

# Added method BinaryNumber.isGreaterThan() in Xcode version. This class is then deprecated.
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

# Already added to Xcode version.
class BitComparator(object):
	"""docstring for BitComparator"""
	def __init__(self):
		super(BitComparator, self).__init__()
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
		# print self.matrix
		
	def compare(self, bit1, bit2):
		basisVector = kron(bit1.getBasisVector(), kron(bit2.getBasisVector(), kron(ket0(), ket0()) ))
		# print 'BitComparator.compare: bit1.state:', bit1.state, 'bit1.basisVector:', bit1.basisVector, ' bit2.state:', bit2.state, ' bit2.basisVector:', bit2.basisVector
		# print 'BitComparator.compare: initial basisVector:', basisVector
		basisVector = dot(self.matrix, basisVector).tolist()
		indexOf1 = basisVector.index(1)
		# print 'BitComparator.compare: resulting basisVector:', basisVector

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

class DiffusionOperator(object):
 	"""docstring for DiffusionOperator"""
 	def __init__(self, dimension):
 		super(DiffusionOperator, self).__init__()
 		# rows = []
 		# for rowIndex in xrange(0,dimension):
 		# 	rows.append(ones(dimension,dtype=int).tolist())
 		# print 'DiffusionOperator: init: matrix:', str(matrix(rows))
 		projector = 2.0 * ((1.0/dimension) * ones((dimension,dimension),dtype=int))
 		self.matrix = projector - identity(dimension, dtype=float)
 		factor = 2.0 * (1.0/dimension)
 		print 'DiffusionOperator: init: factor:', factor
 		print 'DiffusionOperator: init: projector:', projector
 		print 'DiffusionOperator: init: self.matrix:', self.matrix

 	def inverseAboutAverage(self, solutionRegister):
 		amplitudeList = []

 		# Construct a list of amplitudes of all verificationPhaseRegister in solutionRegister
 		solutionCount = 0
 		for verificationPhaseRegister in solutionRegister:
 			if verificationPhaseRegister.getAmplitude() < 0:
 				solutionCount = solutionCount + 1
 			amplitudeList.append(verificationPhaseRegister.getAmplitude())
 		# Print the constructed amplitude list	
 		print 'DiffusionOperator: before inverseAboutAverage: amplitudeList:', amplitudeList

 		# Perform inversion about the mean operation on the list of amplitudes
 		numberOfIteration = int((pi/4.0) * sqrt(len(solutionRegister)/solutionCount))
 		for iterationCount in xrange(0,numberOfIteration):
			amplitudeList = dot(self.matrix, amplitudeList)
		# Print the updated amplitude list
		print 'DiffusionOperator: after inverseAboutAverage: amplitudeList:', amplitudeList
		
		# Update the amplitude of each verificationPhaseRegister in the solutionRegister
		for (index, verificationPhaseRegister) in enumerate(solutionRegister):
 			solutionRegister[index].setAmplitude(amplitudeList[index])
 		
 		# Print the probability for each verificationPhaseRegister
 		probabilities = []
 		totalProbability = 0.0
 		for amplitude in amplitudeList:
 			probabilities.append(amplitude * amplitude)
 			totalProbability = totalProbability + (amplitude * amplitude)
 		print 'DiffusionOperator: Probabilities:', probabilities
 		print 'DiffusionOperator: Total probability:', totalProbability

class ULocOperator(object):
	"""Unitary operator for fetching location of first occurrence of a symbol text[textIndex] in pattern pattern."""
	def __init__(self, text, pattern, symbolSet, locationSet):
		super(ULocOperator, self).__init__()
		textLength = len(text) # should have referred to global variable N instead
		patternLength = len(pattern) # should have referred to global variable M instead
		# matrixDimension = textLength * patternLength # (N - M + 1) * M
		matrixDimension = pow(2, registerBitCount) * pow(2, registerBitCount) # (N - M + 1) * M
		self.matrix = zeros((matrixDimension, matrixDimension), dtype=int)
		for textIndex in xrange(0, textLength):
			if text[textIndex] in symbolSet:
				# print "text["+str(textIndex)+"]: ", text[textIndex], " available in pattern_alphabet and is at index " + str(symbolSet.index(text[textIndex])) + "."
				symbolSetIndex = symbolSet.index(text[textIndex]) # index of symbol text[textIndex] in list pattern_alphabet
				locationSetIndex = locationSet[symbolSetIndex] # since pattern_alphabet and pattern_symbol_locations have same number of elements then the location of first occurrence in pattern of symbol text[textIndex] is at location pattern_symbol_locations[symbolSetIndex]
				self.matrix[textIndex * textLength][(textIndex * textLength) + locationSetIndex] = 1
				self.matrix[(textIndex * textLength) + locationSetIndex][textIndex * textLength] = 1
				for patternIndex in xrange(0, patternLength):
					if patternIndex != 0 and patternIndex != locationSetIndex:
						self.matrix[(textIndex * textLength) + patternIndex][(textIndex * textLength) + patternIndex] = 1
				for index in xrange(patternLength, pow(2,registerBitCount)):
					self.matrix[(textIndex * textLength) + index][(textIndex * textLength) + index] = 1
			else:
				# print "text["+str(textIndex)+"]: ", text[textIndex], " not available in P_sym."
				for index in xrange(0, pow(2,registerBitCount)):
					self.matrix[(textIndex * textLength) + index][(textIndex * textLength) + index] = 1 # just an identity operator

		# print 'ULocOperator:'
		# print self.matrix

	# Apply operator to basis vector kron(indexRegisterBasisVector, startingLocationRegisterBasisVector)
	def identifyLocationOfFirstOccurrence(self, basisVector):
		# print 'ULocOperator: before identifyLocationOfFirstOccurrence:', basisVector
		resultingBasisVector = dot(self.matrix, basisVector).tolist()
		indexOf1 = resultingBasisVector.index(1)
		if indexOf1 % 2 == 0:
			return ket0()
		else:
			return ket1()

		# if indexOf1 in resultingBasisVector[0:(pow(2,registerBitCount)*2)]:
		# 	if indexOf1 in resultingBasisVector[0:pow(2,registerBitCount)]:
		# 		print '1range[0:',(pow(2,registerBitCount)-1),']'
		# 		return ket0()
		# 	else:
		# 		print '2range[', pow(2,registerBitCount), ':', (pow(2,registerBitCount)*2-1),']'
		# 		return ket1()
		# else:
		# 	if indexOf1 in resultingBasisVector[(pow(2,registerBitCount)*2):(pow(2,registerBitCount)*3)]:
		# 		print '3range[', (pow(2,registerBitCount)*2), ':', (pow(2,registerBitCount)*3-1), ']'
		# 		return ket0()
		# 	else:
		# 		print '4range[', (pow(2,registerBitCount)*3), ':', (pow(2,registerBitCount)*4-1),']'
		# 		return ket1()
		# print 'ULocOperator: after identifyLocationOfFirstOccurrence:', resultingBasisVector

# ==================================================================== METHODS ====================================================================

def prepareIndexRegister(index, dimension):
	return convertToBinaryStringFromDecimalNumber(index, dimension)

def prepareStartingLocationRegister(location, dimension):
	return convertToBinaryStringFromDecimalNumber(location, dimension)

def prepareSymbolRegister(symbol):
	binaryCode = symbol.getBinaryCode()
	bits = []
	for bit in binaryCode:
		bits.append(Bit(bit))
	return SymbolRegister(symbol, bits)

def prepareZerosRegister(bitCount):
	return zeros(bitCount, dtype=int).tolist()

def prepareThresholdRegister(threshold, dimension):
	return convertToBinaryNumberFromDecimalNumber(threshold, dimension)

# Convert a decimal number into a binary string (list of 0s and 1s)
def convertToBinaryStringFromDecimalNumber(decimalNumber, dimension):
	binaryString = []
	for x in xrange(0, dimension):
		i = dimension - x - 1
		if decimalNumber >= pow(2,i):
			binaryString.append(1)
			decimalNumber = decimalNumber - pow(2,i)
		else:
			binaryString.append(0)
	binaryString.reverse()
	return binaryString

def convertToBinaryStringFromBasisVector(basisVector, dimension):
	indexOf1 = basisVector.index(1)
	return convertToBinaryStringFromDecimalNumber(indexOf1, dimension)

# Convert a BinaryNumber object into a decimal number
def convertToBasisVectorFromBinaryNumber(binaryNumber):
	decimalNumber = 0
	# print 'convertToBasisVectorFromBinaryNumber: binaryNumber:',binaryNumber.toString()
	for bitIndex in xrange(0,registerBitCount):
		decimalNumber = decimalNumber + pow(2,bitIndex)*binaryNumber.value[bitIndex].state
	return convertToBasisVectorFromDecimalNumber(decimalNumber, pow(2,registerBitCount))

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

def displayFinalProbabilities(solutionRegister):
	amplitudeList = []
	# Construct a list of amplitudes of all verificationPhaseRegister in solutionRegister
	for verificationPhaseRegister in solutionRegister:
		amplitudeList.append(verificationPhaseRegister.getAmplitude())

	plotter.title('Probability of Occurence for Convolute')
	plotter.ylabel('Probability')
	plotter.xlabel('Basis state index')
	x = range(len(solutionRegister))
	width = 0.15
	plotter.bar(x, amplitudeList, width, color="black", alpha=0.75)
	plotter.grid(True)
	plotter.show()

def zeroBits(bitCount):
	bits = []
	for i in xrange(0,bitCount):
		bits.append(Bit(0))
	return bits

def ket0():
	return array([1,0])

def ket1():
	return array([0,1])

# Construct list of distinct symbols in pattern pattern
def constructP_Sym(pattern):
# Given pattern, setify pattern to get all distinct symbols in pattern convert the set into a list by
# listifying it. Sortification is not necessary actually.
	return list(set(pattern))

# Construct a map with location of first occurrence of symbol in pattern_alphabet as key and the symbol value
def constructP_Loc(symbolSet):
# there are many ways to construct pattern_symbol_locations using Python's list comprehension feature but
# using map we explicitly say that there is a one-to-one and onto mapping between 
# pattern_alphabet and pattern_symbol_locations

	# a convenience constructor for map where first parameter
	# is method for determining location of first occurrence as 
	# key and second parameter is method for the value
	return map(getFirstLoc,symbolSet) 

# Get the location of first occurrence of symbol x in pattern pattern
def getFirstLoc(x):
	return pattern.index(x)

# ==================================================================== MAIN ====================================================================

text = 'acgc'
 # + 'aaaaaaaaaa' + 'aaaaaaaaaa' + 'aaaaaaaaaa' + 'aaaaaaaaaa' + 'aaaaaaaaaa' + 'aaaaaaaaaa' + 'aaaaaaaaaa' + 'aaaaaaaaaa' + 'aaaaaaaaaa' + 'a'
pattern = 'gc'
N = len(text)
M = len(pattern)
threshold = 1
logging.info('Input text: '+text)
logging.info('Input pattern: '+pattern)


# FILTERING PHASE ==================================================
# Construct set of distinct symbol in pattern, pattern_alphabet
logging.info('Constructing set pattern_alphabet.')
P_Sym = constructP_Sym(pattern)
print 'pattern_alphabet: ', P_Sym, ", q=", str(len(P_Sym))

# Construct set of locations of first occurrence in pattern of
# all symbols in set pattern_alphabet
logging.info('Constructing set pattern_symbol_locations.')
P_Loc = constructP_Loc(P_Sym)
print 'pattern_symbol_locations: ', P_Loc

# Construct operator for identifying location of first occurrence of each symbol text[i] in pattern pattern
locationOperator = ULocOperator(text, pattern, P_Sym, P_Loc) # TODO: Prepare ULocOperator in the initialization phase instead

logging.info('Preparing initial superposition state.')
# Prepare initial superposition state as a list.
searchSpaceRegister = []
for index in xrange(0, N):
	searchSpaceRegister.append(FilteringPhaseRegister( sqrt(1.0/N), 
														prepareIndexRegister(index, registerBitCount), # int(ceil(log(N-M+1,2))) + 1 
														prepareStartingLocationRegister(0, registerBitCount) # int(ceil(log(N-M+1,2))) + 1
													)
								)
logging.info('Initial superposition state prepared.') 
print 'FilteringPhaseRegisters: =============================='
for filteringPhaseRegister in searchSpaceRegister:
	print filteringPhaseRegister.toString()
print '======================================================='


logging.info('Identifying location of first occurrence of symbols in text text.')
# Identify location of first occurrence pattern_symbol_locations(i) for each symbol text[i] in text text
for filteringPhaseRegister in searchSpaceRegister:
	filteringPhaseRegister.identifyLocationOfFirstOccurrence()
logging.info('Location of first occurrence of symbols in text text identified.')
print 'FilteringPhaseRegisters: =============================='
for filteringPhaseRegister in searchSpaceRegister:
	print filteringPhaseRegister.toString()
print '======================================================='

logging.info('Identifying possible starting locations in text text.')
# Compute for possible starting location (i - pattern_symbol_locations(i))
for filteringPhaseRegister in searchSpaceRegister:
	filteringPhaseRegister.computeStartingLocation()
logging.info('Possible starting locations in text text identified.')
# print 'FilteringPhaseRegisters: =============================='
# for filteringPhaseRegister in searchSpaceRegister:
# 	print filteringPhaseRegister.toString()
# print '======================================================='
#===== EDITED AS OF JANUARY 22, 2016 =====

# VERIFICATION PHASE ===============================================
# Create a list of M-length substrings of text.
# substrings = []
# for startingIndex in xrange(0, N-M+1):
# 	substrings.append(text[startingIndex : startingIndex + M])
# print substrings


# solutionRegister = []
# for substring in substrings:

# 	substringSymbolRegisters = []
# 	for character in substring:
# 		substringSymbolRegisters.append(prepareSymbolRegister(Symbol(character)))
# 	# for substringSymbolRegister in substringSymbolRegisters:
# 	# 	print substringSymbolRegister.getState(),
# 	# print

# 	patternSymbolRegisters = []
# 	for character in pattern:
# 		patternSymbolRegisters.append(prepareSymbolRegister(Symbol(character)))
# 	# for patternSymbolRegister in patternSymbolRegisters:
# 	# 	print patternSymbolRegister.getState(),
# 	# print


# 	mismatchesIndicatorRegister = prepareZerosRegister(M)
# 	mismatchesCountRegister = prepareZerosRegister(int(ceil(log(M,2))) + 1)
# 	thresholdRegister = prepareThresholdRegister(threshold, int(ceil(log(M,2))) + 1)

# 	verificationPhaseRegister = VerificationPhaseRegister(sqrt(1.0/(N-M+1)), substringSymbolRegisters, patternSymbolRegisters, mismatchesIndicatorRegister, mismatchesCountRegister, thresholdRegister)
# 	print 'Verification phase register:', verificationPhaseRegister.toString()

# 	solutionRegister.append(verificationPhaseRegister)


# for (index,verificationPhaseRegister) in enumerate(solutionRegister):
# 	print '=============================== Candidate solution substring', index, ' ==============================='
# 	# identify mismatches
# 	for substringSymbolRegister, patternSymbolRegister, (index, mismatchesIndicatorRegisterBit) in zip(verificationPhaseRegister.getSubstringRegister(), verificationPhaseRegister.getPatternRegister(), enumerate(verificationPhaseRegister.getMismatchesIndicatorRegister())):
# 		verificationPhaseRegister.getMismatchesIndicatorRegister()[index] = substringSymbolRegister.compare(patternSymbolRegister)
# 	# print 'Verification phase register:', verificationPhaseRegister.toString()

# 	# count number of mismatches
# 	# for verificationPhaseRegister in verificationPhaseRegisters:
# 	verificationPhaseRegister.countMismatches()
# 	# print 'Verification phase register:', verificationPhaseRegister.toString()

# 	verificationPhaseRegister.markIfMismatchesCountIsValid()
# 	print 'Verification phase register:', verificationPhaseRegister.toString()

# # Perform amplitude amplification.
# print '=============================== Amplitude amplification ==============================='
# DiffusionOperator(len(solutionRegister)).inverseAboutAverage(solutionRegister)
# for (index,verificationPhaseRegister) in enumerate(solutionRegister):
# 	print 'Verification phase register:', verificationPhaseRegister.toString()
# print '======================================================================================='

# # Plot into a bar graph the probabilities
# print '=============================== Probabilities ========================================='
# displayFinalProbabilities(solutionRegister)
# print '======================================================================================='

# # Perform measurement.
# # for iterationCount in xrange(0,maxInterationCount):
# # 	solutionState = measureSolutionRegisterState(solutionRegister)
# # 	print 'Solution state:', solutionState















