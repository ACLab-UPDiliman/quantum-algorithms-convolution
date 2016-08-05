from Util import Symbol

# Define here the binary encoding for the symbols of the input strings to be used in the program.
# Define an enum class for all symbols.

class Alphabet(object):
	"""docstring for Alphabet"""
	def __init__(self):
		super(Alphabet, self).__init__()
		self.symbols = [Util.Symbol('a'), Util.Symbol('c'), Util.Symbol('t'), Util.Symbol('g')]
		self.codes = [[0,0],[1,0],[0,1],[1,1]]
		self.codeBitCount = 2
		self.size = len(self.symbols)
	
	def getSymbols():
		return self.symbols

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