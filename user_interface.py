###########################################################
# Author: Jeffrey A. Aborot                               #
# Institution: Algorithms & Complexity Laboratory         #
#              Department of Computer Science             #
#              University of the Philippines Diliman      #
###########################################################

############################## IMPORTS ############################
from numpy import zeros, sum
from math import sqrt
from QConvolute import StringStateEncoding, q_convolute

###################################################################

############################## INPUT ##############################
# Assume an alphabet.
alphabet = ['a', 'c', 't', 'g']

# Ask input text and pattern from user.
text = raw_input("Text: ")
print "T:", text
pattern = raw_input("Pattern: ")
print "P:", pattern
pattern = "".join(reversed(pattern))
L = len(text) + len(pattern) - 1

# Ask number of iterations of Algorithm B from user.
iteration_count = int(raw_input("Number of iterations:"))
print 'Number of iterations:', str(iteration_count)
###################################################################

###################### PROCESS INPUT ##############################
# Create quantum states for binary indicator sequences of text and pattern.
text_quantum_states = []  # Array of quantum states corresponding to binary indicator sequences for text
pattern_quantum_states = []  # Array of quantum states corresponding to binary indicator sequences for pattern
for (symbol_index, symbol) in enumerate(alphabet):
    # Create quantum states only for symbols present in both text and pattern.
    if text.count(symbol) > 0 and pattern.count(symbol) > 0:
        # Create quantum states for text.
        text_quantum_state = zeros(L)
        text_symbol_count = text.count(symbol)
        for (text_index, text_symbol) in enumerate(text):
            if text_symbol == symbol:
                text_quantum_state[text_index] = sqrt(1.0 / text_symbol_count)
        text_quantum_states.append(StringStateEncoding(symbol, text_quantum_state))
        print 'Text quantum state encoding >>', text_quantum_states[symbol_index].to_string()

        # Create quantum states for pattern.
        pattern_quantum_state = zeros(L)
        pattern_symbol_count = pattern.count(symbol)
        for (pattern_index, pattern_symbol) in enumerate(pattern):
            if pattern_symbol == symbol:
                pattern_quantum_state[pattern_index] = sqrt(1.0 / pattern_symbol_count)
        pattern_quantum_states.append(StringStateEncoding(symbol, pattern_quantum_state))
        print 'Pattern quantum state >>', pattern_quantum_states[symbol_index].to_string()
###################################################################


############################## PROCESS #############################
# For each pair of binary indicator sequences T_symbol and P_symbol
def save_probability_graph_to_file(symbol, quantum_convolution):
    # write probability graph of quantum convolution to file
    probabilities = [(pow(value.real, 2) + pow(value.imag, 2)) for value in quantum_convolution]
    print 'Probabilities for symbol', symbol, ':', probabilities
    print 'Sum of probabilities:', sum(probabilities)


for i in xrange(0, len(text_quantum_states)):
    # compute for quantum convolution of quantum states corresponding to pairs of binary indicator sequences
    quantum_convolution = q_convolute(text_quantum_states[i].get_quantum_state_encoding(),
                                      pattern_quantum_states[i].get_quantum_state_encoding())
    save_probability_graph_to_file(text_quantum_states[i].get_symbol(), quantum_convolution)
###################################################################

############################## OUTPUT ##############################
# Output table of results.
