###########################################################
# Author: Jeffrey A. Aborot                               #
# Institution: Algorithms & Complexity Laboratory         #
#              Department of Computer Science             #
#              University of the Philippines Diliman      #
###########################################################

############################## IMPORTS ############################
from numpy import zeros, sum, arange, random, unique
from math import sqrt
from matplotlib import pyplot
from QConvolute import StringStateEncoding, q_convolute
import os
import csv


###################################################################

############################## FUNCTIONS ##########################
def generate_random_string(length=8, symbols=['a', 'c', 't', 'g']):
    return ''.join([random.choice(symbols) for n in xrange(length)])


# For each pair of binary indicator sequences T_symbol and P_symbol
def save_probability_graph_to_file(symbol, probability_distribution):
    # create bar graph of probabilities of indices
    width = 0.10
    subplot = pyplot.subplot(211)
    subplot.bar([left for left in range(0, len(probability_distribution), 1)], probability_distribution, width, color="blue",
                alpha=1)
    subplot.legend(['Pr(x)'])
    subplot.set_title('Probability distribution for symbol ' + symbol)
    subplot.set_ylabel('Pr(x)')
    subplot.set_xlabel('index')
    subplot.set_xticks(arange(len(probability_distribution)))
    pyplot.gcf().set_figheight(15)
    pyplot.gcf().set_figwidth(15)

    # save the bar graph to file
    # extract the directory and filename from path
    path = 'images/' + symbol + '-probabilities'
    ext = 'png'
    should_close = True
    verbose = True
    directory = os.path.split(path)[0]
    filename = "%s.%s" % (os.path.split(path)[1], ext)

    # make directory to user home directory
    if directory == '':
        directory = '.'

    # create directory if not yet existing
    if not os.path.exists(directory):
        os.mkdir(directory)

    # set final path to save to
    save_path = os.path.join(directory, filename)

    if verbose:
        print "Saving file to ...", save_path,

    # save the plot to file
    pyplot.savefig(save_path)

    if should_close:
        pyplot.close()

    if verbose:
        print '. Done'


def compute_probability_distribution(convolution):
    # compute for probabilities of the convolution
    probabilities = [(pow(value.real, 2) + pow(value.imag, 2)) for value in convolution]
    # print 'Probabilities for symbol', symbol, ':', probabilities
    # print 'Sum of probabilities:', sum(probabilities)

    return probabilities


def is_solution(index, _text, _pattern, distance):
    print 'Comparing text[', index, ':', (index + len(_pattern)), ']=', _text[index: index + len(
        _pattern)], ' and pattern=', _pattern
    comparison = [1 if t_symbol == p_symbol else 0 for t_symbol, p_symbol in
                  zip(_text[index: index + len(_pattern)], _pattern)]
    if comparison.count(0) > distance:
        return False
    return True

def write_to_csv(average_of_averages):
    with open('average_of_averages.csv', 'wb') as csvfile:
        csv_writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        csv_writer.writerow(['Text length', 'Average of averages of number of trials'])
        for (index, item) in enumerate(average_of_averages):
            csv_writer.writerow([(10*(index+1)), item])

############################## INPUT ##############################

############ Input alphabet, text and pattern configurations #############

# Assume an alphabet.
alphabet = ['a', 'c', 't', 'g']
# Assume fixed lengths of text i.e. power of 10.
lengths = [10*i for i in range(1, 51)]
print 'Text lengths:', lengths
# Assume pattern length.
pattern_length = 4
# Ask number of iterations of Algorithm B from user.
iteration_count = int(raw_input("Number of iterations per text length: "))
print 'Number of iterations:', iteration_count
# Ask for number of runs.
total_run_count = int(raw_input("Number of total runs: "))
print 'Number of total runs:', total_run_count

run_averages = []
for run in range(total_run_count):

    # array of average number of trials before success per text length
    average_trial_counts = []
    # Assume length of text.
    for text_length in lengths:
        print '================================================================================'
        print 'Current text length:', text_length
        trial_counts = []
        for iteration in range(iteration_count):
            print '--------------------------------------------------------------------------------'
            print '\n\nCurrent iteration count:', iteration, '\n'
            # Generate text and pattern.
            text = generate_random_string(text_length, alphabet)
            print 'text:', text
            # Randomly generate index of pattern in text.
            a = range(text_length - pattern_length + 1)
            # print 'a:', a
            pattern_index = random.choice(a, p = [1.0/len(a) for i in range(len(a))])
            print 'pattern_index:', pattern_index
            pattern = text[pattern_index : pattern_index + pattern_length]
            print 'pattern:', pattern

            #########################################################################

            L = text_length + pattern_length - 1
            # Reverse pattern.
            pattern_reversed = ''.join(reversed(pattern))
            ###################################################################

            ###################### PRE-PROCESS INPUT ##########################
            # Create quantum states for binary indicator sequences of text and pattern.
            text_quantum_states = []  # Array of quantum states corresponding to binary indicator sequences for text
            pattern_quantum_states = []  # Array of quantum states corresponding to binary indicator sequences for pattern
            for (symbol_index, symbol) in enumerate(alphabet):
                # Create quantum states only for symbols present in both text and pattern.
                if text.count(symbol) > 0 and pattern.count(symbol) > 0:
                    # TODO: Implement ESQUID algorithm for synthesizing quantum circuit for each quantum state.
                    # Create quantum states for text.
                    text_quantum_state = zeros(L)
                    text_symbol_count = text.count(symbol)
                    for (text_index, text_symbol) in enumerate(text):
                        if text_symbol == symbol:
                            text_quantum_state[text_index] = sqrt(1.0 / text_symbol_count)
                    text_quantum_states.append(StringStateEncoding(symbol, text_quantum_state))
                    # print 'Text quantum state encoding >>', text_quantum_states[-1].to_string()

                    # Create quantum states for pattern.
                    pattern_quantum_state = zeros(L)
                    pattern_symbol_count = pattern.count(symbol)
                    for (pattern_index, pattern_symbol) in enumerate(pattern_reversed):
                        if pattern_symbol == symbol:
                            pattern_quantum_state[pattern_index] = sqrt(1.0 / pattern_symbol_count)
                    pattern_quantum_states.append(StringStateEncoding(symbol, pattern_quantum_state))
                    # print 'Pattern quantum state >>', pattern_quantum_states[-1].to_string()
            ###################################################################


            ############################## PROCESS (COMPUTE FOR CONVOLUTION) #############################
            convolutions = []
            for i in xrange(0, len(text_quantum_states)):
                # compute for quantum convolution of quantum states corresponding to pairs of binary indicator sequences
                quantum_convolution = q_convolute(text_quantum_states[i].get_quantum_state_encoding(),
                                                  pattern_quantum_states[i].get_quantum_state_encoding())
                convolutions.append(quantum_convolution)

            no_solution = True
            trial_count = 0
            while (no_solution):
                trial_count += 1
                output_indices = []
                # for i in xrange(0, len(text_quantum_states)):
                for convolution in convolutions:
                    # compute for quantum convolution of quantum states corresponding to pairs of binary indicator sequences
                    # quantum_convolution = q_convolute(text_quantum_states[i].get_quantum_state_encoding(),
                    #                                   pattern_quantum_states[i].get_quantum_state_encoding())

                    # compute for probability distribution from quantum convolution
                    probability_distribution = compute_probability_distribution(convolution)

                    # save probability distribution to file
                    # save_probability_graph_to_file(text_quantum_states[i].get_symbol(), probability_distribution)

                    # sample the probability distribution
                    output_index = random.choice(arange(0, L), p=probability_distribution)
                    if (output_index > (pattern_length - 2) and output_index < text_length):
                        output_indices.append(output_index)

                # verify the output indices if they are solution
                verification_result = [True if is_solution(index - (pattern_length - 1), text, pattern, distance=1) else False for index in unique(output_indices)]
                print 'Output indices: ', output_indices
                print 'Adjusted indices: ', [i - (pattern_length - 1) for i in output_indices]
                print 'Verification result:', verification_result
                print 'verification_result.count(True):', verification_result.count(True)

                # stop iterating the algorithm if a solution has already been found
                if verification_result.count(True) > 0:
                    solution_index = verification_result.index(True)
                    no_solution = False
                    print 'Trial count:', trial_count
                    print 'Solution index:', (output_indices[solution_index] - (pattern_length - 1))
                    trial_counts.append(trial_count)

            print '--------------------------------------------------------------------------------'

        # compute for average number trials before success
        print 'Trial counts:', trial_counts
        average_trial_counts.append(float(sum(trial_counts)) / iteration_count)
        # write_to_csv(average_trial_counts)

        print '================================================================================'
        print 'Average number of trials:', average_trial_counts
    run_averages.append(average_trial_counts)
average_of_averages = [sum([item[j] for item in run_averages])/total_run_count for j in range(len(lengths))]
print 'Average of averages:', average_of_averages
write_to_csv(average_of_averages)
##############################################################################################


############################## OUTPUT ##############################
# Output table of results.
