###########################################################
# Author: Jeffrey A. Aborot                               #
# Institution: Algorithms & Complexity Laboratory         #
#              Department of Computer Science             #
#              University of the Philippines Diliman      #
###########################################################

from numpy import zeros, empty, dot
from math import sqrt
import matplotlib.pyplot as plotter

from QOperators import qft, iqft, euclid_norm


def construct_v(a):
    """a function for constructing the unitary matrix operator V with diagonal elements as the normalized values of the qft of quantum state corresponding to a binary indicator sequence of the pattern

    :param a: an array corresponding to vector resulting from qft of a quantum state corresponding to a binary indicator sequence of the pattern
    :type a: ndarray of floats (n-dimensional array where n=1)
    :return: a unitary matrix with the normalized values of vector a as diagonal elements
    :rtype: ndarray
    """
    # TODO: construct operator V as tensor product of a decomposition of elementary unitary operators
    dim = a.shape[0]
    v = zeros(shape=(dim, dim), dtype=complex)
    for i in xrange(0, dim):
        if a[i] != 0.0 + 0.0j:
            norm = euclid_norm(a[i])
            v[i][i] = complex(a[i].real / norm, a[i].imag / norm)
        else:
            v[i][i] = complex(1.0, 0)
    return v


def constructV_Inv(a):
    dim = a.shape[0]
    b = zeros(shape=(dim, dim), dtype=complex)
    for i in xrange(0, dim):
        if a[i][i] != 0.0 + 0.0j:
            b[i][i] = a[i][i].conjugate()
    return b


def getProbabilities(a):
    dim = a.shape[0]
    b = empty(dim, dtype=float)
    for i in xrange(0, dim):
        b[i] = a[i].real ** 2 + a[i].imag ** 2

    return b


def getTotalProbabilities(a):
    dim = a.shape[0]
    total = 0.0
    for i in xrange(0, dim):
        total = total + a[i]

    return total


def getDFTFromQFT(a, strln):
    dim = a.shape[0]
    factor = 1.0 / (sqrt(1.0 / dim) * sqrt(1.0 / strln))
    print 'DFT factor = ' + str(factor)
    b = empty(dim, dtype=complex)
    for i in xrange(0, dim):
        b[i] = complex(factor * a[i].real, factor * a[i].imag)

    return b


def getDFTT_dot_DFTP_from_VQFTP(v, vqft_p, dft_t, Pstrln):
    dim = vqft_p.shape[0]
    qft_P_factor = (1.0 / (sqrt(1.0 / dim) * sqrt(1.0 / Pstrln)))
    b = empty(dim, dtype=complex)
    for i in xrange(0, dim):
        factor = ((1.0 / v[i][i]) * dft_t[i]) * qft_P_factor
        b[i] = complex(factor * vqft_p[i].real, factor * vqft_p[i].imag)

    return b


def getIDFT_from_IQFT(v, iqft_v, dft_t, dft_p, Pstrln):
    dim = iqft_v.shape[0]
    qft_P_factor = 1.0 / (sqrt(1.0 / dim) * sqrt(1.0 / Pstrln))
    iqft_factor = dim / (sqrt(1.0 / dim))
    print 'qft_P_factor = ' + str(qft_P_factor)
    print 'iqft_factor = ' + str(iqft_factor)
    b = empty(dim, dtype=complex)
    for j in xrange(0, dim):
        factor = 0.0
        for i in xrange(0, dim):
            factor = factor + iqft_factor * (1.0 / v[i][i]) * dft_t[i] * qft_P_factor
        b[j] = complex(factor * iqft_v[j].real, factor * iqft_v[j].imag)

    return b


def prettyprint(a):
    dim = a.shape[0]
    string = ''
    for i in xrange(0, dim):
        string = string + '|' + str(i) + '>: ' + str(a[i]) + '\n'

    return string


def q_convolute(text_state_encoding, pattern_state_encoding):
    """
    compute the quantum convolution of input quantum states corresponding to two real or complex sequences
    :param text_state_encoding: array of floats corresponding to a quantum state which encodes a binary indicator sequence of input text
    :param pattern_state_encoding: array of floats corresponding to a quantum state which encodes a binary indicator sequence of input pattern
    :type text_state_encoding: array of floats
    :type pattern_state_encoding: array of floats
    :return: an ndarray of floats corresponding to the quantum convolution of input quantum states
    :rtype: ndarray
    """
    qft_text = qft(text_state_encoding)
    qft_pattern = qft(pattern_state_encoding)
    v = construct_v(qft_text)
    v_qft_pattern = dot(v, qft_pattern)
    quantum_convolution = iqft(v_qft_pattern)

    return quantum_convolution


class StringStateEncoding(object):
    """A class for holding metadata of quantum state encoding of a string with respect to a symbol."""

    def __init__(self, symbol, quantum_state_encoding):
        # super(StringStateEncoding, self).__init__()
        self.symbol = symbol
        self.quantum_state_encoding = quantum_state_encoding

    def get_symbol(self):
        return self.symbol

    def get_quantum_state_encoding(self):
        return self.quantum_state_encoding

    def to_string(self):
        return 'Symbol:' + self.symbol + ' Quantum state encoding:' + str(self.quantum_state_encoding)
