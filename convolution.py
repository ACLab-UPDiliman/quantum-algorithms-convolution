from numpy import empty, double, convolve
from math import cos, sin, pi, sqrt
from matplotlib import pyplot
import os


def dft(x):
    X = empty(len(x), dtype=complex)
    for i in xrange(0, len(x)): # 0 .. N-1
        print "--------------------"
        cycles_per_sample = double(i) / len(x) # 0/N, 1/N, ..., N-1/N
        print "i:",i," cycles per sample:",i,"/",len(x),"=",cycles_per_sample
        cos_correlation = 0
        sin_correlation = 0
        for (index,value) in enumerate(x): # values of input x
            cos_correlation = cos_correlation + (value * cos(2 * pi * index * cycles_per_sample))
            sin_correlation = sin_correlation + (value * sin(2 * pi * index * cycles_per_sample))
            print "index:",index," value:",value," cos:",cos_correlation," sin:",sin_correlation
        X[i] = complex(cos_correlation, sin_correlation)
        print "--------------------"
    X /= sqrt(len(x))
    return X


def idft(X):
    x = empty(len(X), dtype=complex)
    for i in xrange(0, len(X)): # 0 .. N-1
        print "--------------------"
        cycles_per_sample = double(i) / len(X) # 0/N, 1/N, ..., N-1/N
        print "i:",i," cycles per sample:",i,"/",len(X),"=",cycles_per_sample
        cos_correlation = 0
        sin_correlation = 0
        for (index,value) in enumerate(X): # values of input x
            cos_correlation = cos_correlation + ((value.real * cos(2 * pi * index * cycles_per_sample)) + (value.imag * sin(2 * pi * index * cycles_per_sample)))
            sin_correlation = sin_correlation - ((value.real * sin(2 * pi * index * cycles_per_sample)) + (value.imag * cos(2 * pi * index * cycles_per_sample)))
            print "index:",index," value:",value," cos:",cos_correlation," sin:",sin_correlation
        x[i] = complex(cos_correlation, sin_correlation)
        print "--------------------"
    x /= sqrt(len(X))
    return x


def normalize(a):
    a_norm = empty(len(a),dtype=complex)
    for (index,complex_number) in enumerate(a):
        norm = sqrt(complex_number.real**2 + complex_number.imag**2)
        a_norm[index] = complex(complex_number.real/norm, complex_number.imag/norm)
    return a_norm


def save_file(path, ext, close, verbose):
    # extract the directory and filename from path
    directory = os.path.split(path)[0]
    filename = "%s.%s" % (os.path.split(path)[1], ext)

    # make directory to user home directory
    if directory == '':
        directory = '.'

    # create directory if not yet existing
    if not os.path.exists(directory):
        os.mkdir(directory)

    # set final path to save to
    savepath = os.path.join(directory,filename)

    if verbose:
        print "Saving file to ...",savepath,

    # save the plot to file
    pyplot.savefig(savepath)

    if close:
        pyplot.close()

    if verbose:
        print '. Done'


text = ['a','a','a','b','a','a']
x = [1/sqrt(2), 0, 0, 0, 1/sqrt(2), 0, 0] # 1,0,0,0,1,0,0
print "x:",x," len:",len(x)
X = dft(x)
print "X:",X
X_normalized = normalize(X)
print "normalized X:",X_normalized


# plot magnitude for X and X_normalized
X_magnitudes = empty(len(X),dtype=double)
X_normalized_magnitudes = empty(len(X_normalized),dtype=double)
for index in xrange(0,len(X_magnitudes)):
    X_magnitudes[index] = sqrt(X[index].real**2 + X[index].imag**2)
    X_normalized_magnitudes[index] = sqrt(X_normalized[index].real**2 + X_normalized[index].imag**2)
width = 0.10
subplot = pyplot.subplot(211)
rects1 = subplot.bar([left-0.10 for left in xrange(0,len(X_magnitudes))], X_magnitudes, width, color="blue", alpha=1)
rects2 = subplot.bar([left for left in xrange(0,len(X_normalized_magnitudes))], X_normalized_magnitudes, width, color="blue", alpha=0.5)
subplot.legend((rects1[0],rects2[0]),('X','X normalized'))
subplot.set_title('magnitudes: X vs X normalized')
subplot.set_ylabel('magnitude')
subplot.set_xlabel('index')
pyplot.grid(True)
#pyplot.show()
#saveFile('Magnitudes for X vs X normalized','png',True,True)


# plot cosinusoid and sinusoid correlation values for X and X_normalized
X_cosinusoid_correlations = empty(len(X),dtype=double)
X_normalized_cosinusoid_correlations = empty(len(X_normalized),dtype=double)
X_sinusoid_correlations = empty(len(X),dtype=double)
X_normalized_sinusoid_correlations = empty(len(X_normalized),dtype=double)
for index in xrange(0,len(X),1):
    X_cosinusoid_correlations[index] = X[index].real
    X_normalized_cosinusoid_correlations[index] = X_normalized[index].real
    X_sinusoid_correlations[index] = X[index].imag
    X_normalized_sinusoid_correlations[index] = X_normalized[index].imag
width = 0.10
subplot = pyplot.subplot(212)
rects1 = subplot.bar([left-0.2 for left in range(0,len(X),1)], X_cosinusoid_correlations, width, color="blue", alpha=1)
rects2 = subplot.bar([left-0.1 for left in range(0,len(X),1)], X_normalized_cosinusoid_correlations, width, color="blue", alpha=0.5)
rects3 = subplot.bar([left+0.1 for left in range(0,len(X),1)], X_sinusoid_correlations, width, color="red", alpha=1)
rects4 = subplot.bar([left+0.2 for left in range(0,len(X),1)], X_normalized_sinusoid_correlations, width, color="red", alpha=0.5)
subplot.legend((rects1[0],rects2[0],rects3[0],rects4[0]),('X Cos','X normed Cos','X Sin','X normed Sin'))
subplot.set_title('correlations: X vs X normalized')
subplot.set_ylabel('correlation value')
subplot.set_xlabel('index')
#subplot.axis([-1,8,0,1])
pyplot.gcf().set_figheight(15)
pyplot.gcf().set_figwidth(15)
#pyplot.show()
save_file('images/Correlation values for X vs X normalized', 'png', True, True)


# we compute the original sequence x from IDFT of normalized X
# we expect to get a sequence different from x as effect of normalization
x_normalized = idft(X_normalized)
print 'x_computed:',x_normalized
x_normalized_cosinusoid_correlations = empty(len(x_normalized), dtype=double)
for (index,value) in enumerate(x_normalized):
    x_normalized_cosinusoid_correlations[index] = value.real
x_normalized_sinusoid_correlations = empty(len(x_normalized), dtype=double)
for (index,value) in enumerate(x_normalized):
    x_normalized_sinusoid_correlations[index] = value.imag
width = 0.10
subplot = pyplot.subplot(111)
rects1 = subplot.bar([left-0.20 for left in xrange(0,len(x),1)], x, width, color="blue", alpha=1)
rects2 = subplot.bar([left-0.10 for left in range(0,len(x_normalized_cosinusoid_correlations),1)], x_normalized_cosinusoid_correlations, width, color="green", alpha=1)
rects3 = subplot.bar([left for left in range(0,len(x_normalized_sinusoid_correlations),1)], x_normalized_sinusoid_correlations, width, color="red", alpha=1)
subplot.legend((rects1[0],rects2[0],rects3[0]),('x','x norm cos','x norm sin'))
pyplot.title('x vs x norm')
pyplot.ylabel('value')
pyplot.xlabel('index')
#pyplot.show()
save_file('images/x vs x normalized', 'png', True, True)


# We check the effect of DFT to a pattern.
# - The assumption is that even though
#   X undergoes a normalization process, the cosinusoidal and/or sinusoidal frequencies
#   which have high positive correlation with the input sequence x (in frequency domain) will still
#   have high positive correlation even after the normalization process.
# - Even though some frequencies in the frequency domain will have higher correlation
#   values (+ or -) as result of normalization, if the input pattern sequence y do not
#   have high cosinusoidal and sinusoidal correlation values in those frequencies, the
#   result of X point-wise-dot-product-multiplication Y in those frequencies will still be lower as
#   compared to that of frequencies high-correlation-value-in-x-but-normalized and
#   high-correlation-value-in-y
pattern = ['a','b','a']
y = [1/sqrt(2),0,1/sqrt(2),0,0,0,0] # 1,0,1
Y = dft(y)

# Plot of magnitude spectrum of DFT of y.
Y_cosinusoid_correlations = empty(len(Y),dtype=double)
Y_sinusoid_correlations = empty(len(Y),dtype=double)
for (index,value) in enumerate(Y):
    Y_cosinusoid_correlations[index] = value.real
    Y_sinusoid_correlations[index] = value.imag
width = 0.10
subplot = pyplot.subplot(111)
rects1 = subplot.bar([left-0.05 for left in range(0,len(Y),1)],Y_cosinusoid_correlations,width, color="blue", alpha=1)
rects2 = subplot.bar([left+0.05 for left in range(0,len(Y),1)],Y_sinusoid_correlations,width, color="red", alpha=1)
subplot.legend((rects1[0],rects2[0]),('Y cos','Y sin'))
pyplot.title('Y = DFT(y)')
pyplot.ylabel('correlation value')
pyplot.xlabel('index')
#pyplot.show()
save_file('images/Y', 'png', True, True)


# We check the correlation between the cosinusoidal and sinusoidal correlation values
# of x and y expressed in frequency domain X and Y. We do this by applying point-wise
# multiplication to each frequency in X and Y.
Z = empty(len(X),dtype=complex)
for index in xrange(0,len(X)):
    Z[index] = X[index] * Y[index]
Z_cosinusoid_correlations = empty(len(Z),dtype=double)
Z_sinusoid_correlations = empty(len(Z),dtype=double)
for (index,value) in enumerate(Z):
    Z_cosinusoid_correlations[index] = value.real
    Z_sinusoid_correlations[index] = value.imag

Z_normalized = empty(len(X),dtype=complex)
for index in xrange(0,len(X_normalized)):
    Z_normalized[index] = X_normalized[index] * Y[index]
Z_normalized_cosinusoid_correlations = empty(len(Z),dtype=double)
Z_normalized_sinusoid_correlations = empty(len(Z),dtype=double)
for (index,value) in enumerate(Z_normalized):
    Z_normalized_cosinusoid_correlations[index] = value.real
    Z_normalized_sinusoid_correlations[index] = value.imag


subplot = pyplot.subplot(111)
rects1 = subplot.bar([left-0.20 for left in range(0,len(Z),1)],Z_cosinusoid_correlations,width, color="blue", alpha=1)
rects2 = subplot.bar([left-0.10 for left in range(0,len(Z),1)],Z_sinusoid_correlations,width, color="red", alpha=1)
rects3 = subplot.bar([left+0.10 for left in range(0,len(Z),1)],Z_normalized_cosinusoid_correlations,width, color="blue", alpha=.5)
rects4 = subplot.bar([left+0.20 for left in range(0,len(Z),1)],Z_normalized_sinusoid_correlations,width, color="red", alpha=.5)
subplot.legend((rects1[0],rects2[0],rects3[0],rects4[0]),('Z cos','Z sin', 'Z cos normed','Z sin normed'))
pyplot.title('(X dot Y) VS (X normalized dot Y)')
pyplot.ylabel('correlation value')
pyplot.xlabel('index')
#pyplot.show()
save_file('images/Correlation values for Z vs Z normalized', 'png', True, True)


# We compare the magnitude of of each X_k with the their normalized counterpart.
Z_magnitude = empty(len(Z),dtype=double)
Z_normalized_magnitude = empty(len(Z),dtype=double)
Y_magnitude = empty(len(Y),dtype=double)
for index in xrange(0,len(Z)):
  Z_magnitude[index] = sqrt(Z[index].real**2 + Z[index].imag**2)
  Z_normalized_magnitude[index] = sqrt(Z_normalized[index].real**2 + Z_normalized[index].imag**2)
  Y_magnitude[index] = sqrt(Y[index].real**2 + Y[index].imag**2)
subplot = pyplot.subplot(111)
rects1 = subplot.bar([left-0.10 for left in range(0,len(Z),1)],Z_magnitude,width, color="blue", alpha=1)
rects2 = subplot.bar([left for left in range(0,len(Z),1)],Z_normalized_magnitude,width, color="blue", alpha=.5)
rects3 = subplot.bar([left+0.10 for left in range(0,len(Z),1)],Y_magnitude,width, color="red", alpha=.5)
subplot.legend((rects1[0],rects2[0],rects3[0]),('|Z_k|','|Z_k normalized|','|Y_k|'))
pyplot.title('|Z_k| VS |Z_k normalized| VS |Y_k|')
pyplot.ylabel('magnitude')
pyplot.xlabel('index')
#pyplot.show()
save_file('images/Magnitudes for Z vs Z normalized', 'png', True, True)


# We determine the sequence z with cosinusoidal and sinusoidal correlations corresponding
# to what we got from the dot multiplication of cosinusoidal and sinusoidal correlations in
# X and Y. The result of the dot product of the cosinusoidal and sinusoidal correlations in
# the frequency domain will determine the shape of the wave corresponding to sample sequence z.

# plot data set
z = idft(Z)
print  'z:',z
z_cosinusoid_correlations = empty(len(z),dtype=double)
z_sinusoid_correlations = empty(len(z),dtype=double)
for (index,value) in enumerate(z):
    z_cosinusoid_correlations[index] = value.real
    z_sinusoid_correlations[index] = value.imag

z_normalized = idft(Z_normalized)
print 'z normalized:',z_normalized
z_normalized_cosinusoid_correlations = empty(len(z),dtype=double)
z_normalized_sinusoid_correlations = empty(len(z),dtype=double)
for (index,value) in enumerate(z_normalized):
    z_normalized_cosinusoid_correlations[index] = value.real
    z_normalized_sinusoid_correlations[index] = value.imag

subplot = pyplot.subplot(111)
# plot legend
rects1 = subplot.bar([left-0.30 for left in range(0,len(z),1)],z_cosinusoid_correlations,width, color="blue", alpha=1)
rects2 = subplot.bar([left-0.20 for left in range(0,len(z),1)],z_sinusoid_correlations,width, color="red", alpha=1)
rects3 = subplot.bar([left-0.10 for left in range(0,len(z),1)],(z_cosinusoid_correlations**2),width, color="green", alpha=1)
rects4 = subplot.bar([left+0.10 for left in range(0,len(z),1)],z_normalized_cosinusoid_correlations,width, color="blue", alpha=.5)
rects5 = subplot.bar([left+0.20 for left in range(0,len(z),1)],z_normalized_sinusoid_correlations,width, color="red", alpha=.5)
rects6 = subplot.bar([left+0.30 for left in range(0,len(z),1)],(z_normalized_cosinusoid_correlations**2),width, color="green", alpha=.5)
subplot.legend((rects1[0],rects2[0],rects3[0],rects4[0],rects5[0],rects6[0]),('z Cos','z Sin','z Prob', 'z Cos normed', 'z Sin normed', 'z Prob normed'))
subplot.axis([-1,8,0,1])

# plot labels
pyplot.title('z VS z normalized')
pyplot.ylabel('correlation value')
pyplot.xlabel('index')
#pyplot.show()
save_file('images/z vs z normalized', 'png', True, True)

# probabilities
z_probability = empty(len(z),dtype=double)
z_normalized_probability = empty(len(z),dtype=double)
for (index,value) in enumerate(z):
    z_probability[index] = z[index].real**2
    z_normalized_probability[index] = z_normalized[index].real**2
print 'z probabilities:',z_probability
print 'normalized z probabilities:',z_normalized_probability

v = convolve(x,y,mode='full')
print 'v:',v

# sample change with tag commit
