import scipy.io
import scipy.signal
import matplotlib.pyplot
from matplotlib.pyplot import plot, hold, show
import numpy
import math
import pylab

fs = 30000.
nyq = fs/2.
cutoff = 5000.
filterorder = 5
lpf_b,lpf_a = scipy.signal.filter_design.butter(filterorder,cutoff/nyq)

n = 10001
hpf_b = scipy.signal.firwin(n, cutoff = (20. / fs) , window = "hanning", pass_zero=False)

signals = []
filteredsignals = []
filteredsignals1 = []
filteredsignals2 = []
r = range(6,11)

for p in r:
	print(p)
	fullv = scipy.io.loadmat('/media/psf/HOME/Desktop/G001/Expt2.p' + p.__str__() + 'FullV.mat');
	signal = fullv['FullV']['V'];
	signal = signal[0][0][0][0:100000];
	signals.append(signal);


averagesignal = numpy.mean(numpy.transpose(signals),1)

for p in range(signals.__len__()):
	print(p)
	#filteredsignals.append(numpy.array(signals[p]) - numpy.array(averagesignal))
	filteredsignals1.append(scipy.signal.lfilter(lpf_b,lpf_a,numpy.array(signals[p]) - numpy.array(averagesignal)))
	filteredsignals2.append(scipy.signal.lfilter(hpf_b, [1.0], (signals[p] - numpy.array(averagesignal))))
	filteredsignals.append(scipy.signal.lfilter(hpf_b, [1.0], scipy.signal.lfilter(lpf_b,lpf_a, numpy.array(signals[p]) - numpy.array(averagesignal))))



pxx, freqs = pylab.psd(signals[4], Fs=30000 )
decibels = 10*numpy.log10(pxx)

pxx2, freqs2 = pylab.psd(filteredsignals2[4], Fs=30000 )
decibels2 = 10*numpy.log10(pxx2)

pylab.subplot(2,1,1)
hold(True)
pylab.plot(range(0,1000),signals[4][10000:11000],'o-') # plot a few samples
pylab.plot(range(0,1000),filteredsignals2[4][10000:11000],'o-') # plot a few samples

pylab.xlabel('time (sec)')
pylab.ylabel('Volts')
pylab.subplot(2,1,2)
pylab.plot(freqs, decibels, 'o-')
hold(True)
pylab.plot(freqs2, decibels2, 'o-')
pylab.xlabel('frequency')
pylab.ylabel('Power (decibels)')
pylab.savefig('example4.png',dpi=72)
pylab.show()

#plot(signal[1:100000],'b')
hold(True)
plot(filteredsignals[4][1:100000], 'g--')
plot(filteredsignals1[4][1:100000], 'r')
hold(True)
plot(filteredsignals2[4][1:100000], 'm')
#matplotlib.pyplot.plot(filteredsignal[1:10000])
#matplotlib.pyplot.plot((signal-filteredsignal)[1:10000])
matplotlib.pyplot.show()

