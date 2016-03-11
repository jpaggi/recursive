from scipy.signal import butter, lfilter, freqz, detrend
import matplotlib.pyplot as plt
import sys


def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y


for line in sys.stdin:
	chrom, start, end, offsets, ratchets, strand, expression = line.strip().split('\t')
	expression = map(int, expression.split(','))
	expression = detrend(expression)
	print '\t'.join([chrom, start, end, offsets, ratchets, strand, ','.join(map(str, expression))])