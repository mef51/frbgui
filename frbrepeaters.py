"""
This file implements the following pipeline for processing dynamic spectra of frbs
using existing and available frb libraries whenever possible.
Pipeline:
	1. read [.fil, .fits, .ar] --> .npy
		* .ar files must be converted to PSRFITS. `psrconv` can be used for example.
	2. --> masking, burst finding/windowing, splitting, simple noise removal
	3. --> dm variations, processburst, parameter normalizing
	4. --> save, model values, figures
"""
import numpy as np
import matplotlib.pyplot as plt
import pypulse, your
import driftrate, driftlaw

def loadpsrfits(filename):
	ar = pypulse.Archive(filename, prepare=True)
	print('loaded')
	wfall, subfall = ar.getData(), None
	try:
		subfall = driftrate.subsample(wfall, 32, wfall.shape[1]//8)
	except ValueError:
		subfall = driftrate.subsample(wfall, 32, wfall.shape[1]//2)

	sstored = filename.split('.')[0]+'_sub'+'.npy'
	stored = filename.split('.')[0]+'.npy'
	# np.save(sstored, subfall)
	# np.save(stored, wfall)

	ts = np.nanmean(subfall, axis=0)
	pkidx = np.nanargmax(ts)
	# ar.getTbin(), ar.getTimeUnit()
	print("done loading")
	return subfall, pkidx, wfall

if __name__ == '__main__':
	filename = 'data/oostrum2020/R1_frb121102/R1_B07.rf'
	# subfall, pkidx, wfall = loadpsrfits(filename)
	subfall = np.load(filename.split('.')[0]+'_sub.npy')
	wfall   = np.load(filename.split('.')[0]+'.npy')
	pkidx   = np.nanargmax(np.nanmean(subfall, axis=0))

	ts  = np.nanmean(wfall, axis=0)
	snr = np.max(ts)/np.std(ts)

	test_sub = driftrate.subsample(wfall, 32, wfall.shape[1]//12)
	print(wfall.shape, test_sub.shape)
	# plt.imshow(subfall, origin='lower', aspect='auto')
	# plt.show()
