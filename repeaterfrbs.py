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
	stored = filename.split('.')[0]+'.npy'
	try:
		subfall = np.load(stored)
	except FileNotFoundError:
		ar = pypulse.Archive(filename, prepare=True)
		wfall = ar.getData()
		try:
			subfall = driftrate.subsample(wfall, 32, wfall.shape[1]//8)
		except ValueError:
			subfall = driftrate.subsample(wfall, 32, wfall.shape[1]//2)
		np.save(stored, subfall)

	ts = np.nanmean(subfall, axis=0)
	pkidx = np.nanargmax(ts)
	# ar.getTbin(), ar.getTimeUnit()

	return subfall, pkidx

if __name__ == '__main__':
	loadpsrfits('data/oostrum2020/R1_frb121102/R1_B01.rf')

