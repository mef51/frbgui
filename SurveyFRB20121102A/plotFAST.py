#!/usr/bin/env python
import numpy as np
#import pyfits
import astropy.io.fits as pyfits
import os
import datetime
import time
import sys
from decimal import Decimal
import driftrate
secperday = 3600 * 24

if (len(sys.argv)<4):
  print( 'too few input parameters!')
  print( 'example:')
  print( 'python fitsio_cut_freq_time_splitpol.py startn endn startchan endchan FAST.fits')
  sys.exit()

startn=int(sys.argv[1])
endn=int(sys.argv[2])
startfreq=int(sys.argv[3])
endfreq=int(sys.argv[4])
filename=sys.argv[5]


hdulist = pyfits.open(filename)
print( 'hdu list length', len(hdulist) )
hdu0 = hdulist[0]
hdu1 = hdulist[1]
#hdu2 = hdulist[2]
data1 = hdu1.data['data']
#data2 = hdu2.data
header1 = hdu1.header
print( hdu0.header)
print( hdu1.header)
print( data1.shape)

fchannel = hdulist['SUBINT'].data[0]['DAT_FREQ']
fch1 = fchannel[0]
obsfreq = hdu0.header['OBSFREQ']
obsnchan = hdu0.header['OBSNCHAN']
obsbw = hdu0.header['OBSBW']
fmin = obsfreq - obsbw/2.
fmax = obsfreq + obsbw/2.
nf = obsnchan
df = hdu1.header['CHAN_BW']
tsamp = hdu1.header['TBIN']
nsubint = hdu1.header['NAXIS2']
samppersubint = int(hdu1.header['NSBLK'])
nsamp = nsubint * samppersubint
sourename = hdu0.header['SRC_NAME']
ra = hdu0.header['RA']
dec = hdu0.header['DEC']
#tstart = Decimal("%d.%d" % (hdu0.header['STT_IMJD'], hdu0.header['STT_SMJD']))
subintoffset = hdu1.header['NSUBOFFS']
tstart = "%.13f" % (Decimal(hdu0.header['STT_IMJD']) + Decimal(hdu0.header['STT_SMJD'] + tsamp * samppersubint * subintoffset )/secperday )
nbits = hdu0.header['BITPIX']
header = hdu0.header + hdu1.header
dtype = ''

print( 'freq %s MHz, nchan %d, bw %s MHz' % ( obsfreq, obsnchan, obsbw))
print( 'MJD:', tstart)
print( 'fch1, df', fch1, df )

print( 'data.shape:', data1.shape)

a,b,c,d,e = data1.shape
d = (endfreq - startfreq)
if c > 1:
    data = data1[startn:endn,:,0,startfreq:endfreq,:].squeeze().reshape((-1,d))
else:
    data = data1.squeeze().reshape((-1,d))
l, m = data.shape
print( data.shape)
#data = data.reshape(l/64, 64, d).sum(axis=1)
#data -= data.mean(axis=0).transpose().astype(np.uint64)

from pylab import *
extents, correxts = driftrate.getExtents(data.T, df=df, dt=tsamp*1000, lowest_freq=1000)
imshow(data.T, aspect='auto',cmap=get_cmap("hot"),origin="lower", extent=extents)
xlabel("Time (ms)")
ylabel("Frequency (MHz)")
colorbar()
#plot(data.sum(axis=0))
#plot(data.sum(axis=1))
show()
