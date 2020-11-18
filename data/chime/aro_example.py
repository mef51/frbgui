import matplotlib.pyplot as plt
import numpy as np

def get_stokes(data):
	X = data[:,:,0]
	Y = data[:,:,1]
	I = abs(X) ** 2 + abs(Y) ** 2
	Q = abs(X) ** 2 - abs(Y) ** 2
	U = 2 * np.real(X * np.conj(Y))
	V = -2 * np.imag(X * np.conj(Y))
	return I, Q, U, V

data = np.load("aro_SGR1935+2154_20200428_baseband.npz")
print(data.files)
cv = data["V"]

# complex voltages are shaped (nt, nf, npol)
nt, nf, _ = cv.shape
tstart = np.datetime64(str(data["start_time"]))
tstop = np.datetime64(str(data["stop_time"]))
dt = (tstop - tstart) / nt
dm = data["DM"]
freq = np.linspace(800, 400, 1024, endpoint=False)[::-1]
I, Q, U, V = get_stokes(cv)

# change shape to (nf, nt) with the bottom frequency at index 0
intensity = np.flipud(I.T)

# self-calibrate data
for ii in range(nf):
	chan = intensity[ii,:]
	if np.nansum(chan) == 0.:
		continue
	mean = np.nanmean(chan)
	chan[:] = chan[:] / mean
	chan[:] = chan[:] - 1
	var = np.nanvar(chan)
	chan[:] = chan[:] / var

# downsampling factors
ds = 384
sub_factor = 4

# downsample if necessary
if ds > 1:
	new_num_spectra = int(nt / ds)
	num_to_trim = nt % ds
	if num_to_trim > 0:
		intensity = intensity[:,:-num_to_trim]
	intensity = np.array(np.column_stack(
		[np.mean(intensities, axis=1) for intensities \
		 in np.hsplit(intensity, new_num_spectra)]))
	nf, nt = intensity.shape

# subband if necessary
if sub_factor > 1:
	intensity = np.nanmean(intensity.reshape(-1, sub_factor, intensity.shape[1]), axis=1)
	freq = np.nanmean(freq.reshape(-1, sub_factor), axis=1)
time_s = np.arange(tstart, tstop - dt * ds, dt * ds)
variance = np.nanvar(intensity, axis=1)

# zap outlier channels
intensity[variance > 0.004, ...] = 0.

# plot waterfall
print('timeres', dt*ds)
print(freq)
np.savez('aro_detection', intensity=intensity)
plt.imshow(intensity, origin="lower", interpolation="nearest", aspect="auto")
plt.savefig("aro_wfall.png")
