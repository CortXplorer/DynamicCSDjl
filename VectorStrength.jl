# from https://journals.physiology.org/doi/pdf/10.1152/jn.01109.2007 :
# The vector strength and mean phase were obtained by expressing spike times relative to the phase of the modulator, representing each spike as a unit vector with orientation given by that phase, and computing the sum of the unit vectors. The vector strength was given by the length of the resultant vector divided by the number of vectors. It could range from 0 (no phase locking) to 1 (all spikes at identical phase); spike probability exactly following the sine modulator waveform would yield a vector strength of 0.5. The mean phase was given by the orientation of the resultant vector. The mean phase lag tended to increase linearly with increasing modulator frequency, as shown in RESULTS.


# prep the data which will fluxuate
peaklat   = Int.(ceil.(rand(10)*100)) # list of fake latencies
clickfreq = 2
sr        = 1000

# generate the AM wave used and get instantaneous phase
wave = getAMwave(clickfreq,sr)
wavephase = angle.(signal.hilbert(wave)) # Compute phase of wave

# get the phase at each time point of signal peak amplitude
orientation = wavephase[peaklat]
magnitudes = ones(length(orientation))

# break vector into components i and j to compute sums, averages, and then the resultant vector 
i = cos.(orientation) # get component i 
j = sin.(orientation)
sumi = sum(i) # 
sumj = sum(j)
avgi = sumi/10 # average of component i
avgj = sumj/10 # average of component j

vstrength = sqrt((avgi^2)+(avgj^2)) # find the vector based on i and j components (pythagoreon)
meanphase = atan(avgj/avgi) # find the mean phase in radian, 
# test1 = acos(sumi/vstrength) -- There's an issue here I'm not quite sure how to address. What of the signs? If I do atan with two negative components, then asin and acos give me a different radian value because vstrength will *always* be positive. Should I take the absolute values to avoid mismatches or just always take atan so that the signs of the averages speak for themselves? Probably the latter.

# let's plot that as a simple histogram 
histogram(r,xlims=(-π,π),xticks=-π:4/π:π,label="moderator phase at peak latency")

function getAMwave(clickfreq=2,sr=1000)

    t  = 0:1/sr:1 # time course
    wave = sin.((2*π*clickfreq*t).+3π/2)
    return wave

end