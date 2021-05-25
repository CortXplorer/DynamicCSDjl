# from https://journals.physiology.org/doi/pdf/10.1152/jn.01109.2007 :
# The vector strength and mean phase were obtained by expressing spike times relative to the phase of the modulator, representing each spike as a unit vector with orientation given by that phase, and computing the sum of the unit vectors. The vector strength was given by the length of the resultant vector divided by the number of vectors. It could range from 0 (no phase locking) to 1 (all spikes at identical phase); spike probability exactly following the sine modulator waveform would yield a vector strength of 0.5. The mean phase was given by the orientation of the resultant vector. The mean phase lag tended to increase linearly with increasing modulator frequency.

using DSP, FFTW
using CSV, DataFrames

home    = @__DIR__
data    = joinpath(home,"Data")
func    = joinpath(home,"functions")
# group   = joinpath(home,"groups")
# figs    = joinpath(home,"figs")
# spect   = joinpath(data,"Spectral")

# include(joinpath(group,"callGroup.jl"))
include(joinpath(func,"VSfunc.jl"))

# unchanging parameters:
sr        = 1000 # sampling rate

# open the dataframe in the workspace
PeakDatST = CSV.File(joinpath(data,"AVRECPeakAMST.csv")) |> DataFrame

# determine which click frequency 
clickfreq = 2    # stimulus frequency

# sort which parameters we're pulling
KIC = PeakDatST[PeakDatST[:,:Group] .== "KIC",:]
KIC = KIC[KIC[:,:Animal] .== "KIC02",:]
KIC = KIC[KIC[:,:Layer] .== "All",:]
KIC = KIC[KIC[:,:Measurement] .== "preAM_1",:]
KIC = KIC[KIC[:,:ClickFreq] .== clickfreq,:]

# our final list of latencies:
KIC = filter(row -> ! isnan(row.PeakLat), KIC)
peaklat = Int.(KIC[:,:PeakLat])

# generate the AM wave used and get instantaneous phase
wave = getAMwave(clickfreq,sr) # currently below
wavephase = angle.(hilbert(wave)) # Compute phase of wave

# get the phase at each time point of signal peak amplitude
orientation = wavephase[peaklat]

### Rayleigh Test of circular uniformity ### 
# _____________________________________________________________________________________
# break vector into components i and j to compute sums, averages, and then the resultant vector 
i = cos.(orientation) # get component i 
j = sin.(orientation)
sumi = sum(i) # 
sumj = sum(j)
avgi = sumi/length(i) # average of component i
avgj = sumj/length(j) # average of component j

vstrength = sqrt((avgi^2)+(avgj^2)) # find the vector based on i and j components (pythagoreon)

# zscore and subesquest p value:
Z = vstrength^2 * length(orientation)
p = exp(-Z) 

# let's plot that as a simple histogram 
# histogram(orientation,xlims=(-π,π),xticks=-π:4/π:π,label="moderator phase at peak latency")
# _____________________________________________________________________________________

meanphase = atan(abs(avgj)/abs(avgi)) # find the mean phase in radian; 
# this gives the corresponding angle upper right quadrant. we need to translate it to the correct quadrant based on the polar graph and signs of i and j
if avgi > 0 && avgj > 0     # upper right target
    meanphase = 0 + meanphase
elseif avgi > 0 && avgj < 0 # lower right target
    meanphase = 2π - meanphase
elseif avgi < 0 && avgj < 0 # lower left target
    meanphase = π + meanphase
elseif avgi < 0 && avgj > 0 # upper left target
    meanphase = π - meanphase
end

# return vstrength, Z, p, meanphase

# Fuck yes. I got this.
