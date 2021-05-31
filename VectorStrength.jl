# from https://journals.physiology.org/doi/pdf/10.1152/jn.01109.2007 :
# The vector strength and mean phase were obtained by expressing spike times relative to the phase of the modulator, representing each spike as a unit vector with orientation given by that phase, and computing the sum of the unit vectors. The vector strength was given by the length of the resultant vector divided by the number of vectors. It could range from 0 (no phase locking) to 1 (all spikes at identical phase); spike probability exactly following the sine modulator waveform would yield a vector strength of 0.5. The mean phase was given by the orientation of the resultant vector. The mean phase lag tended to increase linearly with increasing modulator frequency.

using DSP, FFTW
using CSV, DataFrames

home    = @__DIR__
data    = joinpath(home,"Data")
func    = joinpath(home,"functions")
# figs    = joinpath(home,"figs")

include(joinpath(func,"VSfunc.jl"))

# parameters:
sr   = 1000 # sampling rate
type = ["AM" "CL"]
clickfreq = [2 5 10 20 40]   # stimulus frequency

for curtype = type
    for curCF = clickfreq
        VSstats(data, sr, curtype, curCF)
    end
end

