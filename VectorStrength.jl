# from https://journals.physiology.org/doi/pdf/10.1152/jn.01109.2007 :
# The vector strength and mean phase were obtained by expressing spike times relative to the phase of the modulator, representing each spike as a unit vector with orientation given by that phase, and computing the sum of the unit vectors. The vector strength was given by the length of the resultant vector divided by the number of vectors. It could range from 0 (no phase locking) to 1 (all spikes at identical phase); spike probability exactly following the sine modulator waveform would yield a vector strength of 0.5. The mean phase was given by the orientation of the resultant vector. The mean phase lag tended to increase linearly with increasing modulator frequency.

using DSP, FFTW
using CSV, DataFrames
using StatsPlots

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


type = "CL" # "AM" "CL"
clickfreq = 40 # 2 5 10 20 40 
layer = "IV" # "All" "I_II" "IV" "VI"

FileName = type * "_" * string(clickfreq) * "_VectorStrength_Stats.csv" 
statout  = CSV.File(joinpath(data,"VSoutput",FileName)) |> DataFrame

CurRun = statout[statout[:,:Layer] .== layer,:]

# vstrengthplot = @df CurRun groupedboxplot(:Measurement,:VectorStrenght, group = :Group, markerstrokewidth=0) 
# savename = FileName[1:end-4] * ".png"
# savefig(vstrengthplot, joinpath(data,"VSoutput",savename))

CurRun = CurRun[CurRun[:,:P] .< 0.001,:]

# mphaseplot = @df CurRun groupedviolin(:Measurement,:MeanPhase, group = :Group) #BIMODAL!!!
# savename = "Phase" * FileName[1:end-4] * ".png"
# savefig(mphaseplot, joinpath(data,"VSoutput",savename))

vstrengthplot = @df CurRun groupedboxplot(:Measurement,:VectorStrenght, group = :Group, markerstrokewidth=0) 
savename = "Corrected" * FileName[1:end-4] * ".png"
savefig(vstrengthplot, joinpath(data,"VSoutput",savename))

# scatter(statout[:,:P])
# histogram(statout[:,:VectorStrenght])
# histogram(statout[:,:MeanPhase])
# histogram(statout[:,:Z])
