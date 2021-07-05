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

# curtype = "CL" # "AM" "CL"
# curCF = 10 # 2 5 10 20 40 
layer = ["All" "I_II" "IV" "V" "VI"] # curlay = "IV" 

for curtype = type
    for curCF = clickfreq
        for curlay = layer
            FileName = curtype * "_" * string(curCF) * "_VectorStrength_Stats.csv" 
            statout  = CSV.File(joinpath(data,"VSoutput",FileName)) |> DataFrame

            # basic vector strength, no consideration for significance
            CurRun = statout[statout[:,:Layer] .== curlay,:]
            CurRun = filter(row -> ! isnan(row.VectorStrenght), CurRun)

            vstrengthplot = @df CurRun groupedboxplot(:Measurement,:VectorStrenght, group = :Group, markerstrokewidth=0, ylims = (0,1), legend = false) 
            savename = FileName[1:end-4] * "_" * curlay * ".png"
            savefig(vstrengthplot, joinpath(data,"VSoutput",savename))

            #rule out anything non-significant
            CurRun = CurRun[CurRun[:,:P] .< 0.001,:]

            vstrengthplot = @df CurRun groupedboxplot(:Measurement,:VectorStrenght, group = :Group, markerstrokewidth=0, ylims = (0,1), legend = false) 
            savename = "Corrected" * FileName[1:end-4] * "_" * curlay * ".png"
            savefig(vstrengthplot, joinpath(data,"VSoutput",savename))

            mphaseplot = @df CurRun groupedviolin(:Measurement,:MeanPhase, group = :Group, ylims = (0,6.5), legend = false) #BIMODAL!!!
            savename = "Phase" * FileName[1:end-4] * "_" * curlay * ".png"
            savefig(mphaseplot, joinpath(data,"VSoutput",savename))
        end
    end
end

