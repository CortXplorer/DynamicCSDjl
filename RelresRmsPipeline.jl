### This pipeline takes extracted features from relres trace data in the form of csv (output from matlab) and produces graphs and plots for all t tests and cohen's d tests

## input:   home/Data/RELRESrms**.csv
## output:  home/figs/Relres1Rms, .../CohensDPlots -&- home/Data/RelresStats

# # do for all packages: # #
# type ']' to enter pkg> 
# then in pkg> type 'add Plots' 
# then you can backspace to return to julia> and use 'using Plots'
using Plots, OhMyREPL
using CSV, DataFrames
using StatsPlots
using HypothesisTests, EffectSizes
# using Infiltrator # only used when adding breakpoints

# # make sure Julia is in the directory you would like to use: # #
# type '@__DIR__' to print out the directory julia is currently pointing to
home    = @__DIR__
func    = joinpath(home,"functions")
figs    = joinpath(home,"figs")
data    = joinpath(home,"Data")
include(joinpath(func,"RelresRmsStats.jl")) # contains all stats functions and plot function

savetype = ".pdf" # choose how all figures are saved, default ".pdf"
stimtype = ["CL" "AM"] # CL = click train and AM = amplitude modulation
freqtype = ["5Hz" "10Hz"] # "2Hz" "5Hz" "10Hz" "20Hz" "40Hz"# frequency of train or modulation

# This loop runs through stim type (click and AM), and frequency (2 - 40 Hz). It analyzes groups KIT, KIC, and KIV and can handle changes in group size if animals are added. Analysis is run for trial average and single-trial (ST) data
for iTyp = 1:length(stimtype)
    # Load in data from matlab table csv file which contains 2 and 5 hz peak amp and latency. 
    if stimtype[iTyp] == "CL"
        RmsDatST = CSV.File(joinpath(data,"RELRESrmsCLST.csv")) |> DataFrame # single trial
        println("Click Trains")
    elseif stimtype[iTyp] == "AM"
        RmsDatST = CSV.File(joinpath(data,"RELRESrmsAMST.csv")) |> DataFrame # single trial
        println("Amplitude Modulation")
    end

    for iFrq = 1:length(freqtype)
        println("At frequency: " * freqtype[iFrq])
        StimHz = RmsDatST[RmsDatST[:,:ClickFreq] .== parse(Int,freqtype[iFrq][begin:end-2]),:]

        ## RMS response difference
        peaks = ["1st"] # "2nd" "3rd" "4th" "5th" "6th" "7th" "8th" "9th" "10th"
        for ipeak = 1:length(peaks)

            if ipeak <= StimHz.ClickFreq[1] # cut this to amount of detection windows
                # pull out current response window (i.e. 1st peak)
                Stat = StimHz[StimHz[:,:OrderofClick] .== ipeak,:]
                Relres1Rms(figs,Stat,peaks[ipeak],freqtype[iFrq],savetype,stimtype[iTyp], "ST")
                Rms1_Between(data,Stat,peaks[ipeak],freqtype[iFrq],stimtype[iTyp],"ST")
                Rms1_Within(data,Stat,peaks[ipeak],freqtype[iFrq],stimtype[iTyp],"ST")
            end
        end
    end
end

include(joinpath(func,"CohensProg.jl")) # contains cohen's d plotting function (stats output already has cohen's d results per row)

CohensProgRelres(figs, data, freqtype, stimtype, savetype)

# Run RAnova.R before or after this