### This pipeline takes extracted features from avrec and layer trace data in the form of csv (output from matlab) and produces graphs and plots for all t tests and cohen's d tests as well as a few extra plots to explore the data

## input:   home/Data/AVRECPeak**.csv
## output:  home/figs/AM_TemporalFlow, .../Avrec1Peak, .../AvrecPeakPlots_againstLayer, .../AvrecPeakPlots_againstMeasurement, .../AvrecPeakRatio, .../AvrecScatter, .../CohensDPlots -&- home/Data/AvrecPeakStats

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
include(joinpath(func,"AvrecPeakPlots.jl")) # contains all plotting functions 
include(joinpath(func,"AvrecPeakStats.jl")) # contains all stats functions

savetype = ".png" # choose how all figures are saved, default ".pdf"
stimtype = ["CL" "AM"] # CL = click train and AM = amplitude modulation
freqtype = ["2Hz" "5Hz" "10Hz" "20Hz" "40Hz"] # frequency of train or modulation

# This loop runs through stim type (click and AM), and frequency (2 - 40 Hz). It analyzes groups KIT, KIC, and KIV and can handle changes in group size if animals are added. Analysis is run for trial average and single-trial (ST) data
for iTyp = 1:length(stimtype)
    # Load in data from matlab table csv file which contains 2 and 5 hz peak amp and latency. 
    if stimtype[iTyp] == "CL"
        PeakDatTA = CSV.File(joinpath(data,"AVRECPeakCL.csv")) |> DataFrame # trial average
        PeakDatST = CSV.File(joinpath(data,"AVRECPeakCLST.csv")) |> DataFrame # single trial
        println("Click Trains")
    elseif stimtype[iTyp] == "AM"
        PeakDatTA = CSV.File(joinpath(data,"AVRECPeakAM.csv")) |> DataFrame # trial average
        PeakDatST = CSV.File(joinpath(data,"AVRECPeakAMST.csv")) |> DataFrame # single trial
        println("Amplitude Modulation")
    end

    # seperate by group
    KIC = PeakDatTA[PeakDatTA[:,:Group] .== "KIC",:]
    KIT = PeakDatTA[PeakDatTA[:,:Group] .== "KIT",:]
    KIV = PeakDatTA[PeakDatTA[:,:Group] .== "KIV",:]

    for iFrq = 1:length(freqtype)
        println("At frequency: " * freqtype[iFrq])
        # further seperate by stimulus
        KICfreq = KIC[KIC[:,:ClickFreq] .== parse(Int,freqtype[iFrq][begin:end-2]),:]
        KITfreq = KIT[KIT[:,:ClickFreq] .== parse(Int,freqtype[iFrq][begin:end-2]),:]
        KIVfreq = KIV[KIV[:,:ClickFreq] .== parse(Int,freqtype[iFrq][begin:end-2]),:]

        ## Box Plots First ### 
        # Output: figures in folder AvrecPeakPlots_againstMeasurement/_againstLayer of Peak Amplitude over measurement/layer per layer/measurement
        println("Peak vs Measurement")
        AvrecPeakvsMeas(figs,KICfreq,"KIC",freqtype[iFrq],savetype,stimtype[iTyp])
        AvrecPeakvsMeas(figs,KITfreq,"KIT",freqtype[iFrq],savetype,stimtype[iTyp])
        AvrecPeakvsMeas(figs,KIVfreq,"KIV",freqtype[iFrq],savetype,stimtype[iTyp])

        println("Peak vs Layer")
        AvrecPeakvsLay(figs,KICfreq,"KIC",freqtype[iFrq],savetype)
        AvrecPeakvsLay(figs,KITfreq,"KIT",freqtype[iFrq],savetype)
        AvrecPeakvsLay(figs,KIVfreq,"KIV",freqtype[iFrq],savetype)

        ### Now Stats ###
        #-----------------------------------------------------------------------
        TAorST, StatsTitle = ["TA" "ST"], ["Average Trial Stats" "Single Trial Stats"]
        for iStat = 1:length(TAorST)
            println(StatsTitle[iStat])
            # seperate just stimulus presentation from full table
            if TAorST[iStat] == "TA"
                StimHz = PeakDatTA[PeakDatTA[:,:ClickFreq] .== parse(Int,freqtype[iFrq][begin:end-2]),:]
            elseif TAorST[iStat] == "ST"
                StimHz = PeakDatST[PeakDatST[:,:ClickFreq] .== parse(Int,freqtype[iFrq][begin:end-2]),:]
            end

            ## Peak Amp/Lat/RMS response difference
            peaks = ["1st" "2nd" "3rd" "4th" "5th" "6th" "7th" "8th" "9th" "10th"] # limited to first 10 now, no reason to go beyond to 20 and 40 yet
            for ipeak = 1:length(peaks)

                if ipeak <= StimHz.ClickFreq[1] # cut this to amount of detection windows
                    # pull out current response window (i.e. 1st peak)
                    Stat = StimHz[StimHz[:,:OrderofClick] .== ipeak,:]
                    Avrec1Peak(figs,Stat,peaks[ipeak],freqtype[iFrq],savetype,stimtype[iTyp],TAorST[iStat])
                    Peak1_Between(data,Stat,peaks[ipeak],freqtype[iFrq],stimtype[iTyp],TAorST[iStat])
                    Peak1_Within(data,Stat,peaks[ipeak],freqtype[iFrq],stimtype[iTyp],TAorST[iStat])
                end

            end

            ### Peak Amp/RMS Ratio of Last/First response
            # -> J. Heck, in her paper, used EPSP5/EPSP1 to show the difference between groups of the ratio from the last to first stimulus response

            # seperate the 1st and last click
            Stat  = StimHz[StimHz[:,:OrderofClick] .== 1,:]
            Last  = StimHz[StimHz[:,:OrderofClick] .== parse(Int,freqtype[iFrq][begin:end-2]),:]

            # divide the last by first
            RatioAMP  = Last[:,:PeakAmp] ./ Stat[!,:PeakAmp] .* 100
            RatioRMS  = Last[:,:RMS] ./ Stat[!,:RMS] .* 100
            # add this column to the table to keep tags
            Stat.RatioAMP, Stat.RatioRMS = RatioAMP, RatioRMS

            # cheeky plots first;
            AvrecPeakRatio(figs,Stat,freqtype[iFrq],savetype,stimtype[iTyp],TAorST[iStat])

            # And stats for overlay
            PeakRatio_Between(data,Stat,freqtype[iFrq],stimtype[iTyp],TAorST[iStat])
            PeakRatio_Within(data,Stat,freqtype[iFrq],stimtype[iTyp],TAorST[iStat])

            # Scatter Plots for visualization and possibly use for Brown Forsythe stats overlay
            AvrecScatter(figs,StimHz,freqtype[iFrq],savetype,stimtype[iTyp],TAorST[iStat])

        end
    end
end

include(joinpath(func,"CohensProg.jl")) # contains cohen's d plotting function (stats output already has cohen's d results per row)

CohensProg(figs, data, freqtype, stimtype, savetype)

# Run RAnova.R  