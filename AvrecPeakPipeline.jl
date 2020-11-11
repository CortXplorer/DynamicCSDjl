using CSV, DataFrames
using StatsPlots
using HypothesisTests
using Infiltrator

home    = @__DIR__
func    = joinpath(home,"functions")
figs    = joinpath(home,"figs")
data    = joinpath(home,"Data")
include(joinpath(func,"AvrecPeakPlots.jl"))
include(joinpath(func,"AvrecPeakStats.jl"))

savetype = ".png" # choose how all figures are saved, default ".pdf"
stimtype = ["CL" "AM"]
freqtype = ["2Hz" "5Hz" "10Hz" "20Hz" "40Hz"]

# This loop runs through stim type click and am, and frequency 2 - 40 Hz. It analyzes groups KIT, KIC, and KIV but can generally handle changes in group size if animals are added. Data in is taken from Matlab in the form of .csv's 
for iTyp = 1:length(stimtype)
    # Load in data from matlab table csv file which contains 2 and 5 hz peak amp and latency. 
    if stimtype[iTyp] == "CL"
        PeakDatTA = CSV.File("AVRECPeakCL.csv") |> DataFrame # trial average
        PeakDatST = CSV.File("AVRECPeakCLST.csv") |> DataFrame # single trial
        println("Click Trains")
    elseif stimtype[iTyp] == "AM"
        PeakDatTA = CSV.File("AVRECPeakAM.csv") |> DataFrame # trial average
        PeakDatST = CSV.File("AVRECPeakAMST.csv") |> DataFrame # single trial
        println("Amplitude Modulation")
    end

    # seperate by group
    KIC = PeakDatTA[PeakDatTA[!,:Group] .== "KIC",:]
    KIT = PeakDatTA[PeakDatTA[!,:Group] .== "KIT",:]
    KIV = PeakDatTA[PeakDatTA[!,:Group] .== "KIV",:]

    for iFrq = 1:length(freqtype)
        println("At frequency: " * freqtype[iFrq])
        # further seperate by stimulus
        KIC = KIC[KIC[!,:ClickFreq] .== parse(Int,freqtype[iFrq][begin:end-2]),:]
        KIT = KIT[KIT[!,:ClickFreq] .== parse(Int,freqtype[iFrq][begin:end-2]),:]
        KIV = KIV[KIV[!,:ClickFreq] .== parse(Int,freqtype[iFrq][begin:end-2]),:]

        ## Box Plots First ### 
        # Output: figures in folder AvrecPeakPlots_againstMeasurement/_againstLayer of Peak Amplitude over measurement/layer per layer/measurement
        println("Peak vs Measurement")
        AvrecPeakvsMeas(figs,KIC,"KIC",freqtype[iFrq],savetype,stimtype[iTyp])
        AvrecPeakvsMeas(figs,KIT,"KIT",freqtype[iFrq],savetype,stimtype[iTyp])
        AvrecPeakvsMeas(figs,KIV,"KIV",freqtype[iFrq],savetype,stimtype[iTyp])

        println("Peak vs Layer")
        AvrecPeakvsLay(figs,KIC,"KIC",freqtype[iFrq],savetype)
        AvrecPeakvsLay(figs,KIT,"KIT",freqtype[iFrq],savetype)
        AvrecPeakvsLay(figs,KIV,"KIV",freqtype[iFrq],savetype)

        ### Now Stats ###
        #-----------------------------------------------------------------------
        TAorAT,StatsTitle = ["TA" "ST"], ["Average Trial Stats" "Single Trial Stats"]
        for iStat = 1:length(TAorST)
            println(StatsTitle[iStat])
            # seperate just stimulus presentation from full table
            if TAorAT[iStat] == "TA"
                StimHz = PeakDatTA[PeakDatTA[!,:ClickFreq] .== parse(Int,freqtype[iFrq][begin:end-2]),:]
            elseif TAorAT[iStat] == "ST"
                StimHz = PeakDatST[PeakDatST[!,:ClickFreq] .== parse(Int,freqtype[iFrq][begin:end-2]),:]
            end

            # remove group KIV for now until it is a proper size
            StimHz = StimHz[StimHz[!,:Group] .!= "KIV",:]

            ## Peak Amp/Lat/RMS response difference
            peaks = ["1st" "2nd" "3rd" "4th" "5th" "6th" "7th" "8th" "9th" "10th"]
            for ipeak = 1:length(peaks)

                if ipeak <= StimHz.ClickFreq[1] # cut this to amount of detection windows
                    Stat = StimHz[StimHz[!,:OrderofClick] .== ipeak,:]
                    Avrec1Peak(figs,Stat,peaks[ipeak],freqtype[iFrq],savetype,stimtype[iTyp],TAorAT[iStat])
                    Peak1_Between(data,Stat,peaks[ipeak],freqtype[iFrq],stimtype[iTyp],TAorAT[iStat])
                    Peak1_Within(data,Stat,peaks[ipeak],freqtype[iFrq],stimtype[iTyp],TAorAT[iStat])
                end

            end

            ### Peak Amp/RMS Ratio of Last/First response
            # -> J. Heck, in her paper, used EPSP5/EPSP1 to show the difference between groups of the ratio from the last to first stimulus response

            # seperate the 1st and last click
            Stat  = StimHz[Stim2Hz[!,:OrderofClick] .== 1,:]
            Last  = StimHz[Stim2Hz[!,:OrderofClick] .== parse(Int,freqtype[iFrq][begin:end-2]),:]

            # divide the last by first
            RatioAMP  = Last[!,:PeakAmp] ./ Stat[!,:PeakAmp] .* 100
            RatioRMS  = Last[!,:RMS] ./ Stat[!,:RMS] .* 100
            # add this column to the table to keep tags
            Stat.RatioAMP, Stat.RatioRMS = RatioAMP, RatioRMS

            # cheeky plots first;
            AvrecPeakRatio(figs,Stat,freqtype[iFrq],savetype,stimtype[iTyp],TAorAT[iStat])

            # And stats for overlay
            PeakRatio_Between(data,Stat,freqtype[iFrq],stimtype[iTyp],TAorAT[iStat])
            PeakRatio_Within(data,Stat,freqtype[iFrq],stimtype[iTyp],TAorAT[iStat])

            # Scatter Plots for visualization and possibly use for Brown Forsythe stats overlay
            AvrecScatter(figs,StimHz,freqtype[iFrq],savetype,stimtype[iTyp],TAorAT[iStat])

        end
    end
end