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

savetype = ".pdf" # choose how all figures are saved, default ".pdf"
stimtype = ["CL" "AM"]

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
    # further seperate by stimulus
    KIC2, KIC5, KIC10 = KIC[KIC[!,:ClickFreq] .== 2,:], KIC[KIC[!,:ClickFreq] .== 5,:], KIC[KIC[!,:ClickFreq] .== 10,:]
    KIT2, KIT5, KIT10 = KIT[KIT[!,:ClickFreq] .== 2,:], KIT[KIT[!,:ClickFreq] .== 5,:], KIT[KIT[!,:ClickFreq] .== 10,:]
    KIV2, KIV5, KIV10 = KIV[KIV[!,:ClickFreq] .== 2,:], KIV[KIV[!,:ClickFreq] .== 5,:], KIV[KIV[!,:ClickFreq] .== 10,:]

    ## Box Plots First ### 
    # Output: figures in folder AvrecPeakPlots_againstMeasurement/_againstLayer of Peak Amplitude over measurement/layer per layer/measurement
    AvrecPeakvsMeas(figs,KIC2,"KIC","2Hz",savetype,stimtype[iTyp])
    AvrecPeakvsMeas(figs,KIC5,"KIC","5Hz",savetype,stimtype[iTyp])
    AvrecPeakvsMeas(figs,KIC10,"KIC","10Hz",savetype,stimtype[iTyp])

    AvrecPeakvsMeas(figs,KIT2,"KIT","2Hz",savetype,stimtype[iTyp])
    AvrecPeakvsMeas(figs,KIT5,"KIT","5Hz",savetype,stimtype[iTyp])
    AvrecPeakvsMeas(figs,KIT10,"KIT","10Hz",savetype,stimtype[iTyp])

    AvrecPeakvsMeas(figs,KIV2,"KIV","2Hz",savetype,stimtype[iTyp])
    AvrecPeakvsMeas(figs,KIV5,"KIV","5Hz",savetype,stimtype[iTyp])
    AvrecPeakvsMeas(figs,KIV10,"KIV","10Hz",savetype,stimtype[iTyp])

    AvrecPeakvsLay(figs,KIC2,"KIC","2Hz",savetype)
    AvrecPeakvsLay(figs,KIC5,"KIC","5Hz",savetype)
    AvrecPeakvsLay(figs,KIC10,"KIC","10Hz",savetype)

    AvrecPeakvsLay(figs,KIT2,"KIT","2Hz",savetype)
    AvrecPeakvsLay(figs,KIT5,"KIT","5Hz",savetype)
    AvrecPeakvsLay(figs,KIT10,"KIT","10Hz",savetype)

    AvrecPeakvsLay(figs,KIV2,"KIV","2Hz",savetype)
    AvrecPeakvsLay(figs,KIV5,"KIV","5Hz",savetype)
    AvrecPeakvsLay(figs,KIV10,"KIV","10Hz",savetype)

    ### Now Stats ###
    # Average Trials -----------------------------------------------------------------------

    # seperate just stimulus presentation from full table
    Stim2Hz = PeakDatTA[PeakDatTA[!,:ClickFreq] .== 2,:]
    Stim5Hz = PeakDatTA[PeakDatTA[!,:ClickFreq] .== 5,:]
    Stim10Hz = PeakDatTA[PeakDatTA[!,:ClickFreq] .== 10,:]

    ## Peak Amp/Lat/RMS response difference
    peaks = ["1st" "2nd" "3rd" "4th" "5th" "6th" "7th" "8th" "9th" "10th"]
    for ipeak = 1:length(peaks)
        whichpeak = peaks[ipeak]

        if ipeak <= 2
            Stat2 = Stim2Hz[Stim2Hz[!,:OrderofClick] .== ipeak,:]
            Avrec1Peak(figs,Stat2,whichpeak,"2Hz",savetype,stimtype[iTyp],"TA")
            Peak1_Between(data,Stat2,whichpeak,"2Hz",stimtype[iTyp],"TA")
            Peak1_Within(data,Stat2,whichpeak,"2Hz",stimtype[iTyp],"TA")
        end
        if ipeak <= 5
            Stat5 = Stim5Hz[Stim5Hz[!,:OrderofClick] .== ipeak,:]
            Avrec1Peak(figs,Stat5,whichpeak,"5Hz",savetype,stimtype[iTyp],"TA")
            Peak1_Between(data,Stat5,whichpeak,"5Hz",stimtype[iTyp],"TA")
            Peak1_Within(data,Stat5,whichpeak,"5Hz",stimtype[iTyp],"TA")
        end
        Stat10 = Stim10Hz[Stim10Hz[!,:OrderofClick] .== ipeak,:]
        Avrec1Peak(figs,Stat10,whichpeak,"10Hz",savetype,stimtype[iTyp],"TA")
        Peak1_Between(data,Stat10,whichpeak,"10Hz",stimtype[iTyp],"TA")
        Peak1_Within(data,Stat10,whichpeak,"10Hz",stimtype[iTyp],"TA")
    end

    ### Peak Amp/RMS Ratio of Last/First response
    # -> J. Heck, in her paper, used EPSP5/EPSP1 to show the difference between groups of the ratio from the last to first stimulus response

    # seperate the 1st and last click
    Stat2  = Stim2Hz[Stim2Hz[!,:OrderofClick] .== 1,:]
    Last2  = Stim2Hz[Stim2Hz[!,:OrderofClick] .== 2,:]
    Stat5  = Stim5Hz[Stim5Hz[!,:OrderofClick] .== 1,:]
    Last5  = Stim5Hz[Stim5Hz[!,:OrderofClick] .== 5,:]
    Stat10 = Stim10Hz[Stim10Hz[!,:OrderofClick] .== 1,:]
    Last10 = Stim10Hz[Stim10Hz[!,:OrderofClick] .== 10,:]
    # divide the last by first
    Ratio2AMP  = Last2[!,:PeakAmp] ./ Stat2[!,:PeakAmp]
    Ratio5AMP  = Last5[!,:PeakAmp] ./ Stat5[!,:PeakAmp]
    Ratio10AMP = Last10[!,:PeakAmp] ./ Stat10[!,:PeakAmp]
    Ratio2RMS  = Last2[!,:RMS] ./ Stat2[!,:RMS]
    Ratio5RMS  = Last5[!,:RMS] ./ Stat5[!,:RMS]
    Ratio10RMS = Last10[!,:RMS] ./ Stat10[!,:RMS]
    # add this column to the table to keep tags
    Stat2.RatioAMP, Stat2.RatioRMS = Ratio2AMP, Ratio2RMS
    Stat5.RatioAMP, Stat5.RatioRMS = Ratio5AMP, Ratio5RMS
    Stat10.RatioAMP, Stat10.RatioRMS = Ratio10AMP, Ratio10RMS

    # cheeky plots first;
    AvrecPeakRatio(figs,Stat2,"2Hz",savetype,stimtype[iTyp],"TA")
    AvrecPeakRatio(figs,Stat5,"5Hz",savetype,stimtype[iTyp],"TA")
    AvrecPeakRatio(figs,Stat10,"10Hz",savetype,stimtype[iTyp],"TA")
    # And stats for overlay
    PeakRatio_Between(data,Stat2,"2Hz",stimtype[iTyp],"TA")
    PeakRatio_Within(data,Stat2,"2Hz",stimtype[iTyp],"TA")
    PeakRatio_Between(data,Stat5,"5Hz",stimtype[iTyp],"TA")
    PeakRatio_Within(data,Stat5,"5Hz",stimtype[iTyp],"TA")
    PeakRatio_Between(data,Stat10,"10Hz",stimtype[iTyp],"TA")
    PeakRatio_Within(data,Stat10,"10Hz",stimtype[iTyp],"TA")

    # Scatter Plots for visualization and possibly use for Brown Forsythe stats overlay
    AvrecScatter(figs,Stim2Hz,"2Hz",savetype,stimtype[iTyp],"TA")
    AvrecScatter(figs,Stim5Hz,"5Hz",savetype,stimtype[iTyp],"TA")
    AvrecScatter(figs,Stim10Hz,"10Hz",savetype,stimtype[iTyp],"TA")

    # Single Trials ------------------------------------------------------------------------

    # seperate just stimulus presentation from full table
    Stim2Hz  = PeakDatST[PeakDatST[!,:ClickFreq] .== 2,:]
    Stim5Hz  = PeakDatST[PeakDatST[!,:ClickFreq] .== 5,:]
    Stim10Hz = PeakDatST[PeakDatST[!,:ClickFreq] .== 10,:]

    ## Peak Amp response difference
    peaks = ["1st" "2nd" "3rd" "4th" "5th" "6th" "7th" "8th" "9th" "10th"]
    for ipeak = 1:length(peaks)
        whichpeak = peaks[ipeak]

        if ipeak <= 2
            Stat2 = Stim2Hz[Stim2Hz[!,:OrderofClick] .== ipeak,:]
            Avrec1Peak(figs,Stat2,whichpeak,"2Hz",savetype,stimtype[iTyp],"ST")
            Peak1_Between(data,Stat2,whichpeak,"2Hz",stimtype[iTyp],"ST")
            Peak1_Within(data,Stat2,whichpeak,"2Hz",stimtype[iTyp],"ST") # note - # trials not the same so not a one but two sample t test
        end
        if ipeak <= 5
            Stat5 = Stim5Hz[Stim5Hz[!,:OrderofClick] .== ipeak,:]
            Avrec1Peak(figs,Stat5,whichpeak,"5Hz",savetype,stimtype[iTyp],"ST")
            Peak1_Between(data,Stat5,whichpeak,"5Hz",stimtype[iTyp],"ST")
            Peak1_Within(data,Stat5,whichpeak,"5Hz",stimtype[iTyp],"ST") # note - # trials not the same so not a one but two sample t test
        end
        Stat10 = Stim10Hz[Stim10Hz[!,:OrderofClick] .== ipeak,:]
        Avrec1Peak(figs,Stat10,whichpeak,"10Hz",savetype,stimtype[iTyp],"ST")
        Peak1_Between(data,Stat10,whichpeak,"10Hz",stimtype[iTyp],"ST")
        Peak1_Within(data,Stat10,whichpeak,"10Hz",stimtype[iTyp],"ST")
    end
    # Scatter Plots for visualization and possibly use for Brown Forsythe stats overlay
    AvrecScatter(figs,Stim2Hz,"2Hz",savetype,stimtype[iTyp],"ST")
    AvrecScatter(figs,Stim5Hz,"5Hz",savetype,stimtype[iTyp],"ST")
    AvrecScatter(figs,Stim10Hz,"10Hz",savetype,stimtype[iTyp],"ST")

    ### Peak Ratio of Last/First response

    # seperate the 1st and last click
    Stat2  = Stim2Hz[Stim2Hz[!,:OrderofClick] .== 1,:]
    Last2  = Stim2Hz[Stim2Hz[!,:OrderofClick] .== 2,:]
    Stat5  = Stim5Hz[Stim5Hz[!,:OrderofClick] .== 1,:]
    Last5  = Stim5Hz[Stim5Hz[!,:OrderofClick] .== 5,:]
    Stat10 = Stim10Hz[Stim10Hz[!,:OrderofClick] .== 1,:]
    Last10 = Stim10Hz[Stim10Hz[!,:OrderofClick] .== 10,:]
    # divide the last by first
    Ratio2AMP = Last2[!,:PeakAmp] ./ Stat2[!,:PeakAmp] 
    Ratio5AMP = Last5[!,:PeakAmp] ./ Stat5[!,:PeakAmp]
    Ratio10AMP = Last10[!,:PeakAmp] ./ Stat10[!,:PeakAmp]
    Ratio2RMS = Last2[!,:RMS] ./ Stat2[!,:RMS]
    Ratio5RMS = Last5[!,:RMS] ./ Stat5[!,:RMS]
    Ratio10RMS = Last10[!,:RMS] ./ Stat10[!,:RMS]
    # add this column to the table to keep tags
    Stat2.RatioAMP, Stat2.RatioRMS = Ratio2AMP, Ratio2RMS
    Stat5.RatioAMP, Stat5.RatioRMS = Ratio5AMP, Ratio5RMS
    Stat10.RatioAMP, Stat10.RatioRMS = Ratio10AMP, Ratio10RMS

    AvrecPeakRatio(figs,Stat2,"2Hz",savetype,stimtype[iTyp],"ST")
    AvrecPeakRatio(figs,Stat5,"5Hz",savetype,stimtype[iTyp],"ST")
    AvrecPeakRatio(figs,Stat10,"10Hz",savetype,stimtype[iTyp],"ST")

    PeakRatio_Between(data,Stat2,"2Hz",stimtype[iTyp],"ST")
    PeakRatio_Within(data,Stat2,"2Hz",stimtype[iTyp],"ST")
    PeakRatio_Between(data,Stat5,"5Hz",stimtype[iTyp],"ST")
    PeakRatio_Within(data,Stat5,"5Hz",stimtype[iTyp],"ST")
    PeakRatio_Between(data,Stat10,"10Hz",stimtype[iTyp],"ST")
    PeakRatio_Within(data,Stat10,"10Hz",stimtype[iTyp],"ST")
end