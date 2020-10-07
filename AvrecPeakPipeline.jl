using CSV, DataFrames
using StatsPlots
using HypothesisTests

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
        PeakDatSTRatio = CSV.File("AVRECPeakCLST.csv") |> DataFrame # single trial, no skipped peaks
        println("Click Trains")
    elseif stimtype[iTyp] == "AM"
        PeakDatTA = CSV.File("AVRECPeakAM.csv") |> DataFrame # trial average
        PeakDatST = CSV.File("AVRECPeakAMST.csv") |> DataFrame # single trial
        PeakDatSTRatio = CSV.File("AVRECPeakAMST.csv") |> DataFrame # single trial, no skipped peaks
        println("Amplitude Modulation")
    end

    filter(row -> ! isnan(row.PeakAmp), PeakDatTA)
    filter(row -> ! isnan(row.PeakAmp), PeakDatST)

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
    # Average Trials ---------------------------------------------------------------------------

    # seperate just stimulus presentation from full table
    Stim2Hz = PeakDatTA[PeakDatTA[!,:ClickFreq] .== 2,:]
    Stim5Hz = PeakDatTA[PeakDatTA[!,:ClickFreq] .== 5,:]

    ## Peak Amp/Lat/RMS response difference
    peaks = ["First" "Second" "Third" "Fourth" "Fifth"]
    for ipeak = 1:length(peaks)
        whichpeak = peaks[ipeak]

        if ipeak <= 2
            Stat2  = Stim2Hz[Stim2Hz[!,:OrderofClick] .== ipeak,:]
            Avrec1Peak(figs,Stat2,whichpeak,"2Hz",savetype,"TA")
            Peak1_Between(data,Stat2,whichpeak,"2Hz","TA")
            Peak1_Within(data,Stat2,whichpeak,"2Hz","TA")
        end

        Stat5  = Stim5Hz[Stim5Hz[!,:OrderofClick] .== ipeak,:]
        Avrec1Peak(figs,Stat5,whichpeak,"5Hz",savetype,"TA")
        Peak1_Between(data,Stat5,whichpeak,"5Hz","TA")
        Peak1_Within(data,Stat5,whichpeak,"5Hz","TA")
    end

    ### Peak Amp/RMS Ratio of Last/First response
    # -> J. Heck, in her paper, used EPSP5/EPSP1 to show the difference between groups of the ratio from the last to first stimulus response

    # seperate the 1st and last click
    Stat2  = Stim2Hz[Stim2Hz[!,:OrderofClick] .== 1,:]
    Last2  = Stim2Hz[Stim2Hz[!,:OrderofClick] .== 2,:]
    Stat5  = Stim5Hz[Stim5Hz[!,:OrderofClick] .== 1,:]
    Last5  = Stim5Hz[Stim5Hz[!,:OrderofClick] .== 5,:]
    # divide the last by first
    Ratio2AMP = Last2[!,:PeakAmp] ./ Stat2[!,:PeakAmp]
    Ratio5AMP = Last5[!,:PeakAmp] ./ Stat5[!,:PeakAmp]
    Ratio2RMS = Last2[!,:RMS] ./ Stat2[!,:RMS]
    Ratio5RMS = Last5[!,:RMS] ./ Stat5[!,:RMS]
    # add this column to the table to keep tags
    Stat2.RatioAMP, Stat2.RatioRMS = Ratio2AMP, Ratio2RMS
    Stat5.RatioAMP, Stat5.RatioRMS = Ratio5AMP, Ratio5RMS

    # cheeky plots first;
    AvrecPeakRatio(figs,Stat2,Stat5)
    # And stats for overlay
    PeakRatio_Between(data,Stat2,"2Hz")
    PeakRatio_Within(data,Stat2,"2Hz")
    PeakRatio_Between(data,Stat5,"5Hz")
    PeakRatio_Within(data,Stat5,"5Hz")

    # Scatter Plots for visualization and possibly use for Brown Forsythe stats overlay
    AvrecScatter(figs,Stim2Hz,"2Hz",savetype,"TA")
    AvrecScatter(figs,Stim5Hz,"5Hz",savetype,"TA")


    # Single Trials ----------------------------------------------------------------------------

    # seperate just stimulus presentation from full table
    Stim2Hz = PeakDatST[PeakDatST[!,:ClickFreq] .== 2,:]
    Stim5Hz = PeakDatST[PeakDatST[!,:ClickFreq] .== 5,:]

    ## Peak Amp response difference
    peaks = ["First" "Second" "Third" "Fourth" "Fifth"]
    for ipeak = 1:length(peaks)
        whichpeak = peaks[ipeak]

        if ipeak <= 2
            Stat2  = Stim2Hz[Stim2Hz[!,:OrderofClick] .== ipeak,:]
            Avrec1Peak(figs,Stat2,whichpeak,"2Hz",savetype,"ST")
            Peak1_Between(data,Stat2,whichpeak,"2Hz","ST")
            Peak1_Within(data,Stat2,whichpeak,"2Hz","ST")
        end

        Stat5  = Stim5Hz[Stim5Hz[!,:OrderofClick] .== ipeak,:]
        Avrec1Peak(figs,Stat5,whichpeak,"5Hz",savetype,"ST")
        Peak1_Between(data,Stat5,whichpeak,"5Hz","ST")
        Peak1_Within(data,Stat5,whichpeak,"5Hz","ST")
    end

    # Scatter Plots for visualization and possibly use for Brown Forsythe stats overlay
    AvrecScatter(figs,Stim2Hz,"2Hz",savetype,"ST")
    AvrecScatter(figs,Stim5Hz,"5Hz",savetype,"ST")

    ### Peak Ratio of Last/First response
    Stim2Hz = PeakDatSTRatio[PeakDatSTRatio[!,:ClickFreq] .== 2,:]
    Stim5Hz = PeakDatSTRatio[PeakDatSTRatio[!,:ClickFreq] .== 5,:]

    # seperate the 1st and last click
    Stat2  = Stim2Hz[Stim2Hz[!,:OrderofClick] .== 1,:]
    Last2  = Stim2Hz[Stim2Hz[!,:OrderofClick] .== 2,:]
    Stat5  = Stim5Hz[Stim5Hz[!,:OrderofClick] .== 1,:]
    Last5  = Stim5Hz[Stim5Hz[!,:OrderofClick] .== 5,:]
    # divide the last by first
    Ratio2AMP = Last2[!,:PeakAmp] ./ Stat2[!,:PeakAmp]
    Ratio5AMP = Last5[!,:PeakAmp] ./ Stat5[!,:PeakAmp]
    Ratio2RMS = Last2[!,:RMS] ./ Stat2[!,:RMS]
    Ratio5RMS = Last5[!,:RMS] ./ Stat5[!,:RMS]
    # add this column to the table to keep tags
    Stat2.RatioAMP, Stat2.RatioRMS = Ratio2AMP, Ratio2RMS
    Stat5.RatioAMP, Stat5.RatioRMS = Ratio5AMP, Ratio5RMS

    AvrecPeakRatio(figs,Stat2,Stat5, savetype,"ST")
    PeakRatio_Between(data,Stat2,"2Hz","ST")
    PeakRatio_Within(data,Stat2,"2Hz","ST")
    PeakRatio_Between(data,Stat5,"5Hz","ST")
    PeakRatio_Within(data,Stat5,"5Hz","ST")
end