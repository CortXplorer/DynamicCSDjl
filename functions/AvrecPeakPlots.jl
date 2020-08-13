function AvrecPeakvsLay(figs,Tab2,Tab5,GroupName)
    # Input: folder path figs, table for 2 hz and 5 hz, current group being plotted
    # Output: figures in folder AvrecPeakPlots_againstLayer of Peak Amplitude over layer per measurement condition

    foldername = "AvrecPeakPlots_againstLayer"
    if !isdir(joinpath(figs,foldername))
        mkdir(joinpath(figs,foldername))
    end

    MeasList = unique(Tab2[!,:Measurement])

    for iMeas = 1:length(MeasList)
        ### peak amp by layer per measurement ###
        Tab2_Sort = Tab2[Tab2[!,:Measurement] .== MeasList[iMeas],:]
        Tab5_Sort = Tab5[Tab5[!,:Measurement] .== MeasList[iMeas],:]

        Title = GroupName * " PeakAmp of " * MeasList[iMeas] * " at 2 Hz"
        avrecplot = @df Tab2_Sort groupedboxplot(:Layer, :PeakAmp, group = :OrderofClick, bar_position = :dodge, lab= ["Peak 1" "Peak 2"], title=Title, xlab = "Layer", ylab = "Peak Amplitude");

        name = joinpath(figs,foldername,Title) * ".png"
        savefig(avrecplot, name);

        Title = GroupName * " PeakAmp of " * MeasList[iMeas] * " at 5 Hz"
        avrecplot = @df Tab5_Sort groupedboxplot(:Layer, :PeakAmp, group = :OrderofClick, bar_position = :dodge, lab= ["Peak 1" "Peak 2" "Peak 3" "Peak 4" "Peak 5"], title=Title, xlab = "Layer", ylab = "Peak Amplitude");

        name = joinpath(figs,foldername,Title) * ".png"
        savefig(avrecplot, name);
    end
end

function AvrecPeakvsMeas(figs,Tab2,Tab5,GroupName)
    # Input: folder path figs, table for 2 hz and 5 hz, current group being plotted
    # Output: figures in folder AvrecPeakPlots_againstMeasurement of Peak Amplitude over measurement per layer

    foldername = "AvrecPeakPlots_againstMeasurement"
    if !isdir(joinpath(figs,foldername))
        mkdir(joinpath(figs,foldername))
    end

    LayList = unique(Tab2[!,:Layer])

    for iLay = 1:length(LayList)
        ### peak amp by measurement per Layer ###
        Tab2_Sort = Tab2[Tab2[!,:Layer] .== LayList[iLay],:]
        Tab5_Sort = Tab5[Tab5[!,:Layer] .== LayList[iLay],:]

        Title = GroupName * " PeakAmp of " * LayList[iLay] * " at 2 Hz"
        avrecplot = @df Tab2_Sort groupedboxplot(:Measurement, :PeakAmp, group = :OrderofClick, bar_position = :dodge, lab= ["Peak 1" "Peak 2"], title=Title, xlab = "Measurement", ylab = "Peak Amplitude");

        name = joinpath(figs,foldername,Title) * ".png"
        savefig(avrecplot, name);

        Title = GroupName * " PeakAmp of " * LayList[iLay] * " at 5 Hz"
        avrecplot = @df Tab5_Sort groupedboxplot(:Measurement, :PeakAmp, group = :OrderofClick, bar_position = :dodge, lab= ["Peak 1" "Peak 2" "Peak 3" "Peak 4" "Peak 5"], title=Title, xlab = "Measurement", ylab = "Peak Amplitude");

        name = joinpath(figs,foldername,Title) * ".png"
        savefig(avrecplot, name);

        Title = GroupName * " Peak Latency of " * LayList[iLay] * " at 2 Hz"
        avrecplot = @df Tab2_Sort groupedboxplot(:Measurement, :PeakLat, group = :OrderofClick, bar_position = :dodge, lab= ["Peak 1" "Peak 2"], title=Title, xlab = "Measurement", ylab = "Peak Amplitude");

        name = joinpath(figs,foldername,Title) * ".png"
        savefig(avrecplot, name);

        Title = GroupName * " Peak Latency of " * LayList[iLay] * " at 5 Hz"
        avrecplot = @df Tab5_Sort groupedboxplot(:Measurement, :PeakLat, group = :OrderofClick, bar_position = :dodge, lab= ["Peak 1" "Peak 2" "Peak 3" "Peak 4" "Peak 5"], title=Title, xlab = "Measurement", ylab = "Peak Amplitude");

        name = joinpath(figs,foldername,Title) * ".png"
        savefig(avrecplot, name);

        Title = GroupName * " RMS of " * LayList[iLay] * " at 2 Hz"
        avrecplot = @df Tab2_Sort groupedboxplot(:Measurement, :RMS, group = :OrderofClick, bar_position = :dodge, lab= ["Peak 1" "Peak 2"], title=Title, xlab = "Measurement", ylab = "Peak Amplitude");

        name = joinpath(figs,foldername,Title) * ".png"
        savefig(avrecplot, name);

        Title = GroupName * " RMS of " * LayList[iLay] * " at 5 Hz"
        avrecplot = @df Tab5_Sort groupedboxplot(:Measurement, :RMS, group = :OrderofClick, bar_position = :dodge, lab= ["Peak 1" "Peak 2" "Peak 3" "Peak 4" "Peak 5"], title=Title, xlab = "Measurement", ylab = "Peak Amplitude");

        name = joinpath(figs,foldername,Title) * ".png"
        savefig(avrecplot, name);
    end
end


function AvrecPeakRatio(figs,Tab2,Tab5,trialtype="TA")
    # Input: folder path figs, table for ratio of 2 hz and 5 hz, all groups being plotted
    # Output: figures in folder AvrecPeakRatio of the level of synaptic depression shown as a function of the last divided by the first response peak amplitude
 
    foldername = "AvrecPeakRatio"
    if !isdir(joinpath(figs,foldername))
        mkdir(joinpath(figs,foldername))
    end

    LayList = unique(Tab2[!,:Layer])

    for iLay = 1:length(LayList)
        ### peak amp by measurement per Layer ###
        Tab2_Sort = Tab2[Tab2[!,:Layer] .== LayList[iLay],:]
        Tab5_Sort = Tab5[Tab5[!,:Layer] .== LayList[iLay],:]

        Title = "Synaptic dep ratio of " * LayList[iLay] * " at 2 Hz "
        ratioplot = @df Tab2_Sort groupedboxplot(:Measurement, :Ratio, group = :Group, bar_position = :dodge, lab= ["Control" "Treated" "Virus Control"], title=Title, xlab = "Measurement", ylab = "Ratio of last to first Peak Amp")

        name = joinpath(figs,foldername,Title) * trialtype * ".png"
        savefig(ratioplot, name);

        Title = "Synaptic dep ratio of " * LayList[iLay] * " at 5 Hz "
        avrecplot = @df Tab5_Sort groupedboxplot(:Measurement, :Ratio, group = :Group, bar_position = :dodge, lab= ["Control" "Treated" "Virus Control"], title=Title, xlab = "Measurement", ylab = "Ratio of last to first Peak Amp")

        name = joinpath(figs,foldername,Title) * trialtype * ".png"
        savefig(avrecplot, name);
    end
end

function Avrec1Peak(figs,Tab,whichpeak="First",whichstim="2Hz",trialtype="TA")
    # Input: folder path figs, table for peak amp/lat of 2 hz and 5 hz, all groups being plotted
    # Output: figures in folder Avrec1stPeak of the first peak response and latency for all groups

    foldername = "Avrec1Peak"
    if !isdir(joinpath(figs,foldername))
        mkdir(joinpath(figs,foldername))
    end

    LayList = unique(Tab[!,:Layer])

    for iLay = 1:length(LayList)
        ### peak amp by measurement per Layer ###
        Tab_Sort = Tab[Tab[!,:Layer] .== LayList[iLay],:]

        Title = whichpeak * " peak amplitude of " * LayList[iLay] * " at " * whichstim * " "
        ratioplot = @df Tab_Sort groupedboxplot(:Measurement, :PeakAmp, group = :Group, bar_position = :dodge, lab= ["Control" "Treated" "Virus Control"], title=Title, xlab = "Measurement", ylab = "Ratio of last to first Peak Amp")

        name = joinpath(figs,foldername,Title) * trialtype * ".png"
        savefig(ratioplot, name);

        Title = whichpeak * " peak latency of " * LayList[iLay] * " at " * whichstim * " "
        ratioplot = @df Tab_Sort groupedboxplot(:Measurement, :PeakLat, group = :Group, bar_position = :dodge, lab= ["Control" "Treated" "Virus Control"], title=Title, xlab = "Measurement", ylab = "Ratio of last to first Peak Amp")

        name = joinpath(figs,foldername,Title) * trialtype * ".png"
        savefig(ratioplot, name);

        Title = whichpeak * " RMS of " * LayList[iLay] * " at " * whichstim * " "
        ratioplot = @df Tab_Sort groupedboxplot(:Measurement, :RMS, group = :Group, bar_position = :dodge, lab= ["Control" "Treated" "Virus Control"], title=Title, xlab = "Measurement", ylab = "Ratio of last to first Peak Amp")

        name = joinpath(figs,foldername,Title) * trialtype * ".png"
        savefig(ratioplot, name);

    end
end