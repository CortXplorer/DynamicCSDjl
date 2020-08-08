function AvrecPeakvsLay(figs,Tab2,Tab5,GroupName)

    foldername = "AvrecPeakPlots_againstLayer"
    if !isdir(joinpath(figs,foldername))
        mkdir(joinpath(figs,foldername))
    end

    MeasList = unique(Tab2[!,:Measurement])

    for iMeas = 1:length(MeasList)
        ### peak amp by layer per measurement ###
        Tab2_PreCL = Tab2[Tab2[!,:Measurement] .== MeasList[iMeas],:]
        Tab5_PreCL = Tab5[Tab5[!,:Measurement] .== MeasList[iMeas],:]

        Title = GroupName * " PeakAmp of " * MeasList[iMeas] * " at 2 Hz"
        avrecplot = @df Tab2_PreCL groupedboxplot(:Layer, :PeakAmp, group = :OrderofClick, bar_position = :dodge, lab= ["Peak 1" "Peak 2"], title=Title, xlab = "Layer", ylab = "Peak Amplitude");

        name = joinpath(figs,foldername,Title) * ".pdf"
        savefig(avrecplot, name)

        Title = GroupName * " PeakAmp of " * MeasList[iMeas] * " at 5 Hz"
        @df Tab5_PreCL groupedboxplot(:Layer, :PeakAmp, group = :OrderofClick, bar_position = :dodge, lab= ["Peak 1" "Peak 2" "Peak 3" "Peak 4" "Peak 5"], title=Title, xlab = "Layer", ylab = "Peak Amplitude");

        name = joinpath(figs,foldername,Title) * ".pdf"
        savefig(avrecplot, name)
    end

end

function AvrecPeakvsMeas(figs,Tab2,Tab5,GroupName)

    foldername = "AvrecPeakPlots_againstMeasurement"
    if !isdir(joinpath(figs,foldername))
        mkdir(joinpath(figs,foldername))
    end

    LayList = unique(Tab2[!,:Layer])

    for iLay = 1:length(LayList)
        ### peak amp by layer per Layer ###
        Tab2_PreCL = Tab2[Tab2[!,:Layer] .== LayList[iLay],:]
        Tab5_PreCL = Tab5[Tab5[!,:Layer] .== LayList[iLay],:]

        Title = GroupName * " PeakAmp of " * LayList[iLay] * " at 2 Hz"
        avrecplot = @df Tab2_PreCL groupedboxplot(:Measurement, :PeakAmp, group = :OrderofClick, bar_position = :dodge, lab= ["Peak 1" "Peak 2"], title=Title, xlab = "Layer", ylab = "Peak Amplitude");

        name = joinpath(figs,foldername,Title) * ".pdf"
        savefig(avrecplot, name)

        Title = GroupName * " PeakAmp of " * LayList[iLay] * " at 5 Hz"
        @df Tab5_PreCL groupedboxplot(:Measurement, :PeakAmp, group = :OrderofClick, bar_position = :dodge, lab= ["Peak 1" "Peak 2" "Peak 3" "Peak 4" "Peak 5"], title=Title, xlab = "Layer", ylab = "Peak Amplitude");

        name = joinpath(figs,foldername,Title) * ".pdf"
        savefig(avrecplot, name)
    end

end