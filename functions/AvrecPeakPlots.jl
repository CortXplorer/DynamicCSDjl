function AvrecPeakvsLay(figs,Tab,GroupName,whichstim="2Hz",savetype=".pdf")
    # Input: folder path figs, table for 2 hz and 5 hz, current group being plotted
    # Output: figures in folder AvrecPeakPlots_againstLayer of Peak Amplitude over layer per measurement condition

    foldername = "AvrecPeakPlots_againstLayer"
    if !isdir(joinpath(figs,foldername))
        mkdir(joinpath(figs,foldername))
    end

    MeasList = unique(Tab[!,:Measurement])

    if whichstim == "2Hz"
        labels = ["Peak 1" "Peak 2"]
    elseif whichstim == "5Hz"
        labels = ["Peak 1" "Peak 2" "Peak 3" "Peak 4" "Peak 5"]
    elseif whichstim == "10Hz"
        labels = ["Peak 1" "Peak 2" "Peak 3" "Peak 4" "Peak 5" "Peak 6" "Peak 7" "Peak 8" "Peak 9" "Peak 10"]
    end

    for iMeas = 1:length(MeasList)
        ### peak amp by layer per measurement ###
        Tab_Sort = Tab[Tab[!,:Measurement] .== MeasList[iMeas],:]
        # take out nan peak amp rows
        Tab_Sort = filter(row -> ! isnan(row.PeakAmp), Tab_Sort)

        Title = GroupName * " PeakAmp of " * MeasList[iMeas] * " at " * whichstim
        avrecplot = @df Tab_Sort groupedboxplot(:Layer, :PeakAmp, group = :OrderofClick, bar_position = :dodge, lab= labels, title=Title, xlab = "Layer", ylab = "Peak Amplitude");

        name = joinpath(figs,foldername,Title) * savetype
        savefig(avrecplot, name);

    end
end

function AvrecPeakvsMeas(figs,Tab,GroupName,whichstim="2Hz",savetype=".pdf",stimtype="CL")
    # Input: folder path figs, table for 2 hz and 5 hz, current group being plotted
    # Output: figures in folder AvrecPeakPlots_againstMeasurement of Peak Amplitude over measurement per layer

    foldername = "AvrecPeakPlots_againstMeasurement"
    if !isdir(joinpath(figs,foldername))
        mkdir(joinpath(figs,foldername))
    end

    LayList = unique(Tab[!,:Layer])

    if whichstim == "2Hz"
        labels = ["Peak 1" "Peak 2"]
    elseif whichstim == "5Hz"
        labels = ["Peak 1" "Peak 2" "Peak 3" "Peak 4" "Peak 5"]
    elseif whichstim == "10Hz"
        labels = ["Peak 1" "Peak 2" "Peak 3" "Peak 4" "Peak 5" "Peak 6" "Peak 7" "Peak 8" "Peak 9" "Peak 10"]
    end

    for iLay = 1:length(LayList)
        ### peak amp by measurement per Layer ###
        Tab_Sort    = Tab[Tab[!,:Layer] .== LayList[iLay],:]
        Tab_Sortamp = filter(row -> ! isnan(row.PeakAmp), Tab_Sort)

        Title = GroupName * " PeakAmp of " * LayList[iLay] * " at " * whichstim * " " * stimtype
        avrecplot = @df Tab_Sortamp groupedboxplot(:Measurement, :PeakAmp, group = :OrderofClick, bar_position = :dodge, lab= labels, title=Title, xlab = "Measurement", ylab = "Peak Amplitude");

        name = joinpath(figs,foldername,Title) * savetype
        savefig(avrecplot, name);

        name = joinpath(figs,foldername,Title) * savetype
        savefig(avrecplot, name);

        Title = GroupName * " Peak Latency of " * LayList[iLay] * " at " * whichstim * " " * stimtype
        avrecplot = @df Tab_Sortamp groupedboxplot(:Measurement, :PeakLat, group = :OrderofClick, bar_position = :dodge, lab= labels, title=Title, xlab = "Measurement", ylab = "Peak Amplitude");

        name = joinpath(figs,foldername,Title) * savetype
        savefig(avrecplot, name);

        Title = GroupName * " RMS of " * LayList[iLay] * " at " * whichstim * " " * stimtype
        avrecplot = @df Tab_Sort groupedboxplot(:Measurement, :RMS, group = :OrderofClick, bar_position = :dodge, lab= labels, title=Title, xlab = "Measurement", ylab = "Peak Amplitude");

        name = joinpath(figs,foldername,Title) * savetype
        savefig(avrecplot, name);
    end
end


function AvrecPeakRatio(figs,Tab,whichstim="2Hz",savetype=".pdf",stimtype="CL",trialtype="TA")
    # Input: folder path figs, table for ratio of 2 hz and 5 hz, all groups being plotted
    # Output: figures in folder AvrecPeakRatio of the level of synaptic depression shown as a function of the last divided by the first response peak amplitude
 
    foldername = "AvrecPeakRatio"
    if !isdir(joinpath(figs,foldername))
        mkdir(joinpath(figs,foldername))
    end

    LayList = unique(Tab[!,:Layer])

    for iLay = 1:length(LayList)
        ### peak amp by measurement per Layer ###
        Tab_Sort = Tab[Tab[!,:Layer] .== LayList[iLay],:]
        Tab_Sortamp = filter(row -> ! isnan(row.RatioAMP), Tab_Sort)
        Tab_Sortrms = filter(row -> ! isnan(row.RatioRMS), Tab_Sort)

        Title = "Synaptic dep ratio of " * LayList[iLay] * " at " * whichstim * " Peak Amplitude " * stimtype
        ratioplot = @df Tab_Sortamp groupedboxplot(:Measurement, :RatioAMP, group = :Group, bar_position = :dodge, lab= ["Control" "Treated" "Virus Control"], title=Title, xlab = "Measurement", ylab = "Ratio of last to first Peak Amp")

        name = joinpath(figs,foldername,Title) * " " * trialtype * savetype
        savefig(ratioplot, name);

        Title = "Synaptic dep ratio of " * LayList[iLay] * " at " * whichstim * " RMS " * stimtype
        ratioplot = @df Tab_Sortrms groupedboxplot(:Measurement, :RatioRMS, group = :Group, bar_position = :dodge, lab= ["Control" "Treated" "Virus Control"], title=Title, xlab = "Measurement", ylab = "Ratio of last to first Peak RMS")

        name = joinpath(figs,foldername,Title) * " " * trialtype * savetype
        savefig(ratioplot, name);
    end
end

function Avrec1Peak(figs,Tab,whichpeak="1st",whichstim="2Hz",savetype=".pdf",stimtype="CL",trialtype="TA")
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
        Tab_Sortamp = filter(row -> ! isnan(row.PeakAmp), Tab_Sort)
        Tab_Sortrms = filter(row -> ! isnan(row.RMS), Tab_Sort)

        Title = whichpeak * " peak amplitude of " * LayList[iLay] * " at " * whichstim * " " * stimtype * " "
        ratioplot = @df Tab_Sortamp groupedboxplot(:Measurement, :PeakAmp, group = :Group, bar_position = :dodge, lab= ["Control" "Treated" "Virus Control"], title=Title, xlab = "Measurement", ylab = "Ratio of last to first Peak Amp")

        name = joinpath(figs,foldername,Title) * trialtype * savetype
        savefig(ratioplot, name);

        Title = whichpeak * " peak latency of " * LayList[iLay] * " at " * whichstim * " " * stimtype * " "
        ratioplot = @df Tab_Sortamp groupedboxplot(:Measurement, :PeakLat, group = :Group, bar_position = :dodge, lab= ["Control" "Treated" "Virus Control"], title=Title, xlab = "Measurement", ylab = "Ratio of last to first Peak Amp")

        name = joinpath(figs,foldername,Title) * trialtype * savetype
        savefig(ratioplot, name);

        Title = whichpeak * " RMS of " * LayList[iLay] * " at " * whichstim * " " * stimtype * " "
        ratioplot = @df Tab_Sortrms groupedboxplot(:Measurement, :RMS, group = :Group, bar_position = :dodge, lab= ["Control" "Treated" "Virus Control"], title=Title, xlab = "Measurement", ylab = "Ratio of last to first Peak Amp")

        name = joinpath(figs,foldername,Title) * trialtype * savetype
        savefig(ratioplot, name);

    end
end

function AvrecScatter(figs,Scat,whichstim="2Hz",savetype=".pdf",stimtype="CL",trialtype="TA")

    foldername = "AvrecScatter"
    if !isdir(joinpath(figs,foldername))
        mkdir(joinpath(figs,foldername))
    end

    MeasList = unique(Scat[!,:Measurement])
    LayList  = unique(Scat[!,:Layer])
    for iMeas = 1:length(MeasList)
        ### per measurement ###
        Scat_Meas = Scat[Scat[!,:Measurement] .== MeasList[iMeas],:]
        
        for iLay = 1:length(LayList)
            Scat_Lay = Scat_Meas[Scat_Meas[!,:Layer] .== LayList[iLay],:]
            Scat_Lay = filter(row -> ! isnan(row.PeakAmp), Scat_Lay)
            # edit peak latency time to be after each stim time 
            if whichstim == "2Hz"
                Scat_Lay[Scat_Lay[!,:OrderofClick].==2,:PeakLat] = Scat_Lay[Scat_Lay[!,:OrderofClick] .== 2,:PeakLat] .+ 500
            elseif whichstim == "5Hz"
                Scat_Lay[Scat_Lay[!,:OrderofClick].==2,:PeakLat] = Scat_Lay[Scat_Lay[!,:OrderofClick] .== 2,:PeakLat] .+ 200
                Scat_Lay[Scat_Lay[!,:OrderofClick].==3,:PeakLat] = Scat_Lay[Scat_Lay[!,:OrderofClick] .== 3,:PeakLat] .+ 400
                Scat_Lay[Scat_Lay[!,:OrderofClick].==4,:PeakLat] = Scat_Lay[Scat_Lay[!,:OrderofClick] .== 4,:PeakLat] .+ 600
                Scat_Lay[Scat_Lay[!,:OrderofClick].==5,:PeakLat] = Scat_Lay[Scat_Lay[!,:OrderofClick] .== 5,:PeakLat] .+ 800
            elseif whichstim == "10Hz"
                Scat_Lay[Scat_Lay[!,:OrderofClick].==2,:PeakLat] = Scat_Lay[Scat_Lay[!,:OrderofClick] .== 2,:PeakLat] .+ 200
                Scat_Lay[Scat_Lay[!,:OrderofClick].==3,:PeakLat] = Scat_Lay[Scat_Lay[!,:OrderofClick] .== 3,:PeakLat] .+ 300
                Scat_Lay[Scat_Lay[!,:OrderofClick].==4,:PeakLat] = Scat_Lay[Scat_Lay[!,:OrderofClick] .== 4,:PeakLat] .+ 400
                Scat_Lay[Scat_Lay[!,:OrderofClick].==5,:PeakLat] = Scat_Lay[Scat_Lay[!,:OrderofClick] .== 5,:PeakLat] .+ 500
                Scat_Lay[Scat_Lay[!,:OrderofClick].==3,:PeakLat] = Scat_Lay[Scat_Lay[!,:OrderofClick] .== 3,:PeakLat] .+ 600
                Scat_Lay[Scat_Lay[!,:OrderofClick].==4,:PeakLat] = Scat_Lay[Scat_Lay[!,:OrderofClick] .== 4,:PeakLat] .+ 700
                Scat_Lay[Scat_Lay[!,:OrderofClick].==5,:PeakLat] = Scat_Lay[Scat_Lay[!,:OrderofClick] .== 5,:PeakLat] .+ 800
                Scat_Lay[Scat_Lay[!,:OrderofClick].==3,:PeakLat] = Scat_Lay[Scat_Lay[!,:OrderofClick] .== 3,:PeakLat] .+ 900
                Scat_Lay[Scat_Lay[!,:OrderofClick].==4,:PeakLat] = Scat_Lay[Scat_Lay[!,:OrderofClick] .== 4,:PeakLat] .+ 1000
                Scat_Lay[Scat_Lay[!,:OrderofClick].==5,:PeakLat] = Scat_Lay[Scat_Lay[!,:OrderofClick] .== 5,:PeakLat] .+ 1100
            end

            Title = "PeakAmp against Latency " * LayList[iLay] * " " * MeasList[iMeas] * " at " * whichstim * " " * stimtype * " " * trialtype
            scatterplot = @df Scat_Lay scatter(:PeakLat, :PeakAmp, group = :Group, markersize=3, markerstrokewidth=0, markerstrokealpha=0, markerstrokecolor = :tab10)
        
            name = joinpath(figs,foldername,Title) * savetype
            savefig(scatterplot, name);
            
        end # layer
    end # measurement
end # function