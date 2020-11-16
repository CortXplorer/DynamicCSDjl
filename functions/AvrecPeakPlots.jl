function AvrecPeakvsLay(figs,Tab,GroupName,whichstim="2Hz",savetype=".pdf")
    # Input: folder path figs, table for 2 hz and 5 hz, current group being plotted
    # Output: figures in folder AvrecPeakPlots_againstLayer of Peak Amplitude over layer per measurement condition

    foldername = "AvrecPeakPlots_againstLayer"
    if !isdir(joinpath(figs,foldername))
        mkdir(joinpath(figs,foldername))
    end

    MeasList = unique(Tab[:,:Measurement])

    for iMeas = 1:length(MeasList)
        ### peak amp by layer per measurement ###
        Tab_Sort = Tab[Tab[:,:Measurement] .== MeasList[iMeas],:]
        # take out nan peak amp rows
        Tab_Sort = filter(row -> ! isnan(row.PeakAmp), Tab_Sort)

        Title = GroupName * " PeakAmp of " * MeasList[iMeas] * " at " * whichstim
        avrecplot = @df Tab_Sort groupedboxplot(:Layer, :PeakAmp, group = :OrderofClick, bar_position = :dodge, lab= [1:parse(Int,whichstim[begin:end-2])...], title=Title, xlab = "Layer", ylab = "Peak Amplitude [mV/mm²]");

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

    LayList = unique(Tab[:,:Layer])

    for iLay = 1:length(LayList)
        ### peak amp by measurement per Layer ###
        Tab_Sort    = Tab[Tab[:,:Layer] .== LayList[iLay],:]
        Tab_Sortamp = filter(row -> ! isnan(row.PeakAmp), Tab_Sort)
        Tab_Sortrms = filter(row -> ! isnan(row.RMS), Tab_Sort)

        Title = GroupName * " PeakAmp of " * LayList[iLay] * " at " * whichstim * " " * stimtype
        avrecplot = @df Tab_Sortamp groupedboxplot(:Measurement, :PeakAmp, group = :OrderofClick, bar_position = :dodge, legend=false, title=Title, xlab = "Measurement", ylab = "Peak Amplitude [mV/mm²]"); # lab = [1:parse(Int,whichstim[begin:end-2])...] for legend 

        name = joinpath(figs,foldername,Title) * savetype
        savefig(avrecplot, name);

        Title = GroupName * " Peak Latency of " * LayList[iLay] * " at " * whichstim * " " * stimtype
        avrecplot = @df Tab_Sortamp groupedboxplot(:Measurement, :PeakLat, group = :OrderofClick, bar_position = :dodge, legend=false, title=Title, xlab = "Measurement", ylab = "Peak Latency [ms]");

        name = joinpath(figs,foldername,Title) * savetype
        savefig(avrecplot, name);

        Title = GroupName * " RMS of " * LayList[iLay] * " at " * whichstim * " " * stimtype
        avrecplot = @df Tab_Sortrms groupedboxplot(:Measurement, :RMS, group = :OrderofClick, bar_position = :dodge, legend=false, title=Title, xlab = "Measurement", ylab = "Root Mean Square [mV/mm²]");

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

    LayList = unique(Tab[:,:Layer])

    for iLay = 1:length(LayList)
        ### peak amp by measurement per Layer ###
        Tab_Sort = Tab[Tab[:,:Layer] .== LayList[iLay],:]
        Tab_Sortamp = filter(row -> ! isnan(row.RatioAMP), Tab_Sort)
        Tab_Sortrms = filter(row -> ! isnan(row.RatioRMS), Tab_Sort)

        Title = "Synaptic dep ratio of " * LayList[iLay] * " at " * whichstim * " Peak Amplitude " * stimtype
        ratioplot = @df Tab_Sortamp groupedboxplot(:Measurement, :RatioAMP, group = :Group, bar_position = :dodge, lab= ["Control" "Treated" "Virus Control"], title=Title, xlab = "Measurement", ylab = "Ratio of last to first Peak Amplitude [%]")

        if maximum(Tab_Sortamp.RatioAMP) > 500
            plot!(ylims=(0,500))
        end

        name = joinpath(figs,foldername,Title) * " " * trialtype * savetype
        savefig(ratioplot, name);

        Title = "Synaptic dep ratio of " * LayList[iLay] * " at " * whichstim * " RMS " * stimtype
        ratioplot = @df Tab_Sortrms groupedboxplot(:Measurement, :RatioRMS, group = :Group, bar_position = :dodge, lab= ["Control" "Treated" "Virus Control"], title=Title, xlab = "Measurement", ylab = "Ratio of last to first Peak RMS [%]")

        if maximum(Tab_Sortrms.RatioRMS) > 500
            plot!(ylims=(0,500))
        end

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

    LayList = unique(Tab[:,:Layer])

    for iLay = 1:length(LayList)
        ### peak amp by measurement per Layer ###
        Tab_Sort = Tab[Tab[:,:Layer] .== LayList[iLay],:]
        Tab_Sortamp = filter(row -> ! isnan(row.PeakAmp), Tab_Sort)
        Tab_Sortrms = filter(row -> ! isnan(row.RMS), Tab_Sort)

        Title = whichpeak * " peak amplitude of " * LayList[iLay] * " at " * whichstim * " " * stimtype * " "
        ratioplot = @df Tab_Sortamp groupedboxplot(:Measurement, :PeakAmp, group = :Group, bar_position = :dodge, lab= ["Control" "Treated" "Virus Control"], title=Title, xlab = "Measurement", ylab = "Peak Amplitude [mV/mm²]")

        name = joinpath(figs,foldername,Title) * trialtype * savetype
        savefig(ratioplot, name);

        Title = whichpeak * " peak latency of " * LayList[iLay] * " at " * whichstim * " " * stimtype * " "
        ratioplot = @df Tab_Sortamp groupedboxplot(:Measurement, :PeakLat, group = :Group, bar_position = :dodge, lab= ["Control" "Treated" "Virus Control"], title=Title, xlab = "Measurement", ylab = "Peak Latency [ms]")

        name = joinpath(figs,foldername,Title) * trialtype * savetype
        savefig(ratioplot, name);

        Title = whichpeak * " RMS of " * LayList[iLay] * " at " * whichstim * " " * stimtype * " "
        ratioplot = @df Tab_Sortrms groupedboxplot(:Measurement, :RMS, group = :Group, bar_position = :dodge, lab= ["Control" "Treated" "Virus Control"], title=Title, xlab = "Measurement", ylab = "RMS [mV/mm²]")

        name = joinpath(figs,foldername,Title) * trialtype * savetype
        savefig(ratioplot, name);

    end
end

function AvrecScatter(figs,Scat,whichstim="2Hz",savetype=".pdf",stimtype="CL",trialtype="TA")

    foldername = "AvrecScatter"
    if !isdir(joinpath(figs,foldername))
        mkdir(joinpath(figs,foldername))
    end

    MeasList = unique(Scat[:,:Measurement])
    LayList  = unique(Scat[:,:Layer])
    for iMeas = 1:length(MeasList)
        ### per measurement ###
        Scat_Meas = Scat[Scat[:,:Measurement] .== MeasList[iMeas],:]
        
        for iLay = 1:length(LayList)
            Scat_Lay = Scat_Meas[Scat_Meas[:,:Layer] .== LayList[iLay],:]
            Scat_Lay = filter(row -> ! isnan(row.PeakAmp), Scat_Lay)
            
            # correct peak latency time to be after each stim start time 
            startwin = Int.([0:1000/parse(Int,whichstim[begin:end-2]):1000...])

            for iWin = 1:length(startwin) - 1
                Scat_Lay[Scat_Lay[:,:OrderofClick].== iWin,:PeakLat] = Scat_Lay[Scat_Lay[:,:OrderofClick] .== iWin,:PeakLat] .+ startwin[iWin]
            end

            Title = "PeakAmp against Latency " * LayList[iLay] * " " * MeasList[iMeas] * " at " * whichstim * " " * stimtype * " " * trialtype
            scatterplot = @df Scat_Lay scatter(:PeakLat, :PeakAmp, group = :Group, markersize=3, markerstrokewidth=0, markerstrokealpha=0, markerstrokecolor = :tab10)
        
            name = joinpath(figs,foldername,Title) * savetype
            savefig(scatterplot, name);
            
        end # layer
    end # measurement
end # function