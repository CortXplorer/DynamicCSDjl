using CSV, DataFrames
using StatsPlots
using HypothesisTests
using Infiltrator

home    = @__DIR__
func    = joinpath(home,"functions")
figs    = joinpath(home,"figs")
data    = joinpath(home,"Data")

function AM_Temporal(figs,Tab,GroupName,whichstim="2Hz",savetype=".pdf",stimtype="AM")
    foldername = "AM_TemporalFlow"
    if !isdir(joinpath(figs,foldername))
        mkdir(joinpath(figs,foldername))
    end

    LayList = unique(Tab[!,:Layer])

    for iLay = 1:length(LayList)
        ### peak amp by measurement per Layer ###
        Tab_Sort    = Tab[Tab[!,:Layer] .== LayList[iLay],:]
        Tab_Sortrms = filter(row -> ! isnan(row.RMS), Tab_Sort)

        Title = GroupName * " RMS of " * LayList[iLay] * " at " * whichstim * " " * stimtype
        avrecplot = @df Tab_Sortrms groupedboxplot(:Measurement, :RMS, group = :Which100ms, bar_position = :dodge, legend=false, title=Title, xlab = "Measurement", ylab = "Root Mean Square [mV/mm²]");

        name = joinpath(figs,foldername,Title) * savetype
        savefig(avrecplot, name);
    end
end

AMDat = CSV.File("AM_RMS100msx10.csv") |> DataFrame # trial average

savetype = ".png"

KIC = AMDat[AMDat[!,:Group] .== "KIC",:]
KIT = AMDat[AMDat[!,:Group] .== "KIT",:]
# further seperate by stimulus
KIC2, KIC5, KIC10 = KIC[KIC[!,:ClickFreq] .== 2,:], KIC[KIC[!,:ClickFreq] .== 5,:], KIC[KIC[!,:ClickFreq] .== 10,:]
KIT2, KIT5, KIT10 = KIT[KIT[!,:ClickFreq] .== 2,:], KIT[KIT[!,:ClickFreq] .== 5,:], KIT[KIT[!,:ClickFreq] .== 10,:]

AM_Temporal(figs,KIC2,"KIC","2Hz",savetype)
AM_Temporal(figs,KIC5,"KIC","5Hz",savetype)
AM_Temporal(figs,KIC10,"KIC","10Hz",savetype)

AM_Temporal(figs,KIT2,"KIT","2Hz",savetype)
AM_Temporal(figs,KIT5,"KIT","5Hz",savetype)
AM_Temporal(figs,KIT10,"KIT","10Hz",savetype)