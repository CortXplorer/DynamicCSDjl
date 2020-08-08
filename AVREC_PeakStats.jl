using Statistics, DSP
using Colors, StatsPlots
using CSV
using DataFrames

home    = @__DIR__
func    = joinpath(home,"functions")
figs    = joinpath(home,"figs")
include(joinpath(func,"AvrecPeakPlots.jl"))

# Load in data from matlab table csv file which contains 2 and 5 hz peak amp and latency. 
PeakData = CSV.File("AVRECPeakData.csv") |> DataFrame

KIC = PeakData[PeakData[!,:Group] .== "KIC",:]
KIT = PeakData[PeakData[!,:Group] .== "KIT",:]
KIV = PeakData[PeakData[!,:Group] .== "KIV",:]

KIC2, KIC5 = KIC[KIC[!,:ClickFreq] .== 2,:], KIC[KIC[!,:ClickFreq] .== 5,:]
KIT2, KIT5 = KIT[KIT[!,:ClickFreq] .== 2,:], KIT[KIT[!,:ClickFreq] .== 5,:]
KIV2, KIV5 = KIV[KIV[!,:ClickFreq] .== 2,:], KIV[KIV[!,:ClickFreq] .== 5,:]

## Plots First 
AvrecPeakvsMeas(figs,KIC2,KIC5,"KIC")
AvrecPeakvsMeas(figs,KIT2,KIT5,"KIT")
AvrecPeakvsMeas(figs,KIV2,KIV5,"KIV")
AvrecPeakvsLay(figs,KIC2,KIC5,"KIC")
AvrecPeakvsLay(figs,KIT2,KIT5,"KIT")
AvrecPeakvsLay(figs,KIV2,KIV5,"KIV")

## Now Stats