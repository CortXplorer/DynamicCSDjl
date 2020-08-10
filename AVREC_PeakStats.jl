using StatsPlots
using CSV, DataFrames

home    = @__DIR__
func    = joinpath(home,"functions")
figs    = joinpath(home,"figs")
include(joinpath(func,"AvrecPeakPlots.jl"))

# Load in data from matlab table csv file which contains 2 and 5 hz peak amp and latency. 
PeakData = CSV.File("AVRECPeakData.csv") |> DataFrame

# seperate by group
KIC = PeakData[PeakData[!,:Group] .== "KIC",:]
KIT = PeakData[PeakData[!,:Group] .== "KIT",:]
KIV = PeakData[PeakData[!,:Group] .== "KIV",:]
# further seperate by stimulus
KIC2, KIC5 = KIC[KIC[!,:ClickFreq] .== 2,:], KIC[KIC[!,:ClickFreq] .== 5,:]
KIT2, KIT5 = KIT[KIT[!,:ClickFreq] .== 2,:], KIT[KIT[!,:ClickFreq] .== 5,:]
KIV2, KIV5 = KIV[KIV[!,:ClickFreq] .== 2,:], KIV[KIV[!,:ClickFreq] .== 5,:]

## Box Plots First 
AvrecPeakvsMeas(figs,KIC2,KIC5,"KIC")
AvrecPeakvsMeas(figs,KIT2,KIT5,"KIT")
AvrecPeakvsMeas(figs,KIV2,KIV5,"KIV")
AvrecPeakvsLay(figs,KIC2,KIC5,"KIC")
AvrecPeakvsLay(figs,KIT2,KIT5,"KIT")
AvrecPeakvsLay(figs,KIV2,KIV5,"KIV")

## Now Stats
# J. Heck, in her paper, used EPSP5/EPSP1 to show the difference between groups of the ratio from the last to first stimulus response
# seperate just stimulus presentation from full table
Stim2Hz = PeakData[PeakData[!,:ClickFreq] .== 2,:]
Stim5Hz = PeakData[PeakData[!,:ClickFreq] .== 5,:]
# seperate the 1st and last click
First2 = Stim2Hz[Stim2Hz[!,:OrderofClick] .== 1,:]
Last2  = Stim2Hz[Stim2Hz[!,:OrderofClick] .== 2,:]
First5 = Stim5Hz[Stim5Hz[!,:OrderofClick] .== 1,:]
Last5  = Stim5Hz[Stim5Hz[!,:OrderofClick] .== 5,:]
# divide the last by first
Ratio2 = Last2[!,:PeakAmp] ./ First2[!,:PeakAmp]
Ratio5 = Last5[!,:PeakAmp] ./ First5[!,:PeakAmp]
# add this column to the table to keep tags
First2.Ratio = Ratio2
First5.Ratio = Ratio5

AvrecPeakRatio(figs,First2,First5)


