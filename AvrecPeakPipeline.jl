using CSV, DataFrames
using StatsPlots
using HypothesisTests

home    = @__DIR__
func    = joinpath(home,"functions")
figs    = joinpath(home,"figs")
data    = joinpath(home,"Data")
include(joinpath(func,"AvrecPeakPlots.jl"))
include(joinpath(func,"AvrecPeakStats.jl"))

# Load in data from matlab table csv file which contains 2 and 5 hz peak amp and latency. 
PeakDataTA = CSV.File("AVRECPeakData.csv") |> DataFrame # trial average
PeakDataST = CSV.File("AVRECPeakDataST.csv") |> DataFrame # single trial

# seperate by group
KIC = PeakDataTA[PeakDataTA[!,:Group] .== "KIC",:]
KIT = PeakDataTA[PeakDataTA[!,:Group] .== "KIT",:]
KIV = PeakDataTA[PeakDataTA[!,:Group] .== "KIV",:]
# further seperate by stimulus
KIC2, KIC5 = KIC[KIC[!,:ClickFreq] .== 2,:], KIC[KIC[!,:ClickFreq] .== 5,:]
KIT2, KIT5 = KIT[KIT[!,:ClickFreq] .== 2,:], KIT[KIT[!,:ClickFreq] .== 5,:]
KIV2, KIV5 = KIV[KIV[!,:ClickFreq] .== 2,:], KIV[KIV[!,:ClickFreq] .== 5,:]

## Box Plots First ### 
# Output: figures in folder AvrecPeakPlots_againstMeasurement/_againstLayer of Peak Amplitude over measurement/layer per layer/measurement
AvrecPeakvsMeas(figs,KIC2,KIC5,"KIC")
AvrecPeakvsMeas(figs,KIT2,KIT5,"KIT")
AvrecPeakvsMeas(figs,KIV2,KIV5,"KIV")
AvrecPeakvsLay(figs,KIC2,KIC5,"KIC")
AvrecPeakvsLay(figs,KIT2,KIT5,"KIT")
AvrecPeakvsLay(figs,KIV2,KIV5,"KIV")


### Now Stats ###
# Average Trials ---------------------------------------------------------------------------
# seperate just stimulus presentation from full table
Stim2Hz = PeakDataTA[PeakDataTA[!,:ClickFreq] .== 2,:]
Stim5Hz = PeakDataTA[PeakDataTA[!,:ClickFreq] .== 5,:]

## Peak Amp response difference
# let's just check the first peak response amp and latency for now
Stat2  = Stim2Hz[Stim2Hz[!,:OrderofClick] .== 1,:]
Stat5  = Stim5Hz[Stim5Hz[!,:OrderofClick] .== 1,:]
# plot it first for group comparison 
Avrec1Peak(figs,Stat2,Stat5)
# now the stats;
Peak1_Between(data,Stat2,Stat5)
Peak1_Within(data,Stat2,Stat5)

### Peak Ratio of Last/First response
# -> J. Heck, in her paper, used EPSP5/EPSP1 to show the difference between groups of the ratio from the last to first stimulus response

# seperate the 1st and last click
Stat2  = Stim2Hz[Stim2Hz[!,:OrderofClick] .== 1,:]
Last2  = Stim2Hz[Stim2Hz[!,:OrderofClick] .== 2,:]
Stat5  = Stim5Hz[Stim5Hz[!,:OrderofClick] .== 1,:]
Last5  = Stim5Hz[Stim5Hz[!,:OrderofClick] .== 5,:]
# divide the last by first
Ratio2 = Last2[!,:PeakAmp] ./ Stat2[!,:PeakAmp]
Ratio5 = Last5[!,:PeakAmp] ./ Stat5[!,:PeakAmp]
# add this column to the table to keep tags
Stat2.Ratio = Ratio2
Stat5.Ratio = Ratio5

# cheeky plots first;
AvrecPeakRatio(figs,Stat2,Stat5)
# And stats for overlay
PeakRatio_Between(data,Stat2,Stat5)
PeakRatio_Within(data,Stat2,Stat5)


# Single Trials ----------------------------------------------------------------------------

# seperate just stimulus presentation from full table
Stim2Hz = PeakDataST[PeakDataST[!,:ClickFreq] .== 2,:]
Stim5Hz = PeakDataST[PeakDataST[!,:ClickFreq] .== 5,:]

## Peak Amp response difference
Stat2  = Stim2Hz[Stim2Hz[!,:OrderofClick] .== 1,:]
Stat5  = Stim5Hz[Stim5Hz[!,:OrderofClick] .== 1,:]
 
Avrec1stPeak(figs,Stat2,Stat5,"ST")
Peak1_Between(data,Stat2,Stat5,"ST")
Peak1_Within(data,Stat2,Stat5,"ST") # currently failing only in ST mode

### Peak Ratio of Last/First response
# seperate the 1st and last click
Stat2  = Stim2Hz[Stim2Hz[!,:OrderofClick] .== 1,:]
Last2  = Stim2Hz[Stim2Hz[!,:OrderofClick] .== 2,:]
Stat5  = Stim5Hz[Stim5Hz[!,:OrderofClick] .== 1,:]
Last5  = Stim5Hz[Stim5Hz[!,:OrderofClick] .== 5,:]
# divide the last by first
Ratio2 = Last2[!,:PeakAmp] ./ Stat2[!,:PeakAmp]
Ratio5 = Last5[!,:PeakAmp] ./ Stat5[!,:PeakAmp]
# add this column to the table to keep tags
Stat2.Ratio = Ratio2
Stat5.Ratio = Ratio5

AvrecPeakRatio(figs,Stat2,Stat5,"ST")
PeakRatio_Between(data,Stat2,Stat5,"ST")
PeakRatio_Within(data,Stat2,Stat5,"ST") # currently failing only in ST mode