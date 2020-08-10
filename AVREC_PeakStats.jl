using StatsPlots
using CSV, DataFrames
using HypothesisTests

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
# jk, cheeky plots first
AvrecPeakRatio(figs,Stat2,Stat5)
# Setup for simple 2 sample t test from HypothesisTests ->
# 3 between group comparisons across 5 measurement types and 5 layers:
CompList    = ["KIC vs KIT" "KIC vs KIV" "KIT vs KIV"]
LayList     = unique(Stat2[!,:Layer])
MeasList    = unique(Stat2[!,:Measurement])
# containers for the table columns: 
Comparison  = String[]
Layer       = String[]
Measurement = String[]   
P2Hz        = ones(length(CompList) * length(LayList) * length(MeasList)) # number of rows
P5Hz        = ones(length(CompList) * length(LayList) * length(MeasList))
count       = [1]
# loop through the above - create a table output! :) 
for iComp = 1:length(CompList)
    # pull out groups from ratio table
    G1_2Hz = Stat2[Stat2[!,:Group] .== CompList[iComp][1:3] ,:]
    G2_2Hz = Stat2[Stat2[!,:Group] .== CompList[iComp][8:10],:]
    G1_5Hz = Stat5[Stat5[!,:Group] .== CompList[iComp][1:3] ,:]
    G2_5Hz = Stat5[Stat5[!,:Group] .== CompList[iComp][8:10],:]

    for iMeas = 1:length(MeasList)
        # pull out the measurment for those groups
        G1_2Hz_meas = G1_2Hz[G1_2Hz[!,:Measurement] .== MeasList[iMeas],:]
        G2_2Hz_meas = G2_2Hz[G2_2Hz[!,:Measurement] .== MeasList[iMeas],:]
        G1_5Hz_meas = G1_5Hz[G1_5Hz[!,:Measurement] .== MeasList[iMeas],:]
        G2_5Hz_meas = G2_5Hz[G2_5Hz[!,:Measurement] .== MeasList[iMeas],:]

        for iLay = 1:length(LayList)
            # further pull out each layer definition
            G1_2Hz_lay = G1_2Hz_meas[G1_2Hz_meas[!,:Layer] .== LayList[iLay],:]
            G2_2Hz_lay = G2_2Hz_meas[G2_2Hz_meas[!,:Layer] .== LayList[iLay],:]
            G1_5Hz_lay = G1_5Hz_meas[G1_5Hz_meas[!,:Layer] .== LayList[iLay],:]
            G2_5Hz_lay = G2_5Hz_meas[G2_5Hz_meas[!,:Layer] .== LayList[iLay],:]

            # find the p value outcome for this comparison at both stimuli conditions
            P2Hz[count[1]] = pvalue(UnequalVarianceTTest(G1_2Hz_lay[!,:Ratio],G2_2Hz_lay[!,:Ratio]))
            P5Hz[count[1]] = pvalue(UnequalVarianceTTest(G1_5Hz_lay[!,:Ratio],G2_5Hz_lay[!,:Ratio]))
            count[1] = count[1] + 1
            # store the appropriate tags at the same positions in their respective lists
            push!(Comparison,CompList[iComp])
            push!(Measurement,MeasList[iMeas])
            push!(Layer,LayList[iLay])
        end # layer
    end # measurement type
end # comparison of which groups

BetweenGroup = DataFrame(Comparison=Comparison, Measurement=Measurement, Layer=Layer, P2Hz=P2Hz, P5Hz=P5Hz)