using Statistics, DSP
using Colors, StatsPlots
using CSV
using DataFrames

home    = @__DIR__
func    = joinpath(home,"functions")
figs    = joinpath(home,"figs")

foldername = "AvrecPeakPlots"
if !isdir(joinpath(figs,foldername))
    mkdir(joinpath(figs,foldername))
end

# Load in data from matlab table csv file which contains 2 and 5 hz peak amp and latency. 
PeakData = CSV.File("AVRECPeakData.csv") |> DataFrame

KIC = PeakData[PeakData[!,:Group] .== "KIC",:]
KIT = PeakData[PeakData[!,:Group] .== "KIT",:]
KIV = PeakData[PeakData[!,:Group] .== "KIV",:]

KIC2, KIC5 = KIC[KIC[!,:ClickFreq] .== 2,:], KIC[KIC[!,:ClickFreq] .== 5,:]
KIT2, KIT5 = KIT[KIT[!,:ClickFreq] .== 2,:], KIT[KIT[!,:ClickFreq] .== 5,:]
KIV2, KIV5 = KIV[KIV[!,:ClickFreq] .== 2,:], KIV[KIV[!,:ClickFreq] .== 5,:]

# peak amp by layer per measurement
KIC2_PreCL = KIC2[KIC2[!,:Measurement] .== "preCL_1",:]
KIC5_PreCL = KIC5[KIC5[!,:Measurement] .== "preCL_1",:]

Title = "PeakAmp of PreCL 2 Hz"
avrecplot = @df KIC2_PreCL groupedboxplot(:Layer, :PeakAmp, group = :OrderofClick, bar_position = :dodge, lab= ["Peak 1" "Peak 2"], title=Title, xlab = "Layer", ylab = "Peak Amplitude")

name = joinpath(figs,foldername,Title) * ".pdf"
savefig(avrecplot, name)

Title = "PeakAmp of PreCL 5 Hz"
@df KIC5_PreCL groupedboxplot(:Layer, :PeakAmp, group = :OrderofClick, bar_position = :dodge, lab= ["Peak 1" "Peak 2"], title=Title, xlab = "Layer", ylab = "Peak Amplitude")

name = joinpath(figs,foldername,Title) * ".pdf"
savefig(avrecplot, name)