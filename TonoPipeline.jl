using MAT, StatsBase, DSP
using DataFrames, CSV
using Peaks

home    = @__DIR__
groups  = joinpath(home,"groups")
data    = joinpath(home,"Data")
func    = joinpath(home,"functions")

include(joinpath(groups,"callGroup.jl"))
include(joinpath(func,"TonotopyAnalysis.jl"))
include(joinpath(func,"nanmean.jl"))
include(joinpath(func,"ifindpeaks.jl"))

# function takes all tonotopies in current data (specific call by measurement names inside function) and performs feature extraction to compare with click and AM features 
# input:    .../Data/*_Data.mat; callGroup.jl
# output:   .../Data/Tonotopy/TonoRMS.csv
TonotopyAnalysis(data) 

#____________________________________________________________________________________
#____________________________________________________________________________________

using StatsPlots, DataFrames, CSV
TonoRMS = CSV.File(name) |> DataFrame

TonoProgressionAVREC = @df TonoRMS groupedboxplot(:Measurement,:AVRECrms, group = :Group, markerstrokewidth=0, legend = true, xrotation = 45)

TonoRMSIV = filter(row -> ! isnan(row.LIVrms), TonoRMS)

TonoProgressionLIV = @df TonoRMSIV groupedboxplot(:Measurement,:LIVrms, group = :Group, markerstrokewidth=0, legend = true, xrotation = 45)

