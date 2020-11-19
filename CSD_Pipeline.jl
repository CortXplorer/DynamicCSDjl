### Set up
using Statistics, DSP
using Plots, Colors
using MAT, JLD, CSV
using DataFrames

# get to correct directory and label it home
home    = @__DIR__
raw     = joinpath(home,"raw")
func    = joinpath(home,"functions")
figs    = joinpath(home,"figs")
group   = joinpath(home,"groups")
datap   = joinpath(home,"Data")

include(joinpath(func,"csdStruct.jl"))
include(joinpath(func,"Dynamic_CSD.jl"))
include(joinpath(func,"SingleTrialCSD.jl"))
include(joinpath(func,"get_csd.jl")) # used in SingleTrialCSD.jl
include(joinpath(func,"sink_dura.jl"))
include(joinpath(func,"functions.jl"))
include(joinpath(group,"callGroup.jl"))

# determine data to read -- this will be a more complicated process later!
GroupList = ["KIC" "KIT"]
#"Pre" "preAM" "preAMtono" "preCL" "preCLtono" "spPre1" "spPost1" "CL" "CLtono" "spPre2" "spPost2" "AM" "AMtono" "spEnd"
CondName = ["Pre" "CL"]

# loop through groups
for iGr = 1:length(GroupList)
    # generate lists of channels, layers, and measurements for each animal in this group
    animalList,chanList,LIIList,LIVList,LVList,LVIList,CondList = callGroup(GroupList[iGr]) # in groups folder
    Group = GroupList[iGr]
    # loop through each animal (dictionary per animal)
    Animal = Dict()
    for iAn = 1:length(animalList)

        AnimalName = animalList[iAn]

        channels = chanList[iAn]
        LII  = LIIList[iAn]
        LIV  = LIVList[iAn]
        LV   = LVList[iAn]
        LVI  = LVIList[iAn]
        
        csdStruct(raw,figs,datap,Animal,AnimalName,Group,CondList,CondName,channels,LII,LIV,LV,LVI,iAn)

    end
end

# to load: Data  = load("D:\\DynamicCSDjl\\Data\\KIC02_Data.jld") -- NOT WORKING
