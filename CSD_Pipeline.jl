### This scrip runs everything necessary to get a full CSD dataset out of the plexon-mat converted files
# Input:    converted raw data and callGroup function from /mouse_input (e.g. KIC_0006.mat) 
# Output:   jld2 file type with a Dictionary per condition/measurement. Each dictionary contains the average and single-trail LFP, CSD, AVREC, ABSREC, and RELSRES

# Set up
using Statistics, DSP
using Plots, Colors
using MAT, CSV, JLD2, FileIO
using DataFrames, OrderedCollections

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
include(joinpath(group,"callGroupAw.jl"))

# determine data to read -- this may be a more complicated process later
GroupList = ["CIC"]
CondName = ["StartTono" "StartSpont" "preCL" "preAMBF" "preAMoBF" "pre1Spont" "post1Spont" "postCL" "postAMBF" "pre2Spont" "post2Spont" "postAMoBF" "EndTono"] # "StartTono" "StartSpont" "preCL" "preAMBF" "preAMoBF" "pre1Spont" "post1Spont" "postCL" "postAMBF" "pre2Spont" "post2Spont" "postAMoBF" "EndTono" 

# loop through groups
for iGr = GroupList # for testing: iGr = "CIC"
    # generate lists of channels, layers, and measurements for each animal in this group
    animalList,chanList,LIIList,LIVList,LVList,LVIList,CondList = callGroupAw(iGr); # in groups folder
    # loop through each animal 
    Animal = Dict()
    for iAn = 1:length(animalList)

        AnimalName = animalList[iAn]

        channels = chanList[iAn]
        LII  = LIIList[iAn]
        LIV  = LIVList[iAn]
        LV   = LVList[iAn]
        LVI  = LVIList[iAn]
        
        csdStruct(raw,figs,datap,Animal,AnimalName,iGr,CondList,CondName,channels,LII,LIV,LV,LVI,iAn)

    end
end

## How to load back out into workspace:
AnimalName = "CIC01"
varname  = AnimalName * "_Data"
filename = joinpath(datap,AnimalName) * "_Data.jld2"
Datout = load(filename)[varname]
