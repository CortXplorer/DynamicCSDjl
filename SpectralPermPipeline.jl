using MAT, Statistics, DSP, Colors
using JLD2, FileIO, CSV
using Random, DataFrames, StatsPlots
using HypothesisTests

home    = @__DIR__
data    = joinpath(home,"Data")
func    = joinpath(home,"functions")
group   = joinpath(home,"groups")
figs    = joinpath(home,"figs")
spect   = joinpath(data,"Spectral")

include(joinpath(group,"callGroup.jl"))
include(joinpath(func,"CWTfunc.jl"))

# Parameter
params = (
    sampleRate=1000, startTime=-200, timeLimits=[-200 1377], frequencyLimits=[8 100], timeBandWidth=40, 
    stimList=["twoHz","fiveHz","tenHz","twentyHz","fortyHz"], 
    layers=["I_II","IV","V","VI"], 
    osciBands=["alpha","beta low","beta high","gamma low","gamma high"], bandRanges=[[1:8...],[9:15...],[16:23...],[24:35...],[36:40...]])
# real ranges captured by these data: alpha=(8:12), beta_low=(13:18), beta_high=(19:30), gamma_low=(31:60), gamma_high=(61:100)
takepic = 1


## GROUP determination
GroupList = ["KIC" "KIT" "KIV"] 
## Conditions to run
CLList  = ["preCL" "CL" "preAM" "AM"] 
## Conditional picture; takepic == 1 if you do want figure output

# Loop through Groups to use function CWT_Loop which splits by animals, condition, measurement, stim frequency, and layer and outputs a dictionary with a power and phase coherence dataset for each chunk
KIC_WT = []
KIT_WT = []
KIV_WT = []

for iGr = 1:length(GroupList)
    Group = GroupList[iGr]
    animalList,_,LIIList,LIVList,LVList,LVIList,CondList = callGroup(Group); 
    anipar    = (;LIIList,LIVList,LVList,LVIList)
    if Group == "KIC"
        global KIC_WT = CWT_Loop(figs, animalList, CondList, CLList, params, anipar, takepic)
    elseif Group == "KIT"
        global KIT_WT = CWT_Loop(figs, animalList, CondList, CLList, params, anipar, takepic)
    elseif Group == "KIV"
        global KIV_WT = CWT_Loop(figs, animalList, CondList, CLList, params, anipar, takepic)
    else
        error("Group name does not match what is run through this script, please edit names or script")
    end
end


save(joinpath(spect,"KIC_WT.jld2"),"KIC_WT",KIC_WT)
save(joinpath(spect,"KIT_WT.jld2"),"KIT_WT",KIT_WT)
save(joinpath(spect,"KIV_WT.jld2"),"KIV_WT",KIV_WT)

# if you want to start from here:
KIC_WT = load(joinpath(spect,"KIC_WT.jld2"))["KIC_WT"]
KIT_WT = load(joinpath(spect,"KIT_WT.jld2"))["KIT_WT"]
KIV_WT = load(joinpath(spect,"KIV_WT.jld2"))["KIV_WT"]

MeasList = ["preCL_1" "CL_1"  "preAM_1" "AM_1"] 

# generate spectral plots (osci bands sorted by group and then frequency)

include(joinpath(func,"Permfunctions.jl"))
include(joinpath(func,"RunPermuations.jl"))

# this specifies the range of time accross which to permute 
cuttime = 195:295 # 1:1377 is full length, 200:300 is first 100 ms after tone onset
PermBetween(figs,spect,KIT_WT,KIC_WT,KIV_WT,MeasList,params,cuttime)

# within group
#PermWithin(figs,spect,KIT_WT,KIC_WT,KIV_WT,MeasList,params)