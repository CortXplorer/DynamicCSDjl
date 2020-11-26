using MAT, Statistics, DSP, Colors
using JLD2, FileIO
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
include(joinpath(func,"Permfunctions.jl"))
include(joinpath(func,"PowerPermute.jl"))
include(joinpath(func,"PhasePermute.jl"))

# Parameter
params = (
    sampleRate=1000, startTime=-200, timeLimits=[-200 1377], frequencyLimits=[8 100], timeBandWidth=40, 
    stimList=["twoHz","fiveHz","tenHz","twentyHz","fortyHz"], 
    layers=["I_II","IV","V","VI"], 
    osciBands=["alpha","beta low","beta high","gamma low","gamma high"], bandRanges=[[1:8...],[9:15...],[16:23...],[24:35...],[36:40...]])
# real ranges captured by these data: alpha=(8:12), beta_low=(13:18), beta_high=(19:30), gamma_low=(31:60), gamma_high=(61:100)
## GROUP determination
GroupList = ["KIC" "KIT" "KIV"] 
## Conditions to run
CLList  = ["preCL" "CL" "preAM" "AM"] 
## Conditional picture; takepic == 1 if you do want figure output
takepic = 1

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

MeasList = ["preCL_1" "CL_1" "CL_2" "CL_3" "CL_4" "preAM_1" "AM_1" "AM_2" "AM_3"]

# generate spectral plots (osci bands sorted by group and then frequency)

println("Power Permutation Between Groups")
# ADD: Cohen's D output column
PermOut = []
for iMeas = 1:length(MeasList)
    curMeas = MeasList[iMeas]
    for iSti = 1:length(params.stimList)
        curStim = params.stimList[iSti]
        for iLay = 1:length(params.layers)
            curLay = params.layers[iLay]
            # FIX: ERROR: KeyError: key "AM_4" not found (have it check for this key and skip if it's not there)
            println("Power Measurement $curMeas, Stimulus $curStim, Layer $iLay")
            Permout1 = PowerPermBetween(figs,KIT_WT,KIC_WT,"KIT","KIC",curMeas,curStim,curLay,takepic)
            Permout2 = PowerPermBetween(figs,KIT_WT,KIV_WT,"KIT","KIV",curMeas,curStim,curLay,takepic)
            Permout3 = PowerPermBetween(figs,KIC_WT,KIV_WT,"KIC","KIV",curMeas,curStim,curLay,takepic)

            if !isempty(PermOut)
                PermOut = vcat(PermOut,Permout1,Permout2,Permout3)
            else
                PermOut = vcat(Permout1,Permout2,Permout3)
            end
        end # layer
    end # Stim frequency
end # Measurement

save(joinpath(spect,"PowerBetween.jld2"),"PowerBetween",PermOut)
# PowerBetween = load(joinpath(spect,"PowerBetween.jld2"))["PowerBetween"] #yaaaasss

println("Phase Coherence Permutation Between Groups")
PermOut = []
for iMeas = 1:length(MeasList)
    curMeas = MeasList[iMeas]
    for iSti = 1:length(params.stimList)
        curStim = params.stimList[iSti]
        for iLay = 1:length(params.layers)
            curLay = params.layers[iLay]
            # FIX: ERROR: KeyError: key "AM_4" not found (have it check for this key and skip if it's not there)
            println("Phase Measurement $curMeas, Stimulus $curStim, Layer $iLay")
            Permout1 = PhasePermBetween(figs,KIT_WT,KIC_WT,"KIT","KIC",curMeas,curStim,curLay,takepic)
            Permout2 = PhasePermBetween(figs,KIT_WT,KIV_WT,"KIT","KIV",curMeas,curStim,curLay,takepic)
            Permout3 = PhasePermBetween(figs,KIC_WT,KIV_WT,"KIC","KIV",curMeas,curStim,curLay,takepic)

            if !isempty(PermOut)
                PermOut = vcat(PermOut,Permout1,Permout2,Permout3)
            else
                PermOut = vcat(Permout1,Permout2,Permout3)
            end
        end # layer
    end # Stim frequency
end # Measurement

save(joinpath(spect,"PhaseBetween.jld2"),"PhaseBetween",PermOut)
# PhaseBetween = load(joinpath(spect,"PhaseBetween.jld2"))["PhaseBetween"] 

# within group

# within group