# compute scalograms 
using MAT, Statistics, DSP, Colors

home    = @__DIR__
data    = joinpath(home,"Data")
func    = joinpath(home,"functions")
group   = joinpath(home,"groups")
figs    = joinpath(home,"figs")

include(joinpath(group,"callGroup.jl"))
include(joinpath(func,"CWTfunc.jl"))

# Parameters
params = (
    sampleRate=1000, startTime=-200, timeLimits=[-200 1377], frequencyLimits=[8 100], timeBandWidth=40, 
    stimList=["2Hz","5Hz","10Hz","20Hz","40Hz"], 
    layers=["I_II","IV","V","VI"], 
    osciBands=["alpha","beta low","beta high","gamma low","gamma high"], bandRanges=[[1:8...],[9:15...],[16:23...],[24:35...],[36:40...]])
# real ranges captured by these data: alpha=(8:12), beta_low=(13:18), beta_high=(19:30), gamma_low=(31:60), gamma_high=(61:100)
## GROUP determination
GroupList = ["KIC" "KIT"]
## Conditions to run
CLList  = ["preCL" "CL"]
## Conditional picture; takepic == 1 if you do want figure output
takepic = 0

foldername = "Spectral"
if !isdir(joinpath(data,foldername))
    mkdir(joinpath(data,foldername))
end

# Loop through Groups to use function CWT_Loop which splits by animals, condition, measurement, stim frequency, and layer and outputs a dictionary with a power and phase coherence dataset for each chunk
KIC_WT = []
KIT_WT = []
KIV_WT = []

for iGr = 1:length(GroupList)
    Group = GroupList[iGr]
    animalList,_,LIIList,LIVList,LVList,LVIList,CondList = callGroup(Group); 
    anipar    = (;LIIList,LIVList,LVList,LVIList)
    if Group == "KIC"
        KIC_WT = CWT_Loop(figs, animalList, CondList, CLList, params, anipar, takepic)
    elseif Group == "KIT"
        KIT_WT = CWT_Loop(figs, animalList, CondList, CLList, params, anipar, takepic)
    elseif Group == "KIV"
        KIV_WT = CWT_Loop(figs, animalList, CondList, CLList, params, anipar, takepic)
    else
        error("Group name does not match what is run through this script, please edit names or script")
    end
end

# generate spectral plots (osci bands sorted by group and then frequency)

# spectral power permutations
# between group
# within group

# spectral phase coherence permutations
# between group
# within group