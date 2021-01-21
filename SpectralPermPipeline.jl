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
MeasList = ["preCL_1" "CL_1"  "preAM_1" "AM_1"] 

# Generate group WT jld2s
GenGroupWT(spect,figs,params,GroupList,CLList,takepic)

# if you want to start from here with pre-generated group WTs:
KIC_WT = load(joinpath(spect,"KIC_WT.jld2"))["KIC_WT"]
KIT_WT = load(joinpath(spect,"KIT_WT.jld2"))["KIT_WT"]
KIV_WT = load(joinpath(spect,"KIV_WT.jld2"))["KIV_WT"]

include(joinpath(func,"Permfunctions.jl"))
include(joinpath(func,"RunPermuations.jl"))

# this specifies the range of time accross which to permute 
cuttime = 195:295 # 1:1377 is full length, 200:300 is first 100 ms after tone onset
PermBetween(figs,spect,KIT_WT,KIC_WT,KIV_WT,MeasList,params,cuttime)

# within group
#PermWithin(figs,spect,KIT_WT,KIC_WT,KIV_WT,MeasList,params)

# generate spectral plots (osci bands sorted by group and then frequency)