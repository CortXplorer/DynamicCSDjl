# compute scalograms 
using MAT, Statistics, DSP, Colors

home    = @__DIR__
data    = joinpath(home,"Data")
func    = joinpath(home,"functions")
group   = joinpath(home,"groups")
figs    = joinpath(home,"figs")

include(joinpath(group,"callGroup.jl"))

# Initial Parameters
params = (sampleRate=1000, startTime=-200, timeLimits=[-200 1377], frequencyLimits=[5 100], voicesPerOctave=8, timeBandWidth=40, layers=["I_II","IV","V","VI"])

foldername = "Spectral"
if !isdir(joinpath(data,foldername))
    mkdir(joinpath(data,foldername))
end

GroupList = ["KIC"]
Group     = GroupList[1]

# function runCwtCSD(home, Group, params)
animalList,_,LIIList,LIVList,LVList,LVIList,CondList = callGroup(Group); 
anipar    = (;LIIList,LIVList,LVList,LVIList)
stimList  = ["2Hz" "5Hz" "10Hz" "20Hz" "40Hz"]

iAn = 1

# loop through animals iAn
curAn   = animalList[iAn]
anDat   = matread(joinpath(data,(curAn * "_Data.mat")))
anMeas  = anDat["Data"]["measurement"]

CLList  = ["preCL" "CL"]
iCond = 1
# loop through condition list iCL
curCond = CLList[iCond]
MeasList = CondList[CLList[iCond]][iAn]

iMeas = 1
# loop through measurements in the list iMeas
curMeas = CondList[CLList[iCond]][iAn][iMeas]
curCond = curCond * "_" * string(iMeas)
runthis = [curAn*"_"*curMeas] # generate full measurement name

thisind = findall(anDat["Data"]["measurement"] .== runthis)[1][2] #[1][2] for extracting cartesian index
iSti = 1
# loop through stim frequencies
curCSD  = anDat["Data"]["SglTrl_CSD"][thisind][iSti]
iLay = 2
# loop through layers
if iLay == 2 # change call for each layer
curChan = anipar.LIVList[iAn]
end

if length(curChan) > 3 # take center 3 if greater than 3
    centerChan = Int(ceil(length(curChan)/2))
    curChan = curChan[centerChan-1:centerChan+1]
end

meanLayCSD = mean(curCSD[curChan,:,:],dims=3)
ROI = mean(meanLayCSD,dims=1)





# end function

# generate spectral plots (osci bands sorted by group and then frequency)

# spectral power permutations
# between group
# within group

# spectral phase coherence permutations
# between group
# within group