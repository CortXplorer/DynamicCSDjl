# compute scalograms 
using MAT

home    = @__DIR__
data    = joinpath(home,"Data")
func    = joinpath(home,"functions")
group   = joinpath(home,"groups")

include(joinpath(group,"callGroup.jl"))

# Initial Parameters
params = (sampleRate=1000, startTime=-0.2, timeLimits=[-0.2 0.399], frequencyLimits=[5 500], voicesPerOctave=8, timeBandWidth=54, layers=["I_II","IV","V","VI"])

foldername = "Spectral"
if !isdir(joinpath(data,foldername))
    mkdir(joinpath(data,foldername))
end

# function runCwtCSD(GroupFile, params, home)
GroupList = ["KIC"]
Group     = GroupList[1]
animalList,_,LIIList,LIVList,LVList,LVIList,CondList = callGroup(Group); # in groups folder
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
iMeas = 1
# loop through measurements in the list iMeas
curMeas = CondList[CLList[iCond]][iAn][iMeas]
runthis = [curAn*"_"*curMeas] # generate full measurement name

thisind = findall(anDat["Data"]["measurement"] .== runthis)[1][2] #[1][2] for extracting cartesian index
iSti = 1
# loop through stim frequencies
curCSD  = anDat["Data"]["CSD"][thisind][iSti]

# end function

# generate spectral plots (osci bands sorted by group and then frequency)

# spectral power permutations
# between group
# within group

# spectral phase coherence permutations
# between group
# within group