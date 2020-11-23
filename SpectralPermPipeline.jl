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
params = (sampleRate=1000, startTime=-200, timeLimits=[-200 1377], frequencyLimits=[8 100], timeBandWidth=40, stimList=["2Hz","5Hz","10Hz","20Hz","40Hz"], layers=["I_II","IV","V","VI"], osciBands=["alpha","beta low","beta high","gamma low","gamma high"], bandRanges=[[1:8...],[9:15...],[16:23...],[24:35...],[36:40...]])
# real ranges captured by these data: alpha=(8:12), beta_low=(13:18), beta_high=(19:30), gamma_low=(31:60), gamma_high=(61:100)

foldername = "Spectral"
if !isdir(joinpath(data,foldername))
    mkdir(joinpath(data,foldername))
end

GroupList = ["KIC"]
CLList  = ["preCL" "CL"]
Group     = GroupList[1]

# function runCwtCSD(figs,home, Group, params)
animalList,_,LIIList,LIVList,LVList,LVIList,CondList = callGroup(Group); 
anipar    = (;LIIList,LIVList,LVList,LVIList)

# Loop through animals in Group
Animal = Dict()
for iAn = 1:length(animalList)

    # loop through animals iAn
    curAn   = animalList[iAn]
    anDat   = matread(joinpath(data,(curAn * "_Data.mat")));
    anMeas  = anDat["Data"]["measurement"]

    # Loop through the condition list
    Condition = Dict()
    for iCond = 1:length(CLList)
        # loop through condition list iCL
        MeasList = CondList[CLList[iCond]][iAn]

        # Loop through measurements
        Measurement = Dict()
        for iMeas = 1:length(MeasList)

            curMeas = CondList[CLList[iCond]][iAn][iMeas]
            curCond = CLList[iCond] * "_" * string(iMeas)
            runthis = [curAn*"_"*curMeas] # generate full measurement name

            thisind = findall(anDat["Data"]["measurement"] .== runthis)[1][2] #[1][2] for extracting cartesian index

            # Loop through stimulation frequencies
            StimFreq = Dict()
            for iSti = 1:length(params.stimList)

                curCSD  = anDat["Data"]["SglTrl_CSD"][thisind][iSti];
                curStim = params.stimList[iSti]

                # Loop through layers
                Layers = Dict()
                for iLay = 1:length(params.layers)

                    curLay = params.layers[iLay]
                    if curLay == "I_II"
                        curChan = anipar.LIIList[iAn]
                    elseif curLay == "IV"
                        curChan = anipar.LIVList[iAn]
                    elseif curLay == "V"
                        curChan = anipar.LVList[iAn]
                    elseif curLay == "VI"
                        curChan = anipar.LVIList[iAn]
                    end

                    if length(curChan) > 3 # take center 3 if greater than 3
                        centerChan = Int(ceil(length(curChan)/2))
                        curChan = curChan[centerChan-1:centerChan+1]
                    end

                    # Here is the region of interest accross trials
                    ROI = mean(curCSD[curChan,:,:],dims=1)

                    Datpower, Datphsco = CWTanalysis(figs,ROI,params,curAn,curCond,curLay,curStim)
                    Layers[curLay] = Datpower, Datphsco;
                end # Layers
                StimFreq[curStim] = Layers
            end # Stim Frequencies
            Measurement[curCond] = StimFreq
        end # Condition/Measurement
        Condition[CLList[iCond]] = Measurement
    end # Condition List
    Animal[curAn] = Condition
end # Animal
# pull out probably dictionary of all animal data for use in group plots below
# return SPower, PhaseCo
# end function

# generate spectral plots (osci bands sorted by group and then frequency)

# spectral power permutations
# between group
# within group

# spectral phase coherence permutations
# between group
# within group