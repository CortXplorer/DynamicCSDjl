function csdStruct(raw,figs,datap,Animal,AnimalName,Group,CondList,CondName,channels,LII,LIV,LV,LVI,iAn)

# loop through each type of measurement condition (specified by CondName)
Condition = Dict()
for iCond = 1:length(CondName)
    # loop through each measurement within that condition for that animal
    Measurement = Dict()
    for iMeas = 1:length(CondList[CondName[iCond]][iAn])

        measurement = AnimalName * "_" * CondList[CondName[iCond]][iAn][iMeas] * ".mat"
        println("Analyzing measurement: " * measurement[1:end-4])

        csdData,snkData = Dynamic_CSD(measurement,channels,LII,LIV,
            LV,LVI,raw,figs,GroupList[iGr]);

        Measurement[CondList[CondName[iCond]][iAn][iMeas]] = csdData,snkData;
    end
    Condition[CondName[iCond]] = Measurement
end
Animal[AnimalName] = Condition
name = joinpath(datap,AnimalName) * "_Data.jld"
save(name, Animal[AnimalName])

end