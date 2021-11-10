function csdStruct(raw,figs,datap,Animal,AnimalName,Group,CondList,CondName,channels,LII,LIV,LV,LVI,iAn)


# loop through each type of measurement condition (specified by CondName)
Condition = Dict()
for iCond = 1:length(CondName)
    # loop through each measurement within that condition for that animal
    Measurement = OrderedDict()
    for iMeas = 1:length(CondList[CondName[iCond]][iAn])

        if ismissing(CondList[CondName[iCond]][iAn][iMeas])
            println("Skipping missing measurement")
            Measurement[CondList[CondName[iCond]][iAn][iMeas]] = missing
        else 
            measurement = AnimalName * "_" * CondList[CondName[iCond]][iAn][iMeas] * ".mat"
            println("Analyzing measurement: " * measurement[1:end-4])

            csdData,snkData = Dynamic_CSD(measurement,channels,LII,LIV,
                LV,LVI,raw,figs,Group)

            Measurement[CondList[CondName[iCond]][iAn][iMeas]] = csdData,snkData;
        end
    end
    Condition[CondName[iCond]] = Measurement
end
Animal[AnimalName] = Condition
varname  = AnimalName * "_Data"
filename = joinpath(datap,AnimalName) * "_Data.jld2"
save(filename, varname, Animal[AnimalName])

end