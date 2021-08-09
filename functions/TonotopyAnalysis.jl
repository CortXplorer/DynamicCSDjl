function TonotopyAnalysis(data) 
        
    TonoRMS = DataFrame(Group=String[], Animal=String[], Measurement=String[], Trial=Int64[], AVRECrms=Float64[], AVRECamp=Float64[], AVREClat=Float64[], LIIrms=Float64[], LIIamp=Float64[], LIIlat=Float64[], LIVrms=Float64[], LIVamp=Float64[], LIVlat=Float64[], LVrms=Float64[], LVamp=Float64[], LVlat=Float64[], LVIrms=Float64[], LVIamp=Float64[], LVIlat=Float64[])

    datafiles = readdir("Data")                         # get files in data folder
    datafiles = datafiles[occursin.(".mat",datafiles)]  # filter for mat files 

    println("Start")

    for curAn = 1:length(datafiles)                      # loop through animal mat files
        if length(datafiles[curAn]) !== 14               # quick check for correct file
            error("there is a mat file in your directory that is not standard animal data")
        end
        Data = matread(joinpath(data,datafiles[curAn]))  # load cur animal data
        Data = Data["Data"]                              # shortcut first dict heading

        aniName = datafiles[curAn][1:5]  
        animalList,_,LIIList,LIVList,LVList,LVIList,CondList = callGroup(aniName)
        if sum(isequal.(aniName, animalList)) == 0       # quick check if animal is in list
            continue
        end
        Anindx = findfirst(isequal(aniName), animalList)[2] # animal index in lists

        println("Running " * aniName)

        # all tonotopy measurements: note that first four tonotopies are missing currently
        measurements = ["Pre" "preAMtono" "preCLtono" "CLtono" "AMtono"]

        for iMeasty = 1:length(measurements)              # loop through measurement types
            curMeaslist = CondList[measurements[iMeasty]][Anindx] 
            if ismissing(curMeaslist[1])
                continue
            end

            for iMeas = 1:length(curMeaslist)             # loop through measurements
                curMeas = curMeaslist[iMeas]

                # measurement index in lists
                DataMeas = Data["measurement"]
                DataMeas[findall(isempty.(DataMeas))] .= "no" # (string placeholder for missing)
                Meindx = findfirst(contains(curMeas), DataMeas) 
                #get animal BF (granular)
                curBF = findfirst(isequal(Data["BF_IV"][Meindx]), Data["Frqz"][Meindx]) 
                if isnothing(curBF);                        # if no BF (weird), take 4k
                    curBF = CartesianIndex(1, 3); 
                    println("artificial BF of 4000 selected")
                end 
                curCSD = Data["SglTrl_CSD"][Meindx][curBF]

                for iTri = 1:size(curCSD,3)                 # loop through trials weeeee
                    curCSDst = curCSD[:,:,iTri]

                    AVRECtrace = mean(abs.(curCSDst),dims=1)'[200:400]
                    AVRECrms   = rms(AVRECtrace)
                    AVREClat,AVRECamp = ifindpeaks(AVRECtrace,0)

                    LII = curCSDst[LIIList[Anindx],:]       # pull out layer II
                    LII[findall(LII .> 0)] .= NaN           # NaN all source 
                    LIItrace = nanmean(LII)[200:400]        # custom nanmean function
                    LIIrms = rms(filter(!isnan, LIItrace))  # take rms of first 200 ms
                    LIIlat, LIIamp = ifindpeaks(LIItrace,1) # take peak amp and lat

                    LIV = curCSDst[LIVList[Anindx],:]       # pull out layer IV
                    LIV[findall(LIV .> 0)] .= NaN           # NaN all source 
                    LIVtrace = nanmean(LIV)[200:400]        # custom nanmean function
                    LIVrms = rms(filter(!isnan, LIVtrace))  # take rms of first 200 ms
                    LIVlat, LIVamp = ifindpeaks(LIVtrace,1) # take peak amp and lat

                    LV  = curCSDst[LVList[Anindx],:]        # pull out layer IV
                    LV[findall(LIV .> 0)]  .= NaN           # NaN all source 
                    LVtrace = nanmean(LV)[200:400]          # custom nanmean function
                    LVrms = rms(filter(!isnan, LVtrace))    # take rms of first 200 ms
                    LVlat, LVamp = ifindpeaks(LVtrace,1)    # take peak amp and lat

                    LVI = curCSDst[LVIList[Anindx],:]       # pull out layer IV
                    LVI[findall(LVI .> 0)] .= NaN           # NaN all source 
                    LVItrace = nanmean(LVI)[200:400]        # custom nanmean function
                    LVIrms = rms(filter(!isnan, LVItrace))  # take rms of first 200 ms
                    LVIlat, LVIamp = ifindpeaks(LVItrace,1) # take peak amp and lat

                    push!(TonoRMS, [aniName[1:3] aniName measurements[iMeasty]*"_"*string(iMeas) iTri AVRECrms AVRECamp AVREClat LIIrms LIIamp LIIlat LIVrms LIVamp LIVlat LVrms LVamp LVlat LVIrms LVIamp LVIlat])
                end
            end
        end
    end

    if !isdir(joinpath(data,"Tonotopy"))
        mkdir(joinpath(data,"Tonotopy"))
    end

    name = joinpath(data,"Tonotopy","TonoRMS") * ".csv"
    CSV.write(name, TonoRMS)

end