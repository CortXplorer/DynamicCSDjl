function PeakRatio_Between(data,StatTab,whichstim="2Hz",trialtype="TA")
    # Input: folder path data, table for ratio of 2 hz and 5 hz, all groups being tested
    # Output: table in folder Data/AvrecPeakStats which contains the pvalue result for an unequal test of variance (2 sample t test) of the ratio of final to first peak response between each group

    # Setup for simple 2 sample t test from HypothesisTests ->
    # 3 between group comparisons across 5 measurement types and 5 layers:
    CompList    = ["KIC vs KIT" "KIC vs KIV" "KIT vs KIV"]
    LayList     = unique(StatTab[!,:Layer])
    MeasList    = unique(StatTab[!,:Measurement])
    # containers for the table columns: 
    Comparison  = String[]
    Layer       = String[]
    Measurement = String[]   
    PAMP        = ones(length(CompList) * length(LayList) * length(MeasList)) # number of rows
    PRMS        = ones(length(PAMP))
    
    count       = [1]
    # loop through the above - create a table output! :) 
    for iComp = 1:length(CompList)
        # pull out groups from ratio table
        G1 = StatTab[StatTab[!,:Group] .== CompList[iComp][1:3] ,:]
        G2 = StatTab[StatTab[!,:Group] .== CompList[iComp][8:10],:]

        for iMeas = 1:length(MeasList)
            # pull out the measurment for those groups
            G1_meas = G1[G1[!,:Measurement] .== MeasList[iMeas],:]
            G2_meas = G2[G2[!,:Measurement] .== MeasList[iMeas],:]

            for iLay = 1:length(LayList)
                # further pull out each layer definition       
                G1_lay = G1_meas[G1_meas[!,:Layer] .== LayList[iLay],:]
                G2_lay = G2_meas[G2_meas[!,:Layer] .== LayList[iLay],:]

                # find the p value outcome for this comparison at both stimuli conditions
                PAMP[count[1]] = pvalue(UnequalVarianceTTest(G1_lay[!,:RatioAMP],G2_lay[!,:RatioAMP]))
                PRMS[count[1]] = pvalue(UnequalVarianceTTest(G1_lay[!,:RatioRMS],G2_lay[!,:RatioRMS]))

                count[1] = count[1] + 1
                # store the appropriate tags at the same positions in their respective lists
                push!(Comparison,CompList[iComp])
                push!(Measurement,MeasList[iMeas])
                push!(Layer,LayList[iLay])
            end # layer
        end # measurement type
    end # comparison of which groups

    BetweenGroup = DataFrame(Comparison=Comparison, Measurement=Measurement, Layer=Layer, PAMP=PAMP, PRMS=PRMS)

    foldername = "AvrecPeakStats"
    if !isdir(joinpath(data,foldername))
        mkdir(joinpath(data,foldername))
    end
    title = "Ratio_BetweenGroups_" * whichstim * "_" * trialtype
    name = joinpath(data,foldername,title) * ".csv"
    CSV.write(name, BetweenGroup)
end

function PeakRatio_Within(data,StatTab,whichstim="2Hz",trialtype="TA")
    # Input: folder path data, table for ratio of 2 hz and 5 hz, all groups being tested
    # Output: table in folder Data/AvrecPeakStats which contains the pvalue result for an equal test of variance (2 sample t test) of the ratio of final to first peak response between each measurement to the first

    # Setup for simple 2 sample t test from HypothesisTests ->
    # 4 within group comparisons per group across 5 measurement types and 5 layers:
    GroupList    = ["KIC" "KIT" "KIV"]
    LayList     = unique(StatTab[!,:Layer])
    MeasList    = unique(StatTab[!,:Measurement])
    # containers for the table columns: 
    Group       = String[]
    Layer       = String[]
    # Measurement = String[]   
    Prev1_AMP  = ones(length(GroupList) * length(LayList)) # number of rows (more columns this time)
    Prev2_AMP  = ones(length(Prev1_AMP))
    Prev3_AMP  = ones(length(Prev1_AMP))
    Prev4_AMP  = ones(length(Prev1_AMP))

    Prev1_RMS  = ones(length(Prev1_AMP))
    Prev2_RMS  = ones(length(Prev1_AMP))
    Prev3_RMS  = ones(length(Prev1_AMP))
    Prev4_RMS  = ones(length(Prev1_AMP))

    count       = [1]
    # loop through the above - create a table output! :) 
    for iGrp = 1:length(GroupList)
        # pull out groups from ratio table
        Gr_Stat = StatTab[StatTab[!,:Group] .== GroupList[iGrp],:]

        for iLay = 1:length(LayList)
            # pull out the layer one at a time
            Gr_lay = Gr_Stat[Gr_Stat[!,:Layer] .== LayList[iLay],:]

            # further pull out each measurement all together
            Gr_Pre = Gr_lay[Gr_lay[!,:Measurement] .== MeasList[1],:]
            Gr_1   = Gr_lay[Gr_lay[!,:Measurement] .== MeasList[2],:]
            Gr_2   = Gr_lay[Gr_lay[!,:Measurement] .== MeasList[3],:]
            Gr_3   = Gr_lay[Gr_lay[!,:Measurement] .== MeasList[4],:]
            Gr_4   = Gr_lay[Gr_lay[!,:Measurement] .== MeasList[5],:]

            # find the p value outcome for group between measurements
            Prev1_AMP[count[1]]  = pvalue(EqualVarianceTTest(Gr_Pre[!,:RatioAMP], Gr_1[!,:RatioAMP]))
            Prev2_AMP[count[1]]  = pvalue(EqualVarianceTTest(Gr_Pre[!,:RatioAMP], Gr_2[!,:RatioAMP]))
            Prev3_AMP[count[1]]  = pvalue(EqualVarianceTTest(Gr_Pre[!,:RatioAMP], Gr_3[!,:RatioAMP]))
            Prev4_AMP[count[1]]  = pvalue(EqualVarianceTTest(Gr_Pre[!,:RatioAMP], Gr_4[!,:RatioAMP]))

            Prev1_RMS[count[1]]  = pvalue(EqualVarianceTTest(Gr_Pre[!,:RatioRMS], Gr_1[!,:RatioRMS]))
            Prev2_RMS[count[1]]  = pvalue(EqualVarianceTTest(Gr_Pre[!,:RatioRMS], Gr_2[!,:RatioRMS]))
            Prev3_RMS[count[1]]  = pvalue(EqualVarianceTTest(Gr_Pre[!,:RatioRMS], Gr_3[!,:RatioRMS]))
            Prev4_RMS[count[1]]  = pvalue(EqualVarianceTTest(Gr_Pre[!,:RatioRMS], Gr_4[!,:RatioRMS]))

            
            count[1] = count[1] + 1
            # store the appropriate tags at the same positions in their respective lists
            push!(Group,GroupList[iGrp])
            push!(Layer,LayList[iLay])

        end # layer
    end # comparison of which groups

    BetweenGroup = DataFrame(Group=Group, Layer=Layer, Prev1_AMP=Prev1_AMP, Prev2_AMP=Prev2_AMP, Prev3_AMP=Prev3_AMP, Prev4_AMP=Prev4_AMP, Prev1_RMS=Prev1_RMS, Prev2_RMS=Prev2_RMS, Prev3_RMS=Prev3_RMS, Prev4_RMS=Prev4_RMS)

    foldername = "AvrecPeakStats"
    if !isdir(joinpath(data,foldername))
        mkdir(joinpath(data,foldername))
    end
    title = "Ratio_WithinGroups_" * whichstim * "_" * trialtype
    name = joinpath(data,foldername,title) * ".csv"
    CSV.write(name, BetweenGroup)
end


function Peak1_Between(data,StatTab,whichpeak="First",whichstim="2Hz",trialtype="TA")
    # Input: folder path data, table for peak amp/lat of 2 hz or 5 hz, all groups being tested
    # Output: table in folder Data/AvrecPeakStats which contains the pvalue result for an unequal test of variance (2 sample t test) of the peak amp and lat of the first response between each group

    # Setup for simple 2 sample t test from HypothesisTests ->
    # 3 between group comparisons across 5 measurement types and 5 layers:
    CompList    = ["KIC vs KIT" "KIC vs KIV" "KIT vs KIV"]
    LayList     = unique(StatTab[!,:Layer])
    MeasList    = unique(StatTab[!,:Measurement])
    # containers for the table columns: 
    Comparison  = String[]
    Layer       = String[]
    Measurement = String[]   
    PAmp        = ones(length(CompList) * length(LayList) * length(MeasList)) # number of rows
    PLat    = ones(length(PAmp))
    PRMS    = ones(length(PAmp))

    count       = [1]
    # loop through the above - create a table output! :) 
    for iComp = 1:length(CompList)
        # pull out groups from ratio table
        G1 = StatTab[StatTab[!,:Group] .== CompList[iComp][1:3] ,:]
        G2 = StatTab[StatTab[!,:Group] .== CompList[iComp][8:10],:]

        for iMeas = 1:length(MeasList)
            # pull out the measurment for those groups
            G1_meas = G1[G1[!,:Measurement] .== MeasList[iMeas],:]
            G2_meas = G2[G2[!,:Measurement] .== MeasList[iMeas],:]

            for iLay = 1:length(LayList)
                # further pull out each layer definition
                G1_lay = G1_meas[G1_meas[!,:Layer] .== LayList[iLay],:]
                G2_lay = G2_meas[G2_meas[!,:Layer] .== LayList[iLay],:]

                # find the p value outcome for this comparison at both stimuli conditions
                PAmp[count[1]] = pvalue(UnequalVarianceTTest(G1_lay[!,:PeakAmp],G2_lay[!,:PeakAmp]))
                PLat[count[1]] = pvalue(UnequalVarianceTTest(G1_lay[!,:PeakLat],G2_lay[!,:PeakLat]))
                PRMS[count[1]] = pvalue(UnequalVarianceTTest(G1_lay[!,:RMS],G2_lay[!,:RMS]))
                count[1] = count[1] + 1
                # store the appropriate tags at the same positions in their respective lists
                push!(Comparison,CompList[iComp])
                push!(Measurement,MeasList[iMeas])
                push!(Layer,LayList[iLay])
            end # layer
        end # measurement type
    end # comparison of which groups

    BetweenGroup = DataFrame(Comparison=Comparison, Measurement=Measurement, Layer=Layer, PAmp=PAmp, PLat=PLat, PRMS=PRMS)

    foldername = "AvrecPeakStats"
    if !isdir(joinpath(data,foldername))
        mkdir(joinpath(data,foldername))
    end
    title = whichpeak * "_BetweenGroups_" * whichstim * "_" * trialtype
    name = joinpath(data,foldername,title) * ".csv"
    CSV.write(name, BetweenGroup)
end

function Peak1_Within(data,StatTab,whichpeak="First",whichstim="2Hz",trialtype="TA")
    # Input: folder path data, table for peak amp/lat of 2 hz and 5 hz, all groups being tested
    # Output: table in folder Data/AvrecPeakStats which contains the pvalue result for an equal test of variance (2 sample t test) of the peak amp and lat of the first response between each measurement to the first

    # Setup for simple 2 sample t test from HypothesisTests ->
    # 4 within group comparisons per group across 5 measurement types and 5 layers:
    GroupList    = ["KIC" "KIT" "KIV"]
    LayList     = unique(StatTab[!,:Layer])
    MeasList    = unique(StatTab[!,:Measurement])
    # containers for the table columns: 
    Group       = String[]
    Layer       = String[]
    # Measurement = String[]   
    Prev1_Amp  = ones(length(GroupList) * length(LayList)) # number of rows (more columns this time)
    Prev2_Amp  = ones(length(Prev1_Amp))
    Prev3_Amp  = ones(length(Prev1_Amp))
    Prev4_Amp  = ones(length(Prev1_Amp))

    Prev1_Lat  = ones(length(Prev1_Amp))
    Prev2_Lat  = ones(length(Prev1_Amp))
    Prev3_Lat  = ones(length(Prev1_Amp))
    Prev4_Lat  = ones(length(Prev1_Amp))

    Prev1_RMS  = ones(length(Prev1_Amp))
    Prev2_RMS  = ones(length(Prev1_Amp))
    Prev3_RMS  = ones(length(Prev1_Amp))
    Prev4_RMS  = ones(length(Prev1_Amp))

    count       = [1]
    # loop through the above - create a table output! :) 
    for iGrp = 1:length(GroupList)
        # pull out groups from ratio table
        Gr_Stat = StatTab[StatTab[!,:Group] .== GroupList[iGrp],:]

        for iLay = 1:length(LayList)
            # pull out the layer one at a time
            Gr_lay = Gr_Stat[Gr_Stat[!,:Layer] .== LayList[iLay],:]

            # further pull out each measurement all together
            Gr_Pre = Gr_lay[Gr_lay[!,:Measurement] .== MeasList[1],:]
            Gr_1   = Gr_lay[Gr_lay[!,:Measurement] .== MeasList[2],:]
            Gr_2   = Gr_lay[Gr_lay[!,:Measurement] .== MeasList[3],:]
            Gr_3   = Gr_lay[Gr_lay[!,:Measurement] .== MeasList[4],:]
            Gr_4   = Gr_lay[Gr_lay[!,:Measurement] .== MeasList[5],:]

            # find the p value outcome for group between measurements
            Prev1_Amp[count[1]]  = pvalue(EqualVarianceTTest(Gr_Pre[!,:PeakAmp], Gr_1[!,:PeakAmp]))
            Prev2_Amp[count[1]]  = pvalue(EqualVarianceTTest(Gr_Pre[!,:PeakAmp], Gr_2[!,:PeakAmp]))
            Prev3_Amp[count[1]]  = pvalue(EqualVarianceTTest(Gr_Pre[!,:PeakAmp], Gr_3[!,:PeakAmp]))
            Prev4_Amp[count[1]]  = pvalue(EqualVarianceTTest(Gr_Pre[!,:PeakAmp], Gr_4[!,:PeakAmp]))

            Prev1_Lat[count[1]]  = pvalue(EqualVarianceTTest(Gr_Pre[!,:PeakLat], Gr_1[!,:PeakLat]))
            Prev2_Lat[count[1]]  = pvalue(EqualVarianceTTest(Gr_Pre[!,:PeakLat], Gr_2[!,:PeakLat]))
            Prev3_Lat[count[1]]  = pvalue(EqualVarianceTTest(Gr_Pre[!,:PeakLat], Gr_3[!,:PeakLat]))
            Prev4_Lat[count[1]]  = pvalue(EqualVarianceTTest(Gr_Pre[!,:PeakLat], Gr_4[!,:PeakLat]))

            Prev1_RMS[count[1]]  = pvalue(EqualVarianceTTest(Gr_Pre[!,:RMS], Gr_1[!,:RMS]))
            Prev2_RMS[count[1]]  = pvalue(EqualVarianceTTest(Gr_Pre[!,:RMS], Gr_2[!,:RMS]))
            Prev3_RMS[count[1]]  = pvalue(EqualVarianceTTest(Gr_Pre[!,:RMS], Gr_3[!,:RMS]))
            Prev4_RMS[count[1]]  = pvalue(EqualVarianceTTest(Gr_Pre[!,:RMS], Gr_4[!,:RMS]))
            
            count[1] = count[1] + 1
            # store the appropriate tags at the same positions in their respective lists
            push!(Group,GroupList[iGrp])
            push!(Layer,LayList[iLay])

        end # layer
    end # comparison of which groups

    BetweenGroup = DataFrame(Group=Group, Layer=Layer, Prev1_Amp=Prev1_Amp, Prev2_Amp=Prev2_Amp, Prev3_Amp=Prev3_Amp, Prev4_Amp=Prev4_Amp,  Prev1_Lat=Prev1_Lat, Prev2_Lat=Prev2_Lat, Prev3_Lat=Prev3_Lat, Prev4_Lat=Prev4_Lat, Prev1_RMS=Prev1_RMS, Prev2_RMS=Prev2_RMS, Prev3_RMS=Prev3_RMS, Prev4_RMS=Prev4_RMS)

    foldername = "AvrecPeakStats"
    if !isdir(joinpath(data,foldername))
        mkdir(joinpath(data,foldername))
    end
    title = whichpeak * "_WithinGroups_" * whichstim * "_" * trialtype
    name = joinpath(data,foldername,title) * ".csv"
    CSV.write(name, BetweenGroup)
end