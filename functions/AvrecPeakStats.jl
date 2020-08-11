function PeakRatio_Between(data,Stat2,Stat5)
    # Input: folder path data, table for ratio of 2 hz and 5 hz, all groups being tested
    # Output: table in folder Data/AvrecPeakStats which contains the pvalue result for an unequal test of variance (2 sample t test) of the ratio of final to first peak response between each group

    # Setup for simple 2 sample t test from HypothesisTests ->
    # 3 between group comparisons across 5 measurement types and 5 layers:
    CompList    = ["KIC vs KIT" "KIC vs KIV" "KIT vs KIV"]
    LayList     = unique(Stat2[!,:Layer])
    MeasList    = unique(Stat2[!,:Measurement])
    # containers for the table columns: 
    Comparison  = String[]
    Layer       = String[]
    Measurement = String[]   
    P2Hz        = ones(length(CompList) * length(LayList) * length(MeasList)) # number of rows
    P5Hz        = ones(length(CompList) * length(LayList) * length(MeasList))
    count       = [1]
    # loop through the above - create a table output! :) 
    for iComp = 1:length(CompList)
        # pull out groups from ratio table
        G1_2Hz = Stat2[Stat2[!,:Group] .== CompList[iComp][1:3] ,:]
        G2_2Hz = Stat2[Stat2[!,:Group] .== CompList[iComp][8:10],:]
        G1_5Hz = Stat5[Stat5[!,:Group] .== CompList[iComp][1:3] ,:]
        G2_5Hz = Stat5[Stat5[!,:Group] .== CompList[iComp][8:10],:]

        for iMeas = 1:length(MeasList)
            # pull out the measurment for those groups
            G1_2Hz_meas = G1_2Hz[G1_2Hz[!,:Measurement] .== MeasList[iMeas],:]
            G2_2Hz_meas = G2_2Hz[G2_2Hz[!,:Measurement] .== MeasList[iMeas],:]
            G1_5Hz_meas = G1_5Hz[G1_5Hz[!,:Measurement] .== MeasList[iMeas],:]
            G2_5Hz_meas = G2_5Hz[G2_5Hz[!,:Measurement] .== MeasList[iMeas],:]

            for iLay = 1:length(LayList)
                # further pull out each layer definition
                G1_2Hz_lay = G1_2Hz_meas[G1_2Hz_meas[!,:Layer] .== LayList[iLay],:]
                G2_2Hz_lay = G2_2Hz_meas[G2_2Hz_meas[!,:Layer] .== LayList[iLay],:]
                G1_5Hz_lay = G1_5Hz_meas[G1_5Hz_meas[!,:Layer] .== LayList[iLay],:]
                G2_5Hz_lay = G2_5Hz_meas[G2_5Hz_meas[!,:Layer] .== LayList[iLay],:]

                # find the p value outcome for this comparison at both stimuli conditions
                P2Hz[count[1]] = pvalue(UnequalVarianceTTest(G1_2Hz_lay[!,:Ratio],G2_2Hz_lay[!,:Ratio]))
                P5Hz[count[1]] = pvalue(UnequalVarianceTTest(G1_5Hz_lay[!,:Ratio],G2_5Hz_lay[!,:Ratio]))
                count[1] = count[1] + 1
                # store the appropriate tags at the same positions in their respective lists
                push!(Comparison,CompList[iComp])
                push!(Measurement,MeasList[iMeas])
                push!(Layer,LayList[iLay])
            end # layer
        end # measurement type
    end # comparison of which groups

    BetweenGroup = DataFrame(Comparison=Comparison, Measurement=Measurement, Layer=Layer, P2Hz=P2Hz, P5Hz=P5Hz)

    foldername = "AvrecPeakStats"
    if !isdir(joinpath(data,foldername))
        mkdir(joinpath(data,foldername))
    end
    title = "Ratio_BetweenGroups"
    name = joinpath(data,foldername,title) * ".csv"
    CSV.write(name, BetweenGroup)
end

function PeakRatio_Within(data,Stat2,Stat5)
    # Input: folder path data, table for ratio of 2 hz and 5 hz, all groups being tested
    # Output: table in folder Data/AvrecPeakStats which contains the pvalue result for an equal test of variance (2 sample t test) of the ratio of final to first peak response between each measurement to the first

    # Setup for simple 2 sample t test from HypothesisTests ->
    # 4 within group comparisons per group across 5 measurement types and 5 layers:
    GroupList    = ["KIC" "KIT" "KIV"]
    LayList     = unique(Stat2[!,:Layer])
    MeasList    = unique(Stat2[!,:Measurement])
    # containers for the table columns: 
    Group       = String[]
    Layer       = String[]
    # Measurement = String[]   
    P2Hz_Prev1  = ones(length(GroupList) * length(LayList)) # number of rows (more columns this time)
    P2Hz_Prev2  = ones(length(P2Hz_Prev1))
    P2Hz_Prev3  = ones(length(P2Hz_Prev1))
    P2Hz_Prev4  = ones(length(P2Hz_Prev1))
    P5Hz_Prev1  = ones(length(P2Hz_Prev1))
    P5Hz_Prev2  = ones(length(P2Hz_Prev1))
    P5Hz_Prev3  = ones(length(P2Hz_Prev1))
    P5Hz_Prev4  = ones(length(P2Hz_Prev1))

    count       = [1]
    # loop through the above - create a table output! :) 
    for iGrp = 1:length(GroupList)
        # pull out groups from ratio table
        G_2Hz = Stat2[Stat2[!,:Group] .== GroupList[iGrp],:]
        G_5Hz = Stat5[Stat5[!,:Group] .== GroupList[iGrp],:]

        for iLay = 1:length(LayList)
            # pull out the layer one at a time
            G_2Hz_lay = G_2Hz[G_2Hz[!,:Layer] .== LayList[iLay],:]
            G_5Hz_lay = G_5Hz[G_2Hz[!,:Layer] .== LayList[iLay],:]

            # further pull out each measurement all together
            G1_2Hz_Pre = G_2Hz_lay[G_2Hz_lay[!,:Measurement] .== MeasList[1],:]
            G1_2Hz_1   = G_2Hz_lay[G_2Hz_lay[!,:Measurement] .== MeasList[2],:]
            G1_2Hz_2   = G_2Hz_lay[G_2Hz_lay[!,:Measurement] .== MeasList[3],:]
            G1_2Hz_3   = G_2Hz_lay[G_2Hz_lay[!,:Measurement] .== MeasList[4],:]
            G1_2Hz_4   = G_2Hz_lay[G_2Hz_lay[!,:Measurement] .== MeasList[5],:]

            G1_5Hz_Pre = G_5Hz_lay[G_5Hz_lay[!,:Measurement] .== MeasList[1],:]
            G1_5Hz_1   = G_5Hz_lay[G_5Hz_lay[!,:Measurement] .== MeasList[2],:]
            G1_5Hz_2   = G_5Hz_lay[G_5Hz_lay[!,:Measurement] .== MeasList[3],:]
            G1_5Hz_3   = G_5Hz_lay[G_5Hz_lay[!,:Measurement] .== MeasList[4],:]
            G1_5Hz_4   = G_5Hz_lay[G_5Hz_lay[!,:Measurement] .== MeasList[5],:]

            # find the p value outcome for group between measurements
            P2Hz_Prev1[count[1]]  = pvalue(EqualVarianceTTest(G1_2Hz_Pre[!,:Ratio], G1_2Hz_1[!,:Ratio]))
            P2Hz_Prev2[count[1]]  = pvalue(EqualVarianceTTest(G1_2Hz_Pre[!,:Ratio], G1_2Hz_2[!,:Ratio]))
            P2Hz_Prev3[count[1]]  = pvalue(EqualVarianceTTest(G1_2Hz_Pre[!,:Ratio], G1_2Hz_3[!,:Ratio]))
            P2Hz_Prev4[count[1]]  = pvalue(EqualVarianceTTest(G1_2Hz_Pre[!,:Ratio], G1_2Hz_4[!,:Ratio]))

            P5Hz_Prev1[count[1]]  = pvalue(EqualVarianceTTest(G1_5Hz_Pre[!,:Ratio], G1_5Hz_1[!,:Ratio]))
            P5Hz_Prev2[count[1]]  = pvalue(EqualVarianceTTest(G1_5Hz_Pre[!,:Ratio], G1_5Hz_2[!,:Ratio]))
            P5Hz_Prev3[count[1]]  = pvalue(EqualVarianceTTest(G1_5Hz_Pre[!,:Ratio], G1_5Hz_3[!,:Ratio]))
            P5Hz_Prev4[count[1]]  = pvalue(EqualVarianceTTest(G1_5Hz_Pre[!,:Ratio], G1_5Hz_4[!,:Ratio]))
            
            count[1] = count[1] + 1
            # store the appropriate tags at the same positions in their respective lists
            push!(Group,GroupList[iGrp])
            push!(Layer,LayList[iLay])

        end # layer
    end # comparison of which groups

    BetweenGroup = DataFrame(Group=Group, Layer=Layer, P2Hz_Prev1=P2Hz_Prev1, P2Hz_Prev2=P2Hz_Prev2, P2Hz_Prev3=P2Hz_Prev3, P2Hz_Prev4=P2Hz_Prev4, P5Hz_Prev1=P5Hz_Prev1, P5Hz_Prev2=P5Hz_Prev2, P5Hz_Prev3=P5Hz_Prev3, P5Hz_Prev4=P5Hz_Prev4)

    foldername = "AvrecPeakStats"
    if !isdir(joinpath(data,foldername))
        mkdir(joinpath(data,foldername))
    end
    title = "Ratio_WithinGroups"
    name = joinpath(data,foldername,title) * ".csv"
    CSV.write(name, BetweenGroup)
end


function Peak1_Between(data,Stat2,Stat5)
    # Input: folder path data, table for peak amp/lat of 2 hz and 5 hz, all groups being tested
    # Output: table in folder Data/AvrecPeakStats which contains the pvalue result for an unequal test of variance (2 sample t test) of the peak amp and lat of the first response between each group

    # Setup for simple 2 sample t test from HypothesisTests ->
    # 3 between group comparisons across 5 measurement types and 5 layers:
    CompList    = ["KIC vs KIT" "KIC vs KIV" "KIT vs KIV"]
    LayList     = unique(Stat2[!,:Layer])
    MeasList    = unique(Stat2[!,:Measurement])
    # containers for the table columns: 
    Comparison  = String[]
    Layer       = String[]
    Measurement = String[]   
    P2Hz_Amp    = ones(length(CompList) * length(LayList) * length(MeasList)) # number of rows
    P2Hz_Lat    = ones(length(P2Hz_Amp))
    P5Hz_Amp    = ones(length(P2Hz_Amp))
    P5Hz_Lat    = ones(length(P2Hz_Amp))
    count       = [1]
    # loop through the above - create a table output! :) 
    for iComp = 1:length(CompList)
        # pull out groups from ratio table
        G1_2Hz = Stat2[Stat2[!,:Group] .== CompList[iComp][1:3] ,:]
        G2_2Hz = Stat2[Stat2[!,:Group] .== CompList[iComp][8:10],:]
        G1_5Hz = Stat5[Stat5[!,:Group] .== CompList[iComp][1:3] ,:]
        G2_5Hz = Stat5[Stat5[!,:Group] .== CompList[iComp][8:10],:]

        for iMeas = 1:length(MeasList)
            # pull out the measurment for those groups
            G1_2Hz_meas = G1_2Hz[G1_2Hz[!,:Measurement] .== MeasList[iMeas],:]
            G2_2Hz_meas = G2_2Hz[G2_2Hz[!,:Measurement] .== MeasList[iMeas],:]
            G1_5Hz_meas = G1_5Hz[G1_5Hz[!,:Measurement] .== MeasList[iMeas],:]
            G2_5Hz_meas = G2_5Hz[G2_5Hz[!,:Measurement] .== MeasList[iMeas],:]

            for iLay = 1:length(LayList)
                # further pull out each layer definition
                G1_2Hz_lay = G1_2Hz_meas[G1_2Hz_meas[!,:Layer] .== LayList[iLay],:]
                G2_2Hz_lay = G2_2Hz_meas[G2_2Hz_meas[!,:Layer] .== LayList[iLay],:]
                G1_5Hz_lay = G1_5Hz_meas[G1_5Hz_meas[!,:Layer] .== LayList[iLay],:]
                G2_5Hz_lay = G2_5Hz_meas[G2_5Hz_meas[!,:Layer] .== LayList[iLay],:]

                # find the p value outcome for this comparison at both stimuli conditions
                P2Hz_Amp[count[1]] = pvalue(UnequalVarianceTTest(G1_2Hz_lay[!,:PeakAmp],G2_2Hz_lay[!,:PeakAmp]))
                P2Hz_Lat[count[1]] = pvalue(UnequalVarianceTTest(G1_2Hz_lay[!,:PeakLat],G2_2Hz_lay[!,:PeakLat]))
                P5Hz_Amp[count[1]] = pvalue(UnequalVarianceTTest(G1_5Hz_lay[!,:PeakAmp],G2_5Hz_lay[!,:PeakAmp]))
                P5Hz_Lat[count[1]] = pvalue(UnequalVarianceTTest(G1_5Hz_lay[!,:PeakLat],G2_5Hz_lay[!,:PeakLat]))
                count[1] = count[1] + 1
                # store the appropriate tags at the same positions in their respective lists
                push!(Comparison,CompList[iComp])
                push!(Measurement,MeasList[iMeas])
                push!(Layer,LayList[iLay])
            end # layer
        end # measurement type
    end # comparison of which groups

    BetweenGroup = DataFrame(Comparison=Comparison, Measurement=Measurement, Layer=Layer, P2Hz_Amp=P2Hz_Amp, P2Hz_Lat=P2Hz_Lat, P5Hz_Amp=P5Hz_Amp, P5Hz_Lat=P5Hz_Lat)

    foldername = "AvrecPeakStats"
    if !isdir(joinpath(data,foldername))
        mkdir(joinpath(data,foldername))
    end
    title = "Peak1_BetweenGroups"
    name = joinpath(data,foldername,title) * ".csv"
    CSV.write(name, BetweenGroup)
end

function Peak1_Within(data,Stat2,Stat5)
    # Input: folder path data, table for peak amp/lat of 2 hz and 5 hz, all groups being tested
    # Output: table in folder Data/AvrecPeakStats which contains the pvalue result for an equal test of variance (2 sample t test) of the peak amp and lat of the first response between each measurement to the first

    # Setup for simple 2 sample t test from HypothesisTests ->
    # 4 within group comparisons per group across 5 measurement types and 5 layers:
    GroupList    = ["KIC" "KIT" "KIV"]
    LayList     = unique(Stat2[!,:Layer])
    MeasList    = unique(Stat2[!,:Measurement])
    # containers for the table columns: 
    Group       = String[]
    Layer       = String[]
    # Measurement = String[]   
    P2Hz_Prev1_Amp  = ones(length(GroupList) * length(LayList)) # number of rows (more columns this time)
    P2Hz_Prev2_Amp  = ones(length(P2Hz_Prev1))
    P2Hz_Prev3_Amp  = ones(length(P2Hz_Prev1))
    P2Hz_Prev4_Amp  = ones(length(P2Hz_Prev1))
    P5Hz_Prev1_Amp  = ones(length(P2Hz_Prev1))
    P5Hz_Prev2_Amp  = ones(length(P2Hz_Prev1))
    P5Hz_Prev3_Amp  = ones(length(P2Hz_Prev1))
    P5Hz_Prev4_Amp  = ones(length(P2Hz_Prev1))

    P2Hz_Prev1_Lat  = ones(length(P2Hz_Prev1))
    P2Hz_Prev2_Lat  = ones(length(P2Hz_Prev1))
    P2Hz_Prev3_Lat  = ones(length(P2Hz_Prev1))
    P2Hz_Prev4_Lat  = ones(length(P2Hz_Prev1))
    P5Hz_Prev1_Lat  = ones(length(P2Hz_Prev1))
    P5Hz_Prev2_Lat  = ones(length(P2Hz_Prev1))
    P5Hz_Prev3_Lat  = ones(length(P2Hz_Prev1))
    P5Hz_Prev4_Lat  = ones(length(P2Hz_Prev1))

    count       = [1]
    # loop through the above - create a table output! :) 
    for iGrp = 1:length(GroupList)
        # pull out groups from ratio table
        G_2Hz = Stat2[Stat2[!,:Group] .== GroupList[iGrp],:]
        G_5Hz = Stat5[Stat5[!,:Group] .== GroupList[iGrp],:]

        for iLay = 1:length(LayList)
            # pull out the layer one at a time
            G_2Hz_lay = G_2Hz[G_2Hz[!,:Layer] .== LayList[iLay],:]
            G_5Hz_lay = G_5Hz[G_2Hz[!,:Layer] .== LayList[iLay],:]

            # further pull out each measurement all together
            G1_2Hz_Pre = G_2Hz_lay[G_2Hz_lay[!,:Measurement] .== MeasList[1],:]
            G1_2Hz_1   = G_2Hz_lay[G_2Hz_lay[!,:Measurement] .== MeasList[2],:]
            G1_2Hz_2   = G_2Hz_lay[G_2Hz_lay[!,:Measurement] .== MeasList[3],:]
            G1_2Hz_3   = G_2Hz_lay[G_2Hz_lay[!,:Measurement] .== MeasList[4],:]
            G1_2Hz_4   = G_2Hz_lay[G_2Hz_lay[!,:Measurement] .== MeasList[5],:]

            G1_5Hz_Pre = G_5Hz_lay[G_5Hz_lay[!,:Measurement] .== MeasList[1],:]
            G1_5Hz_1   = G_5Hz_lay[G_5Hz_lay[!,:Measurement] .== MeasList[2],:]
            G1_5Hz_2   = G_5Hz_lay[G_5Hz_lay[!,:Measurement] .== MeasList[3],:]
            G1_5Hz_3   = G_5Hz_lay[G_5Hz_lay[!,:Measurement] .== MeasList[4],:]
            G1_5Hz_4   = G_5Hz_lay[G_5Hz_lay[!,:Measurement] .== MeasList[5],:]

            # find the p value outcome for group between measurements
            P2Hz_Prev1_Amp[count[1]]  = pvalue(EqualVarianceTTest(G1_2Hz_Pre[!,:PeakAmp], G1_2Hz_1[!,:PeakAmp]))
            P2Hz_Prev2_Amp[count[1]]  = pvalue(EqualVarianceTTest(G1_2Hz_Pre[!,:PeakAmp], G1_2Hz_2[!,:PeakAmp]))
            P2Hz_Prev3_Amp[count[1]]  = pvalue(EqualVarianceTTest(G1_2Hz_Pre[!,:PeakAmp], G1_2Hz_3[!,:PeakAmp]))
            P2Hz_Prev4_Amp[count[1]]  = pvalue(EqualVarianceTTest(G1_2Hz_Pre[!,:PeakAmp], G1_2Hz_4[!,:PeakAmp]))

            P5Hz_Prev1_Amp[count[1]]  = pvalue(EqualVarianceTTest(G1_5Hz_Pre[!,:PeakAmp], G1_5Hz_1[!,:PeakAmp]))
            P5Hz_Prev2_Amp[count[1]]  = pvalue(EqualVarianceTTest(G1_5Hz_Pre[!,:PeakAmp], G1_5Hz_2[!,:PeakAmp]))
            P5Hz_Prev3_Amp[count[1]]  = pvalue(EqualVarianceTTest(G1_5Hz_Pre[!,:PeakAmp], G1_5Hz_3[!,:PeakAmp]))
            P5Hz_Prev4_Amp[count[1]]  = pvalue(EqualVarianceTTest(G1_5Hz_Pre[!,:PeakAmp], G1_5Hz_4[!,:PeakAmp]))

            P2Hz_Prev1_Lat[count[1]]  = pvalue(EqualVarianceTTest(G1_2Hz_Pre[!,:PeakLat], G1_2Hz_1[!,:PeakLat]))
            P2Hz_Prev2_Lat[count[1]]  = pvalue(EqualVarianceTTest(G1_2Hz_Pre[!,:PeakLat], G1_2Hz_2[!,:PeakLat]))
            P2Hz_Prev3_Lat[count[1]]  = pvalue(EqualVarianceTTest(G1_2Hz_Pre[!,:PeakLat], G1_2Hz_3[!,:PeakLat]))
            P2Hz_Prev4_Lat[count[1]]  = pvalue(EqualVarianceTTest(G1_2Hz_Pre[!,:PeakLat], G1_2Hz_4[!,:PeakLat]))

            P5Hz_Prev1_Lat[count[1]]  = pvalue(EqualVarianceTTest(G1_5Hz_Pre[!,:PeakLat], G1_5Hz_1[!,:PeakLat]))
            P5Hz_Prev2_Lat[count[1]]  = pvalue(EqualVarianceTTest(G1_5Hz_Pre[!,:PeakLat], G1_5Hz_2[!,:PeakLat]))
            P5Hz_Prev3_Lat[count[1]]  = pvalue(EqualVarianceTTest(G1_5Hz_Pre[!,:PeakLat], G1_5Hz_3[!,:PeakLat]))
            P5Hz_Prev4_Lat[count[1]]  = pvalue(EqualVarianceTTest(G1_5Hz_Pre[!,:PeakLat], G1_5Hz_4[!,:PeakLat]))
            
            count[1] = count[1] + 1
            # store the appropriate tags at the same positions in their respective lists
            push!(Group,GroupList[iGrp])
            push!(Layer,LayList[iLay])

        end # layer
    end # comparison of which groups

    BetweenGroup = DataFrame(Group=Group, Layer=Layer, P2Hz_Prev1_Amp=P2Hz_Prev1_Amp, P2Hz_Prev2_Amp=P2Hz_Prev2_Amp, P2Hz_Prev3_Amp=P2Hz_Prev3_Amp, P2Hz_Prev4_Amp=P2Hz_Prev4_Amp, P5Hz_Prev1_Amp=P5Hz_Prev1_Amp, P5Hz_Prev2_Amp=P5Hz_Prev2_Amp, P5Hz_Prev3_Amp=P5Hz_Prev3_Amp, P5Hz_Prev4_Amp=P5Hz_Prev4_Amp, P2Hz_Prev1_Lat=P2Hz_Prev1_Lat, P2Hz_Prev2_Lat=P2Hz_Prev2_Lat, P2Hz_Prev3_Lat=P2Hz_Prev3_Lat, P2Hz_Prev4_Lat=P2Hz_Prev4_Lat, P5Hz_Prev1_Lat=P5Hz_Prev1_Lat, P5Hz_Prev2_Lat=P5Hz_Prev2_Lat, P5Hz_Prev3_Lat=P5Hz_Prev3_Lat, P5Hz_Prev4_Lat=P5Hz_Prev4_Lat)

    foldername = "AvrecPeakStats"
    if !isdir(joinpath(data,foldername))
        mkdir(joinpath(data,foldername))
    end
    title = "Peak1_WithinGroups"
    name = joinpath(data,foldername,title) * ".csv"
    CSV.write(name, BetweenGroup)
end