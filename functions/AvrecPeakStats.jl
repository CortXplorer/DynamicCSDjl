function PeakRatio_Between(data,StatTab,whichstim="2Hz",stimtype="CL",trialtype="TA")
    # Input: folder path data, table for ratio of 2 hz and 5 hz, all groups being tested
    # Output: table in folder Data/AvrecPeakStats which contains the pvalue result for an unequal test of variance (2 sample t test) of the ratio of final to first peak response between each group

    # Setup for simple 2 sample t test from HypothesisTests ->
    # 3 between group comparisons across 5 measurement types and 5 layers:
    CompList    = ["KIC vs KIT" "KIC vs KIV" "KIT vs KIV"]
    LayList     = unique(StatTab[:,:Layer])
    MeasList    = unique(StatTab[:,:Measurement])
    # containers for the table columns: 
    Comparison  = String[]
    Layer       = String[]
    Measurement = String[]   
    PAMP        = ones(length(CompList) * length(LayList) * length(MeasList)) * -1 # number of rows
    PRMS        = ones(length(PAMP)) * -1
    
    count       = [1]
    # loop through the above - create a table output! :) 
    for iComp = 1:length(CompList)
        # pull out groups from ratio table
        G1 = StatTab[StatTab[:,:Group] .== CompList[iComp][1:3] ,:]
        G2 = StatTab[StatTab[:,:Group] .== CompList[iComp][8:10],:]

        for iMeas = 1:length(MeasList)
            # pull out the measurment for those groups
            G1_meas = G1[G1[:,:Measurement] .== MeasList[iMeas],:]
            G2_meas = G2[G2[:,:Measurement] .== MeasList[iMeas],:]

            for iLay = 1:length(LayList)
                # further pull out each layer definition       
                G1_lay = G1_meas[G1_meas[:,:Layer] .== LayList[iLay],:]
                G1_layamp = filter(row -> ! isnan(row.RatioAMP), G1_lay)
                G1_layrms = filter(row -> ! isnan(row.RatioRMS), G1_lay)
                G2_lay = G2_meas[G2_meas[:,:Layer] .== LayList[iLay],:]
                G2_layamp = filter(row -> ! isnan(row.RatioAMP), G2_lay)
                G2_layrms = filter(row -> ! isnan(row.RatioRMS), G2_lay)

                if isempty(G1_lay) || isempty(G2_lay)
                    continue 
                end
                # find the p value outcome for this comparison at both stimuli conditions
                PAMP[count[1]] = pvalue(UnequalVarianceTTest(G1_layamp[:,:RatioAMP],G2_layamp[:,:RatioAMP]))
                PRMS[count[1]] = pvalue(UnequalVarianceTTest(G1_layrms[:,:RatioRMS],G2_layrms[:,:RatioRMS]))

                count[1] = count[1] + 1
                # store the appropriate tags at the same positions in their respective lists
                push!(Comparison,CompList[iComp])
                push!(Measurement,MeasList[iMeas])
                push!(Layer,LayList[iLay])
            end # layer
        end # measurement type
    end # comparison of which groups

    PAMP = PAMP[PAMP .!= -1] # need to remove unused rows that still equal 1
    PRMS = PRMS[PRMS .!= -1]

    BetweenGroup = DataFrame(Comparison=Comparison, Measurement=Measurement, Layer=Layer, PAMP=PAMP, PRMS=PRMS)

    foldername = "AvrecPeakStats"
    if !isdir(joinpath(data,foldername))
        mkdir(joinpath(data,foldername))
    end
    title = "Ratio_BetweenGroups_" * whichstim * "_" * stimtype * "_" * trialtype
    name = joinpath(data,foldername,title) * ".csv"
    CSV.write(name, BetweenGroup)
end

function PeakRatio_Within(data,StatTab,whichstim="2Hz",stimtype="CL",trialtype="TA")
    # Input: folder path data, table for ratio of 2 hz and 5 hz, all groups being tested
    # Output: table in folder Data/AvrecPeakStats which contains the pvalue result for an equal test of variance (2 sample t test) of the ratio of final to first peak response between each measurement to the first

    # Setup for simple 2 sample t test from HypothesisTests ->
    # 4 within group comparisons per group across 5 measurement types and 5 layers:
    GroupList    = ["KIC" "KIT"] #  "KIV"
    LayList     = unique(StatTab[:,:Layer])
    MeasList    = unique(StatTab[:,:Measurement])
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
        Gr_Stat = StatTab[StatTab[:,:Group] .== GroupList[iGrp],:]

        for iLay = 1:length(LayList)
            # pull out the layer one at a time
            Gr_lay = Gr_Stat[Gr_Stat[:,:Layer] .== LayList[iLay],:]

            # further pull out each measurement all together
            Gr_Pre = Gr_lay[Gr_lay[:,:Measurement] .== MeasList[1],:]
            Gr_Pre = filter(row -> ! isnan(row.RatioAMP), Gr_Pre)
            Gr_1   = Gr_lay[Gr_lay[:,:Measurement] .== MeasList[2],:]
            Gr_1 = filter(row -> ! isnan(row.RatioAMP), Gr_1)
            Gr_2   = Gr_lay[Gr_lay[:,:Measurement] .== MeasList[3],:]
            Gr_2 = filter(row -> ! isnan(row.RatioAMP), Gr_2)
            Gr_3   = Gr_lay[Gr_lay[:,:Measurement] .== MeasList[4],:]
            Gr_3 = filter(row -> ! isnan(row.RatioAMP), Gr_3)
            Gr_4   = Gr_lay[Gr_lay[:,:Measurement] .== MeasList[5],:]
            Gr_4 = filter(row -> ! isnan(row.RatioAMP), Gr_4)

            # find the p value outcome for group between measurements
            Prev1_AMP[count[1]]  = pvalue(EqualVarianceTTest(Gr_Pre[:,:RatioAMP], Gr_1[:,:RatioAMP]))
            Prev2_AMP[count[1]]  = pvalue(EqualVarianceTTest(Gr_Pre[:,:RatioAMP], Gr_2[:,:RatioAMP]))
            Prev3_AMP[count[1]]  = pvalue(EqualVarianceTTest(Gr_Pre[:,:RatioAMP], Gr_3[:,:RatioAMP]))
            Prev4_AMP[count[1]]  = pvalue(EqualVarianceTTest(Gr_Pre[:,:RatioAMP], Gr_4[:,:RatioAMP]))

            Prev1_RMS[count[1]]  = pvalue(EqualVarianceTTest(Gr_Pre[:,:RatioRMS], Gr_1[:,:RatioRMS]))
            Prev2_RMS[count[1]]  = pvalue(EqualVarianceTTest(Gr_Pre[:,:RatioRMS], Gr_2[:,:RatioRMS]))
            Prev3_RMS[count[1]]  = pvalue(EqualVarianceTTest(Gr_Pre[:,:RatioRMS], Gr_3[:,:RatioRMS]))
            Prev4_RMS[count[1]]  = pvalue(EqualVarianceTTest(Gr_Pre[:,:RatioRMS], Gr_4[:,:RatioRMS]))

            
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
    title = "Ratio_WithinGroups_" * whichstim * "_" * stimtype * "_" * trialtype
    name = joinpath(data,foldername,title) * ".csv"
    CSV.write(name, BetweenGroup)
end


function Peak1_Between(data,StatTab,whichpeak="First",whichstim="2Hz",stimtype="CL",trialtype="TA")
    # Input: folder path data, table for peak amp/lat of 2 hz or 5 hz, all groups being tested
    # Output: table in folder Data/AvrecPeakStats which contains the pvalue result for an unequal test of variance (2 sample t test) of the peak amp and lat of the first response between each group

    # Setup for simple 2 sample t test from HypothesisTests ->
    # 3 between group comparisons across 5 measurement types and 5 layers:
    CompList    = ["KIC vs KIT" "KIC vs KIV" "KIT vs KIV"]
    LayList     = unique(StatTab[:,:Layer])
    MeasList    = unique(StatTab[:,:Measurement])
    # containers for the table columns: 
    Comparison  = String[]
    Layer       = String[]
    Measurement = String[]   
    PAmp        = ones(length(CompList) * length(LayList) * length(MeasList)) * -1 # number of rows
    PLat    = ones(length(PAmp)) * -1
    PRMS    = ones(length(PAmp)) * -1

    count       = [1]
    # loop through the above - create a table output! :) 
    for iComp = 1:length(CompList)
        # pull out groups from ratio table
        G1 = StatTab[StatTab[:,:Group] .== CompList[iComp][1:3] ,:]
        G2 = StatTab[StatTab[:,:Group] .== CompList[iComp][8:10],:]

        for iMeas = 1:length(MeasList)
            # pull out the measurment for those groups
            G1_meas = G1[G1[:,:Measurement] .== MeasList[iMeas],:]
            G2_meas = G2[G2[:,:Measurement] .== MeasList[iMeas],:]

            for iLay = 1:length(LayList)
                # further pull out each layer definition
                G1_lay = G1_meas[G1_meas[:,:Layer] .== LayList[iLay],:]
                G1_layamp = filter(row -> ! isnan(row.PeakAmp), G1_lay)
                G1_laylat = filter(row -> ! isnan(row.PeakLat), G1_lay)
                G1_layrms = filter(row -> ! isnan(row.RMS), G1_lay)
                G2_lay = G2_meas[G2_meas[:,:Layer] .== LayList[iLay],:]
                G2_layamp = filter(row -> ! isnan(row.PeakAmp), G2_lay)
                G2_laylat = filter(row -> ! isnan(row.PeakLat), G2_lay)
                G2_layrms = filter(row -> ! isnan(row.RMS), G2_lay)

                # find the p value outcome for this comparison at both stimuli conditions
                if isempty(G1_lay) || isempty(G2_lay)
                    continue 
                end
                PAmp[count[1]] = pvalue(UnequalVarianceTTest(G1_layamp[:,:PeakAmp],G2_layamp[:,:PeakAmp]))
                PLat[count[1]] = pvalue(UnequalVarianceTTest(G1_laylat[:,:PeakLat],G2_laylat[:,:PeakLat]))
                PRMS[count[1]] = pvalue(UnequalVarianceTTest(G1_layrms[:,:RMS],G2_layrms[:,:RMS]))
                count[1] = count[1] + 1
                # store the appropriate tags at the same positions in their respective lists
                push!(Comparison,CompList[iComp])
                push!(Measurement,MeasList[iMeas])
                push!(Layer,LayList[iLay])
            end # layer
        end # measurement type
    end # comparison of which groups

    PAmp = PAmp[PAmp .!= -1] # need to remove unused rows that still equal 1
    PLat = PLat[PLat .!= -1] # peak lat can actually be 1 so take if peak amp isn't 1
    PRMS = PRMS[PRMS .!= -1]
    BetweenGroup = DataFrame(Comparison=Comparison, Measurement=Measurement, Layer=Layer, PAmp=PAmp, PLat=PLat, PRMS=PRMS)

    foldername = "AvrecPeakStats"
    if !isdir(joinpath(data,foldername))
        mkdir(joinpath(data,foldername))
    end
    title = whichpeak * "_BetweenGroups_" * whichstim * "_" * stimtype * "_" * trialtype
    name = joinpath(data,foldername,title) * ".csv"
    CSV.write(name, BetweenGroup)
end

function Peak1_Within(data,StatTab,whichpeak="First",whichstim="2Hz",stimtype="CL",trialtype="TA")
    # Input: folder path data, table for peak amp/lat of 2 hz and 5 hz, all groups being tested
    # Output: table in folder Data/AvrecPeakStats which contains the pvalue result for an equal test of variance (2 sample t test) of the peak amp and lat of the first response between each measurement to the first

    # Setup for simple 2 sample t test from HypothesisTests ->
    # 4 within group comparisons per group across 5 measurement types and 5 layers:
    GroupList    = ["KIC" "KIT"] #"KIV"
    LayList     = unique(StatTab[:,:Layer])
    MeasList    = unique(StatTab[:,:Measurement])
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
        Gr_Stat = StatTab[StatTab[:,:Group] .== GroupList[iGrp],:]
        
        for iLay = 1:length(LayList)
            # pull out the layer one at a time
            Gr_lay = Gr_Stat[Gr_Stat[:,:Layer] .== LayList[iLay],:]

            # further pull out each measurement all together
            Gr_Pre = Gr_lay[Gr_lay[:,:Measurement] .== MeasList[1],:]
            GP_amp = filter(row -> ! isnan(row.PeakAmp), Gr_Pre)
            GP_rms = filter(row -> ! isnan(row.RMS), Gr_Pre)
            Gr_1   = Gr_lay[Gr_lay[:,:Measurement] .== MeasList[2],:]
            G1_amp = filter(row -> ! isnan(row.PeakAmp), Gr_1)
            G1_rms = filter(row -> ! isnan(row.RMS), Gr_1)
            Gr_2   = Gr_lay[Gr_lay[:,:Measurement] .== MeasList[3],:]
            G2_amp = filter(row -> ! isnan(row.PeakAmp), Gr_2)
            G2_rms = filter(row -> ! isnan(row.RMS), Gr_2)
            Gr_3   = Gr_lay[Gr_lay[:,:Measurement] .== MeasList[4],:]
            G3_amp = filter(row -> ! isnan(row.PeakAmp), Gr_3)
            G3_rms = filter(row -> ! isnan(row.RMS), Gr_3)
            Gr_4   = Gr_lay[Gr_lay[:,:Measurement] .== MeasList[5],:]
            G4_amp = filter(row -> ! isnan(row.PeakAmp), Gr_4)
            G4_rms = filter(row -> ! isnan(row.RMS), Gr_4)

            # find the p value outcome for group between measurements
            Prev1_Amp[count[1]]  = pvalue(EqualVarianceTTest(GP_amp[:,:PeakAmp], G1_amp[:,:PeakAmp]))
            Prev2_Amp[count[1]]  = pvalue(EqualVarianceTTest(GP_amp[:,:PeakAmp], G2_amp[:,:PeakAmp]))
            Prev3_Amp[count[1]]  = pvalue(EqualVarianceTTest(GP_amp[:,:PeakAmp], G3_amp[:,:PeakAmp]))
            Prev4_Amp[count[1]]  = pvalue(EqualVarianceTTest(GP_amp[:,:PeakAmp], G4_amp[:,:PeakAmp]))

            Prev1_Lat[count[1]]  = pvalue(EqualVarianceTTest(GP_amp[:,:PeakLat], G1_amp[:,:PeakLat]))
            Prev2_Lat[count[1]]  = pvalue(EqualVarianceTTest(GP_amp[:,:PeakLat], G2_amp[:,:PeakLat]))
            Prev3_Lat[count[1]]  = pvalue(EqualVarianceTTest(GP_amp[:,:PeakLat], G3_amp[:,:PeakLat]))
            Prev4_Lat[count[1]]  = pvalue(EqualVarianceTTest(GP_amp[:,:PeakLat], G4_amp[:,:PeakLat]))

            Prev1_RMS[count[1]]  = pvalue(EqualVarianceTTest(GP_rms[:,:RMS], G1_rms[:,:RMS]))
            Prev2_RMS[count[1]]  = pvalue(EqualVarianceTTest(GP_rms[:,:RMS], G2_rms[:,:RMS]))
            Prev3_RMS[count[1]]  = pvalue(EqualVarianceTTest(GP_rms[:,:RMS], G3_rms[:,:RMS]))
            Prev4_RMS[count[1]]  = pvalue(EqualVarianceTTest(GP_rms[:,:RMS], G4_rms[:,:RMS]))
            
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
    title = whichpeak * "_WithinGroups_" * whichstim * "_" * stimtype * "_" * trialtype
    name = joinpath(data,foldername,title) * ".csv"
    CSV.write(name, BetweenGroup)
end