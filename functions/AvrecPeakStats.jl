function PeakRatio_Between(data,Stat,whichstim="2Hz",stimtype="CL",trialtype="TA")
    # Input: folder path data, table for ratio of 2 hz and 5 hz, all groups being tested
    # Output: table in folder Data/AvrecPeakStats which contains the pvalue result for an unequal test of variance (2 sample t test) of the ratio of final to first peak response between each group

    # Setup for simple 2 sample t test from HypothesisTests ->
    # 3 between group comparisons across 5 measurement types and 5 layers:
    CompList    = ["KIC vs KIT" "KIC vs KIV" "KIT vs KIV"]
    LayList     = unique(Stat[:,:Layer])
    MeasList    = unique(Stat[:,:Measurement])
    # containers for the table columns: 
    Comparison  = String[]
    Layer       = String[]
    Measurement = String[]   
    PAMP        = ones(length(CompList) * length(LayList) * length(MeasList)) * -1 # number of rows
    PRMS,CDAMP,CDRMS = ones(length(PAMP))*-1,ones(length(PAMP))*-1,ones(length(PAMP))*-1
    
    count       = [1]
    # loop through the above - create a table output! :) 
    for iComp = 1:length(CompList)
        # pull out groups from ratio table
        G1 = Stat[Stat[:,:Group] .== CompList[iComp][1:3] ,:]
        G2 = Stat[Stat[:,:Group] .== CompList[iComp][8:10],:]

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
                CDAMP[count[1]] = effectsize(CohenD(G1_layamp[:,:RatioAMP],G2_layamp[:,:RatioAMP]))
                CDRMS[count[1]] = effectsize(CohenD(G1_layrms[:,:RatioRMS],G2_layrms[:,:RatioRMS]))

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
    CDAMP,CDRMS = CDAMP[CDAMP .!= -1], CDRMS[CDRMS .!= -1]

    BetweenGroup = DataFrame(Comparison=Comparison, Measurement=Measurement, Layer=Layer, PAMP=PAMP, PRMS=PRMS, CDAMP=CDAMP, CDRMS=CDRMS)

    foldername = "AvrecPeakStats"
    if !isdir(joinpath(data,foldername))
        mkdir(joinpath(data,foldername))
    end
    title = "Ratio_BetweenGroups_" * whichstim * "_" * stimtype * "_" * trialtype
    name = joinpath(data,foldername,title) * ".csv"
    CSV.write(name, BetweenGroup)
end

function PeakRatio_Within(data,Stat,whichstim="2Hz",stimtype="CL",trialtype="TA")
    # Input: folder path data, table for ratio of 2 hz and 5 hz, all groups being tested
    # Output: table in folder Data/AvrecPeakStats which contains the pvalue result for an equal test of variance (2 sample t test) of the ratio of final to first peak response between each measurement to the first

    # Setup for simple 2 sample t test from HypothesisTests ->
    # 4 within group comparisons per group across 5 measurement types and 5 layers:
    GroupList    = ["KIC" "KIT"] #  "KIV"
    LayList     = unique(Stat[:,:Layer])
    MeasList    = unique(Stat[:,:Measurement])
    # containers for the table columns: 
    Group       = String[]
    Layer       = String[]
    # Measurement = String[]   
    Prev1_AMP  = ones(length(GroupList) * length(LayList)) # number of rows (more columns this time)
    Prev2_AMP,Prev3_AMP,Prev4_AMP = ones(length(Prev1_AMP)),ones(length(Prev1_AMP)),ones(length(Prev1_AMP))
    Prev1_RMS,Prev2_RMS,Prev3_RMS,Prev4_RMS = ones(length(Prev1_AMP)),ones(length(Prev1_AMP)),ones(length(Prev1_AMP)),ones(length(Prev1_AMP))

    Prev1_CDAMP,Prev2_CDAMP,Prev3_CDAMP,Prev4_CDAMP = ones(length(Prev1_AMP)),ones(length(Prev1_AMP)),ones(length(Prev1_AMP)),ones(length(Prev1_AMP))
    Prev1_CDRMS,Prev2_CDRMS,Prev3_CDRMS,Prev4_CDRMS = ones(length(Prev1_AMP)),ones(length(Prev1_AMP)),ones(length(Prev1_AMP)),ones(length(Prev1_AMP))

    count       = [1]
    # loop through the above - create a table output! :) 
    for iGrp = 1:length(GroupList)
        # pull out groups from ratio table
        Gr_Stat = Stat[Stat[:,:Group] .== GroupList[iGrp],:]

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
            # Cohen's D
            Prev1_CDAMP[count[1]]  = effectsize(CohenD(Gr_Pre[:,:RatioAMP], Gr_1[:,:RatioAMP]))
            Prev2_CDAMP[count[1]]  = effectsize(CohenD(Gr_Pre[:,:RatioAMP], Gr_2[:,:RatioAMP]))
            Prev3_CDAMP[count[1]]  = effectsize(CohenD(Gr_Pre[:,:RatioAMP], Gr_3[:,:RatioAMP]))
            Prev4_CDAMP[count[1]]  = effectsize(CohenD(Gr_Pre[:,:RatioAMP], Gr_4[:,:RatioAMP]))

            Prev1_CDRMS[count[1]]  = effectsize(CohenD(Gr_Pre[:,:RatioRMS], Gr_1[:,:RatioRMS]))
            Prev2_CDRMS[count[1]]  = effectsize(CohenD(Gr_Pre[:,:RatioRMS], Gr_2[:,:RatioRMS]))
            Prev3_CDRMS[count[1]]  = effectsize(CohenD(Gr_Pre[:,:RatioRMS], Gr_3[:,:RatioRMS]))
            Prev4_CDRMS[count[1]]  = effectsize(CohenD(Gr_Pre[:,:RatioRMS], Gr_4[:,:RatioRMS]))
            
            count[1] = count[1] + 1
            # store the appropriate tags at the same positions in their respective lists
            push!(Group,GroupList[iGrp])
            push!(Layer,LayList[iLay])

        end # layer
    end # comparison of which groups

    BetweenGroup = DataFrame(Group=Group, Layer=Layer, Prev1_AMP=Prev1_AMP, Prev2_AMP=Prev2_AMP, Prev3_AMP=Prev3_AMP, Prev4_AMP=Prev4_AMP, Prev1_RMS=Prev1_RMS, Prev2_RMS=Prev2_RMS, Prev3_RMS=Prev3_RMS, Prev4_RMS=Prev4_RMS, Prev1_CDAMP=Prev1_CDAMP, Prev2_CDAMP=Prev2_CDAMP, Prev3_CDAMP=Prev3_CDAMP, Prev4_CDAMP=Prev4_CDAMP, Prev1_CDRMS=Prev1_CDRMS, Prev2_CDRMS=Prev2_CDRMS, Prev3_CDRMS=Prev3_CDRMS, Prev4_CDRMS=Prev4_CDRMS)

    foldername = "AvrecPeakStats"
    if !isdir(joinpath(data,foldername))
        mkdir(joinpath(data,foldername))
    end
    title = "Ratio_WithinGroups_" * whichstim * "_" * stimtype * "_" * trialtype
    name = joinpath(data,foldername,title) * ".csv"
    CSV.write(name, BetweenGroup)
end


function Peak1_Between(data,Stat,whichpeak="1st",whichstim="2Hz",stimtype="CL",trialtype="TA")
    # Input: folder path data, table for peak amp/lat of 2 hz or 5 hz, all groups being tested
    # Output: table in folder Data/AvrecPeakStats which contains the pvalue result for an unequal test of variance (2 sample t test) of the peak amp and lat of the first response between each group

    # Setup for simple 2 sample t test from HypothesisTests ->
    # 3 between group comparisons across 5 measurement types and 5 layers:
    CompList    = ["KIC vs KIT" "KIC vs KIV" "KIT vs KIV"]
    LayList     = unique(Stat[:,:Layer])
    MeasList    = unique(Stat[:,:Measurement])
    # containers for the table columns: 
    Comparison  = String[]
    Layer       = String[]
    Measurement = String[]   
    PAmp        = ones(length(CompList) * length(LayList) * length(MeasList)) * -1 # number of rows
    PLat, PRMS  = ones(length(PAmp))*-1 , ones(length(PAmp))*-1
    CDAmp,CDLat,CDRMS = ones(length(PAmp))*-1 , ones(length(PAmp))*-1, ones(length(PAmp))*-1

    count       = [1]
    # loop through the above - create a table output! :) 
    for iComp = 1:length(CompList)
        # pull out groups from ratio table
        G1 = Stat[Stat[:,:Group] .== CompList[iComp][1:3] ,:]
        G2 = Stat[Stat[:,:Group] .== CompList[iComp][8:10],:]

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
                if isempty(G1_lay) || isempty(G2_lay) || size(G1_layamp)[1] <= 3 || size(G2_layamp)[1] <= 3
                    continue 
                end

                PAmp[count[1]]  = pvalue(UnequalVarianceTTest(G1_layamp[:,:PeakAmp],G2_layamp[:,:PeakAmp]))
                CDAmp[count[1]] = effectsize(CohenD(G1_layamp[:,:PeakAmp],G2_layamp[:,:PeakAmp]))
                PLat[count[1]]  = pvalue(UnequalVarianceTTest(G1_laylat[:,:PeakLat],G2_laylat[:,:PeakLat]))
                CDLat[count[1]] = effectsize(CohenD(G1_laylat[:,:PeakLat],G2_laylat[:,:PeakLat]))
                PRMS[count[1]]  = pvalue(UnequalVarianceTTest(G1_layrms[:,:RMS],G2_layrms[:,:RMS]))
                CDRMS[count[1]] = effectsize(CohenD(G1_layamp[:,:RMS],G2_layamp[:,:RMS]))
                count[1] = count[1] + 1
                # store the appropriate tags at the same positions in their respective lists
                push!(Comparison,CompList[iComp])
                push!(Measurement,MeasList[iMeas])
                push!(Layer,LayList[iLay])
            end # layer
        end # measurement type
    end # comparison of which groups

    # need to remove unused rows that still equal -1
    PAmp, PLat, PRMS  = PAmp[PAmp .!= -1], PLat[PLat .!= -1], PRMS[PRMS .!= -1]
    CDAmp,CDLat,CDRMS = CDAmp[CDAmp .!= -1], CDLat[CDLat .!= -1], CDRMS[CDRMS .!= -1] 
    BetweenGroup = DataFrame(Comparison=Comparison, Measurement=Measurement, Layer=Layer, PAmp=PAmp, PLat=PLat, PRMS=PRMS, CDAmp=CDAmp, CDLat=CDLat, CDRMS=CDRMS)

    foldername = "AvrecPeakStats"
    if !isdir(joinpath(data,foldername))
        mkdir(joinpath(data,foldername))
    end
    title = whichpeak * "_BetweenGroups_" * whichstim * "_" * stimtype * "_" * trialtype
    name = joinpath(data,foldername,title) * ".csv"
    CSV.write(name, BetweenGroup)
end

function Peak1_Within(data,Stat,whichpeak="1st",whichstim="2Hz",stimtype="CL",trialtype="TA")
    # Input: folder path data, table for peak amp/lat of 2 hz and 5 hz, all groups being tested
    # Output: table in folder Data/AvrecPeakStats which contains the pvalue result for an equal test of variance (2 sample t test) of the peak amp and lat of the first response between each measurement to the first

    # Setup for simple 2 sample t test from HypothesisTests ->
    # 4 within group comparisons per group across 5 measurement types and 5 layers:
    GroupList    = ["KIC" "KIT"] #"KIV"
    LayList     = unique(Stat[:,:Layer])
    MeasList    = unique(Stat[:,:Measurement])
    # containers for the table columns: 
    Group       = String[]
    Layer       = String[]
    # Measurement = String[]   
    Prev1_Amp  = ones(length(GroupList) * length(LayList))*-1 # number of rows 
    Prev2_Amp,Prev3_Amp,Prev4_Amp  = ones(length(Prev1_Amp))*-1,ones(length(Prev1_Amp))*-1,ones(length(Prev1_Amp))*-1
    Prev1_Lat,Prev2_Lat,Prev3_Lat,Prev4_Lat  = ones(length(Prev1_Amp))*-1,ones(length(Prev1_Amp))*-1,ones(length(Prev1_Amp))*-1,ones(length(Prev1_Amp))*-1
    Prev1_RMS,Prev2_RMS,Prev3_RMS,Prev4_RMS  = ones(length(Prev1_Amp))*-1,ones(length(Prev1_Amp))*-1,ones(length(Prev1_Amp))*-1,ones(length(Prev1_Amp))*-1
    
    Prev1_CDAmp,Prev2_CDAmp,Prev3_CDAmp,Prev4_CDAmp  = ones(length(Prev1_Amp))*-1,ones(length(Prev1_Amp))*-1,ones(length(Prev1_Amp))*-1,ones(length(Prev1_Amp))*-1
    Prev1_CDLat,Prev2_CDLat,Prev3_CDLat,Prev4_CDLat  = ones(length(Prev1_Amp))*-1,ones(length(Prev1_Amp))*-1,ones(length(Prev1_Amp))*-1,ones(length(Prev1_Amp))*-1
    Prev1_CDRMS,Prev2_CDRMS,Prev3_CDRMS,Prev4_CDRMS  = ones(length(Prev1_Amp))*-1,ones(length(Prev1_Amp))*-1,ones(length(Prev1_Amp))*-1,ones(length(Prev1_Amp))*-1

    count       = [1]
    # loop through the above - create a table output! :) 
    for iGrp = 1:length(GroupList)
        # pull out groups from ratio table
        Gr_Stat = Stat[Stat[:,:Group] .== GroupList[iGrp],:]
        
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

            # if size(GP_amp)[1] >= 3 || size(G1_amp)[1] >= 3 || size(GP_rms)[1] >= 3 || size(G1_rms)[1] >= 3 
            #     continue
            # end - this is causing it to not take the pvalues for ANY of it but I do need a solution for the Cohen's D output if there is an error there
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
            #Cohen's D
            if isnan(Prev1_Amp[count[1]]) # I hate this solution so much but CohenD throws an error instead of giving a NaN output like the pvalue above
                Prev1_CDAmp[count[1]]  = NaN
                Prev1_CDLat[count[1]]  = NaN
            else
                Prev1_CDAmp[count[1]]  = effectsize(CohenD(GP_amp[:,:PeakAmp], G1_amp[:,:PeakAmp]))
                Prev1_CDLat[count[1]]  = effectsize(CohenD(GP_amp[:,:PeakLat], G1_amp[:,:PeakLat]))
            end
            if isnan(Prev2_Amp[count[1]])
                Prev2_CDAmp[count[1]]  = NaN
                Prev2_CDLat[count[1]]  = NaN
            else
                Prev2_CDAmp[count[1]]  = effectsize(CohenD(GP_amp[:,:PeakAmp], G2_amp[:,:PeakAmp]))
                Prev2_CDLat[count[1]]  = effectsize(CohenD(GP_amp[:,:PeakLat], G2_amp[:,:PeakLat]))
            end
            if isnan(Prev3_Amp[count[1]])
                Prev3_CDAmp[count[1]]  = NaN
                Prev3_CDLat[count[1]]  = NaN
            else
                Prev3_CDAmp[count[1]]  = effectsize(CohenD(GP_amp[:,:PeakAmp], G3_amp[:,:PeakAmp]))
                Prev3_CDLat[count[1]]  = effectsize(CohenD(GP_amp[:,:PeakLat], G3_amp[:,:PeakLat]))
            end
            if isnan(Prev4_Amp[count[1]])
                Prev4_CDAmp[count[1]]  = NaN
                Prev4_CDLat[count[1]]  = NaN 
            else
                Prev4_CDAmp[count[1]]  = effectsize(CohenD(GP_amp[:,:PeakAmp], G4_amp[:,:PeakAmp]))
                Prev4_CDLat[count[1]]  = effectsize(CohenD(GP_amp[:,:PeakLat], G4_amp[:,:PeakLat]))
            end

            if isnan(Prev1_RMS[count[1]])
                Prev1_CDRMS[count[1]]  = NaN
            else
                Prev1_CDRMS[count[1]]  = effectsize(CohenD(GP_rms[:,:RMS], G1_rms[:,:RMS]))
            end
            if isnan(Prev2_RMS[count[1]])
                Prev2_CDRMS[count[1]]  = NaN
            else
                Prev2_CDRMS[count[1]]  = effectsize(CohenD(GP_rms[:,:RMS], G2_rms[:,:RMS]))
            end
            if isnan(Prev3_RMS[count[1]])
                Prev3_CDRMS[count[1]]  = NaN
            else
                Prev3_CDRMS[count[1]]  = effectsize(CohenD(GP_rms[:,:RMS], G3_rms[:,:RMS]))
            end
            if isnan(Prev4_RMS[count[1]])
                Prev4_CDRMS[count[1]]  = NaN
            else
                Prev4_CDRMS[count[1]]  = effectsize(CohenD(GP_rms[:,:RMS], G4_rms[:,:RMS]))
            end

            count[1] = count[1] + 1
            # store the appropriate tags at the same positions in their respective lists
            push!(Group,GroupList[iGrp])
            push!(Layer,LayList[iLay])

        end # layer
    end # comparison of which groups

    Prev1_Amp,Prev2_Amp,Prev3_Amp,Prev4_Amp  = Prev1_Amp[Prev1_Amp .!= -1], Prev2_Amp[Prev2_Amp .!= -1], Prev3_Amp[Prev3_Amp .!= -1], Prev4_Amp[Prev4_Amp .!= -1]
    Prev1_Lat,Prev2_Lat,Prev3_Lat,Prev4_Lat  = Prev1_Lat[Prev1_Lat .!= -1], Prev2_Lat[Prev2_Lat .!= -1], Prev3_Lat[Prev3_Lat .!= -1], Prev4_Lat[Prev4_Lat .!= -1]
    Prev1_RMS,Prev2_RMS,Prev3_RMS,Prev4_RMS  = Prev1_RMS[Prev1_RMS .!= -1], Prev2_RMS[Prev2_RMS .!= -1], Prev3_RMS[Prev3_RMS .!= -1], Prev4_RMS[Prev4_RMS .!= -1]

    Prev1_CDAmp,Prev2_CDAmp,Prev3_CDAmp,Prev4_CDAmp  = Prev1_CDAmp[Prev1_CDAmp .!= -1], Prev2_CDAmp[Prev2_CDAmp .!= -1], Prev3_CDAmp[Prev3_CDAmp .!= -1], Prev4_CDAmp[Prev4_CDAmp .!= -1]
    Prev1_CDLat,Prev2_CDLat,Prev3_CDLat,Prev4_CDLat  = Prev1_CDLat[Prev1_CDLat .!= -1], Prev2_CDLat[Prev2_CDLat .!= -1], Prev3_CDLat[Prev3_CDLat .!= -1], Prev4_CDLat[Prev4_CDLat .!= -1]
    Prev1_CDRMS,Prev2_CDRMS,Prev3_CDRMS,Prev4_CDRMS  = Prev1_CDRMS[Prev1_CDRMS .!= -1], Prev2_CDRMS[Prev2_CDRMS .!= -1], Prev3_CDRMS[Prev3_CDRMS .!= -1], Prev4_CDRMS[Prev4_CDRMS .!= -1]

    BetweenGroup = DataFrame(Group=Group, Layer=Layer, Prev1_Amp=Prev1_Amp, Prev2_Amp=Prev2_Amp, Prev3_Amp=Prev3_Amp, Prev4_Amp=Prev4_Amp,  Prev1_Lat=Prev1_Lat, Prev2_Lat=Prev2_Lat, Prev3_Lat=Prev3_Lat, Prev4_Lat=Prev4_Lat, Prev1_RMS=Prev1_RMS, Prev2_RMS=Prev2_RMS, Prev3_RMS=Prev3_RMS, Prev4_RMS=Prev4_RMS, Prev1_CDAmp=Prev1_CDAmp, Prev2_CDAmp=Prev2_CDAmp, Prev3_CDAmp=Prev3_CDAmp, Prev4_CDAmp=Prev4_CDAmp,  Prev1_CDLat=Prev1_CDLat, Prev2_CDLat=Prev2_CDLat, Prev3_CDLat=Prev3_CDLat, Prev4_CDLat=Prev4_CDLat, Prev1_CDRMS=Prev1_CDRMS, Prev2_CDRMS=Prev2_CDRMS, Prev3_CDRMS=Prev3_CDRMS, Prev4_CDRMS=Prev4_CDRMS)

    foldername = "AvrecPeakStats"
    if !isdir(joinpath(data,foldername))
        mkdir(joinpath(data,foldername))
    end
    title = whichpeak * "_WithinGroups_" * whichstim * "_" * stimtype * "_" * trialtype
    name = joinpath(data,foldername,title) * ".csv"
    CSV.write(name, BetweenGroup)
end