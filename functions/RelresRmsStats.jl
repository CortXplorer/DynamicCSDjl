function Relres1Rms(figs,Tab,whichpeak="1st",whichstim="2Hz",savetype=".pdf",stimtype="CL",trialtype="ST")
    # Input: folder path figs, table for rms, all groups being plotted
    # Output: figures in folder Relres1Rms of the first rms response all groups

    foldername = "Relres1Rms"
    if !isdir(joinpath(figs,foldername))
        mkdir(joinpath(figs,foldername))
    end

    ### by measurement ###
    Tab_rms = filter(row -> ! isnan(row.RMS), Tab)

    Title = whichpeak * " RMS at " * whichstim * " " * stimtype * " "
    boxplot = @df Tab_rms groupedboxplot(:Measurement, :RMS, group = :Group, bar_position = :dodge, lab= ["Control" "Treated" "Virus Control"], title=Title, xlab = "Measurement", ylab = "RMS [mV/mm²]")

    name = joinpath(figs,foldername,Title) * trialtype * savetype
    savefig(boxplot, name);

end

function Rms1_Between(data,Stat,whichpeak="1st",whichstim="2Hz",stimtype="CL",trialtype="ST")
    # Input: folder path data, table for peak amp/lat of 2 hz or 5 hz, all groups being tested
    # Output: table in folder Data/AvrecPeakStats which contains the pvalue result for a welch's t test (2 sample t test) of the peak amp and lat of the first response between each group

    # Setup for simple 2 sample t test from HypothesisTests ->
    # 3 between group comparisons across 5 measurement types and 5 layers:
    CompList    = ["KIC vs KIT" "KIC vs KIV" "KIT vs KIV"]
    MeasList    = unique(Stat[:,:Measurement])
    # containers for the table columns: 
    Comparison  = String[]
    Measurement = String[]   
    PRMS        = ones(length(CompList) * length(MeasList)) * -1 # number of rows
    CDRMS = ones(length(PRMS))*-1 

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

            # remove NaNs
            G1_rms = filter(row -> ! isnan(row.RMS), G1_meas)
            G2_rms = filter(row -> ! isnan(row.RMS), G2_meas)

            # find the p value outcome for this comparison at both stimuli conditions
            if isempty(G1_rms) || isempty(G2_rms) || size(G1_rms)[1] <= 3 || size(G2_rms)[1] <= 3
                continue 
            end

            PRMS[count[1]]  = pvalue(UnequalVarianceTTest(G1_rms[:,:RMS],G2_rms[:,:RMS]))
            CDRMS[count[1]] = effectsize(CohenD(G1_rms[:,:RMS],G2_rms[:,:RMS]))
            count[1] = count[1] + 1
            # store the appropriate tags at the same positions in their respective lists
            push!(Comparison,CompList[iComp])
            push!(Measurement,MeasList[iMeas])
        end # measurement type
    end # comparison of which groups

    # need to remove unused rows that still equal -1
    PRMS  = PRMS[PRMS .!= -1]
    CDRMS = CDRMS[CDRMS .!= -1] 
    BetweenGroup = DataFrame(Comparison=Comparison, Measurement=Measurement, PRMS=PRMS, CDRMS=CDRMS)

    foldername = "RelresStats"
    if !isdir(joinpath(data,foldername))
        mkdir(joinpath(data,foldername))
    end
    title = whichpeak * "_BetweenGroups_" * whichstim * "_" * stimtype * "_" * trialtype
    name = joinpath(data,foldername,title) * ".csv"
    CSV.write(name, BetweenGroup)
end

function Rms1_Within(data,Stat,whichpeak="1st",whichstim="2Hz",stimtype="CL",trialtype="TA")
    # Input: folder path data, table for peak amp/lat of 2 hz and 5 hz, all groups being tested
    # Output: table in folder Data/AvrecPeakStats which contains the pvalue result for a welch's t test (2 sample t test) of the peak amp and lat of the first response between each measurement to the first. Single trials are NOT equal before and after laser so a paired sample t test would be throwing out data

    # Setup for simple 2 sample t test from HypothesisTests ->
    # 4 within group comparisons per group across 5 measurement types and 5 layers:
    GroupList    = ["KIC" "KIT" "KIV"] 
    MeasList    = unique(Stat[:,:Measurement])
    # containers for the table columns: 
    Group       = String[]
    # Measurement = String[]   
    Prev1_RMS = ones(length(GroupList)) *-1
    Prev2_RMS,Prev3_RMS,Prev4_RMS  = ones(length(GroupList))*-1 , ones(length(GroupList))*-1 , ones(length(GroupList))*-1
    
    Prev1_CDRMS,Prev2_CDRMS,Prev3_CDRMS,Prev4_CDRMS  = ones(length(GroupList))*-1,ones(length(GroupList))*-1,ones(length(GroupList))*-1,ones(length(GroupList))*-1

    count = [1]
    # loop through the above - create a table output! :) 
    for iGrp = 1:length(GroupList)
        # pull out groups from ratio table
        Gr_Stat = Stat[Stat[:,:Group] .== GroupList[iGrp],:]

        # further pull out each measurement all together
        Gr_Pre = Gr_Stat[Gr_Stat[:,:Measurement] .== MeasList[1],:]
        GP_rms = filter(row -> ! isnan(row.RMS), Gr_Pre)
        Gr_1   = Gr_Stat[Gr_Stat[:,:Measurement] .== MeasList[2],:]
        G1_rms = filter(row -> ! isnan(row.RMS), Gr_1)
        Gr_2   = Gr_Stat[Gr_Stat[:,:Measurement] .== MeasList[3],:]
        G2_rms = filter(row -> ! isnan(row.RMS), Gr_2)
        Gr_3   = Gr_Stat[Gr_Stat[:,:Measurement] .== MeasList[4],:]
        G3_rms = filter(row -> ! isnan(row.RMS), Gr_3)
        Gr_4   = Gr_Stat[Gr_Stat[:,:Measurement] .== MeasList[5],:]
        G4_rms = filter(row -> ! isnan(row.RMS), Gr_4)

        # if size(GP_amp)[1] >= 3 || size(G1_amp)[1] >= 3 || size(GP_rms)[1] >= 3 || size(G1_rms)[1] >= 3 
        #     continue
        # end - this is causing it to not take the pvalues for ANY of it but I do need a solution for the Cohen's D output if there is an error there
        # find the p value outcome for group between measurements
        if isempty(GP_rms) || isempty(G1_rms) || (size(GP_rms)[1]==1 && size(G1_rms)[1]==1)
            Prev1_RMS[count[1]] = NaN
        else
            Prev1_RMS[count[1]] = pvalue(UnequalVarianceTTest(GP_rms[:,:RMS], G1_rms[:,:RMS]))
        end
        if isempty(GP_rms) || isempty(G2_rms) || (size(GP_rms)[1]==1 && size(G2_rms)[1]==1)
            Prev2_RMS[count[1]] = NaN
        else
            Prev2_RMS[count[1]] = pvalue(UnequalVarianceTTest(GP_rms[:,:RMS], G2_rms[:,:RMS]))
        end
        if isempty(GP_rms) || isempty(G3_rms) || (size(GP_rms)[1]==1 && size(G3_rms)[1]==1)
            Prev3_RMS[count[1]] = NaN
        else
            Prev3_RMS[count[1]] = pvalue(UnequalVarianceTTest(GP_rms[:,:RMS], G3_rms[:,:RMS]))
        end 
        if isempty(GP_rms) || isempty(G4_rms) || (size(GP_rms)[1]==1 && size(G4_rms)[1]==1)
            Prev4_RMS[count[1]] = NaN
        else
            Prev4_RMS[count[1]] = pvalue(UnequalVarianceTTest(GP_rms[:,:RMS], G4_rms[:,:RMS]))
        end
        #Cohen's D
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

    end # comparison of which groups

    Prev1_RMS,Prev2_RMS,Prev3_RMS,Prev4_RMS  = Prev1_RMS[Prev1_RMS .!= -1], Prev2_RMS[Prev2_RMS .!= -1], Prev3_RMS[Prev3_RMS .!= -1], Prev4_RMS[Prev4_RMS .!= -1]

    Prev1_CDRMS,Prev2_CDRMS,Prev3_CDRMS,Prev4_CDRMS  = Prev1_CDRMS[Prev1_CDRMS .!= -1], Prev2_CDRMS[Prev2_CDRMS .!= -1], Prev3_CDRMS[Prev3_CDRMS .!= -1], Prev4_CDRMS[Prev4_CDRMS .!= -1]

    BetweenGroup = DataFrame(Group=Group, Prev1_RMS=Prev1_RMS, Prev2_RMS=Prev2_RMS, Prev3_RMS=Prev3_RMS, Prev4_RMS=Prev4_RMS,Prev1_CDRMS=Prev1_CDRMS, Prev2_CDRMS=Prev2_CDRMS, Prev3_CDRMS=Prev3_CDRMS, Prev4_CDRMS=Prev4_CDRMS)

    foldername = "RelresStats"
    if !isdir(joinpath(data,foldername))
        mkdir(joinpath(data,foldername))
    end
    title = whichpeak * "_WithinGroups_" * whichstim * "_" * stimtype * "_" * trialtype
    name = joinpath(data,foldername,title) * ".csv"
    CSV.write(name, BetweenGroup)
end