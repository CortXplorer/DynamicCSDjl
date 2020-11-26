function makettest_cluster(Group1,Group2,Tthresh,isobs=1)
    # take the group data, run a t-test, convert to clusters, pull out clustermass number 
    grpsize1 = size(Group1,3)
    grpsize2 = size(Group2,3)

    grp1_mean = mean(Group1,dims=3)[:,:]
    grp1_std = std(Group1,dims=3)[:,:]

    grp2_mean = mean(Group2,dims=3)[:,:]
    grp2_std = std(Group2,dims=3)[:,:]

    ## Permutation step 2 - t test 
    ttest = (grp2_mean .- grp1_mean)./sqrt.((grp2_std.^2/grpsize2)+(grp1_std.^2/grpsize1))
    clusters = abs.(ttest)
    clusters[clusters .< Tthresh] .= 0
    # the total count of significant points:
    clustermass = length(clusters[clusters .!= 0])

    ## Ocillatory bands:
    alpha  = clusters[params.bandRanges[1],:]
    alpha  = length(alpha[alpha .!= 0])
    betal  = clusters[params.bandRanges[2],:]
    betal  = length(betal[betal .!= 0])
    betah  = clusters[params.bandRanges[3],:]
    betah  = length(betah[betah .!= 0])
    gammal = clusters[params.bandRanges[4],:]
    gammal = length(gammal[gammal .!= 0])
    gammah = clusters[params.bandRanges[5],:]
    gammah = length(gammah[gammah .!= 0])
    
    if isobs == 1
        difmeans = grp1_mean .- grp2_mean
        return grp1_mean, grp2_mean, difmeans, clusters, clustermass, alpha, betal, betah, gammal, gammah
    else
        return clustermass, alpha, betal, betah, gammal, gammah
    end
end

function makeMWutest_cluster(thisGrp1,thisGrp2,isobs=1)
    # take the group data, run a mann whitney u-test, convert to clusters, pull out clustermass number 

    grp1_mean = mean(thisGrp1,dims=3)[:,:]
    grp2_mean = mean(thisGrp2,dims=3)[:,:]

    ## Permutation step 2 - t test 
    mwutest = Array{Float64}(undef, 40, 1377);
    for i = 1:40
        for j = 1:1377
            mwutest[i,j] = pvalue(MannWhitneyUTest(thisGrp1[i,j,:], thisGrp2[i,j,:]))
        end
    end

    mwutest[mwutest .> 0.05] .= 0
    # the total count of significant points:
    clustermass = length(mwutest[mwutest .!= 0])

    ## Ocillatory bands:
    alpha  = mwutest[params.bandRanges[1],:]
    alpha  = length(alpha[alpha .!= 0])
    betal  = mwutest[params.bandRanges[2],:]
    betal  = length(betal[betal .!= 0])
    betah  = mwutest[params.bandRanges[3],:]
    betah  = length(betah[betah .!= 0])
    gammal = mwutest[params.bandRanges[4],:]
    gammal = length(gammal[gammal .!= 0])
    gammah = mwutest[params.bandRanges[5],:]
    gammah = length(gammah[gammah .!= 0])
    
    if isobs == 1
        difmeans = grp1_mean .- grp2_mean
        return grp1_mean, grp2_mean, difmeans, mwutest, clustermass, alpha, betal, betah, gammal, gammah
    else
        return clustermass, alpha, betal, betah, gammal, gammah
    end
end

function differenceplots(figs, params, grp1_mean, grp2_mean, difmeans, clusters, Group1, Group2, curComp, curMeas,curStim,curLay,typefig="power")
    if typefig == "power"
        CLIM = (-10,10)
        typecall = " Power.png"
    elseif typefig == "phase"
        CLIM = (-0.5,1.5)
        typecall = "_Phase.png"
    end
    frex = exp10.(range(log10(params.frequencyLimits[1]),log10(params.frequencyLimits[2]), length=params.timeBandWidth))
    figTime = [-200:1176...]
    Gr1fig = heatmap(
        figTime, 
        frex,
        grp1_mean,
        levels=40,
        clim=CLIM,
        xlims=(-200,1176),
        yaxis=:log,
        formatter =x->round(Int, x),
        ytick=exp10.(range(log10(params.frequencyLimits[1]),log10(params.frequencyLimits[2]),length=10)),
        title= Group1 * " " * typefig
    );

    Gr2fig = heatmap(
        figTime, 
        frex,
        grp2_mean,
        levels=40,
        clim=CLIM,
        xlims=(-200,1176),
        yaxis=:log,
        formatter =x->round(Int, x),
        ytick=exp10.(range(log10(params.frequencyLimits[1]),log10(params.frequencyLimits[2]),length=10)),
        title= Group2 * " " * typefig
    );

    Diffig = heatmap(
        figTime, 
        frex,
        difmeans,
        levels=40,
        clim=CLIM,
        xlims=(-200,1176),
        yaxis=:log,
        formatter =x->round(Int, x),
        ytick=exp10.(range(log10(params.frequencyLimits[1]),log10(params.frequencyLimits[2]),length=10)),
        title= "Difference of " * typefig
    );

    Clufig = heatmap(
        figTime, 
        frex,
        clusters,
        levels=40,
        xlims=(-200,1176),
        yaxis=:log,
        formatter =x->round(Int, x),
        ytick=exp10.(range(log10(params.frequencyLimits[1]),log10(params.frequencyLimits[2]),length=10)),
        title="Clusters of significance"
    );

    full_plot = plot(Gr1fig,Gr2fig,Diffig,Clufig,size=(900,600));
    name = joinpath(figs,"Spectral",curComp) * "_" * curMeas * "_" * curStim * "_" * curLay * typecall
            savefig(full_plot, name)
end