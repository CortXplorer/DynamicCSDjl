function PowerPermBetween(figs,WTof1,WTof2,Group1,Group2,curMeas,curStim,curLay,takepic=1)
    curComp = Group1 * "v" * Group2
    animals1,_,_,_,_,_,_ = callGroup(Group1); 
    grpsize1 = length(WTof1)
    thisGrp1 = Array{Float64}(undef, 40, 1377, grpsize1); 
    for iwt = 1:grpsize1
        thisGrp1[:,:,iwt] = WTof1[animals1[iwt]][curMeas[begin:end-2]][curMeas][curStim][curLay][1][:,1:1377]
    end

    animals2,_,_,_,_,_,_ = callGroup(Group2); 
    grpsize2 = length(WTof2)
    thisGrp2 = Array{Float64}(undef, 40, 1377, grpsize2) 
    for iwt = 1:grpsize2
        thisGrp2[:,:,iwt] = WTof2[animals2[iwt]][curMeas[begin:end-2]][curMeas][curStim][curLay][1][:,1:1377]
    end

    ## Degrees of Freedom and t Threshold
    dfree = grpsize1 + grpsize2 - 2
    Tchart = [12.71,4.303,3.182,2.776,2.571,2.447,2.365,2.306,2.262,2.228,2.201,2.179,2.160,2.145,2.131,2.120,2.110,2.101,2.093,2.086]
    Tthresh = Tchart[dfree]

    ## Permutation step 1 - observed differences / step 2 - t test and clustermass
    grp1_mean, grp2_mean, difmeans, clusters, clustermass, alpha, betal, betah, gammal, gammah = makettest_cluster(thisGrp1,thisGrp2,Tthresh,1)

    ## Permutation step 3 - do the permute
    nperms = 1000
    # preallocate
    perm_clustermass = Array{Float64}(undef, nperms)
    perm_alpha = Array{Float64}(undef, nperms)
    perm_betal = Array{Float64}(undef, nperms)
    perm_betah = Array{Float64}(undef, nperms)
    perm_gammal = Array{Float64}(undef, nperms)
    perm_gammah = Array{Float64}(undef, nperms)
    # concatonate groups together and initialize new containers
    thisAll = cat(thisGrp2,thisGrp1,dims=3)
    newGrp1 = Array{Float64}(undef, size(thisAll,1), size(thisAll,2), grpsize1);
    newGrp2 = Array{Float64}(undef, size(thisAll,1), size(thisAll,2), grpsize2);

    for iperm = 1:nperms
        # generate random permutation list
        permorder = randperm(size(thisAll,3))
        for iord = 1:length(permorder)
            if iord <= grpsize1
                newGrp1[:,:,iord] = thisAll[:,:,permorder[iord]]
            else
                newGrp2[:,:,iord-grpsize1] = thisAll[:,:,permorder[iord]]
            end
        end
        perm_clustermass[iperm], perm_alpha[iperm], perm_betal[iperm], perm_betah[iperm], perm_gammal[iperm], perm_gammah[iperm] = makettest_cluster(newGrp1,newGrp2,Tthresh,0)
    end

    ## Check significance
    sig_mass = length(findall(perm_clustermass .> clustermass))
    pVal = sig_mass / nperms
    permMean = mean(perm_clustermass)
    permSTD  = std(perm_clustermass)

    sig_massA = length(findall(perm_alpha .> alpha))
    ApVal = sig_massA / nperms
    permAmean = mean(perm_alpha)
    permASTD  = std(perm_alpha)
    sig_massBL = length(findall(perm_betal .> betal))
    BLpVal = sig_massBL / nperms
    permBLmean = mean(perm_betal)
    permBLSTD  = std(perm_betal)
    sig_massBH = length(findall(perm_betah .> betah))
    BHpVal = sig_massBH / nperms
    permBHmean = mean(perm_betah)
    permBHSTD  = std(perm_betah)
    sig_massGL = length(findall(perm_gammal .> gammal))
    GLpVal = sig_massGL / nperms
    permGLmean = mean(perm_gammal)
    permGLSTD  = std(perm_gammal)
    sig_massGH = length(findall(perm_gammah .> gammah))
    GHpVal = sig_massGH / nperms
    permGHmean = mean(perm_gammah)
    permGHSTD  = std(perm_gammah)

    obsvalues = [clustermass, alpha, betal, betah, gammal, gammah]
    ticknames = ["All", "alpha", "beta low", "beta high", "gamma low", "gamma high"]

    if takepic == 1
        # difference figs
        # NOTE, when loop is made, this will have to be updated to take in which layer/stim etc and name the figs dynamically
        differenceplots(figs, params, grp1_mean, grp2_mean, difmeans, clusters, Group1, Group2, curComp, curMeas,curStim,curLay,"power")

        PermPlot = violin(
            [perm_clustermass, perm_alpha, perm_betal, perm_betah, perm_gammal, perm_gammah],
            lab = [pVal ApVal BLpVal BHpVal GLpVal GHpVal],
            formatter =x->round(Int, x),
            xticks = (1:6,ticknames),
        );
        plot!(obsvalues,markershape=:circle,markersize=10);

        name = joinpath(figs,"Spectral",curComp) * "_" * curMeas * "_" * curStim * "_" * curLay * "_Power_PermViolins.png"
        savefig(PermPlot, name)
    end

    PermOut = DataFrame(Comparison = repeat([curComp],6), Condition = repeat([curMeas],6), StimFreq = repeat([curStim],6), Layer = repeat([curLay],6), Osciband = ticknames, PValues = [pVal,ApVal,BLpVal,BHpVal,GLpVal,GHpVal], Means = [permMean,permAmean,permBLmean,permBHmean,permGLmean,permGHmean], STDs = [permSTD,permASTD,permBLSTD,permBHSTD,permGLSTD,permGHSTD])

    return PermOut

end