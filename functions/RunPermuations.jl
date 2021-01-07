function RunPermuations(figs,spect,KIT_WT,KIC_WT,KIV_WT,MeasList,params)

    println("Power Permutation Between Groups")
    # ADD: Cohen's D output column
    PermOut = []
    for iMeas = 1:length(MeasList)
        curMeas = MeasList[iMeas]
        for iSti = 1:length(params.stimList)
            curStim = params.stimList[iSti]
            for iLay = 1:length(params.layers)
                curLay = params.layers[iLay]
                # FIX: ERROR: KeyError: key "AM_4" not found (have it check for this key and skip if it's not there)
                println("Power Measurement $curMeas, Stimulus $curStim, Layer $iLay")
                Permout1 = PowerPermBetween(figs,KIT_WT,KIC_WT,"KIT","KIC",curMeas,curStim,curLay,takepic)
                Permout2 = PowerPermBetween(figs,KIT_WT,KIV_WT,"KIT","KIV",curMeas,curStim,curLay,takepic)
                Permout3 = PowerPermBetween(figs,KIC_WT,KIV_WT,"KIC","KIV",curMeas,curStim,curLay,takepic)

                if !isempty(PermOut)
                    PermOut = vcat(PermOut,Permout1,Permout2,Permout3)
                else
                    PermOut = vcat(Permout1,Permout2,Permout3)
                end
            end # layer
        end # Stim frequency
    end # Measurement

    CSV.write(joinpath(spect,"PowerBetween.csv"),PermOut)
    # PowerBetween = CSV.File(joinpath(spect,"PowerBetween.csv")) |> DataFrame

    println("Phase Coherence Permutation Between Groups")
    PermOut = []
    for iMeas = 1:length(MeasList)
        curMeas = MeasList[iMeas]
        for iSti = 1:length(params.stimList)
            curStim = params.stimList[iSti]
            for iLay = 1:length(params.layers)
                curLay = params.layers[iLay]
                # FIX: ERROR: KeyError: key "AM_4" not found (have it check for this key and skip if it's not there)
                println("Phase Measurement $curMeas, Stimulus $curStim, Layer $iLay")
                Permout1 = PhasePermBetween(figs,KIT_WT,KIC_WT,"KIT","KIC",curMeas,curStim,curLay,takepic)
                Permout2 = PhasePermBetween(figs,KIT_WT,KIV_WT,"KIT","KIV",curMeas,curStim,curLay,takepic)
                Permout3 = PhasePermBetween(figs,KIC_WT,KIV_WT,"KIC","KIV",curMeas,curStim,curLay,takepic)

                if !isempty(PermOut)
                    PermOut = vcat(PermOut,Permout1,Permout2,Permout3)
                else
                    PermOut = vcat(Permout1,Permout2,Permout3)
                end
            end # layer
        end # Stim frequency
    end # Measurement

    CSV.write(joinpath(spect,"PhaseBetween.csv"),PermOut)
    # PhaseBetween = CSV.File(joinpath(spect,"PhaseBetween.csv")) |> DataFrame

end