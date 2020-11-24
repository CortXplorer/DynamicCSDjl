function CWT_Loop(figs, animalList, CondList, CLList, params, anipar, takepic)
    # Loop through animals in Group
    Animal = Dict()
    for iAn = 1:length(animalList)

        # loop through animals iAn
        curAn   = animalList[iAn]
        anDat   = matread(joinpath(data,(curAn * "_Data.mat")));
        anMeas  = anDat["Data"]["measurement"]

        # Loop through the condition list
        Condition = Dict()
        for iCond = 1:length(CLList)
            # loop through condition list iCL
            MeasList = CondList[CLList[iCond]][iAn]

            # Loop through measurements
            Measurement = Dict()
            for iMeas = 1:length(MeasList)

                curMeas = CondList[CLList[iCond]][iAn][iMeas]
                curCond = CLList[iCond] * "_" * string(iMeas)
                runthis = [curAn*"_"*curMeas] # generate full measurement name

                thisind = findall(anDat["Data"]["measurement"] .== runthis)[1][2] #[1][2] for extracting cartesian index
                println("Analyzing $runthis")

                # Loop through stimulation frequencies
                StimFreq = Dict()
                for iSti = 1:length(params.stimList)

                    curCSD  = anDat["Data"]["SglTrl_CSD"][thisind][iSti];
                    curStim = params.stimList[iSti]

                    # Loop through layers
                    Layers = Dict()
                    for iLay = 1:length(params.layers)

                        curLay = params.layers[iLay]
                        if curLay == "I_II"
                            curChan = anipar.LIIList[iAn]
                        elseif curLay == "IV"
                            curChan = anipar.LIVList[iAn]
                        elseif curLay == "V"
                            curChan = anipar.LVList[iAn]
                        elseif curLay == "VI"
                            curChan = anipar.LVIList[iAn]
                        end

                        if length(curChan) > 3 # take center 3 if greater than 3
                            centerChan = Int(ceil(length(curChan)/2))
                            curChan = curChan[centerChan-1:centerChan+1]
                        end

                        # Here is the region of interest accross trials
                        ROI = mean(curCSD[curChan,:,:],dims=1)

                        Datpower, Datphsco = CWTanalysis(figs,ROI,params,curAn,curCond,curLay,curStim)
                        Layers[curLay] = Datpower, Datphsco;
                    end # Layers
                    StimFreq[curStim] = Layers
                end # Stim Frequencies
                Measurement[curCond] = StimFreq
            end # Condition/Measurement
            Condition[CLList[iCond]] = Measurement
        end # Condition List
        Animal[curAn] = Condition
    end # Animal
    return Animal
end

function CWTanalysis(figs,ROI,params,curAn="KIC02",curCond="preCL_1",curLay="IV",curStim="2Hz")
# Adapted from Asim Hassan Dar's code following the Mike X Cohen analyzing neural time series 23.11.20

    time     = -1:(1/params.sampleRate):1
    ROItime  = [-200:size(ROI)[2]-200-1...]

    min_freq = params.frequencyLimits[1]
    max_freq = params.frequencyLimits[2]
    num_frex = params.timeBandWidth

    n_data   = size(ROI)[2]*size(ROI)[3] #data points sampled * trials

    # exp10. is the same as logspace
    frex     = exp10.(range(log10(min_freq), log10(max_freq), length=num_frex))
    s        = exp10.(range(log10(3),log10(10), length=num_frex)) ./ (2 * pi * frex) # scaling factor

    # Initialize
    datpower = Array{Float64}(undef, num_frex, size(ROI)[2])
    datphsco = Array{Float64}(undef, num_frex, size(ROI)[2])

    baseidx  = ([params.startTime+1 params.startTime+100]) .+ params.startTime*-1 # time of baseline (100 away from 0 to avoid most temporal smoothing)

    for fi=1:num_frex
        wavelet = sqrt(1/(s[fi]*sqrt(pi))) * exp.(2*1im*pi*frex[fi].*time) .* exp.(-time.^2 ./ (2*(s[fi]^2)))

        datconv = same_conv(ROI[:], wavelet)
        datconvshaped = reshape(datconv,size(ROI)[2],size(ROI)[3])
        # Average power over trials
        temppower = mean(abs.(datconvshaped).^2, dims=2) # POWER 
        # decibel change baseline normalization
        temppower = 10*log10.(temppower ./ mean(temppower[baseidx[1]:baseidx[2]]))
        datpower[fi,:] = temppower

        tempphsco = abs.(mean((datconvshaped ./ abs.(datconvshaped)), dims=2)) # PHASE COHERENCE
        datphsco[fi,:] = tempphsco
    end

    if takepic == 1
        power_plot = heatmap(
            ROItime,
            frex,
            datpower,
            levels=40,
            clim=(-10,10),
            xlims=(-200,1170),
            yaxis=:log,
            formatter =x->round(Int, x),
            ytick=exp10.(range(log10(min_freq),log10(max_freq),length=10)),
            title="Power of " * curAn * " " * curCond * " layer " * curLay * " " * curStim
        );

        phsco_plot = heatmap(
            ROItime,
            frex,
            datphsco,
            levels=40,
            clim=(0,1),
            xlims=(-200,1170),
            yaxis=:log,
            formatter =x->round(Int, x),
            ytick=exp10.(range(log10(min_freq),log10(max_freq),length=10)),
            title="Phase Co of " * curAn * " " * curCond * " layer " * curLay * " " * curStim
        );

        full_plot = plot(power_plot, phsco_plot, titlefontsize = 10, size=(900,400))

        foldername = "Spectral"
        if !isdir(joinpath(figs,foldername))
            mkdir(joinpath(figs,foldername))
        end

        name = joinpath(figs,foldername,curAn) * "_" * curCond * "_" * curLay * "_" * curStim * "_Scalograms.pdf"
        savefig(full_plot, name)
    end
    return datpower, datphsco
end

function same_conv(signal, kernel)

    idx_start = ((length(kernel) - 1) ÷ 2) + 1
    idx_end = (length(kernel) - 1) ÷ 2
    same_size_convolution = conv(signal, kernel)[idx_start:end-idx_end]

    if iseven(length(kernel))
        @error "Kernel was even, the last value was truncated."
        return same_size_convolution[1:length(signal)]

    end
    return same_size_convolution
end