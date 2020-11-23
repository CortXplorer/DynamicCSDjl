
function CWTanalysis(ROI,params,curAn="KIC02",curCond="preCL_1",curLay="IV")

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

    baseidx  = ([params.startTime+1 0]) .+ params.startTime*-1 # time of baseline 

    for fi=1:num_frex

        wavelet = sqrt(1/(s[fi]*sqrt(pi))) * exp.(2*1im*pi*frex[fi].*time) .* exp.(-time.^2 ./ (2*(s[fi]^2)))

        datconv = same_conv(ROI[:], wavelet)
        datconvshaped = reshape(datconv,size(ROI)[2],size(ROI)[3])
        # Average power over trials (this code performs baseline transform)
        temppower = mean(abs.(datconvshaped).^2, dims=2) # POWER 
        temppower = 10*log10.(temppower ./ mean(temppower[baseidx[1]:baseidx[2]]))
        datpower[fi,:] = temppower

        tempphsco = abs.(mean((datconvshaped ./ abs.(datconvshaped)), dims=2)) # PHASE COHERENCE
        datphsco[fi,:] = tempphsco
    end

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
        title="Power of " * curAn * " " * curCond * " layer " * params.layers[iLay]
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
        title="Phase Co of " * curAn * " " * curCond * " layer " * params.layers[iLay]
    );

    full_plot = plot(power_plot, phsco_plot, titlefontsize = 10, size=(900,400))

    foldername = "Spectral"
    if !isdir(joinpath(figs,foldername))
        mkdir(joinpath(figs,foldername))
    end

    name = joinpath(figs,foldername,curAn) * "_" * curCond * "_" * params.layers[iLay] * "_ScaloPower.pdf"
    savefig(full_plot, name)

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