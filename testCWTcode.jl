
#function CWTanalysis(ROI,params)
ROI
time = -1:(1/params.sampleRate):1
ROItime = [-200:size(ROI)[2]-200-1...]

min_freq = params.frequencyLimits[1]
max_freq = params.frequencyLimits[2]
num_frex = params.timeBandWidth

n_data          = size(ROI)[2]*size(ROI)[3] #data points sampled * trials
# NOTE: for multiple trial data, this will have to be examined, the same_conv function below required that I cut ROI to [1:1377] for a reasonably sized output to variable datconv
n_wavelet       = length(time)
n_convolution   = n_wavelet+n_data-1
n_conv_pow2     = nextpow(2,n_convolution)
half_wavelet    = Int((n_wavelet-1)/2)


# exp10. is the same as logspace
frex    = exp10.(range(log10(min_freq), log10(max_freq), length=num_frex))
s       = exp10.(range(log10(3),log10(10), length=num_frex)) ./ (2 * pi * frex) # scaling factor

# Initialize
datpower = Array{Float64}(undef, num_frex, size(ROI)[2])

baseidx = ([params.startTime+1 0]) .+ params.startTime*-1 # time of baseline 

for fi=1:num_frex

    wavelet = sqrt(1/(s[fi]*sqrt(pi))) * exp.(2*1im*pi*frex[fi].*time) .* exp.(-time.^2 ./ (2*(s[fi]^2)))

    datconv = same_conv(ROI[:], wavelet)
    # Average power over trials (this code performs baseline transform)
    temppower = mean(abs.(reshape(datconv,size(ROI)[2],size(ROI)[3])).^2, dims=2) 
    # NOTE: this reshaping is currently unnecessary on a single vector of data but I'm leaving it for when I re-examine for multiple trials
    temppower = 10*log10.(temppower ./ mean(temppower[baseidx[1]:baseidx[2]]))
    datpower[fi,:] = temppower

end

log_freq = contourf(
    ROItime,
    frex,
    datpower,
    color=:viridis,
    levels=40,
    # clim=(-3,3),
    xlims=(-200,1170),
    yaxis=:log,
    formatter =x->round(Int, x),
    ytick=exp10.(range(log10(min_freq),log10(max_freq),length=6)),
    title="Logarithmic frequency scaling"
);


plot(log_freq, titlefontsize = 10, size=(400,400))



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