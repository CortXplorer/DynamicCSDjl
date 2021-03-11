function visSpectrum(signal, sr = 1000)

    t   = [1:Int(1000/sr):length(signal)...] #time accounting for sr
    dt  = t[2] - t[1]   # sampling interval
    fNQ = 1 / dt / 2    # Nyquest is half sampling rate (fraction: 0.5)
    N   = length(t)     # number of data points 
    T   = length(layIV) # the total duration of data (T = N here)

    x  = hanning(N) .* layIV[:] # multiply data by hanning taper
    xf = fft(x) # compute Fourier transform
    Sxx = (2*dt)^(2/T) * (xf.*conj(xf)) # compute the spectrum
    Sxx = real(Sxx) # ignore complex components
    df = 1/T # define frequency resolution 
    faxis = [0:df:fNQ+df...] # construct freq axis 
    specplot = plot(
        faxis[1:139], 
        10 * log10.(Sxx[1:139]),
        title  = "Spectrum of CSD FFT",
        ylabel = "Power [mV^2/Hz]",
        xlabel = "Frequency [Hz]", 
    )

    return specplot

end