function ifindpeaks(ROI,flip=0)
    
    # layer traces come in with negative peaks
    if flip == 1
        ROI = ROI*-1
    end

    
    if isempty(argmaxima(ROI)) # if ~full row of NaNs then just move on
        peaklat = NaN
        peakamp = NaN
    else
        # no NaNs for peak detection though: turn them to 0
        ROI[isnan.(ROI)] .= 0
        # find the location of the highest peak prominance and take feautres
        pkprom   = findmax(peakproms(argmaxima(ROI), ROI)[2])[2]
        if peakproms(argmaxima(ROI), ROI)[2][pkprom] >= 0.0003
            peaklat = argmaxima(ROI)[pkprom]
            peakamp = ROI[peaklat]
        else # if the peak prominance is not at least 0.0003 (arbitrary, very low thresh)
            peaklat = NaN
            peakamp = NaN
        end
    end

    return peaklat, peakamp

end