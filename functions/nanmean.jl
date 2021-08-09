function nanmean(varin) 
    trace = zeros(size(varin,2))
    for i = 1:size(varin,2)
        trace[i] = mean(filter(!isnan, varin[:,i]))
    end
    return trace
end