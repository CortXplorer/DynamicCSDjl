function AvgCSDfig(figs,datap,chanlength,timeaxis,daysrecord,subjects,frequencies,Group,Condition)

    if isdir(joinpath(figs,"AvgCSD"))
    else
        mkdir(joinpath(figs,"AvgCSD"))
    end

    # preallocate an full group data matrix
    fullGroup = Array{Union{Missing, Float64}}(missing,chanlength,timeaxis,frequencies,daysrecord*subjects) 

    # find all files belonging to animals of this group in the data folder
    AnFiles = readdir("Data")[contains.(readdir("Data"),Group)]
    theserows = 0 # indexing variable for placement in data out matrix

    for iAn = AnFiles # iAn="CIC01_Data.jld2"
        
        filename = joinpath(datap,iAn) 
        Datout = load(filename)[iAn[1:end-5]] # this takes the longest
        curAn = iAn[1:5]
        
        # preallocate an average picture for the animal for each condition
        fullAnimal = Array{Union{Missing, Float64}}(missing,chanlength,timeaxis,frequencies,daysrecord)
        
        # cycle through the days of recording
        for iday = 1:daysrecord # iday=1
            # cycle through the number of stimuli in the trial
            for ifreq = 1:frequencies # ifreq=1
                
                curDay  = collect(keys(Datout[Condition]))[iday]
                thisCSD = Datout[Condition][curDay][1]["AvgCSD"][ifreq]
                fullAnimal[1:size(thisCSD,1),1:size(thisCSD,2),ifreq,iday] = thisCSD 
                if iday == daysrecord
                    fullAnimal = fullAnimal[1:size(thisCSD,1),1:size(thisCSD,2),1:frequencies,1:daysrecord] 
                end
                
            end # frequencies
            
        end # day of recording
        
        # average over days and drop the singleton dimension at the end
        avgAn_preCL = mean(fullAnimal,dims=4)[:,:,:] 
        # plot all current stimuli
        plot_avgcsd(avgAn_preCL,figs,curAn,Condition) # plot func below
        
        fullGroup[1:size(fullAnimal,1),1:size(fullAnimal,2),1:size(fullAnimal,3),(1:daysrecord).+theserows] = fullAnimal 
        theserows = theserows + daysrecord 
    end # animal / entry in data folder

    # cut time length to 1377 ms, the normal amount of time out
    fullGroup = fullGroup[:,1:1377,:,:]
    # find the largest size that shoots back non-missing output
    while ismissing(mean(fullGroup[:,1,1,1]))
        fullGroup = fullGroup[1:size(fullGroup,1)-1,:,:,:]
    end

    avgGr_preCL = mean(fullGroup,dims=4)[:,:,:]
    plot_avgcsd(avgGr_preCL,figs,Group,Condition) # plot func below
end


# plot function -- move to after when making this whole thing a func
function plot_avgcsd(avgCSDin,figs,Name, Cond)
    avgcsd = plot(layout = (size(avgCSDin,3))) 
    
    if ndims(avgCSDin) > 3
        avgCSDin = avgCSDin[:,:,:]
    end

    for icsd = 1:size(avgCSDin,3)
        if icsd < size(avgCSDin,3)
            heatmap!(
                avgcsd[icsd], # call the figure that was opened
                avgCSDin[:,:,icsd], # plot the map
                size    = (1480,1000),
                c       = :jet, # requires Colors.jl
                clims   = (-0.0003,0.0003),
                yflip   = true,
                title   = Name*Cond,
                colorbar= false)
        else 
            heatmap!(
                avgcsd[icsd], # call the figure that was opened
                avgCSDin[:,:,icsd], # plot the map
                size    = (1480,1000),
                c       = :jet, # requires Colors.jl
                clims   = (-0.0003,0.0003),
                yflip   = true,
                title   = Name*Cond)
        end
    end

    name = joinpath(figs,"AvgCSD",Name) * "_" * Cond * "_CSD.pdf"
    savefig(avgcsd, name)
end

