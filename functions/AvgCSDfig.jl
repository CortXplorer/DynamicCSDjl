# AvgCSDfig

using Statistics
using JLD2, FileIO
# get to correct directory and label it home
func    = @__DIR__ # this script is in the functions folder already
home    = func[1:end-10]
figs    = joinpath(home,"figs")
datap   = joinpath(home,"Data")

if isdir(joinpath(figs,"AvgCSD"))
else
    mkdir(joinpath(figs,"AvgCSD"))
end

# plot function -- move to after when making this whole thing a func
function plot_avgcsd(avgCSDin,figs,Name, Cond)
    avgcsd = plot(layout = (size(avgCSDin,3))) 
    
    if ndims(avgCSDin) > 3
        avgCSDin = avgCSDin[:,:,:]
    end

    for icsd = 1:size(avgCSDin,3)
        heatmap!(
            avgcsd[icsd], # call the figure that was opened
            avgCSDin[:,:,icsd], # plot the map
            size    = (1480,1000),
            c       = :jet, # requires Colors.jl
            clims   = (-0.0003,0.0003),
            yflip   = true,
            title   = Name*Cond,
            colorbar= false)
    end
    name = joinpath(figs,"AvgCSD",Name) * "_" * Cond * "_CSD.pdf"
    savefig(avgcsd, name)
end

# best would be to have this information sent into this dynamically but that may be more possible with this whole script turned into a function call
chanlength  = 21   # max number of channels
timeaxis    = 1400 # actual time usually 1377ms
daysrecord  = 7    # days of recording
subjects    = 4    # number of subjects
frequencies = 5    # number of stimuli in trials
# let's say that we would only be working on one group at a time ever
Group       = "CIC"

# current conditions: "preCL"

# preallocate an full group data matrix
fullGroup = Array{Union{Missing, Float64}}(missing,chanlength,timeaxis,frequencies,daysrecord*subjects) 

# find all files belonging to animals of this group in the data folder
AnFiles = readdir("Data")[contains.(readdir("Data"),Group)]
theserows = 0 # indexing variable for placement in data out matrix

for iAn = AnFiles # iAn="CIC01_Data.jld2"

    filename = joinpath(datap,iAn) 
    Datout = load(filename)[iAn[1:end-5]]
    curAn = iAn[1:5]
    
    # preallocate an average picture for the animal for each condition
    fullAnimal = Array{Union{Missing, Float64}}(missing,chanlength,timeaxis,frequencies,daysrecord)

    # cycle through the days of recording
    for iday = 1:daysrecord # iday=1
        # cycle through the number of stimuli in the trial
        for ifreq = 1:numfreq # ifreq=1

            curDay  = collect(keys(Datout["preCL"]))[iday]
            thisCSD = Datout["preCL"][curDay][1]["AvgCSD"][ifreq]
            fullAnimal[1:size(thisCSD,1),1:size(thisCSD,2),ifreq,iday] = thisCSD 
            if iday == daysrecord
                fullAnimal = fullAnimal[1:size(thisCSD,1),1:size(thisCSD,2),1:frequencies,1:daysrecord] 
            end

        end # frequencies

    end # day of recording
    
    # average over days and drop the singleton dimension at the end
    avgAn_preCL = mean(fullAnimal,dims=4)[:,:,:] 
    # plot all current stimuli
    plot_avgcsd(avgAn_preCL,figs,curAn, "preCL")

    fullGroup[1:size(fullAnimal,1),1:size(fullAnimal,2),1:numfreq,(1:daysrecord).+theserows] = fullAnimal 
    theserows = theserows + daysrecord 
end # animal / entry in data folder

# cut time length to 1377 ms, the normal amount of time out
fullGroup = fullGroup[:,1:1377,:]
# find the largest size that shoots back non-missing output
while ismissing(mean(fullGroup[:,1,1]))
    fullGroup = fullGroup[1:size(fullGroup,1)-1,:,:]
end

avgGr_preCL = mean(fullGroup,dims=3)
plot_avgcsd(avgGr_preCL,figs,iGr,"preCL")




