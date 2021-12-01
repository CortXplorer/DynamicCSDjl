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
    avgcsd = plot(layout = (1,1)) # lul
    
    if ndims(avgCSDin) > 2
        avgCSDin = avgCSDin[:,:]
    end
    heatmap!(
        avgcsd, # call the figure that was opened
        avgCSDin, # plot the map
        size    = (1480,1000),
        c       = :jet, # requires Colors.jl
        clims   = (-0.0003,0.0003),
        yflip   = true,
        title   = Name*Cond)
    
    name = joinpath(figs,"AvgCSD",Name) * "_" * Cond * "_CSD.pdf"
    savefig(avgcsd, name)
end

chanlength = 21   # max number of channels
timeaxis   = 1400 # actual time usually 1377ms
daysrecord = 7    # days of recording
subjects   = 4    # number of subjects

GroupList = ["CIC"]

# current conditions: "preCL"

for iGr = GroupList # iGr="CIC"

    # preallocate an average picture for the group for each condition
    Group_preCL = Array{Union{Missing, Float64}}(missing,chanlength,timeaxis,daysrecord*subjects)

    AnFiles = readdir("Data")[contains.(readdir("Data"),iGr)]
    theserows = 0
    for iAn = AnFiles # iAn="CIC01_Data.jld2"
        filename = joinpath(datap,iAn) 
        Datout = load(filename)[iAn[1:end-5]]
        curAn = iAn[1:5]

        # preallocate an average picture for the animal for each condition
        CurAn_preCL = Array{Union{Missing, Float64}}(missing,chanlength,timeaxis,daysrecord)
        for iday = 1:daysrecord # iday=1
            curDay  = collect(keys(Datout["preCL"]))[iday]
            thisCSD = Datout["preCL"][curDay][1]["AvgCSD"][2]
            CurAn_preCL[1:size(thisCSD,1),1:size(thisCSD,2),iday] = thisCSD # the final 2 indicates 5Hz stimulus!
            if iday == daysrecord
                CurAn_preCL = CurAn_preCL[1:size(thisCSD,1),1:size(thisCSD,2),1:daysrecord] 
            end
        end # condition
        avgAn_preCL = mean(CurAn_preCL,dims=3)
        plot_avgcsd(avgAn_preCL,figs,curAn, "preCL")

        Group_preCL[1:size(CurAn_preCL,1),1:size(CurAn_preCL,2),(1:daysrecord).+theserows] = CurAn_preCL 
        theserows = theserows + daysrecord 
    end # animal / entry in data folder

    # cut time length to 1377 ms, the normal amount of time out
    Group_preCL = Group_preCL[:,1:1377,:]
    # find the largest size that shoots back non-missing output
    while ismissing(mean(Group_preCL[:,1,1]))
        Group_preCL = Group_preCL[1:size(Group_preCL,1)-1,:,:]
    end

    avgGr_preCL = mean(Group_preCL,dims=3)
    plot_avgcsd(avgGr_preCL,figs,iGr,"preCL")
end # group



