### RAnova takes extracted features from avrec and layer trace data in the form of csv (output from matlab) and produces Anova output in txt format

## input:   home/Data/AVRECPeak**.csv
## output:  home/Data/AvrecPeakStats/ANOVASummary.txt

# # Set working directory and create necessary file paths: # #
setwd("D:/DynamicCSDjl") ## Change to your! (don't know how to do this dynamically in R)
home = getwd()
datapath = file.path(home,"Data")

sink(file=file.path(datapath,"AvrecPeakStats","ANOVASummary.txt"), type=c("output")) #save all output to txt file!
TAorST   = list("TA", "ST")
stimtype = list("CL", "AM")
layer    = list("All", "I_II", "IV", "V", "VI")
comptype = list("KIT vs KIC", "KIT vs KIV", "KIC vs KIV")

for (iRun in 1:length(TAorST)) {
  
  if (TAorST[iRun] == "TA") {
    print(paste0("============================== Trial Average =============================="))
  } else {
    print(paste0("============================== Single Trial  =============================="))
  } # one more if statement when loading the dataset in in the next step
    
  for (iStim in 1:length(stimtype)) {
    
    print(paste0("=========================== ANOVAs FOR STIM TYPE ", stimtype[iStim], " ==========================="))
    
    # load full dataset in of appropriate stim type
    if (TAorST[iRun] == "TA") {
      peakdata = read.csv(file = file.path(datapath,paste0("AVRECPeak",stimtype[1],".csv")))
    } else {
      peakdata = read.csv(file = file.path(datapath,paste0("AVRECPeak",stimtype[1],"ST.csv")))
    }
    # pull out stim frequencies present to loop through
    stimfreq = unique(peakdata[c("ClickFreq")])
    
    for (iLay in 1:length(layer)) {
      
      print(paste0("=========================== LAYER ", layer[iLay], " ==========================="))
      LDay = subset(peakdata, Layer == layer[iLay])
      
      for (iFreq in 1:nrow(stimfreq)) {
        
        print(paste0("===================== STIMULUS FREQUENCY ", stimfreq[iFreq,1], " Hz ====================="))
        
        # cut down to one stim frequency 
        FDat = subset(LDay, ClickFreq == stimfreq[iFreq,1])
        
        # loop through each response per stimuli 
        for (iOrd in 1:stimfreq[iFreq,1]) {
          
          print(paste0("=============== RESPONSE ", iOrd, " OF ", stimfreq[iFreq,1], " ==============="))
          OFDat = subset(FDat, OrderofClick == iOrd)
          
          for (iCmp in 1:length(comptype)) {
            
            if (comptype[iCmp] == "KIT vs KIC") {
              subDat = subset(OFDat, !(Group == "KIV"))
            } else if (comptype[iCmp] == "KIT vs KIV") {
              subDat = subset(OFDat, !(Group == "KIC"))
            } else if (comptype[iCmp] == "KIC vs KIV") {
              subDat = subset(OFDat, !(Group == "KIT"))
            }
            
            print(paste0("========= PEAK AMPLITUDE ", comptype[iCmp], " =========")) 
            res.aov3 = aov(PeakAmp ~ Group * Measurement, data = subDat)
            print(summary(res.aov3))
            
            print(paste0("========= PEAK LATITUDE ", comptype[iCmp], " =========")) 
            res.aov3 = aov(PeakLat ~ Group * Measurement, data = subDat)
            print(summary(res.aov3))
            
            print(paste0("========= ROOT MEAN SQUARE ", comptype[iCmp], " =========")) 
            res.aov3 = aov(RMS ~ Group * Measurement, data = subDat)
            print(summary(res.aov3))
          @info "lel"
          } # which comparison
        } # order of click
      } # which stim frequency
    } # which layer
  } # which stim type, CL or AM
} # Trial average or single trial
sink() #restore output to console and finish using txt file
