## Set working directory and create necessary file paths
setwd("D:/DynamicCSDjl")
home = getwd()
datapath = file.path(home,"Data")

sink(file=file.path(datapath,"AvrecPeakStats","ANOVASummary.txt"), type=c("output")) #save all output to txt file!
stimtype = list("CL","AM")
comptype = list("KIT vs KIC", "KIT vs KIV", "KIC vs KIV") #,"KIT vs KIV", "KIC vs KIV"

for (iStim in 1:length(stimtype)) {
  print(paste0("=========================== ANOVAs FOR STIM TYPE ", stimtype[iStim], " ==========================="))
  # load full dataset in of appropriate stim type
  peakdata = read.csv(file = file.path(datapath,paste0("AVRECPeak",stimtype[1],".csv")))
  # pull out stim frequencies present to loop through
  stimfreq = unique(peakdata[c("ClickFreq")])
  
  for (iFreq in 1:nrow(stimfreq)) {
    
    print(paste0("===================== STIMULUS FREQUENCY ", stimfreq[iFreq,1], " Hz ====================="))
    
    # cut down to one stim frequency 
    FDat = subset(peakdata, ClickFreq == stimfreq[iFreq,1])
    
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
      
      }
    } # order of click
  } # which stim frequency
} # which stim type, CL or AM

sink() #restore output to console and finish using txt file
