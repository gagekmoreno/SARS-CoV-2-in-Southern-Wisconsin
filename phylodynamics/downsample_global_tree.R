#### Import variables ####
args<-commandArgs()
outdir<-args[6]
regionname<-args[7]

#### Load required packages and install if not already installed ####
requiredPackages = c('ape','treedater','splitstackshape','lubridate','stringr')
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}

#### randomly downsample from time bins ####
filename<-paste0("global",regionname,sep="_")
filename2<-paste0(filename,"filtered_names.txt",sep="_")
filename3<-paste0(outdir,filename2,sep="/")
global<-read.table(file = filename3, header = FALSE)
global$V2<-global$V1
globalsplit<-cSplit(global, "V2", sep = "|")
dfsplit<-globalsplit[,c(1,4)]
dfsplit$V2_3<-as.POSIXct(dfsplit$V2_3)
dfsplit$V2_3<-decimal_date(dfsplit$V2_3)
dfsplit$diff<-ceiling((dfsplit$V2_3-min(dfsplit$V2_3))*365)
dfsplit$week<-floor(dfsplit$diff/7)+1
dfsplit$week<-paste0("week",dfsplit$week)
list_df <- split(dfsplit, dfsplit$week)

weeks<-as.list(unique(dfsplit$week))
len<-length(weeks)

for (i in 1:len) {
  a<-list_df[[i]]
  b<-ifelse(nrow(a)<100,5,
            ifelse(nrow(a)<200,10,50))
  assign(paste0("Random",i),a[sample(nrow(a),b), ])
}

dfs<-paste0("Random",1:len)
randomlist = lapply(dfs, get)
time_downsampled<-bind_rows(randomlist, .id=NULL)

downsampled_names<-as.data.frame(time_downsampled$V1)

#### Write output ####
output1<-paste0("downsampled_global_names",regionname,sep="_")
output2<-paste0(output1,".txt",sep="")
output3<-paste0(outdir,output3,sep="/")
write.table(downsampled_names, file = output3, row.names=FALSE, col.names=FALSE, quote=FALSE)



