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

#### read iqtree-generated maximum-likelihood tree ####
filename<-paste0(outdir,"dates_filtered_CDS_global.tree",sep="/")
tree<-read.tree(file = filename)

#### extract tip dates ####
tips<-as.data.frame(tree$tip.label)
tips$`tree$tip.label`<-gsub('_2020-', '_#2020-', tips$`tree$tip.label`)
tips$`tree$tip.label`<-gsub('_2019-', '_#2019-', tips$`tree$tip.label`)
tipsplit<-cSplit(tips, "tree$tip.label", sep = "#")
tiplabels<-as.data.frame(tree$tip.label)
tipdates<-as.data.frame(cbind(tiplabels,tipsplit[,2]))
colnames(tipdates)<-c("tip","date")
tipdates$date<-as.POSIXct(tipdates$date)
tipdates$date<-decimal_date(tipdates$date)

dates<-as.vector(tipdates$date)
names(dates)<-tipdates$tip

#### generate a rooted time-scaled tree using treedater ####
treedater<-dater(tree, dates, s=29409, omega0 = c(0.001), clock = c("strict"))


#### identify outlier tips and filter tips below q threshold of 0.05 ####
out_analysis<-outlierTips(treedater, alpha = 0.05)

outliers<-subset(out_analysis,out_analysis$q<=0.05)
out_to_filter<-as.vector(outliers$taxon)

tree_nooutliers<-drop.tip(treedater, out_to_filter)


#### determine cophenetic distances ####
dist<-cophenetic.phylo(treedater)
dist<-as.data.frame(dist)

REGION_names<-read.table(file = "PATH/TO/REGION_names.txt")
colnames(REGION_names)<-c("tip")
REGION_names$tip<-gsub('/', '_', REGION_names$tip)
REGION_names$tip<-gsub('\\|', '_', REGION_names$tip)
row.names(REGION_names)<-REGION_names$tip

dist2<-subset(dist, rownames(dist) %in% rownames(REGION_names))

dist3<-1-dist2
dist3[dist3 == 1]<-0
closest<-as.data.frame(colnames(dist3)[max.col(dist3,ties.method="first")])
row.names(closest)<-row.names(dist3)
colnames(closest)<-c("tip")


#### Select time-stratified tips for exogenous sequences ####
timedtips<-as.data.frame(tree_nooutliers$tip.label)
timedtips$`tree_nooutliers$tip.label`<-gsub('_2020-', '_#2020-', timedtips$`tree_nooutliers$tip.label`)
timedtips$`tree_nooutliers$tip.label`<-gsub('_2019-', '_#2019-', timedtips$`tree_nooutliers$tip.label`)
timedtipsplit<-cSplit(timedtips, "tree_nooutliers$tip.label", sep = "#")
timedtiplabels<-as.data.frame(tree_nooutliers$tip.label)
timedtipdates<-as.data.frame(cbind(timedtiplabels,timedtipsplit[,2]))
colnames(timedtipdates)<-c("tip","date")
timedtipdates$date<-as.POSIXct(timedtipdates$date)
timedtipdates$date<-decimal_date(timedtipdates$date)

`%notin%` <- Negate(`%in%`)
timed_noREGION<-subset(timedtipdates, timedtipdates$tip %notin% REGION_names$tip)

timed_noREGION$diff<-ceiling((timed_noREGION$date-min(timed_noREGION$date))*365)
timed_noREGION$week<-floor(timed_noREGION$diff/7)+1
timed_noREGION$week<-paste0("week",timed_noREGION$week,sep="")

list_timed_noREGION <- split(timed_noREGION, timed_noREGION$week)

weeks<-as.list(unique(timed_noREGION$week))
len<-length(weeks)

for (i in 1:len) {
  a<-list_timed_noREGION[[i]]
  b<-ifelse(nrow(a)<5,nrow(a),
            ifelse(nrow(a)<20,5,10))
  assign(paste0("Random",i),a[sample(nrow(a),b), ])
}

dfs<-paste0("Random",1:len)
randomlist = lapply(dfs, get)
timestrat_timed<-bind_rows(randomlist, .id=NULL)
timestrat_names<-as.data.frame(timestrat_timed$tip)
colnames(timestrat_names)<-c("tip")

timestrat_closest<-merge(timestrat_names, closest, by="tip", all = TRUE)
timestrat_closest<-merge(timestrat_closest, REGION_names, by="tip", all = TRUE)
timestrat_closest$tip<-gsub("_", "/", timestrat_closest$tip)
timestrat_closest$tip<-gsub("hCoV-19_", "hCoV-19/", timestrat_closest$tip)
timestrat_closest$tip<-gsub("/EPI/ISL/", "|EPI_ISL_", timestrat_closest$tip)
timestrat_closest$tip<-gsub("/2020-", "|2020-", timestrat_closest$tip)
timestrat_closest$tip<-gsub("/2019-", "|2019-", timestrat_closest$tip)
timestrat_closest_unique<-as.data.frame(unique(timestrat_closest$tip))

#### Write names of exogenous sequences ####
output1<-paste0(outdir,"timestrat_closest_names.txt",sep="/")
write.table(timestrat_closest_unique, file = output1, row.names=FALSE, col.names=FALSE, quote=FALSE)
