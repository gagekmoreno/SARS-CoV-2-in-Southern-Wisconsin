library(ggplot2) 
library(lubridate)
library(dplyr)
library(RcppRoll)

traj1 <- read.table('PATH/TO/trajectory1.traj', header = TRUE)
traj2 <- read.table('PATH/TO/trajectory2.traj', header = TRUE)

### 10% burn-in ##
traj1 <- tail(traj1,-(floor(0.1*(nrow(traj1)/200))*200))
traj2 <- tail(traj2,-(floor(0.1*(nrow(traj2)/200))*200))

### Combine trajectory tables ###
traj1$Sample<-paste0(traj1$Sample,"_48d", sep="")
traj2$Sample<-paste0(traj2$Sample,"_48e", sep="")
traj1_ids<-as.data.frame(sort(unique(traj1$Sample)))
traj2_ids<-as.data.frame(sort(unique(traj2$Sample)))
traj1_keepids<-sample_frac(traj1_ids, 0.5)
traj2_keepids<-sample_frac(traj2_ids, 0.5)
colnames(traj1_keepids)<-c("Sample")

traj1_keep<-traj1 %>% filter(Sample %in% traj1_keepids$Sample)
traj2_keep<-traj2 %>% filter(Sample %in% traj2_keepids$Sample)

traj_merged<-rbind(traj1_keep, traj2_keep)

traj_merged$newI <- traj_merged$infections-dplyr::lag(traj_merged$infections)
traj_merged$newI[is.na(traj_merged$newI)]<-0
traj_merged$newI<-ifelse(traj_merged$newI<0,0,traj_merged$newI)

dfs <- split( traj_merged, traj_merged$Sample )
taxis <- dfs[[1]]$t 

qs <- c( .5, .025, .975 )


# infectious (prevalence)
Il = do.call( cbind, lapply( dfs, '[[', 'Il' ))
Ih = do.call( cbind, lapply( dfs, '[[', 'Ih' ))
I = Il + Ih 
t(apply( I, MAR=1, FUN= function(x) quantile(x, qs ))) -> Ici 

# new infections
newI <- do.call(cbind, lapply(dfs, '[[', 'newI'))
NIci <- t(apply( newI, MAR=1, FUN= function(x) quantile(x, qs )))

# cases 
cases <- do.call( cbind, lapply( dfs, '[[', 'infections' ))
t(apply( cases, MAR=1, FUN=function(x) quantile(x,qs))) -> casesci 


#exog 
exog <- do.call( cbind, lapply( dfs, '[[', 'exog' ))
t(apply( exog, MAR=1, FUN=function(x) quantile(x, qs )	)) -> exogci 


#exog 
E <- do.call( cbind, lapply( dfs, '[[', 'E' ))
t(apply( E, MAR=1, FUN=function(x) quantile(x, qs )	)) -> Eci 

#### Plot cumulative incidence ####
pldf <- data.frame( Date = ( date_decimal( taxis ) ) , reported=FALSE )
pldf$`Cumulative infections` = casesci[,1]
pldf$`2.5%` = casesci[,2]
pldf$`97.5%` = casesci[,3] 
pldf <- pldf[ with( pldf, Date > as.Date('2020-02-17') & Date < as.Date('2020-04-29') ) , ]
dates<-pldf$Date

testing<-read.csv(file="PATH/TO/positivetests.csv")
testing$Date<-date_decimal(testing$Date)
testing_sub <- testing[ with( testing, Date > as.Date('2020-02-17') & Date < as.Date('2020-04-25') ) , ]
pldf2<-merge(pldf,testing_sub,by='Date',all = TRUE)

pl <- ggplot( pldf2 ) + 
  geom_path( aes(x = Date, y = `Cumulative infections` , group = !reported), lwd=0.5) + 
  geom_ribbon( aes(x = Date, ymin=`2.5%`, ymax=`97.5%`, group = !reported) , alpha = .25 ) +
  geom_path( aes(x = Date, y = `2.5%` , group = !reported), lwd=0.25, alpha = 0.5) +
  geom_path( aes(x = Date, y = `97.5%` , group = !reported), lwd=0.25, alpha = 0.5) +
  theme_classic() + 
  xlab('') + ylab ('Cumulative estimated infections') + 
  scale_y_continuous(breaks = c(1, 10, 100, 1000, 10000, 100000), limits = c(1,100000), trans = 'log10', labels=scales::label_number(accuracy = 1, big.mark=",")) +
  annotation_logticks(base=10, sides="l") +
  scale_x_datetime(date_breaks = "1 week",) +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  geom_point(aes(x=Date, y=DC), color = "red", alpha=0.5, size = 0.7) +
  geom_vline(xintercept = as.numeric(pldf2$Date[c(64)]), linetype='dashed', alpha=0.4, size=0.25)
pl

ggsave(pl, file = 'Cumulative_Incidence_merged.pdf', width = 3.6, height = 3.25, dpi = 600)
#ggsave(pl, file = 'Cumulative_Incidence_merged.svg', width = 3.6, height = 3.25, dpi = 600)

cumulativecases<-as.data.frame(cbind(taxis,casesci))
colnames(cumulativecases)[1]<-c("date")
write.csv(cumulativecases,file='casesci_DaneCounty_all_phydynmerged.csv',row.names = FALSE)

#### Plot daily new infections ####
newdi <- data.frame( Date = ( date_decimal( taxis ) ) , reported=FALSE )
newdi$`Daily new infections` = NIci[,1]
newdi$`2.5%` = NIci[,2]
newdi$`97.5%` = NIci[,3]
newdi <- newdi[ with( newdi, Date > as.Date('2020-02-26') & Date < as.Date('2020-04-19') ) , ]
dates<-newdi$Date

pl2 <- ggplot( newdi ) + 
  geom_path( aes(x = Date, y = `Daily new infections` , group = !reported), lwd=0.5) + 
  geom_ribbon( aes(x = Date, ymin=`2.5%`, ymax=`97.5%`, group = !reported) , alpha = 0.25 ) +
  geom_path( aes(x = Date, y = `2.5%` , group = !reported), lwd=0.25, alpha = 0.5) +
  geom_path( aes(x = Date, y = `97.5%` , group = !reported), lwd=0.25, alpha = 0.5) +
  theme_classic() + 
  xlab('') + ylab ('Incidence') + 
  scale_y_continuous(breaks = c(0, 1, 10, 100, 1000, 10000, 10000), limits=c(1,10000), trans = 'log10') +
  scale_x_datetime(date_breaks = "1 week",) +
  annotation_logticks(base=10, sides="l") +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  geom_vline(xintercept = as.numeric(pldf2$Date[c(57)]), linetype='dashed', alpha=0.4, size=0.2)
pl2

ggsave(pl2, file = 'Incidence_merged.pdf', width = 3.6, height = 3.25, dpi = 600)
#ggsave(pl2, file = 'Incidence_merged.svg', width = 3.6, height = 3.25, dpi = 600)


#### Plot prevalence over time ####
prev <- data.frame( Date = ( date_decimal( taxis ) ) , reported=FALSE )
prev$`Prevalence` = Ici[,1]
prev$`2.5%` = Ici[,2]
prev$`97.5%` = Ici[,3]
prev <- prev[ with( prev, Date > as.Date('2020-02-26') & Date < as.Date('2020-04-29') ) , ]
dates<-prev$Date
pl3 <- ggplot( prev ) + 
  geom_path( aes(x = Date, y = `Prevalence` , group = !reported), lwd=0.5) + 
  geom_ribbon( aes(x = Date, ymin=`2.5%`, ymax=`97.5%`, group = !reported) , alpha = .25 ) +
  geom_path( aes(x = Date, y = `2.5%` , group = !reported), lwd=0.25, alpha = 0.5) +
  geom_path( aes(x = Date, y = `97.5%` , group = !reported), lwd=0.25, alpha = 0.5) +
  theme_classic() + 
  xlab('') + ylab ('Prevalence') + 
  scale_y_continuous(breaks = c(0, 1, 10, 100, 1000, 10000, 100000), trans = 'log10', labels=scales::label_number(accuracy = 1, big.mark=",")) +
  annotation_logticks(base=10, sides="l") +
  scale_x_datetime(date_breaks = "1 week",) +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  geom_vline(xintercept = as.numeric(pldf2$Date[c(57)]), linetype='dashed', alpha=0.4, size=0.2)
pl3

ggsave(pl3, file = 'Prevalence_merged.pdf', width = 3.6, height = 3.25, dpi = 600)


### Compute R0 post-break from alpha and R0 ###
log<-read.table(file = "/PATH/TO/combinedlog.log", header=TRUE)

log$R02<-log$seir.R0*log$seir.a

print(summary(log$R02))
print(quantile(log$R02, c(0.025, 0.975)))
hist(log$R02)





