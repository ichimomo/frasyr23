library(tidyverse)
datafile <- read.csv("~/GoogleFRAdrive/2-kei-issues/type2_data.csv")

Stocks <-unique(datafile$Stock)

ccdata.stock<-list()
for(i in 1:length(unique(datafile$Stock))){
  ccdata.stock[[i]] <- datafile %>% filter(Stock==Stocks[i])
}
ABCs<-list()
n.catch<-5
for(i in 1:length(Stocks)){
  if(max(ccdata.stock[[i]]$Year)-min(ccdata.stock[[i]]$Year)<=4) next
  cpuetmp <-ccdata.stock[[i]]$CPUE
  cpuetmp <- cpuetmp[-na.omit(cpuetmp)]
  if(length(cpuetmp) <=4 ) next
  ccdata<-data.frame(year=ccdata.stock[[i]]$Year,cpue=ccdata.stock[[i]]$CPUE,catch=ccdata.stock[[i]]$Catch)
  filename<-paste0("~/Desktop/2kei-stocks/",Stocks[i],"_default_bt5year.png")
  resabc2 <-calc_abc2(ccdata,BTyear=(max(ccdata$year)-4),tune.par = c(0.5,0.4,0.4),summary_abc = F)
  #resabc2 <-calc_abc2(ccdata,summary_abc = F)
  graph_abc2 <-plot_abc2_fixHC_seqOut(resabc2)
  ABCs[[i]]<-graph_abc2[[1]]
  ABCs[[i]]$stock <- rep(Stocks[i],nrow(ABCs[[i]]))
  ABCs[[i]]$tunepar <- rep(str_c("0.5-0.4-0.4"),nrow(ABCs[[i]]))
  ABCs[[i]]$ABCdeviation <- (ABCs[[i]]$ABC-ABCs[[i]]$ABC[1])/ABCs[[i]]$ABC[1]
  ori.catch <- ccdata$catch
  l.catch <- length(ori.catch)
  Catch5yr <- mean(ori.catch[(l.catch-n.catch+1):l.catch],na.rm = TRUE)
  ABCs[[i]]$Catch5yr <- rep(Catch5yr,nrow(ABCs[[i]]))
  ABCs[[i]]$Catch5yrdeviation <- (ABCs[[i]]$ABC-Catch5yr)/Catch5yr
  #names(ABCs[[i]])<-Stocks[i]
  #ggsave(width=420,height=150,dpi=200,units="mm", graph_abc2[[3]],file=filename)
}

#ABCs1st<-ABCs
#ABCs2nd<-ABCs
ABCsdefault<-ABCs

# for(i in 1:length(Stocks)){
#   if(is.null(ABCs[[i]])) next
#
# }
