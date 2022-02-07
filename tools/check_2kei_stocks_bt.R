library(tidyverse)
datafile <- read.csv("~/GoogleFRAdrive/2-kei-issues/type2_data.csv")

Stocks <-unique(datafile$Stock)

ccdata.stock<-list()
for(i in 1:length(unique(datafile$Stock))){
  ccdata.stock[[i]] <- datafile %>% filter(Stock==Stocks[i])
}
ABCs<-list()
names(ccdata.stock[[1]])
for(i in 1:length(Stocks)){
  if(max(ccdata.stock[[i]]$Year)-min(ccdata.stock[[i]]$Year)<=4) next
  cpuetmp <-ccdata.stock[[i]]$CPUE
  cpuetmp <- cpuetmp[-na.omit(cpuetmp)]
  if(length(cpuetmp) <=4 ) next
  ccdata<-data.frame(year=ccdata.stock[[i]]$Year,cpue=ccdata.stock[[i]]$CPUE,catch=ccdata.stock[[i]]$Catch)
  filename<-paste0("~/Desktop/2kei-stocks/",Stocks[i],"_2ndbest_bt5year.png")
  resabc2 <-calc_abc2(ccdata,BTyear=(max(ccdata$year)-4),tune.par = c(0.4,0.4,0.8),summary_abc = F)
  #resabc2 <-calc_abc2(ccdata,summary_abc = F)
  graph_abc2 <-plot_abc2_fixHC_seqOut(resabc2,outABCs = T)
  ABCs[[i]]<-graph_abc2[[1]]
  #names(ABCs[[i]])<-Stocks[i]
  ggsave(width=420,height=150,dpi=200,units="mm", graph_abc2[[3]],file=filename)
}
