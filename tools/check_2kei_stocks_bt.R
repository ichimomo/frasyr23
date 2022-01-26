library(tidyverse)
datafile <- read.csv("~/GoogleFRAdrive/2-kei-issues/type2_data.csv")

Stocks <-unique(datafile$Stock)

ccdata.stock<-list()
for(i in 1:length(unique(datafile$Stock))){
  ccdata.stock[[i]] <- datafile %>% filter(Stock==Stocks[i])
}

names(ccdata.stock[[1]])
for(i in 1:length(Stocks)){
  if(max(ccdata.stock[[i]]$Year)-min(ccdata.stock[[i]]$Year)<=4) next
  cpuetmp <-ccdata.stock[[i]]$CPUE
  cpuetmp <- cpuetmp[-na.omit(cpuetmp)]
  if(length(cpuetmp) <=4 ) next
  ccdata<-data.frame(year=ccdata.stock[[i]]$Year,cpue=ccdata.stock[[i]]$CPUE,catch=ccdata.stock[[i]]$Catch)
  filename<-paste0("~/Desktop/",Stocks[i],".png")
  resabc2 <-calc_abc2(ccdata,BTyear=(max(ccdata$year)-4),summary_abc = F)
  #resabc2 <-calc_abc2(ccdata,summary_abc = F)
  graph_abc2 <-plot_abc2(resabc2,BThcr = T,detABC = 0)
  ggsave(width=420,height=150,dpi=200,units="mm", graph_abc2[[2]],file=filename)
}
