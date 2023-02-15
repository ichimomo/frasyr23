source("./tools/generate-testdata/abc_t23.r")

#data_akaをつかう。
data("data_aka")

n.catch <- 5
year <- data_aka$year
catch <- data_aka$catch
cpue <- data_aka$cpue

ccTF<-!is.na(catch)*!is.na(cpue)
ccavail<-TRUE
availy<-0
while(ccavail){
  availy<-availy+1
  ccavail<-!as.logical(ccTF[availy])
}

ABC<-c()
for(i in length(ccTF):(availy+n.catch-1)){
  aka_abc2 <- abc_t23(catch=catch,cpue=cpue)
  ABC <- c(ABC,aka_abc2$ABC)
  catch <- catch[-(length(catch))]
  cpue <- cpue[-(length(cpue))]
}
ABCyear<-c(year[length(ccTF)]:year[availy+n.catch-1])
aka_retro2<-data.frame(year=rev(ABCyear),abc=rev(ABC))
save(aka_retro2,file="./inst/extdata/res_aka_retro2.rda")
