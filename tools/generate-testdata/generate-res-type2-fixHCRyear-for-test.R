source("./tools/generate-testdata/abc_t23.r")

#data_akaをつかう。
data("data_aka")
#head(data_aka)

# セッティングパラメータ
n.catch <- 5
BT <- 0.8
PL <- 0.7
PB <- 0
AAV <- 0.4
tune.par <- c(0.5,0.5,0.4)
beta <- 1.0

BTyear<-2010

delta1 <- tune.par[1]   # velocity to go to BT
delta2 <- tune.par[2]   # correction factor when D <= BL
delta3 <- tune.par[3]   # tuning parameter for updating BT
BT <- BT      # Btarget
BL <- PL*BT      # Blimit
BB <- PB*BT      # Bban
BRP <- c(BT, BL, BB)

n<-nrow(data_aka)
l.catch <- length(data_aka$catch)
mean.catch <- mean(data_aka$catch[(l.catch-n.catch+1):l.catch],na.rm = TRUE)

ccdata_forBt <- data_aka[which(data_aka$year <= BTyear),]
cpue <- ccdata_forBt$cpue
ori.cpue <- cpue
cpue <- cpue[!is.na(ccdata_forBt$cpue)]
target.cpue <- data_aka$cpue[nrow(data_aka$cpue)]

cD <- pnorm(target.cpue,mean(cpue),sd(cpue))
D <- pnorm(scale(cpue),0,1)
mD <- attributes(D)$'scaled:center'
sD <- attributes(D)$'scaled:scale'
icum.cpue <- function(x) sD*qnorm(x,0,1)+mD   # inverse function from D to CPUE

ifelse(delta3 > 0, AAV <- aav.f(data_aka$cpue), AAV <- 0)

if(cD <= BB) alpha <- 0
if(BB < cD & cD < BL){
  k <- delta1 + delta2* exp(delta3*log(AAV^2+1)) * (BL-cD)/(cD-BB)
  alpha <- exp(k*(cD-BT))
}
if(cD >= BL) alpha <- exp(delta1*(cD-BT))
alpha<-alpha*beta

ABC <- mean.catch * alpha
Obs_BRP <- c(icum.cpue(BT), icum.cpue(BL), icum.cpue(BB))
names(BRP) <- names(Obs_BRP) <- c("Target","Limit","Ban")
Current_Status <- c(D[n],cpue[n])
names(Current_Status) <- c("Level","CPUE")

aka_abc2_bt2010 <- list(BRP=BRP,Obs_BRP=Obs_BRP,Current_Status=Current_Status,AAV=AAV,tune.par=tune.par,ABC=ABC)

save(aka_abc2_bt2010,file="./inst/extdata/res_aka_abc2_bt2010.rda")

