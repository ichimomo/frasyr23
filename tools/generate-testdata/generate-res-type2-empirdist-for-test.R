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

# empir.distのsimple.empir or not
simple <- FALSE

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
cpue <- data_aka$cpue[!is.na(data_aka$cpue)]

# empir.dist=T, simple.empir=F
cum.cpue.empir <- ecdf(cpue) # cumulative empirical distribution
cD <- cum.cpue.empir(cpue[n])
D <- cum.cpue.empir(cpue)

ifelse(delta3 > 0, AAV <- aav.f(data_aka$cpue), AAV <- 0)

cum.cpue <- ecdf(cpue)
cpue.order <- sort(unique(cpue))
cpue.prob <-cum.cpue(cpue.order)

trans_empir_prob<-function(prob.seq,cpue.prob){ #CPUEの累積確率分布を使ってあるCPUEの確率を求める
  empir.prob <- prob.seq
  n <- length(cpue.prob)
  for(i in 1:length(prob.seq)){
    loopclose <-FALSE
    j<-1
    while(j < n){
      if(prob.seq[i] < min(cpue.prob)){
        empir.prob[i] <- 0
      }else if(prob.seq[i] >= cpue.prob[j] & prob.seq[i] < cpue.prob[j+1]) {
        empir.prob[i] <-cpue.prob[j]
      }else if(prob.seq[i] == cpue.prob[n]) {
        empir.prob[i] <- cpue.prob[n]
        loopclose <-TRUE
      }
      j <- j+1
    }
  }
  return(empir.prob)
}
cD <- trans_empir_prob(cD,cpue.prob) #データ最終年のCPUEを経験分布をつかって確率へ変換
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

aka_abc2_empir <- list(BRP=BRP,Obs_BRP=Obs_BRP,Current_Status=Current_Status,AAV=AAV,tune.par=tune.par,ABC=ABC)

save(aka_abc2_empir,file="./inst/extdata/res_aka_abc2_empir.rda")

# empir.distのsimple.empir=Tの場合

cD <- (cpue[n]-min(cpue))/(max(cpue)-min(cpue))
D <- (cpue-min(cpue))/(max(cpue)-min(cpue))

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

aka_abc2_empir_simple <- list(BRP=BRP,Obs_BRP=Obs_BRP,Current_Status=Current_Status,AAV=AAV,tune.par=tune.par,ABC=ABC)

save(aka_abc2_empir_simple,file="./inst/extdata/res_aka_abc2_empir_simple.rda")
