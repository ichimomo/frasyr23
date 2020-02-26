#
# New ABC Program for Type 2 and 3 Fishery
#

abc_t23 <- function(
  catch,   # catch timeseries data
  cpue=NULL,   # cpue timeseries data
  BT=0.8,   # initial target level
  tune.par = c(0.5,0.4,0.4), #  tuning parameters: (beta, delta, lambda)
  PL=0.7,   #  BL = PL*BT
  PB=0.0,   #  BB = PB*BT
  catch.only=FALSE,    # catch only method
  m=1,    # modification coefficient
  default=TRUE,    # use the default setting
  max.change=NULL,   # allowed maximum change
  n.catch=5   #  period for averaging the past catches
){
  #
  # C[t+1] = C[t]*exp(k*(D-BT))
  #

  if(is.null(cpue))catch.only<-TRUE

  if (default) if (catch.only) {BT=0.1;PL=6;PB=10;tune.par=c(2.0,1.0,0.0);n.catch=3} else {BT=0.8;PL=0.7;PB=0.0;tune.par=c(0.5,0.4,0.4);n.catch=5}

  beta <- tune.par[1]   # velocity to go to BT
  delta <- tune.par[2]   # correction factor when D <= BL
  lambda <- tune.par[3]   # tuning parameter for updating BT

  n <- length(catch)   # the number of catch data

  cum.cpue <- function(x) pnorm(scale(x),0,1) # cumulative normal distribution

  if (catch.only){
    max.cat <- max(catch,na.rm=TRUE)    # max of catch
    D <- catch/max.cat
    cD <- mean(D[(n-n.catch+1):n],na.rm=TRUE)            # current catch level

    AAV <- NA

    BT <- BT      # Btarget
    BL <- PL*BT      # Blimit
    BB <- PB*BT      # Bban

    BRP <- c(BT, BL, BB)

    k <- ifelse(cD < BB, -(beta+(cD >= BL)*delta*(BL-cD)/(cD-BB)), -Inf)     #  calculation of k
    abc <- m*ifelse(cD < BB, mean(catch[(n-n.catch+1):n],na.rm=TRUE)*exp(k*(cD-BT)), 0)     # calculation of ABC

    Obs_BRP <- max.cat*c(BT, BL, BB)
    Current_Status <- c(cD,mean(catch[(n-n.catch+1):n],na.rm=TRUE))
    names(Current_Status) <- c("Level","Catch")
  } else {
    D <- cum.cpue(as.numeric(cpue))              # cumulative probability of cpue
    mD <-  attributes(D)$'scaled:center'         # mean of cpue
    sD <- attributes(D)$'scaled:scale'           # standard deviation of cpue
    cD <- D[n]                                   # final depletion

    icum.cpue <- function(x) sD*qnorm(x,0,1)+mD   # inverse function from D to CPUE

    BT <- BT      # Btarget
    BL <- PL*BT      # Blimit
    BB <- PB*BT      # Bban

    BRP <- c(BT, BL, BB)

    if (lambda > 0) AAV <- aav.f(cpue) else AAV <- 0

    k <- ifelse(cD > BB, beta+(cD <= BL)*delta*exp(lambda*log(AAV^2+1))*(BL-cD)/(cD-BB), Inf)    #  calculation of k
    abc <- m*ifelse(cD > BB & cpue[n] > 0, mean(catch[(n-n.catch+1):n],na.rm=TRUE)*exp(k*(cD-BT)), 0)    # calculation of ABC

    Obs_BRP <- c(icum.cpue(BT), icum.cpue(BL), icum.cpue(BB))
    Current_Status <- c(D[n],cpue[n])
    names(Current_Status) <- c("Level","CPUE")
  }

    if (!is.null(max.change)) {
      if (abc > catch[n]*(1+max.change)) abc <- catch[n]*(1+max.change)
      if (abc < catch[n]*(1-max.change)) abc <- catch[n]*(1-max.change)
    }

    names(BRP) <- names(Obs_BRP) <- c("Target","Limit","Ban")

    output <- list(BRP=BRP,Obs_BRP=Obs_BRP,Current_Status=Current_Status,catch.only=catch.only,AAV=AAV,tune.par=tune.par,k=k,m=m,ABC=abc)

  return(output)
}

aav.f <- function(x){
  n <- length(x)
  aav <- 2*abs(x[2:n]-x[1:(n-1)])/(x[2:n]+x[1:(n-1)])
  aav <- ifelse(aav < Inf, aav, NA)
  mean(aav,na.rm=TRUE)
}
