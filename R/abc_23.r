#' 2,3系のABC計算をするためのパッケージ
#'
#' @import dplyr
#' @import tibble
#' @import stringr
#' @import magrittr
#' @import ggplot2
#' 
#'　
#' 

"_PACKAGE"

#' 2系のABC計算用の関数
#' 
#' @param ccdata year,cpue,catchのデータフレーム
#' @param BT  initial target level
#' @param PL  BL = PL*BT
#' @param PB  BB = PB*BT
#' @param tune.par  tuning parameters: (beta, delta, lambda)   
#' @param AAV if this value is given, given value is used
#' @param n.catch  period for averaging the past catches
#'
#'
#' @examples
#' library(frasyr23)
#' catch <- c(10,5,10,3,2,1,3)
#' cpue <- c(11,10,2,3,2,5,2)
#' example_data <- data.frame(year=1:7,cpue=cpue,catch=catch)
#'
#' # 2系
#' example_abc2 <- calc_abc2(example_data)
#' example_abc2$ABC
#' # $ABC=1.894191
#' graph_abc2 <- plot_abc2(example_abc2)
#'
#' @export
#'

calc_abc2 <- function(
  ccdata,   # data frame for year, CPUE, and catch
  BT=0.8,   # initial target level
  PL=0.7,   #  BL = PL*BT
  PB=0.0,   #  BB = PB*BT
  tune.par = c(0.5,0.4,0.4), #  tuning parameters: (beta, delta, lambda)   
  AAV=NULL,
  n.catch=5   #  period for averaging the past catches
){
    argname <- ls() # 引数をとっておいて再現できるようにする
    arglist <- lapply(argname,function(xx) eval(parse(text=xx)))
    names(arglist) <- argname

    cpue <- ccdata$cpue
    cpue <- cpue[!is.na(cpue)]
    catch <- ccdata$catch
    catch <- catch[!is.na(catch)]
  
    beta <- tune.par[1]   # velocity to go to BT
    delta <- tune.par[2]   # correction factor when D <= BL
    lambda <- tune.par[3]   # tuning parameter for updating BT
    BT <- BT      # Btarget
    BL <- PL*BT      # Blimit
    BB <- PB*BT      # Bban
    BRP <- c(BT, BL, BB)
    
    n <- length(catch)   # the number of catch data

    cum.cpue <- function(x) pnorm(scale(x),0,1) # cumulative normal distribution
    cum.cpue2 <- function(x) pnorm(x,mean(x),sd(x)) # cumulative normal distribution    
    mean.catch <- mean(catch[(n-n.catch+1):n])

    D <- cum.cpue(as.numeric(cpue))              # cumulative probability of cpue
    mD <- attributes(D)$'scaled:center'         # mean of cpue
    sD <- attributes(D)$'scaled:scale'           # standard deviation of cpue
    cD <- D[n]                                   # final depletion
    
    icum.cpue <- function(x) sD*qnorm(x,0,1)+mD   # inverse function from D to CPUE
  

    if(is.null(AAV)) if (lambda > 0) AAV <- aav.f(cpue) else AAV <- 0

#    k <- ifelse(cD > BB, beta+(cD <= BL)*delta*exp(lambda*log(AAV^2+1))*(BL-cD)/(cD-BB), Inf)    #  calculation of k
#    ABC <- ifelse(cD > BB & cpue[n] > 0, mean.catch*exp(k*(cD-BT)), 0)    # calculation of ABC
    #    alpha <- exp(k*(cD-BT))
    alpha <- type2_func(cD,cpue[n],BT=BT,PL=PL,PB=PB,AAV=AAV,tune.par=tune.par)
    ABC <- mean.catch * alpha
    
    Obs_BRP <- c(icum.cpue(BT), icum.cpue(BL), icum.cpue(BB))
    Obs_percent <- icum.cpue(seq(from=0.1,to=0.9,by=0.1))
    Current_Status <- c(D[n],cpue[n])
    names(Current_Status) <- c("Level","CPUE")

    names(BRP) <- names(Obs_BRP) <- c("Target","Limit","Ban")
    
    output <- list(BRP=BRP,Obs_BRP=Obs_BRP,Current_Status=Current_Status,
                   AAV=AAV,tune.par=tune.par,ABC=ABC,arglist=arglist,mean.catch=mean.catch,Obs_percent=Obs_percent,
                   alpha=alpha)

  return(output)
}

type2_func <- function(cD,cpue.n,BT=0.8,PL=0.7,PB=0,AAV=0.4,tune.par=c(0.5,0.5,0.4)){
    beta <- tune.par[1]   # velocity to go to BT
    delta <- tune.par[2]   # correction factor when D <= BL
    lambda <- tune.par[3]   # tuning parameter for updating BT
    BT <- BT      # Btarget
    BL <- PL*BT      # Blimit
    BB <- PB*BT      # Bban
    
    if(cD <= BB) alpha <- 0
    if(BB < cD & cD < BL){
        k <- beta + delta* exp(lambda*log(AAV^2+1)) * (BL-cD)/(cD-BB)
        alpha <- exp(k*(cD-BT))
    }
    if(cD > BL) alpha <- exp(beta*(cD-BT))
    return(alpha)
    # cpue.nは必要か？
    #    k <- ifelse(cD > BB, beta+(cD <= BL)*delta*exp(lambda*log(AAV^2+1))*(BL-cD)/(cD-BB), Inf)    #  calculation of k
    #    ifelse(cD > BB & cpue.n > 0, exp(k*(cD-BT)), 0)    # calculation of ABC    
}

type2_func_wrapper <- function(DL,...){
    purrr::map_dbl(DL,type2_func,...)    
}


#' 3系のABC計算用の関数
#' 
#' @param catch catch timeseries data
#' @param BT  initial target level
#' @param PL  BL = PL*BT
#' @param PB  BB = PB*BT
#' @param tune.par  tuning parameters: (beta, delta)   
#' @param n.catch  period for averaging the past catches
#'
#' @examples
#' 
#' library(frasyr23)
#' 
#' catch <- c(10,5,10,3,2,1,3)
#' cpue <- c(11,10,2,3,2,5,2)
#' example_data <- data.frame(year=1:7,cpue=cpue,catch=catch)
#'
#' example_abc3 <- calc_abc3(example_data)
#' example_abc3$ABC
#' # [1] 2.658756
#' graph_abc3 <- plot_abc3(example_abc3)
#'
#' @export
#'

calc_abc3 <- function(
  ccdata,   # data frame for year, CPUE, and catch
  BT=0.1,   # initial target level
  PL=2,   #  BL = PL*BT
  PB=10,   #  BB = PB*BT
  tune.par = c(1.5,2), #  tuning parameters: (beta, delta, lambda)   
  n.catch=5   #  period for averaging the past catches
){
    argname <- ls() # 引数をとっておいて再現できるようにする
    arglist <- lapply(argname,function(xx) eval(parse(text=xx)))
    names(arglist) <- argname

    catch <- ccdata$catch
    catch <- catch[!is.na(catch)]
  
    beta <- tune.par[1]   # velocity to go to BT
    delta <- tune.par[2]   # correction factor when D <= BL

    n <- length(catch)   # the number of catch data
    mean.catch <- mean(catch[(n-n.catch+1):n])
  
    max.cat <- max(catch,na.rm=TRUE)    # max of catch
    D <- catch/max.cat
    cD <- D[n]            # current catch level

    BT <- BT      # Btarget
    BL <- PL*BT      # Blimit
    BB <- PB*BT      # Bban
    BRP <- c(BT, BL, BB)
 
    #    k <- ifelse(cD < BB, -(beta+(cD >= BL)*delta*(BL-cD)/(cD-BB)), -Inf)     #  calculation of k
    #    alpha <- exp(k*(cD-BT))
    alpha <- type3_func(cD,BT=BT,PL=PL,PB=PB,tune.par=tune.par)
    ABC <- mean.catch * alpha    
    
#    ABC <- ifelse(cD < BB, mean.catch*alpha, 0)     # calculation of ABC
  
    Obs_BRP <- max.cat*c(BT, BL, BB)
    Current_Status <- c(D[n],catch[n])
    names(Current_Status) <- c("Level","Catch")
    names(BRP) <- names(Obs_BRP) <- c("Target","Limit","Ban")
    
    output <- list(BRP=BRP,Obs_BRP=Obs_BRP,Current_Status=Current_Status,
                   tune.par=tune.par,ABC=ABC,arglist=arglist,mean.catch=mean.catch,
                   alpha=alpha)

  return(output)
}

type3_func <- function(cD,BT=0.1,PL=2,PB=10,tune.par=c(0.5,0.5)){
    beta <- tune.par[1]   # velocity to go to BT
    delta <- tune.par[2]   # correction factor when D <= BL
    BT <- BT      # Btarget
    BL <- PL*BT      # Blimit
    BB <- PB*BT      # Bban
    #    k <- ifelse(cD < BB, -(beta+(cD >= BL)*delta*(BL-cD)/(cD-BB)), -Inf)     #  calculation of k
    if(cD<BL) k <- beta
    if(BL <= cD & cD < BB) k <- beta+delta*(BL-cD)/(cD-BB)
    if(cD >= BB) k <- Inf
    alpha <- exp(-k*(cD-BT))
    return(alpha)
}

type3_func_wrapper <- function(DL,...){
    purrr::map_dbl(DL,type3_func,...)    
}


aav.f <- function(x){
  n <- length(x)
  aav <- 2*abs(x[2:n]-x[1:(n-1)])/(x[2:n]+x[1:(n-1)])
  aav <- ifelse(aav < Inf, aav, NA)
  mean(aav,na.rm=TRUE)
}

diag.plot <- function(dat,res,lwd=3,cex=1.5,legend.location="topleft",main=""){
  Year <- dat$Year
  Catch <- dat$Catch
  CPUE <- dat$CPUE

  n <- length(Catch)   # the number of catch data

  if (res$catch.only) {
    D <- Catch/max(Catch,na.rm=TRUE)
    Level.name <- "Harvest Level"
  } else {
    cum.cpue <- function(x) pnorm(scale(x),0,1) # cumulative normal distribution
    D <- cum.cpue(as.numeric(CPUE))
    Level.name <- "Depletion Level"
  }
  
  Catch.change <- c(Catch[2:n]/Catch[1:(n-1)],res$ABC/Catch[n])
  D.change <- c(D[-1],NA)
  
  plot.range <- range(c(0,1),Catch.change, D.change,na.rm=TRUE)
  
  year <- c(Year[-1],Year[n]+1)
  
  plot(year,Catch.change,type="b",xlab="Year",ylab="",lwd=lwd,col=c(rep("blue",n-1),"brown"),ylim=plot.range,pch=c(rep(16,n-1),18),cex=cex,main=main,lty=2)
  lines(year,D.change,col="orange",lwd=lwd,lty=1)
  abline(h=res$BRP,lty=1:3,col=c("green","yellow","red"),lwd=lwd)
  legend(legend.location,c(expression(C[t]/C[t-1]),Level.name),col=c("blue","orange"),lwd=lwd,lty=c(1,2),cex=cex)
}


#' 2系のABC計算結果のをプロットするための関数
#' 
#' @param res calc_abc2の返り値
#'
#' @export
#' 

plot_abc2 <- function(res){
    # plot
    ccdata <- res$arglist$ccdata    
    years <- ccdata$year
    last.year <- rev(years)[1]
    data_catch <- tibble(year=c((last.year-res$arglist$n.catch+1):last.year,last.year+2),
                         catch=c(rep(res$mean.catch,res$arglist$n.catch),res$ABC),
                         type=c(rep(str_c(res$arglist$n.catch,"年平均"),5),"ABC"))
    data_BRP <- tibble(BRP=names(res$BRP),value_obs=res$Obs_BRP,
                       value_ratio=res$BRP)
    
    data_percent <- tibble(x=rep(min(years),9),
                               y=res$Obs_percent,
                               label=str_c(seq(from=0.1,to=0.9,by=0.1)*100,"%"))        
    (g.cpue <- ccdata %>% ggplot() +
         geom_hline(yintercept=res$Obs_percent,color="gray",linetype=2)+
         geom_text(data=data_percent,aes(x=x,y=y,label=label))+
         geom_hline(data=data_BRP,mapping=aes(yintercept=value_obs,color=BRP))+
         geom_path(aes(x=year,y=cpue))+
         theme_bw()+ylab("資源量指数")+
         ylim(0,NA)+
         theme(legend.position="top"))

    BT <- res$arglist$BT
    PL <- res$arglist$PL
    PB <- res$arglist$PB
    tune.par <- res$arglist$tune.par
    
    (g.hcr <- ggplot(data=data.frame(X=c(0,1.2)), aes(x=X)) +
         stat_function(fun=type2_func_wrapper,
                       args=list(BT=BT,PL=0,PB=PB,tune.par=tune.par,AAV=res$AAV),color="gray")+
         stat_function(fun=type2_func_wrapper,
                       args=list(BT=BT,PL=PL,PB=PB,tune.par=tune.par,AAV=res$AAV),color="black")+         
         geom_point(aes(x=res$Current_Status[1],y=res$alpha),color=2)+
         geom_vline(data=data_BRP,mapping=aes(xintercept=value_ratio,color=BRP))+
         theme_bw()+
         theme(legend.position="top"))

    (g.catch <- ccdata %>% ggplot() +
         geom_path(data=data_catch,mapping=aes(x=year,y=catch,color=type),lwd=3)+
         geom_point(data=data_catch,mapping=aes(x=year,y=catch,color=type),lwd=3)+     
         geom_path(aes(x=year,y=catch))+
         theme_bw()+ylab("漁獲量")+
         ylim(0,NA)+     
         theme(legend.position="top"))

    graph.all <- gridExtra::grid.arrange(g.cpue,g.hcr,g.catch,ncol=3)
    graph.all
}


#' 3系のABC計算結果のをプロットするための関数
#' 
#' @param res calc_abc3の返り値
#'
#' @export
#' 

plot_abc3 <- function(res){
    # plot
    ccdata <- res$arglist$ccdata    
    years <- ccdata$year
    last.year <- rev(years)[1]
    data_catch <- tibble(year=c((last.year-res$arglist$n.catch+1):last.year,last.year+2),
                         catch=c(rep(res$mean.catch,res$arglist$n.catch),res$ABC),
                         type=c(rep(str_c(res$arglist$n.catch,"年平均"),5),"ABC"))
    
    data_BRP <- tibble(BRP=names(res$BRP),value_obs=res$Obs_BRP,
                       value_ratio=res$BRP)
    BT <- res$arglist$BT
    PL <- res$arglist$PL
    PB <- res$arglist$PB
    tune.par <- res$arglist$tune.par
    
    (g.hcr <- ggplot(data=data.frame(X=c(0,1.2)), aes(x=X)) +
         stat_function(fun=type3_func_wrapper,
                       args=list(BT=BT,PL=10,PB=PB,tune.par=tune.par),color="gray")+
         stat_function(fun=type3_func_wrapper,
                       args=list(BT=BT,PL=PL,PB=PB,tune.par=tune.par),color="black")+         
         geom_point(aes(x=res$Current_Status[1],y=res$alpha),color=2)+
         geom_vline(data=data_BRP,mapping=aes(xintercept=value_ratio,color=BRP))+
         theme_bw()+
         theme(legend.position="top"))

    (g.catch <- ccdata %>% ggplot() +
         geom_path(data=data_catch,mapping=aes(x=year,y=catch,color=type),lwd=3)+
         geom_point(data=data_catch,mapping=aes(x=year,y=catch,color=type),lwd=3)+     
         geom_path(aes(x=year,y=catch))+
         theme_bw()+ylab("漁獲量")+
         geom_hline(data=data_BRP,mapping=aes(yintercept=value_obs,color=BRP))+         
         ylim(0,NA)+     
         theme(legend.position="top"))

    graph.all <- gridExtra::grid.arrange(g.hcr,g.catch,ncol=2)
    graph.all
}


#' 2系・3系のABC計算関数．岡村さん作成のプロトタイプ．チェック用に使う．
#'
#' @export
#'

abc_t23_proto1 <- function(
  catch,   # catch timeseries data
  cpue=NULL,   # cpue timeseries data
  BT=0.8,   # initial target level
  tune.par = c(0.5,0.4,0.4), #  tuning parameters: (beta, delta, lambda) 
  PL=0.7,   #  BL = PL*BT
  PB=0.0,   #  BB = PB*BT
  catch.only=FALSE,    # catch only method
  default=TRUE,
  n.catch=5   #  period for averaging the past catches
){
  #
  # C[t+1] = C[t]*exp(k*(D-BT))
  # 
  
  if (default & (sum(cpue,na.rm=TRUE)==0 | catch.only)) {catch.only <- TRUE; tune.par <- c(1.5,2,0); BT <- 0.1; PL <- 2; PB <- 10}
  
  beta <- tune.par[1]   # velocity to go to BT
  delta <- tune.par[2]   # correction factor when D <= BL
  lambda <- tune.par[3]   # tuning parameter for updating BT

  n <- length(catch)   # the number of catch data

  cum.cpue <- function(x) pnorm(scale(x),0,1) # cumulative normal distribution
  
  if (catch.only){
    max.cat <- max(catch,na.rm=TRUE)    # max of catch
    D <- catch/max.cat
    cD <- D[n]            # current catch level

    AAV <- NA

    BT <- BT      # Btarget
    BL <- PL*BT      # Blimit
    BB <- PB*BT      # Bban

    BRP <- c(BT, BL, BB)
 
    k <- ifelse(cD < BB, -(beta+(cD >= BL)*delta*(BL-cD)/(cD-BB)), -Inf)     #  calculation of k    
    abc <- ifelse(cD < BB, mean(catch[(n-n.catch+1):n])*exp(k*(cD-BT)), 0)     # calculation of ABC
  
    Obs_BRP <- max.cat*c(BT, BL, BB)
    Current_Status <- c(D[n],catch[n])
    names(Current_Status) <- c("Level","Catch")
  } else {
    D <- cum.cpue(as.numeric(cpue))              # cumulative probability of cpue
    mD <- attributes(D)$'scaled:center'         # mean of cpue
    sD <- attributes(D)$'scaled:scale'           # standard deviation of cpue
    cD <- D[n]                                   # final depletion
    
    icum.cpue <- function(x) sD*qnorm(x,0,1)+mD   # inverse function from D to CPUE
  
    BT <- BT      # Btarget
    BL <- PL*BT      # Blimit
    BB <- PB*BT      # Bban

    BRP <- c(BT, BL, BB)

    if (lambda > 0) AAV <- aav.f(cpue) else AAV <- 0

    k <- ifelse(cD > BB, beta+(cD <= BL)*delta*exp(lambda*log(AAV^2+1))*(BL-cD)/(cD-BB), Inf)    #  calculation of k
    abc <- ifelse(cD > BB & cpue[n] > 0, mean(catch[(n-n.catch+1):n],na.rm=TRUE)*exp(k*(cD-BT)), 0)    # calculation of ABC
    
    Obs_BRP <- c(icum.cpue(BT), icum.cpue(BL), icum.cpue(BB))
    Current_Status <- c(D[n],cpue[n])
    names(Current_Status) <- c("Level","CPUE")
  }

    names(BRP) <- names(Obs_BRP) <- c("Target","Limit","Ban")
    
    output <- list(BRP=BRP,Obs_BRP=Obs_BRP,Current_Status=Current_Status,catch.only=catch.only,AAV=AAV,tune.par=tune.par,k=k,ABC=abc)

  return(output)
}


