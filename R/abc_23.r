#' 2,3系のABC計算をするためのパッケージ
#'
#' @import dplyr
#' @import tibble
#' @import stringr
#' @import magrittr
#' @import ggplot2
#' @import graphics
#' @import stats
#'
#'　
#'

col.BRP <- c("#00533E","#edb918","#C73C2E")

"_PACKAGE"

#' 2系のABC計算用の関数
#'
#' @param ccdata year,cpue,catchのデータフレーム（ラベルは小文字でyear, cpue, catchとする）
#' @param BT  目標水準のレベル（デフォルト=0.8)
#' @param PL  目標水準の何割を限界水準にするか（BL = PL*BT、デフォルトは0.7）
#' @param PB  目標水準の何割を禁漁水準にするか（BB = PB*BT、デフォルトは0）
#' @param tune.par  c(delta1, delta2, delta3) に対応するチューニングパラメータ。デフォルトの値はc(0.5,0.4,0.4)
#' @param beta  保守的なABCのためのパラメータ(デフォルトは1)
#' @param AAV default:"auto", "auto"の場合、CPUEのAAVを内部で計算した値を使う。明示的に数値を与える場合は、その数字が使われる
#' @param n.catch 過去の漁獲量を平均する年数（デフォルトの値は5）
#'
#'
#' @examples
#' library(frasyr23)
#' catch <- c(15,20,13,14,11,10,5,10,3,2,1,3)
#' cpue <- c(10,9,8,4,8,11,10,2,3,2,5,2)
#' example_data <- data.frame(year=2001:2012,cpue=cpue,catch=catch)
#'
#' # 2系
#' example_abc2 <- calc_abc2(example_data,beta=1)
#' graph_abc2 <- plot_abc2(example_abc2)
#'
#' @export
#'

calc_abc2 <- function(
  ccdata,   # data frame for year, CPUE, and catch
  BT=0.8,   # initial target level
  PL=0.7,   #  BL = PL*BT
  PB=0.0,   #  BB = PB*BT
  tune.par = c(0.5,0.4,0.4), #  tuning parameters: (delta1, delta2, delta3)
  AAV="auto", #
  n.catch=5,   #  period for averaging the past catches
  beta = 1.0
){
    argname <- ls() # 引数をとっておいて再現できるようにする
    arglist <- lapply(argname,function(xx) eval(parse(text=xx)))
    names(arglist) <- argname

    cpue <- ccdata$cpue
    ori.cpue <- cpue
    cpue <- cpue[!is.na(ccdata$cpue)]
    ori.catch <- ccdata$catch
    catch <- ccdata$catch
    catch <- catch[!is.na(ccdata$cpue)]
    catch.na.warning <- FALSE

    delta1 <- tune.par[1]   # velocity to go to BT
    delta2 <- tune.par[2]   # correction factor when D <= BL
    delta3 <- tune.par[3]   # tuning parameter for updating BT
    BT <- BT      # Btarget
    BL <- PL*BT      # Blimit
    BB <- PB*BT      # Bban
    BRP <- c(BT, BL, BB)

    n <- length(catch)   # the number of catch data
    l.catch <- length(ori.catch)
    l.cpue <- length(ori.cpue)

    cum.cpue <- function(x) pnorm(scale(x),0,1) # cumulative normal distribution
    cum.cpue2 <- function(x) pnorm(x,mean(x),sd(x)) # cumulative normal distribution
    mean.catch <- mean(ori.catch[(l.catch-n.catch+1):l.catch],na.rm = TRUE)
    for(i in 0:(n.catch-1)){
      if(is.na (ori.catch[l.catch-i])){
        catch.na.warning <- TRUE
      }
    }

    D <- cum.cpue(as.numeric(cpue))              # cumulative probability of cpue
    mD <- attributes(D)$'scaled:center'         # mean of cpue
    sD <- attributes(D)$'scaled:scale'           # standard deviation of cpue
    cD <- D[n]                                   # final depletion

    icum.cpue <- function(x) sD*qnorm(x,0,1)+mD   # inverse function from D to CPUE

    if (delta3 > 0){
        if(AAV=="auto"){
            AAV <- aav.f(ori.cpue)
        }
        else{
            AAV <- AAV
        }}
    else{
            AAV <- 0
        }

#    k <- ifelse(cD > BB, delta1+(cD <= BL)*delta2*exp(delta3*log(AAV^2+1))*(BL-cD)/(cD-BB), Inf)    #  calculation of k
#    ABC <- ifelse(cD > BB & cpue[n] > 0, mean.catch*exp(k*(cD-BT)), 0)    # calculation of ABC
    #    alpha <- exp(k*(cD-BT))
    alpha <- type2_func(cD,cpue[n],BT=BT,PL=PL,PB=PB,AAV=AAV,tune.par=tune.par,beta)
    ABC <- mean.catch * alpha

    Obs_BRP <- c(icum.cpue(BT), icum.cpue(BL), icum.cpue(BB))
    Obs_percent <- icum.cpue(c(0.05,seq(from=0.1,to=0.9,by=0.1),0.95))
    Current_Status <- c(D[n],cpue[n])
    names(Current_Status) <- c("Level","CPUE")

    names(BRP) <- names(Obs_BRP) <- c("Target","Limit","Ban")

    cat("---------------------\n")
    cat(stringr::str_c("Target CPUE value and Level: ",round(Obs_BRP[1],2)," and ", round(BRP[1],2) ,"\n",
                       "Limit CPUE value and Level: ",round(Obs_BRP[2],2)," and ", round(BRP[2],2) ,"\n",
                       "Last year's CPUE value and Level: ",round(cpue[n],3)," and ",
                       round(D[n],3),"\n",
                       "AAV of CPUE: ",round(AAV,3),"\n",
                       "alpha: ",round(alpha,3),"\n",
                       "Average catch: ",round(mean.catch,3),"\n",
                       "ABC in ",max(ccdata$year,na.rm=T)+2,": ",round(ABC,3),"\n"))
    cat("---------------------\n")
if(isTRUE(catch.na.warning))cat("Warning! Recent n.catch year data contains NA.")

    output <- list(BRP=BRP,Obs_BRP=Obs_BRP,Current_Status=Current_Status,
                   AAV=AAV,tune.par=tune.par,ABC=ABC,arglist=arglist,
                   mean.catch=mean.catch,Obs_percent=Obs_percent,
                   D=D,
                   alpha=alpha,beta=beta)

  return(output)
}

type2_func <- function(cD,cpue.n,BT=0.8,PL=0.7,PB=0,AAV=0.4,tune.par=c(0.5,0.5,0.4),beta=1.0){
    delta1 <- tune.par[1]   # velocity to go to BT
    delta2 <- tune.par[2]   # correction factor when D <= BL
    delta3 <- tune.par[3]   # tuning parameter for updating BT
    BT <- BT      # Btarget
    BL <- PL*BT      # Blimit
    BB <- PB*BT      # Bban

    if(cD <= BB) alpha <- 0
    if(BB < cD & cD < BL){
        k <- delta1 + delta2* exp(delta3*log(AAV^2+1)) * (BL-cD)/(cD-BB)
        alpha <- exp(k*(cD-BT))
    }
    if(cD > BL) alpha <- exp(delta1*(cD-BT))
    return(alpha*beta)
    # cpue.nは必要か？
    #    k <- ifelse(cD > BB, delta1+(cD <= BL)*delta2*exp(delta3*log(AAV^2+1))*(BL-cD)/(cD-BB), Inf)    #  calculation of k
    #    ifelse(cD > BB & cpue.n > 0, exp(k*(cD-BT)), 0)    # calculation of ABC
}

type2_func_wrapper <- function(DL,type=NULL,...){
    if(type=="%") DL <- DL/100
    purrr::map_dbl(DL,type2_func,...)
}


#' 3系のABC計算用の関数
#'
#' @param catch catch timeseries data
#' @param BT  initial target level
#' @param PL  BL = PL*BT
#' @param PB  BB = PB*BT
#' @param tune.par  tuning parameters: (gamma1, gamma2)
#' @param n.catch  period for averaging the past catches (default value = 3)
#'
#' @examples
#'
#' library(frasyr23)
#'
#' catch <- c(15,20,13,14,11,10,5,10,3,2,1,3)
#' cpue <- c(10,9,8,4,8,11,10,2,3,2,5,2)
#' example_data <- data.frame(year=2001:2012,cpue=cpue,catch=catch)
#'
#' example_abc3 <- calc_abc3(example_data)
#' graph_abc3 <- plot_abc3(example_abc3)
#'
#' @export
#'

calc_abc3 <- function(
  ccdata,   # data frame for year, CPUE, and catch
  BT=0.1,   # initial target level
  PL=6,   #  BL = PL*BT
  PB=10,   #  BB = PB*BT
  tune.par = c(2,1), #  tuning parameters: (gamma1, gamma2, gamma3)
#  BT=0.05,   # initial target level
#  PL=8,   #  BL = PL*BT
#  PB=20,   #  BB = PB*BT
#  tune.par = c(3.5,3.5), #  tuning parameters: (gamma1, gamma2, gamma3)
  n.catch=3   #  period for averaging the past catches
){
    argname <- ls() # 引数をとっておいて再現できるようにする
    arglist <- lapply(argname,function(xx) eval(parse(text=xx)))
    names(arglist) <- argname

    catch <- ccdata$catch
    #catch <- catch[!is.na(catch)]
    catch.na.warning <- FALSE

    gamma1 <- tune.par[1]   # velocity to go to BT
    gamma2 <- tune.par[2]   # correction factor when D <= BL

    n <- length(catch)   # the number of catch data
    mean.catch <- mean(catch[(n-n.catch+1):n],na.rm=TRUE)

    for(i in 0:(n.catch-1)){
      if(is.na (catch[n-i])){
        catch.na.warning <- TRUE
      }
    }

    max.cat <- max(catch,na.rm=TRUE)    # max of catch
    D <- catch/max.cat
    cD <- mean(D[(n-n.catch+1):n],na.rm=TRUE)            # current catch level

    BT <- BT      # Btarget
    BL <- PL*BT      # Blimit
    BB <- PB*BT      # Bban
    BRP <- c(BT, BL, BB)

    #    k <- ifelse(cD < BB, -(gamma1+(cD >= BL)*gamma2*(BL-cD)/(cD-BB)), -Inf)     #  calculation of k
    #    alpha <- exp(k*(cD-BT))
    alpha <- type3_func(cD,BT=BT,PL=PL,PB=PB,tune.par=tune.par)
    ABC <- mean.catch * alpha

#    ABC <- ifelse(cD < BB, mean.catch*alpha, 0)     # calculation of ABC

    Obs_BRP <- max.cat*c(BT, BL, BB)
    #Current_Status <- c(D[n],catch[n])
    Current_Status <- c(cD,mean(catch[(n-n.catch+1):n],na.rm=TRUE))
    names(Current_Status) <- c("Level","Catch")
    names(BRP) <- names(Obs_BRP) <- c("Target","Limit","Ban")

    cat("---------------------\n")
    cat(stringr::str_c("Target catch value and Level: ",round(Obs_BRP[1],2)," and ", round(BRP[1],2) ,"\n",
                       "Limit catch value and Level: ",round(Obs_BRP[2],2)," and ", round(BRP[2],2) ,"\n",
                       "Last year's catch value and Level: ",round(catch[n],3)," and ",
                       round(D[n],3),"\n",
                       "alpha: ",round(alpha,3),"\n",
                       "Average catch: ",round(mean.catch,3),"\n",
                       "ABC in ",max(ccdata$year,na.rm=T)+2,": ",round(ABC,3),"\n"))
    cat("---------------------\n")
    if(isTRUE(catch.na.warning))cat("Warning! Recent n.catch year data contains NA.")

    output <- list(BRP=BRP,Obs_BRP=Obs_BRP,Current_Status=Current_Status,
                   tune.par=tune.par,ABC=ABC,arglist=arglist,mean.catch=mean.catch,
                   alpha=alpha)

  return(output)
}

type3_func <- function(cD,BT=0.1,PL=6,PB=10,tune.par=c(2.0,1.0)){
                       #BT=0.05,PL=8,PB=20,tune.par=c(3.5,3.5)){
    gamma1 <- tune.par[1]   # velocity to go to BT
    gamma2 <- tune.par[2]   # correction factor when D <= BL
    BT <- BT      # Btarget
    BL <- PL*BT      # Blimit
    BB <- PB*BT      # Bban
    #    k <- ifelse(cD < BB, -(gamma1+(cD >= BL)*gamma2*(BL-cD)/(cD-BB)), -Inf)     #  calculation of k
    if(cD<BL) k <- gamma1
    if(BL <= cD & cD < BB) k <- gamma1+gamma2*(BL-cD)/(cD-BB)
    if(cD >= BB) k <- Inf
    alpha <- exp(-k*(cD-BT))
    return(alpha)
}

type3_func_wrapper <- function(DL,type=NULL,...){
    if(type=="%") DL <- DL/100
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

plot_abc2 <- function(res,stock.name=NULL){
    # plot
    ccdata <- res$arglist$ccdata
    n.catch <- res$arglist$n.catch
    years <- ccdata$year
    last.year <- rev(years)[1]
    data_catch <- tibble(year=c((last.year-res$arglist$n.catch+1):last.year,last.year+2),
                         catch=c(rep(res$mean.catch,res$arglist$n.catch),res$ABC),
                         type=c(rep(str_c(res$arglist$n.catch,"年平均漁獲量"),n.catch),"ABC"))
    data_BRP <- tibble(BRP=names(res$BRP),value_obs=res$Obs_BRP,
                       value_ratio=res$BRP)

    data_percent <- tibble(x=rep(max(years)+2,11),
                               y=res$Obs_percent,
                               label=str_c(c(0.05,seq(from=0.1,to=0.9,by=0.1),0.95)*100,"%"))
    font_MAC <- "HiraKakuProN-W3"#"Japan1GothicBBB"#
    legend.labels <-c("目標水準案","限界水準案","禁漁水準案")
    legend.labels2 <-c(str_c(res$arglist$n.catch,"年平均漁獲量"),"算定漁獲量")
    g.cpue <- ccdata %>% ggplot() +
         geom_hline(yintercept=res$Obs_percent,color="gray",linetype=2)+
         geom_text(data=data_percent,aes(x=x,y=y,label=label))+
         geom_text(aes(x=max(years),y=min(data_percent$y)*0.85,label="(資源量水準)"),size=4)+
        geom_hline(data=data_BRP,mapping=aes(yintercept=value_obs,color=BRP))+
      scale_color_manual(name="",values=rev(c(col.BRP)),labels=rev(c(legend.labels)))+
      geom_path(aes(x=year,y=cpue),size=1)+
         theme_bw()+ylab("資源量指標値")+xlab("漁期年")+
        ylim(0,NA)+theme_custom()+
         ggtitle("資源量指標値のトレンド")+
         theme(legend.position="top")

    if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){
      g.cpue <- ccdata %>% ggplot() +
        geom_hline(yintercept=res$Obs_percent,color="gray",linetype=2)+
        geom_text(data=data_percent,aes(x=x,y=y,label=label))+
        geom_text(aes(x=max(years),y=min(data_percent$y)*0.85,family=font_MAC,label="(資源量水準)"),size=4)+
        geom_hline(data=data_BRP,mapping=aes(yintercept=value_obs,color=BRP))+
        scale_color_manual(name="",values=rev(c(col.BRP)),labels=rev(c(legend.labels)))+
        geom_path(aes(x=year,y=cpue),size=1)+
        theme_bw()+ylab("資源量指標値")+xlab("漁期年")+
        ylim(0,NA)+theme_custom()+
        ggtitle("資源量指標値のトレンド")+
        theme(legend.position="top") +
        theme(text = element_text(family = font_MAC))

    }

    BT <- res$arglist$BT
    PL <- res$arglist$PL
    PB <- res$arglist$PB
    tune.par <- res$arglist$tune.par
    beta <- res$arglist$beta

    g.hcr <- ggplot(data=data.frame(X=c(0,120)), aes(x=X)) +
         stat_function(fun=type2_func_wrapper,
                       args=list(BT=BT,PL=0,PB=PB,tune.par=tune.par,beta=beta,AAV=res$AAV,type="%"),
                       color="gray")+
         stat_function(fun=type2_func_wrapper,
                       args=list(BT=BT,PL=PL,PB=PB,tune.par=tune.par,beta=beta,AAV=res$AAV,type="%"),
                       color="black",size=1)+
         geom_point(aes(x=res$Current_Status[1]*100,y=res$alpha),color=2,size=2)+
        geom_vline(data=data_BRP,mapping=aes(xintercept=value_ratio*100,color=BRP))+
      scale_color_manual(name="",values=rev(c(col.BRP)),labels=rev(c(legend.labels)))+
      theme_bw()+theme_custom()+
         ggtitle("漁獲管理規則")+
         xlab("資源量水準(%)")+ylab(str_c("α (ABC年の漁獲量の削減率)"))+
         theme(legend.position="top")

    if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){
      g.hcr <- ggplot(data=data.frame(X=c(0,120)), aes(x=X)) +
        stat_function(fun=type2_func_wrapper,
                      args=list(BT=BT,PL=0,PB=PB,tune.par=tune.par,beta=beta,AAV=res$AAV,type="%"),
                      color="gray")+
        stat_function(fun=type2_func_wrapper,
                      args=list(BT=BT,PL=PL,PB=PB,tune.par=tune.par,beta=beta,AAV=res$AAV,type="%"),
                      color="black",size=1)+
        geom_point(aes(x=res$Current_Status[1]*100,y=res$alpha),color="red",size=2)+
        geom_vline(data=data_BRP,mapping=aes(xintercept=value_ratio*100,color=BRP))+
        scale_color_manual(name="",values=rev(c(col.BRP)),labels=rev(c(legend.labels)))+
        theme_bw()+theme_custom()+
        ggtitle("漁獲管理規則")+
        xlab("資源量水準(%)")+ylab(str_c("α (ABC年の漁獲量の削減率)"))+
        theme(legend.position="top") +
        theme(text = element_text(family = font_MAC))

    }

    g.catch <- ccdata %>% ggplot() +
        geom_path(data=data_catch,mapping=aes(x=year,y=catch,color=type),lwd=2)+
        geom_point(data=data_catch,mapping=aes(x=year,y=catch,color=type),lwd=3)+
        scale_color_manual(name="",values=c(1,2),labels=legend.labels2)+
        # geom_point(data=dplyr::filter(data_catch,type=="ABC"),
#                    mapping=aes(x=year,y=catch),lwd=2,color=1)+
#         geom_line(data=dplyr::filter(data_catch,type!="ABC"),
#                    mapping=aes(x=year,y=catch),lwd=2,color="gray")+
         geom_path(aes(x=year,y=catch),size=1)+
        ylab("漁獲量")+xlab("漁期年")+
        ggtitle("漁獲量のトレンドとABC")+
        ylim(0,NA)+ theme_custom()+
        theme(legend.position="top")

    if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){
      g.catch <- ccdata %>% ggplot() +
        geom_path(data=data_catch,mapping=aes(x=year,y=catch,color=type),lwd=2)+
        geom_point(data=data_catch,mapping=aes(x=year,y=catch,color=type),lwd=3)+
        scale_color_manual(name="",values=c("black","red"),labels=legend.labels2)+
        # geom_point(data=dplyr::filter(data_catch,type=="ABC"),
        #                    mapping=aes(x=year,y=catch),lwd=2,color=1)+
        #         geom_line(data=dplyr::filter(data_catch,type!="ABC"),
        #                    mapping=aes(x=year,y=catch),lwd=2,color="gray")+
        geom_path(aes(x=year,y=catch),size=1)+
        ylab("漁獲量")+xlab("漁期年")+
        ggtitle("漁獲量のトレンドとABC")+
        ylim(0,NA)+ theme_custom()+
        theme(legend.position="top") +
        theme(text = element_text(family = font_MAC))
    }

    graph.component <- list(g.cpue,g.hcr,g.catch)
    graph.combined <- gridExtra::grid.arrange(g.cpue,g.hcr,g.catch,ncol=3,top=stock.name)
    return(list(graph.component=graph.component,graph.combined=graph.combined))
}


#' 3系のABC計算結果のをプロットするための関数
#'
#' @param res calc_abc3の返り値
#'
#' @export
#'

plot_abc3 <- function(res,stock.name=NULL){
    # plot
    ccdata <- res$arglist$ccdata
    n.catch <- res$arglist$n.catch
    years <- ccdata$year
    last.year <- rev(years)[1]
    data_catch <- tibble(year=c((last.year-res$arglist$n.catch+1):last.year,last.year+2),
                         catch=c(rep(res$mean.catch,res$arglist$n.catch),res$ABC),
                         type=c(rep(str_c(res$arglist$n.catch,"年平均"),n.catch),"ABC"))

    data_BRP <- tibble(BRP=names(res$BRP),value_obs=res$Obs_BRP,
                       value_ratio=res$BRP)
    data_percent <- tibble(x=rep(min(years),10),
                           y=max(ccdata$catch)*1:10/10,
                           label=str_c(1:10/10*100,"%"))
    font_MAC <- "HiraKakuProN-W3"#"Japan1GothicBBB"#
    legend.labels <-c("目標水準案","限界水準案","禁漁水準案")
    legend.labels2 <-c(str_c(res$arglist$n.catch,"年平均漁獲量"),"算定漁獲量",rev(c(legend.labels)))

    BT <- res$arglist$BT
    PL <- res$arglist$PL
    PB <- res$arglist$PB
    tune.par <- res$arglist$tune.par

    g.hcr <- plot_hcr3(res)

    ## (g.hcr <- ggplot(data=data.frame(X=c(0,120)), aes(x=X)) +
    ##      stat_function(fun=type3_func_wrapper,
    ##                    args=list(BT=BT,PL=PL,PB=PB,tune.par=tune.par,type="%")
    ##                   ,color="black",size=1)+
    ##      geom_point(aes(x=res$Current_Status[1]*100,y=res$alpha),color=2,size=1)+
    ## geom_vline(data=data_BRP,mapping=aes(xintercept=value_ratio*100,color=BRP))+
    ## scale_color_manual(values=rev(col.BRP))+
    ##      theme_bw()+theme_custom()+
    ## xlab("漁獲量水準 (漁獲量/最大漁獲量, %)")+ylab(str_c("alpha (ABC年の漁獲量の削減率)"))+
    ## ggtitle("漁獲管理規則")+
    ##      theme(legend.position="top"))

    (g.catch <- ccdata %>% ggplot() +
        geom_path(data=data_catch,mapping=aes(x=year,y=catch,color=type),lwd=2)+
        geom_point(data=data_catch,mapping=aes(x=year,y=catch,color=type),lwd=3)+
#         geom_point(data=dplyr::filter(data_catch,type=="ABC"),
#                    mapping=aes(x=year,y=catch),lwd=2,color=1)+
#         geom_line(data=dplyr::filter(data_catch,type!="ABC"),
#                    mapping=aes(x=year,y=catch),lwd=2,color="gray")+
         geom_path(aes(x=year,y=catch),size=1)+
         geom_text(data=data_percent,aes(x=x,y=y,label=label))+
         geom_text(aes(x=min(ccdata$year)+2,y=min(data_percent$y)*0.85,label="(漁獲量水準)"),size=4)+
         theme_bw()+ylab("漁獲量")+xlab("漁期年")+theme_custom()+
    geom_hline(data=data_BRP,mapping=aes(yintercept=value_obs,color=BRP))+
        scale_color_manual(name="",values=c(1,2,rev(col.BRP)),labels=legend.labels2)+
        ylim(0,NA)+xlim(min(ccdata$year)-1,NA)+
         ggtitle("漁獲量のトレンドとABC")+
         theme(legend.position="top"))

    if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){
      col.set <- c("#000000","#FF0000",rev(col.BRP))
      g.catch <- ccdata %>% ggplot() +
         theme_bw(base_family = font_MAC)+
         geom_path(data=data_catch,mapping=aes(x=year,y=catch,color=type),lwd=2)+
         geom_point(data=data_catch,mapping=aes(x=year,y=catch,color=type),lwd=3)+
        #          geom_point(data=dplyr::filter(data_catch,type=="ABC"),
        #                      mapping=aes(x=year,y=catch),lwd=2,color="red")+
        #          geom_line(data=dplyr::filter(data_catch,type!="ABC"),
        #                      mapping=aes(x=year,y=catch),lwd=3,color="black")+
         geom_path(aes(x=year,y=catch),size=1)+
         ylab("漁獲量")+xlab("漁期年")+theme_custom()+geom_text(data=data_percent,aes(x=x,y=y,label=label),family = font_MAC)+
         geom_text(aes(x=min(ccdata$year)+2,y=min(data_percent$y)*0.85,family=font_MAC,label="(漁獲量水準)"),size=4)+
         geom_hline(data=data_BRP,mapping=aes(yintercept=value_obs,color=BRP))+
        scale_color_manual(name="",values=col.set,labels=legend.labels2)+
ylim(0,NA)+xlim(min(ccdata$year)-1,NA)+
         ggtitle("漁獲量のトレンドとABC")+
         theme(text = element_text(family = font_MAC))+
         theme(legend.position="top")
    }

    graph.component <- list(g.hcr,g.catch)
    graph.combined <- gridExtra::grid.arrange(g.hcr,g.catch,ncol=2,top=stock.name)
    return(list(graph.component=graph.component,graph.combined=graph.combined))
}

#' 3系のHCRを比較するための関数
#'
#' @param res.list calc_abc3の返り値のリスト
#'
#' @export
#'

plot_hcr3 <- function(res.list,stock.name=NULL){
  font_MAC <- "HiraKakuProN-W3"#"Japan1GothicBBB"#
  legend.labels <-c("目標水準案","限界水準案","禁漁水準案")
    if("arglist"%in%names(res.list)) res.list <- list(res.list)

    (g.hcr <- ggplot(data=data.frame(X=c(0,100)), aes(x=X)) +
         theme_bw()+theme_custom()+
         xlab("漁獲量水準 (漁獲量/最大漁獲量, %)")+ylab(str_c("α (ABC年の漁獲量の削減率)"))+
         ggtitle("漁獲管理規則")+
         theme(legend.position="top"))

    if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){
      (g.hcr <- ggplot(data=data.frame(X=c(0,100)), aes(x=X)) +
         theme_bw(base_family = font_MAC)+theme_custom()+
         xlab("漁獲量水準 (漁獲量/最大漁獲量, %)")+ylab(str_c("α (ABC年の漁獲量の削減率)"))+
         ggtitle("漁獲管理規則")+
         theme(legend.position="top")+
         theme(text = element_text(family = font_MAC)))
    }

    for(i in 1:length(res.list)){
        res <- res.list[[i]]

        data_BRP <- tibble(BRP=names(res$BRP),value_obs=res$Obs_BRP,
                           value_ratio=res$BRP)

        BT <- res$arglist$BT
        PL <- res$arglist$PL
        PB <- res$arglist$PB
        tune.par <- res$arglist$tune.par

        (g.hcr <- g.hcr +
             stat_function(fun=type3_func_wrapper,
                           args=list(BT=BT,PL=PL,PB=PB,tune.par=tune.par,type="%"),
                           color="black",size=1,linetype=i)+
             geom_point(aes(x=res$Current_Status[1]*100,y=res$alpha),color=2,size=1)+
             geom_vline(data=data_BRP,mapping=aes(xintercept=value_ratio*100,color=BRP),
                        linetype=i)+
            scale_color_manual(name="",values=rev(c(col.BRP)),labels=rev(c(legend.labels))))
        if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){
        (g.hcr <- g.hcr +
            stat_function(fun=type3_func_wrapper,
                          args=list(BT=BT,PL=PL,PB=PB,tune.par=tune.par,type="%"),
                          color="black",size=1,linetype=i)+
            geom_point(aes(x=res$Current_Status[1]*100,y=res$alpha),color="red",size=1)+
            geom_vline(data=data_BRP,mapping=aes(xintercept=value_ratio*100,color=BRP),
                       linetype=i)+
           scale_color_manual(name="",values=rev(c(col.BRP)),labels=rev(c(legend.labels))))}


    }

    return(g.hcr)
}

#' 2系のHCRを比較するための関数
#'
#' @param res.list calc_abc2の返り値のリスト
#'
#' @export
#'

plot_hcr2 <- function(res.list,stock.name=NULL){
  font_MAC <- "HiraKakuProN-W3"#"Japan1GothicBBB"#
  legend.labels <-c("目標水準案","限界水準案","禁漁水準案")
    if("arglist"%in%names(res.list)) res.list <- list(res.list)

    g.hcr <- ggplot(data=data.frame(X=c(0,120)), aes(x=X)) +
        theme_bw()+theme_custom()+
        ggtitle("漁獲管理規則")+
        xlab("資源量水準(%)")+ylab(str_c("α (ABC年の漁獲量の削減率)"))+
        theme(legend.position="top")

    if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){
      g.hcr <- ggplot(data=data.frame(X=c(0,120)), aes(x=X)) +
        theme_bw(base_family = font_MAC)+theme_custom()+
        ggtitle("漁獲管理規則")+
        xlab("資源量水準(%)")+ylab(str_c("α (ABC年の漁獲量の削減率)"))+
        theme(legend.position="top")+
        theme(text = element_text(family = font_MAC))
    }


    for(i in 1:length(res.list)){
        res <- res.list[[i]]

        data_BRP <- tibble(BRP=names(res$BRP),value_obs=res$Obs_BRP,
                           value_ratio=res$BRP)

        BT <- res$arglist$BT
        PL <- res$arglist$PL
        PB <- res$arglist$PB
        tune.par <- res$arglist$tune.par
        beta <- res$arglist$beta
        g.hcr <- g.hcr +
#            stat_function(fun=type2_func_wrapper,
#                          args=list(BT=BT,PL=0,PB=PB,tune.par=tune.par,AAV=res$AAV,type="%"),
#                       color="gray")+
         stat_function(fun=type2_func_wrapper,
                       args=list(BT=BT,PL=PL,PB=PB,tune.par=tune.par,beta=beta,AAV=res$AAV,type="%"),
                       color="black",size=1,linetype=i)+
            geom_point(aes(x=res$Current_Status[1]*100,y=res$alpha),color=2,size=2)+
            geom_vline(data=data_BRP,mapping=aes(xintercept=value_ratio*100,color=BRP))+
            scale_color_manual(name="",values=rev(c(col.BRP)),label=legend.labels)
            }
    g.hcr
}

#' 2系・3系のABC計算関数．岡村さん作成のプロトタイプ．チェック用に使う．
#'
#' @export
#'

abc_t23_proto1 <- function(
  catch,   # catch timeseries data
  cpue=NULL,   # cpue timeseries data
  BT=0.8,   # initial target level
  tune.par = c(0.5,0.4,0.4), #  tuning parameters: (delta1, delta2, delta3)
  PL=0.7,   #  BL = PL*BT
  PB=0.0,   #  BB = PB*BT
  catch.only=FALSE,    # catch only method
  default=TRUE,
  n.catch=5   #  period for averaging the past catches
){
  #
  # C[t+1] = C[t]*exp(k*(D-BT))
  #

  if (default & (sum(cpue,na.rm=TRUE)==0 | catch.only)) {catch.only <- TRUE; tune.par <- c(2.0,1,0); BT <- 0.1; PL <- 6; PB <- 10; n.catch=3}

  delta1 <- tune.par[1]   # velocity to go to BT
  delta2 <- tune.par[2]   # correction factor when D <= BL
  delta3 <- tune.par[3]   # tuning parameter for updating BT

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

    k <- ifelse(cD < BB, -(delta1+(cD >= BL)*delta2*(BL-cD)/(cD-BB)), -Inf)     #  calculation of k
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

    if (delta3 > 0) AAV <- aav.f(cpue) else AAV <- 0

    k <- ifelse(cD > BB, delta1+(cD <= BL)*delta2*exp(delta3*log(AAV^2+1))*(BL-cD)/(cD-BB), Inf)    #  calculation of k
    abc <- ifelse(cD > BB & cpue[n] > 0, mean(catch[(n-n.catch+1):n],na.rm=TRUE)*exp(k*(cD-BT)), 0)    # calculation of ABC

    Obs_BRP <- c(icum.cpue(BT), icum.cpue(BL), icum.cpue(BB))
    Current_Status <- c(D[n],cpue[n])
    names(Current_Status) <- c("Level","CPUE")
  }

    names(BRP) <- names(Obs_BRP) <- c("Target","Limit","Ban")

    output <- list(BRP=BRP,Obs_BRP=Obs_BRP,Current_Status=Current_Status,catch.only=catch.only,AAV=AAV,tune.par=tune.par,k=k,ABC=abc)

  return(output)
}


#' グラフ出力用のテーマ
#'
#' @export
#'
#'

theme_custom <- function(){
    theme_bw(base_size=12) +
    theme(panel.grid = element_blank(),
          axis.text.x=element_text(size=11,color="black"),
          axis.text.y=element_text(size=11,color="black"),
          axis.line.x=element_line(size= 0.3528),
          axis.line.y=element_line(size= 0.3528))
}
