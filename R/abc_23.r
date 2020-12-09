#' 2,3系のABC計算をするためのパッケージ
#'
#' @import dplyr
#' @import tibble
#' @import stringr
#' @import magrittr
#' @import ggplot2
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
#' @param n.cpue 過去の資源量指標値を平均する年数（デフォルトの値は3）,期間中にnaを含む場合はrm.na
#' @param smooth.cpue ABC算出時に平均化CPUEを使うか（デフォルトはFALSE）
#' @param empir.dist CPUEの分布に経験分布を用いる（デフォルトはFALSEで正規分布を仮定）
#'
#' @examples
#' library(frasyr23)
#' catch <- c(15,20,13,14,11,10,5,10,3,2,1,3)
#' cpue <- c(10,9,8,4,8,11,10,2,3,2,5,2)
#' example_data <- data.frame(year=2001:2012,cpue=cpue,catch=catch,D2alpha=0.1)
#'
#' # 2系
#' example_abc2 <- calc_abc2(example_data,beta=1)
#' graph_abc2 <- plot_abc2(example_abc2,fishseason=0,detABC=1,abc4=FALSE,fillarea=FALSE)
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
  n.cpue=3,   #  period for averaging the past cpues
  smooth.cpue = FALSE,  # option using smoothed cpue
  empir.dist = FALSE,   # option for cupe dist
  beta = 1.0,
  D2alpha = NULL,
  summary_abc = TRUE # 浜辺加筆（'20/07/10）
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
    cum.cpue3 <- function(y,x) pnorm(y,mean(x),sd(x)) # cumulative normal distribution
    cum.cpue4 <- ecdf(cpue) # cumulative empirical distribution

    mean.catch <- mean(ori.catch[(l.catch-n.catch+1):l.catch],na.rm = TRUE)
    mean.cpue <- mean(ori.cpue[(l.cpue-n.cpue+1):l.cpue],na.rm = TRUE)

    for(i in 0:(n.catch-1)){
      if(is.na (ori.catch[l.catch-i])){
        catch.na.warning <- TRUE
      }
    }

    D <- cum.cpue(as.numeric(cpue))              # cumulative probability of cpue
    mD <- attributes(D)$'scaled:center'         # mean of cpue
    sD <- attributes(D)$'scaled:scale'           # standard deviation of cpue
    cD <- D[n]                                   # final depletion
    if(smooth.cpue==TRUE) cD <- cum.cpue3(mean.cpue,cpue)
    if(empir.dist==TRUE){
      cD <- cum.cpue4(cpue[n])
      if(smooth.cpue==TRUE) cD <- mean(cum.cpue4(cpue[n:n-n.cpue+1]))
    }

    icum.cpue <- function(x) sD*qnorm(x,0,1)+mD   # inverse function from D to CPUE
    if(empir.dist==TRUE) icum.cpue <- function(x) as.numeric(quantile(cpue,x))   # inverse function from empirical dist D to CPUE

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
    if(smooth.cpue) alpha <- type2_func(cD,mean.cpue,BT=BT,PL=PL,PB=PB,AAV=AAV,tune.par=tune.par,beta)

    if(is.null(D2alpha)){
      alphafromD <- NULL
    }else{
      alphafromD <- type2_func(D2alpha,cpue[n],BT=BT,PL=PL,PB=PB,AAV=AAV,tune.par=tune.par,beta)
    }
    alphafromD01 <- type2_func(0.1,cpue[n],BT=BT,PL=PL,PB=PB,AAV=AAV,tune.par=tune.par,beta)
    alphafromD005 <- type2_func(0.05,cpue[n],BT=BT,PL=PL,PB=PB,AAV=AAV,tune.par=tune.par,beta)

    ABC <- mean.catch * alpha

    Obs_BRP <- c(icum.cpue(BT), icum.cpue(BL), icum.cpue(BB))
    Obs_percent <- icum.cpue(c(0.05,seq(from=0.1,to=0.9,by=0.1),0.95))
    Obs_percent_even <- icum.cpue(c(0.05,seq(from=0.2,to=0.8,by=0.2),0.95))
    Current_Status <- c(D[n],cpue[n])
    names(Current_Status) <- c("Level","CPUE")
    Recent_Status <- c(cD,mean.cpue)
    names(Recent_Status) <- c("Level","CPUE")

    names(BRP) <- names(Obs_BRP) <- c("Target","Limit","Ban")

    if(summary_abc){ # summary_abc=Tなら以下の結果を自動で書く
    cat("---------------------\n")
    cat(stringr::str_c("Target CPUE value and Level: ",round(Obs_BRP[1],2)," and ", round(BRP[1],2) ,"\n",
                       "Limit CPUE value and Level: ",round(Obs_BRP[2],2)," and ", round(BRP[2],2) ,"\n",
                       "Histrical low CPUE value and Level: ",round(min(cpue),3)," and ", round(min(D),3), "  (", ccdata[ccdata$cpue==min(cpue),]$year, ")", "\n"))
    if(smooth.cpue==FALSE) cat(stringr::str_c("Last year's CPUE value and Level: ",round(cpue[n],3)," and ",round(D[n],3),"\n"))
    else cat(stringr::str_c("Recent ", n.cpue, " year's average CPUE value and Level: ",round(mean.cpue,3)," and ",round(cD,3),"\n"))
            cat(stringr::str_c("AAV of CPUE: ",round(AAV,3),"\n",
                       "alpha: ",round(alpha,3),"\n",
                       "Average catch: ",round(mean.catch,3),"\n",
                       "ABC in ",max(ccdata$year,na.rm=T)+2,": ",round(ABC,3),"\n",
                       "CPUE Level and alpha: 0.1  and  ",round(alphafromD01,3),"\n",
                       "CPUE Level and alpha: 0.05 and  ",round(alphafromD005,3),"\n"))
    if(!is.null(D2alpha)) cat("alpha at CPUE Level=",round(D2alpha,3),": ",round(alphafromD,3),"\n")
    cat("---------------------\n")
    if(isTRUE(catch.na.warning))cat("Warning! Recent n.catch year data contains NA.")
    }

    output <- list(BRP=BRP,Obs_BRP=Obs_BRP,Current_Status=Current_Status,
                   AAV=AAV,tune.par=tune.par,ABC=ABC,arglist=arglist,
                   mean.catch=mean.catch,Obs_percent=Obs_percent,Obs_percent_even=Obs_percent_even,
                   D=D,
                   alpha=alpha,beta=beta,D2alpha=alphafromD)

    if(smooth.cpue==TRUE) output <- list(BRP=BRP,Obs_BRP=Obs_BRP,Current_Status=Recent_Status,
                   AAV=AAV,tune.par=tune.par,ABC=ABC,arglist=arglist,
                   mean.catch=mean.catch,Obs_percent=Obs_percent,Obs_percent_even=Obs_percent_even,
                   D=D,
                   alpha=alpha,beta=beta,D2alpha=alphafromD)
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
#' @param fishseason  X軸のラベルを変更（0なら年、1なら漁期年)
#' @param detABC  次漁期の漁獲量の凡例表記を変更（0ならABC、1なら算定漁獲量、2なら予測値）
#' @param abc4  北海道東部の跨り資源で使用する図を描画（TRUEなら使用、デフォルトはFALSE）
#' @param fillarea  資源量指標値の図にkobeプロットに似た色を塗る（TRUEなら塗る、デフォルトはFALSE）
#' @param cpueunit  資源量指標値の縦軸見出しに追記したい指標値の単位（例えば"（トン/網）"のように指定する）
#' @param leftalign  資源量指標値の時系列の長さが漁獲量に比べて短い時、データが無い範囲の空間を削除する（TRUEなら使用、デフォルトはFALSE）
#'
#' @export
#'

plot_abc2 <- function(res, stock.name=NULL, fishseason=0, detABC=2, abc4=FALSE, fillarea=FALSE, cpueunit="", RP=TRUE, leftalign=FALSE, proposal=TRUE){
    # abc4は北海道東部海域の「跨り資源」で資源量指標値の平均水準・過去最低値を描画する際に使用する。その際、calc_abc2の引数BTは0.5に設定すること。

    # 漁期年/年設定 ----
    ifelse(fishseason==1, year.axis.label <- "漁期年", year.axis.label <- "年")

    # plot
    ccdata <<- res$arglist$ccdata
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
    data_percent_even <- tibble(x=rep(max(years)+2,6),
                           y=res$Obs_percent_even,
                           label=str_c(c(0.05,seq(from=0.2,to=0.8,by=0.2),0.95)*100,"%"))
    font_MAC <- "HiraginoSans-W3"#"Japan1GothicBBB"#

    if(proposal==TRUE){
      legend.labels <-c("目標管理基準値（目標水準）案","限界管理基準値（限界水準）案","禁漁水準案")
    }else{
      legend.labels <-c("目標管理基準値（目標水準）","限界管理基準値（限界水準）","禁漁水準")
    }
    linetype.set <- c("dashed","longdash","solid")
    legend.labels2 <-c(str_c(res$arglist$n.catch,"年平均漁獲量"),"ABC")
    legend.labels2.1 <-c(str_c(res$arglist$n.catch,"年平均漁獲量"),"算定漁獲量")
    legend.labels2.2 <-c(str_c(res$arglist$n.catch,"年平均漁獲量"),paste(max(years)+2,"年",gsub("年","",year.axis.label),"の予測値",sep=""))
    col.BRP.hcr <- col.BRP
    data_BRP_hcr <- tibble(BRP=names(res$BRP),value_obs=res$Obs_BRP, value_ratio=res$BRP)

    # PB=0の時の禁漁水準削除設定 ----
    if(res$BRP[3] == 0) {
      if(proposal==TRUE){
        legend.labels <- c("目標管理基準値（目標水準）案","限界管理基準値（限界水準）案")
      }else{
        legend.labels <- c("目標管理基準値（目標水準）","限界管理基準値（限界水準）")
      }
      linetype.set <- c("22","41")
      if(abc4==TRUE){
        col.BRP <- c("blue","red")
      }else{
        col.BRP <- c("#00533E","#edb918")
      }
      data_BRP2 <- data_BRP
      data_BRP <- tibble(BRP=names(res$BRP[-3]),value_obs=res$Obs_BRP[-3],value_ratio=res$BRP[-3])
    }else{
      if(abc4==TRUE){
        col.BRP <- c("blue","red","orange")
      }else{
        col.BRP <- c("#00533E","#edb918","#C73C2E")
      }
    }

    # ABC決定可能/不可能設定 ----
    if(detABC==1){
      g.catch.title <- "漁獲量のトレンドと算定漁獲量"
      g.catch.abcpoint <- "算定漁獲量"
      legend.labels2 <- legend.labels2.1
    }else if(detABC==2){
      g.catch.title <- "漁獲量の推移と予測値"
      g.catch.abcpoint <- "予測値"
      legend.labels2 <- legend.labels2.2
    }else{
      g.catch.title <- "漁獲量のトレンドとABC"
      g.catch.abcpoint <- "ABC"
    }

    #資源量指標値のトレンド ----

    if(fillarea==TRUE){
      #colfill <- c("olivedrab2", "khaki1", "khaki2", "indianred1")
      colfill <- c("olivedrab2", "khaki1", "white", "white")
    }else{
      colfill <- c("white", "white", "white", "white")
    }

    #minyearを追加しポリゴンをコントロール。最後にxlimで制御
    if(leftalign==TRUE){
      minyears <- min(ccdata[!is.na(ccdata$cpue),]$year)
    }else{
      minyears <- min(years)-2
    }
    g.cpue <- ccdata %>% ggplot() +
      geom_polygon(data=tibble(x=c(minyears,max(years)+4,max(years)+4,minyears), y=c(data_BRP2$value_obs[1],data_BRP2$value_obs[1],max(ccdata$cpue,na.rm=T)*1.05,max(ccdata$cpue,na.rm=T)*1.05)), aes(x=x,y=y), fill=colfill[1]) +
      geom_polygon(data=tibble(x=c(minyears,max(years)+4,max(years)+4,minyears), y=c(data_BRP2$value_obs[2],data_BRP2$value_obs[2],data_BRP2$value_obs[1],data_BRP2$value_obs[1])), aes(x=x,y=y), fill=colfill[2]) +
      geom_polygon(data=tibble(x=c(minyears,max(years)+4,max(years)+4,minyears), y=c(data_BRP2$value_obs[3],data_BRP2$value_obs[3],data_BRP2$value_obs[2],data_BRP2$value_obs[2])), aes(x=x,y=y), fill=colfill[3]) +
      geom_polygon(data=tibble(x=c(minyears,max(years)+4,max(years)+4,minyears), y=c(0,0,data_BRP2$value_obs[3],data_BRP2$value_obs[3])), aes(x=x,y=y), fill=colfill[4]) +
      geom_hline(yintercept=res$Obs_percent_even,color="gray",linetype=2)+
      geom_text(data=data_percent_even, aes(x=x,y=y*1.05,label=label))+
      geom_text(aes(x=max(years)-1,y=min(data_percent_even$y)*0.75,label="(資源量水準)"),size=4)
    if(RP==TRUE){
     g.cpue <- g.cpue +
      geom_hline(data=data_BRP, mapping=aes(yintercept=value_obs, color=rev(col.BRP), linetype=rev(linetype.set)), size = 0.9*1.5)+
      #ggrepel::geom_label_repel(data=data_BRP, mapping=aes(x=min(years)+0.5, y=value_obs+0.5, label=legend.labels), box.padding=0.5, nudge_x=1)+
      scale_linetype_manual(name="", values=rev(c(linetype.set)), labels=rev(c(legend.labels))) +
      scale_color_manual(name="",values=rev(c(col.BRP)),labels=rev(c(legend.labels)))
     }
    g.cpue <- g.cpue +
      geom_path(data=ccdata, aes(x=year,y=cpue),size=1)+
      theme_bw()+ylab(paste("資源量指標値",cpueunit))+xlab(year.axis.label)+
      ylim(0,max(ccdata$cpue,na.rm=T)*1.05)+theme_custom()
    if(RP!=TRUE){
     g.cpue <- g.cpue +
      geom_point(mapping=aes(x=rev(year)[1], y=rev(ccdata$cpue)[1], color="red"),size=4, show.legend =TRUE)+
      scale_color_manual(name="",values="red",labels="直近年の資源量指標値")
     }
    g.cpue <- g.cpue +
      ggtitle("") + theme(legend.position="top", legend.spacing=unit(0.25,'lines'), legend.key.width = unit(2.0, 'lines'),legend.justification=c(1,0))
    if(leftalign==TRUE){
     g.cpue <- g.cpue + xlim(minyears, max(ccdata[!is.na(ccdata$cpue),]$year)+4)
     }

    if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){ # plot 設定 for mac----
      g.cpue <- ccdata %>% ggplot() +
        geom_polygon(data=tibble(x=c(minyears,max(years)+4,max(years)+4,minyears), y=c(data_BRP2$value_obs[1],data_BRP2$value_obs[1],max(ccdata$cpue,na.rm=T)*1.05,max(ccdata$cpue,na.rm=T)*1.05)), aes(x=x,y=y), fill=colfill[1]) +
        geom_polygon(data=tibble(x=c(minyears,max(years)+4,max(years)+4,minyears), y=c(data_BRP2$value_obs[2],data_BRP2$value_obs[2],data_BRP2$value_obs[1],data_BRP2$value_obs[1])), aes(x=x,y=y), fill=colfill[2]) +
        geom_polygon(data=tibble(x=c(minyears,max(years)+4,max(years)+4,minyears), y=c(data_BRP2$value_obs[3],data_BRP2$value_obs[3],data_BRP2$value_obs[2],data_BRP2$value_obs[2])), aes(x=x,y=y), fill=colfill[3]) +
        geom_polygon(data=tibble(x=c(minyears,max(years)+4,max(years)+4,minyears), y=c(0,0,data_BRP2$value_obs[3],data_BRP2$value_obs[3])), aes(x=x,y=y), fill=colfill[4]) +
        geom_hline(yintercept=res$Obs_percent_even,color="gray",linetype=2)+
        geom_text(data=data_percent_even,aes(x=x,y=y*1.05,label=label))+
        geom_text(aes(x=max(years)-1,y=min(data_percent_even$y)*0.75,family=font_MAC,label="(資源量水準)"),size=4)+
     if(RP==TRUE){
      g.cpue <- g.cpue +
        geom_hline(data=data_BRP, mapping=aes(yintercept=value_obs, color=rev(col.BRP), linetype=rev(linetype.set)), size = 0.9*1.5)+
        #ggrepel::geom_label_repel(data=data_BRP, mapping=aes(x=min(years)+0.5, y=value_obs+0.5, label=legend.labels,family = font_MAC), box.padding=0.5, nudge_x=1)+
        scale_linetype_manual(name="", values=rev(c(linetype.set)), labels=rev(c(legend.labels))) +
        scale_color_manual(name="",values=rev(c(col.BRP)),labels=rev(c(legend.labels)))
       }
      g.cpue <- g.cpue +
        geom_path(aes(x=year,y=cpue),size=1)+
        theme_bw()+ylab(paste("資源量指標値",cpueunit))+xlab(year.axis.label)+
        ylim(0,max(ccdata$cpue,na.rm=T)*1.05)+theme_custom()
     if(RP!=TRUE){
      g.cpue <- g.cpue +
        geom_point(mapping=aes(x=rev(year)[1], y=rev(ccdata$cpue)[1], color="red"),size=4, show.legend =TRUE)+
        scale_color_manual(name="",values="red",labels="直近年の資源量指標値")
      }
     g.cpue <- g.cpue +
        ggtitle("")+
        theme(legend.position="top",legend.justification = c(1,0), legend.spacing=unit(0.25,'lines'), legend.key.width = unit(2.0, 'lines')) +
        theme(text = element_text(family = font_MAC))
     if(leftalign==TRUE){
      g.cpue <- g.cpue + xlim(minyears,max(ccdata[!is.na(ccdata$cpue),]$year)+4)
      }
    }

    if(isTRUE(abc4)){
     hanrei_label <- rev(c(paste(min(ccdata[!is.na(ccdata$cpue),]$year),"～",max(ccdata[!is.na(ccdata$cpue),]$year),"年", gsub("年","",year.axis.label), "の平均水準",sep=""),"過去最低値"))  ##OS200702
     g.cpue4 <- ccdata %>% ggplot() +
      geom_hline(yintercept=res$Obs_percent_even,color="gray",linetype=2)+
      geom_text(data=data_percent_even,aes(x=x,y=y*1.05,label=label))+
      geom_text(aes(x=max(years)-1,y=min(data_percent_even$y)*0.75,label="(指標値の水準)"),size=4)+
      geom_hline(data=data_BRP, mapping=aes(yintercept=value_obs[1], color=col.BRP[2], linetype ="twodash"), size = 0.9*1.5, show.legend =TRUE)+
      geom_hline(mapping=aes(yintercept=min(cpue, na.rm=TRUE), color=col.BRP[1], linetype ="longdash"), size = 0.9*1.5, show.legend =TRUE)+
      #ggrepel::geom_label_repel(mapping=aes(x=c(min(years, na.rm=TRUE)+0.5,min(years, na.rm=TRUE)+0.5), y=c(min(cpue, na.rm=TRUE),data_BRP$value_obs[1]), label=rev(c("平均水準","過去最低値"))),
      #                          box.padding=0.5, nudge_x=1)+
      scale_linetype_manual(name="", values=(c("twodash","longdash")), labels=hanrei_label) +
      scale_color_manual(name="",values=rev(c(col.BRP)), labels=hanrei_label)+
      geom_path(aes(x=year,y=cpue),linetype=1,size=1)+
      theme_bw()+ylab(paste("資源量指標値",cpueunit))+xlab(year.axis.label)+
      ylim(0,max(ccdata$cpue,na.rm=T)*1.05)+theme_custom()+
      ggtitle("")+
      theme(legend.position="top",legend.justification = c(1,0), legend.key.width = unit(5, 'lines'))
     if(leftalign==TRUE){
      g.cpue4 <- g.cpue4 + xlim(minyears,max(ccdata[!is.na(ccdata$cpue),]$year)+4)
     }

     if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){ # plot 設定 for mac----
       g.cpue4 <- ccdata %>% ggplot() +
         geom_hline(yintercept=res$Obs_percent_even,color="gray",linetype=2)+
         geom_text(data=data_percent_even,aes(x=x,y=y*1.05,label=label))+
         geom_text(aes(x=max(years)-1,y=min(data_percent_even$y)*0.75,family=font_MAC,label="(指標値の水準)"),size=4)+
         geom_hline(data=data_BRP, mapping=aes(yintercept=value_obs[1], color=col.BRP[2]), size = 0.9*1.5, linetype = 2)+
         geom_hline(mapping=aes(yintercept=min(cpue, na.rm=TRUE), color=col.BRP[1]), size = 0.9*2, linetype = 4)+
         #ggrepel::geom_label_repel(mapping=aes(x=c(min(years, na.rm=TRUE)+0.5,min(years, na.rm=TRUE)+0.5), y=c(min(cpue, na.rm=TRUE),data_BRP$value_obs[1]), label=rev(c("平均水準","過去最低値"))),
         #                          box.padding=0.5, nudge_x=1)+
         scale_color_manual(name="",values=rev(c(col.BRP)),labels=rev(c(paste(min(ccdata[!is.na(ccdata$cpue),]$year),"～",max(ccdata[!is.na(ccdata$cpue),]$year),"の平均水準",sep=""),"過去最低値")))+
         geom_path(aes(x=year,y=cpue),linetype=1,size=1)+
         theme_bw()+ylab(paste("資源量指標値",cpueunit))+xlab(year.axis.label)+
         ggtitle("資源量指標値のトレンド")+
         ylim(0,max(ccdata$cpue,na.rm=T)*1.05)+theme_custom()+
         ggtitle("")+
         theme(legend.position="top",legend.justification = c(1,0))+
         theme(text = element_text(family = font_MAC))
       if(leftalign==TRUE){
         g.cpue4 <- g.cpue4 + xlim(minyears,max(ccdata[!is.na(ccdata$cpue),]$year)+4)
       }
     }
    }


    BT <- res$arglist$BT
    PL <- res$arglist$PL
    PB <- res$arglist$PB
    tune.par <- res$arglist$tune.par
    beta <- res$arglist$beta

    #漁獲管理規則案 HCR ----
    g.hcr <- ggplot(data=data.frame(X=c(0,120)), aes(x=X)) +
          #stat_function(fun=type2_func_wrapper,
          #        args=list(BT=BT,PL=0,PB=PB,tune.par=tune.par,beta=beta,AAV=res$AAV,type="%"),
          #              color="gray")+
          stat_function(fun=type2_func_wrapper,
                        args=list(BT=BT,PL=PL,PB=PB,tune.par=tune.par,beta=beta,AAV=res$AAV,type="%"),
                        color="black",size=1)+
          geom_point(aes(x=res$Current_Status[1]*100,y=res$alpha),color=2,size=4)+
          geom_vline(data=data_BRP,mapping=aes(xintercept=value_ratio*100,color=BRP), size = 0.9*1.5, linetype = linetype.set)+
          ggrepel::geom_label_repel(data=data_BRP,
                                    mapping=aes(x=value_ratio*100, y=c(0.5,0.4), label=legend.labels),
                                    box.padding=0.5)+ #, nudge_y=1
          scale_color_manual(name="",values=rev(c(col.BRP)), guide=FALSE)+#,labels=rev(c(legend.labels)))+
          theme_bw()+theme_custom()+
          ggtitle("")+
          xlab("資源量水準(%)")+ylab(str_c("漁獲量を増減させる係数"))+
          theme(legend.position="top",legend.justification = c(1,0))

    if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){ # plot 設定 for mac----
      g.hcr <- ggplot(data=data.frame(X=c(0,120)), aes(x=X)) +
          #stat_function(fun=type2_func_wrapper,
          #        args=list(BT=BT,PL=0,PB=PB,tune.par=tune.par,beta=beta,AAV=res$AAV,type="%"),
          #                  color="gray")+
          stat_function(fun=type2_func_wrapper,
                            args=list(BT=BT,PL=PL,PB=PB,tune.par=tune.par,beta=beta,AAV=res$AAV,type="%"),
                            color="black",size=1)+
          geom_point(aes(x=res$Current_Status[1]*100,y=res$alpha),color="red",size=4)+
          geom_vline(data=data_BRP,mapping=aes(xintercept=value_ratio*100,color=BRP), size = 0.9*1.5, linetype = linetype.set)+
        ggrepel::geom_label_repel(data=data_BRP,
                                  mapping=aes(x=value_ratio*100, y=1.1, label=legend.labels,family = font_MAC),
                                  box.padding=0.5, nudge_y=1)+
          scale_color_manual(name="",values=rev(c(col.BRP)), guide=FALSE )+ #,labels=rev(c(legend.labels)))+
          theme_bw()+theme_custom()+
          ggtitle("")+
          xlab("資源量水準(%)")+ylab(str_c("漁獲量を増減させる係数"))+
          theme(legend.position="top",legend.justification = c(1,0)) +
          theme(text = element_text(family = font_MAC))
    }

    #漁獲量のトレンドとABC/算定漁獲量 ----
    g.catch <- ccdata %>% ggplot() +
      geom_path(data=data_catch,mapping=aes(x=year,y=catch,color=type),lwd=2)+
      geom_point(data=data_catch,mapping=aes(x=year,y=catch,color=type),lwd=3)+
      #ggrepel::geom_label_repel(data=data_catch,
      #                          mapping=aes(x=max(year)-5, y=catch, label=legend.labels2),
      #                          box.padding=0.5, nudge_y=1)+
      scale_color_manual(name="",values=c(1,2),labels=legend.labels2)+
      # geom_point(data=dplyr::filter(data_catch,type=="ABC"),
      #                    mapping=aes(x=year,y=catch),lwd=2,color=1)+
      #         geom_line(data=dplyr::filter(data_catch,type!="ABC"),
      #                    mapping=aes(x=year,y=catch),lwd=2,color="gray")+
      geom_path(aes(x=year,y=catch),size=1)+
      ylab("漁獲量（トン）")+xlab(year.axis.label)+
      ggtitle("")+
      ylim(0,NA)+ theme_custom()+
      theme(legend.position="top",legend.justification = c(1,0))

    if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){# plot 設定 for mac
      g.catch <- ccdata %>% ggplot() +
        geom_path(data=data_catch,mapping=aes(x=year,y=catch,color=type),lwd=2)+
        geom_point(data=data_catch,mapping=aes(x=year,y=catch,color=type),lwd=3)+
        #ggrepel::geom_label_repel(data=data_catch,
        #                          mapping=aes(x=max(year)-5,y=catch[1],label=legend.labels2,family=font_MAC),
        #                          box.padding=0.5, nudge_y=1)+
        scale_color_manual(name="",values=c("black","red"),labels=legend.labels2)+
        # geom_point(data=dplyr::filter(data_catch,type=="ABC"),
        #                    mapping=aes(x=year,y=catch),lwd=2,color=1)+
        #         geom_line(data=dplyr::filter(data_catch,type!="ABC"),
        #                    mapping=aes(x=year,y=catch),lwd=2,color="gray")+
        geom_path(aes(x=year,y=catch),size=1)+
        ylab("漁獲量（トン）")+xlab(year.axis.label)+
        ggtitle("")+
        ylim(0,NA)+ theme_custom()+
        theme(legend.position="top",legend.justification = c(1,0)) +
        theme(text = element_text(family = font_MAC))
    }

    # 出力設定 ----
    if(isTRUE(abc4)){
      graph.component <- list(g.cpue4,g.cpue,g.hcr,g.catch)
      graph.combined <- gridExtra::grid.arrange(g.cpue4,g.cpue,g.hcr,g.catch,ncol=2,top=stock.name)
      return(list(graph.component=graph.component,graph.combined=graph.combined))
    }else{
      graph.component <- list(g.cpue,g.hcr,g.catch)
      graph.combined <- gridExtra::grid.arrange(g.cpue,g.hcr,g.catch,ncol=3,top=stock.name)
      return(list(graph.component=graph.component,graph.combined=graph.combined))
    }
}


#' 3系のABC計算結果のをプロットするための関数
#'
#' @param res calc_abc3の返り値
#'
#' @export
#'

plot_abc3 <- function(res,stock.name=NULL,fishseason=0,detABC=0,proposal=TRUE){
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
    font_MAC <- "HiraginoSans-W3"#"Japan1GothicBBB"#
    if(proposal==TRUE){
      legend.labels <-c("目標水準案","限界水準案","禁漁水準案")
    }else{
      legend.labels <-c("目標水準","限界水準","禁漁水準")
    }
    linetype.set <- c("22","41","solid")
    legend.labels2 <-c(str_c(res$arglist$n.catch,"年平均漁獲量"),"ABC",rev(c(legend.labels)))
    legend.labels2.1 <-c(str_c(res$arglist$n.catch,"年平均漁獲量"),"算定漁獲量",rev(c(legend.labels)))

    #漁期年/年の設定 ----
    ifelse(fishseason==1, year.axis.label <- "漁期年",year.axis.label <- "年")
    #ABC決定可能/不可能の設定 ----
    if(detABC==1){
      g.catch.title <- ""
      g.catch.abcpoint <- "算定漁獲量"
      legend.labels2 <- legend.labels2.1
    }else{
      g.catch.title <- ""
      g.catch.abcpoint <- "ABC"
    }

    BT <- res$arglist$BT
    PL <- res$arglist$PL
    PB <- res$arglist$PB
    tune.par <- res$arglist$tune.par

    # 漁獲管理規則案HCR ----
    g.hcr <- plot_hcr3(res)

    ## (g.hcr <- ggplot(data=data.frame(X=c(0,120)), aes(x=X)) +
    ##      stat_function(fun=type3_func_wrapper,
    ##                    args=list(BT=BT,PL=PL,PB=PB,tune.par=tune.par,type="%")
    ##                   ,color="black",size=1)+
    ##      geom_point(aes(x=res$Current_Status[1]*100,y=res$alpha),color=2,size=1)+
    ## geom_vline(data=data_BRP,mapping=aes(xintercept=value_ratio*100,color=BRP))+
    ## scale_color_manual(values=rev(col.BRP))+
    ##      theme_bw()+theme_custom()+
    ## xlab("漁獲量水準 (漁獲量/最大漁獲量, %)")+ylab(str_c("alpha (漁獲量の削減率)"))+
    ## ggtitle("漁獲管理規則案")+
    ##      theme(legend.position="top",legend.justification = c(1,0)))

    # 漁獲量トレンドとABC/算定漁獲量 ----
    (g.catch <- ccdata %>% ggplot() +
        geom_path(data=data_catch,mapping=aes(x=year,y=catch,color=type),lwd=2)+
        geom_point(data=data_catch,mapping=aes(x=year,y=catch,color=type),lwd=3)+
        #         geom_point(data=dplyr::filter(data_catch,type=="ABC"),
        #                    mapping=aes(x=year,y=catch),lwd=2,color=1)+
        #         geom_line(data=dplyr::filter(data_catch,type!="ABC"),
        #                    mapping=aes(x=year,y=catch),lwd=2,color="gray")+
        geom_path(aes(x=year,y=catch),size=1)+
        geom_text(data=data_percent,aes(x=x,y=y,label=label))+
        geom_text(aes(x=min(ccdata$year)+2,y=min(data_percent$y)*0.75,label="(漁獲量水準)"),size=4)+
        theme_bw()+ylab("漁獲量（トン）")+xlab(year.axis.label)+theme_custom()+
        geom_hline(data=data_BRP,mapping=aes(yintercept=value_obs,color=BRP), size = 0.9*2, linetype = linetype.set)+
        scale_color_manual(name="",values=c(1,2,rev(col.BRP)),labels=legend.labels2)+
        ylim(0,NA)+xlim(min(ccdata$year)-1,NA)+
        ggtitle(g.catch.title)+
        theme(legend.position="top",legend.justification = c(1,0)))

    if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){ # plot 設定 for mac
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
        ylab("漁獲量（トン）")+xlab(year.axis.label)+theme_custom()+geom_text(data=data_percent,aes(x=x,y=y,label=label),family = font_MAC)+
        geom_text(aes(x=min(ccdata$year)+2,y=min(data_percent$y)*0.75,family=font_MAC,label="(漁獲量水準)"),size=4)+
        geom_hline(data=data_BRP,mapping=aes(yintercept=value_obs,color=BRP), size = 0.9*2, linetype = linetype.set)+
        scale_color_manual(name="",values=col.set,labels=legend.labels2)+
        ylim(0,NA)+xlim(min(ccdata$year)-1,NA)+
        ggtitle(g.catch.title)+
        theme(text = element_text(family = font_MAC))+
        theme(legend.position="top",legend.justification = c(1,0))
    }

    #出力設定 ----
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

plot_hcr3 <- function(res.list,stock.name=NULL,proposal=TRUE){
  font_MAC <- "HiraginoSans-W3"#"Japan1GothicBBB"#
  if(proposal==TRUE){
    legend.labels.hcr <-c("目標水準案","限界水準案","禁漁水準案")
  }else{
    legend.labels.hcr <-c("目標水準","限界水準","禁漁水準")
  }
  linetype.set <- c("22","41","solid")
  if("arglist"%in%names(res.list)) res.list <- list(res.list)

    (g.hcr <- ggplot(data=data.frame(X=c(0,100)), aes(x=X)) +
       theme_bw()+theme_custom()+
       xlab("漁獲量水準 (漁獲量/最大漁獲量, %)")+ylab(str_c("漁獲量を増減させる係数"))+
       ggtitle("")+
       theme(legend.position="top",legend.justification = c(1,0)))

    if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){ # plot setting for mac----
      (g.hcr <- ggplot(data=data.frame(X=c(0,100)), aes(x=X)) +
         theme_bw(base_family = font_MAC)+theme_custom()+
         xlab("漁獲量水準 (漁獲量/最大漁獲量, %)")+ylab(str_c("漁獲量を増減させる係数"))+
         ggtitle("")+
         theme(legend.position="top",legend.justification = c(1,0))+
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
          geom_vline(data=data_BRP,mapping=aes(xintercept=value_ratio*100,color=BRP), size = 0.9*2, linetype = linetype.set)+
          ggrepel::geom_label_repel(data=data_BRP,
                                    mapping=aes(x=value_ratio*100, y=1.1, label=legend.labels.hcr),
                                    box.padding=0.5, nudge_y=1)+
          scale_color_manual(name="",values=rev(c(col.BRP)),guide=FALSE)) #labels=rev(c(legend.labels.hcr))))
      if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){ # plot setting for mac----
        (g.hcr <- g.hcr +
           stat_function(fun=type3_func_wrapper,
                         args=list(BT=BT,PL=PL,PB=PB,tune.par=tune.par,type="%"),
                         color="black",size=1,linetype=i)+
           geom_point(aes(x=res$Current_Status[1]*100,y=res$alpha),color="red",size=1)+
           geom_vline(data=data_BRP,mapping=aes(xintercept=value_ratio*100,color=BRP), size = 0.9*2, linetype = linetype.set)+
           ggrepel::geom_label_repel(data=data_BRP,
                                     mapping=aes(x=value_ratio*100, y=1.1, label=legend.labels.hcr, family=font_MAC),
                                     box.padding=0.5, nudge_y=1)+
           scale_color_manual(name="",values=rev(c(col.BRP)),guide=FALSE)) #labels=rev(c(legend.labels.hcr))))}
      }
    }

  return(g.hcr)
}

#' 2系のHCRを比較するための関数
#'
#' @param res.list calc_abc2の返り値のリスト
#'
#' @export
#'

plot_hcr2 <- function(res.list,stock.name=NULL,proposal=TRUE){
  font_MAC <- "HiraginoSans-W3"#"Japan1GothicBBB"#
  if(proposal==TRUE){
    legend.labels.hcr <-c("目標管理基準値（目標水準）案","限界管理基準値（限界水準）案","禁漁水準案")
  }else{
    legend.labels.hcr <-c("目標管理基準値（目標水準）","限界管理基準値（限界水準）","禁漁水準")
  }
  linetype.set <- c("22","41","solid")
  if("arglist"%in%names(res.list)) res.list <- list(res.list)

      g.hcr <- ggplot(data=data.frame(X=c(0,120)), aes(x=X)) +
        theme_bw()+theme_custom()+
        ggtitle("")+
        xlab("資源量水準(%)")+ylab(str_c("漁獲量を増減させる係数"))+
        theme(legend.position="top",legend.justification = c(1,0))
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
          geom_point(aes(x=res$Current_Status[1]*100,y=res$alpha),color=2,size=4)
        }
      g.hcr <- g.hcr + geom_vline(data=data_BRP,mapping=aes(xintercept=value_ratio*100,color=BRP), size = 0.9, linetype = linetype.set)+
          ggrepel::geom_label_repel(data=data_BRP,
                                    mapping=aes(x=value_ratio*100, y=c(1.1,1.0,0.9), label=legend.labels.hcr),
                                    box.padding=0.5)+
          scale_color_manual(name="",values=rev(c(col.BRP)),guide=FALSE) #label=rev(legend.labels.hcr))


      if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){
        g.hcr <- ggplot(data=data.frame(X=c(0,120)), aes(x=X)) +
          theme_bw(base_family = font_MAC)+theme_custom()+
          ggtitle("")+
          xlab("資源量水準(%)")+ylab(str_c("漁獲量を増減させる係数"))+
          theme(legend.position="top",legend.justification = c(1,0))+
          theme(text = element_text(family = font_MAC))

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
          geom_point(aes(x=res$Current_Status[1]*100,y=res$alpha),color=2,size=4)
      }
      g.hcr <- g.hcr + geom_vline(data=data_BRP,mapping=aes(xintercept=value_ratio*100,color=BRP), size = 0.9, linetype = linetype.set)+
                       ggrepel::geom_label_repel(data=data_BRP,
                       mapping=aes(x=value_ratio*100, y=c(1.1,1.0,0.9), label=legend.labels.hcr,family=font_MAC),
                       box.padding=0.5)+
                       scale_color_manual(name="",values=rev(c(col.BRP)),guide=FALSE) #label=rev(legend.labels.hcr))
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
