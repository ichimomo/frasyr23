#' 2,3系のABC計算をするためのパッケージ
#'
#' @import dplyr
#' @import tibble
#' @import stringr
#' @import magrittr
#' @import ggrepel
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
#' @param smooth.dist CPUEの分布に平滑化（n.cpue年移動平均）したCPUEをつかって正規分布に当てる
#' @param empir.dist CPUEの分布に経験分布を用いる（デフォルトはFALSEで正規分布を仮定）
#' @param simple.empir CPUEの分布に経験分布を用いるときに旧2系的に分布を仮定する（デフォルトはFALSE）empir.dist==Tとした上で追加。
#' @param BTyear 管理目標水準を計算するときのCPUE時系列最終年を手動で決める（デフォルトはNULLでccdata$cpueの最終年）
#' @param timelag0 入力データ最終年の翌年のABCとして算出できるケース（デフォルトはFALSE）
#' @param resp 変動緩和措置をとるとき、前年漁獲量を1としたときの指数
#' @examples
#' library(frasyr23)
#' catch <- c(15,20,13,14,11,10,5,10,3,2,1,3)
#' cpue <- c(10,9,8,4,8,11,10,2,3,2,5,2)
#' example_data <- data.frame(year=2001:2012,cpue=cpue,catch=catch)
#'
#' # 2系
#' example_abc2 <- calc_abc2(example_data,beta=1)
#' graph_abc2 <- plot_abc2(example_abc2,fishseason=0,detABC=1,x_right_space=5,abc4=FALSE,fillarea=FALSE)
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
    smooth.dist = FALSE,  # option for cpue dist using smoothed cpue
    empir.dist = FALSE,   # option for cpue dist
    simple.empir = FALSE, # option for empirical cpue dist
    beta = 1.0,
    D2alpha = NULL,
    BTyear = NULL,
    timelag0 = FALSE,
    resp = NULL,
    summary_abc = TRUE # 浜辺加筆（'20/07/10）
){
  argname <- ls() # 引数をとっておいて再現できるようにする
  arglist <- lapply(argname,function(xx) eval(parse(text=xx)))
  names(arglist) <- argname


  if(is.null(BTyear)){
    cpue <- ccdata$cpue
    ori.cpue <- cpue
    cpue <- cpue[!is.na(ccdata$cpue)]

  } else{
    if(BTyear > max(ccdata$year)) stop("BTyear year must be set less than max(ccdata$year)!")
    if(BTyear < min(ccdata$year)) stop("BTyear year must be set larger than min(ccdata$year)!")
    ccdata_fixedBT <- ccdata[which(ccdata$year <= BTyear),]

    cpue <- ccdata_fixedBT$cpue
    ori.cpue <- cpue
    cpue <- cpue[!is.na(ccdata_fixedBT$cpue)]
    target.cpue <- ccdata$cpue[nrow(ccdata)]
  }

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

  smoothed.cpue <- c()
  for(i in n.cpue:l.cpue){
    smoothed.cpue <- cbind(smoothed.cpue,mean(ori.cpue[(i-n.cpue+1):i],na.rm = TRUE))
  }

  cum.cpue <- function(x) pnorm(scale(x),0,1) # cumulative normal distribution
  cum.cpue2 <- function(x) pnorm(x,mean(x),sd(x)) # cumulative normal distribution
  cum.cpue3 <- function(y,x) pnorm(y,mean(x),sd(x)) # cumulative normal distribution
  cum.cpue4 <- ecdf(cpue) # cumulative empirical distribution

  mean.catch <- mean(ori.catch[(l.catch-n.catch+1):l.catch],na.rm = TRUE)

  mean.cpue <- mean(ori.cpue[(l.cpue-n.cpue+1):l.cpue],na.rm = TRUE) # this mean.cpue calculates BT based on BTyear
  mean.cpue.current <- mean(ccdata$cpue[(length(ccdata$cpue)-n.cpue+1):length(ccdata$cpue)],na.rm = TRUE) # this mean.cpue is the most recent cpue (max(ccdata$year))


  for(i in 0:(n.catch-1)){
    if(is.na (ori.catch[l.catch-i])){
      catch.na.warning <- TRUE
    }
  }


  D <- cum.cpue(as.numeric(cpue))              # cumulative probability of cpue
  mD <- attributes(D)$'scaled:center'         # mean of cpue
  sD <- attributes(D)$'scaled:scale'           # standard deviation of cpue
  if(is.null(BTyear)){
    cD <- D[n]                                   # final depletion
  }else{
    cD <- pnorm((target.cpue-mD)/sD,0,1)
    #(CPUE - mD)/sD = qnorm(x,0,1)
    # final depletion
  }


  if(smooth.dist==TRUE){
    D <- cum.cpue(as.numeric(smoothed.cpue))              # cumulative probability of cpue
    mD <- attributes(D)$'scaled:center'         # mean of cpue
    sD <- attributes(D)$'scaled:scale'           # standard deviation of cpue
    if(is.null(BTyear)){
      cD <- D[length(smoothed.cpue)]                                   # final depletion
    }else{
      cD <- pnorm((mean.cpue.current-mD)/sD,0,1)
    }
  }

  if(is.null(BTyear)){
    if(smooth.cpue==TRUE) cD <- cum.cpue3(mean.cpue,cpue)
    if(empir.dist==TRUE){
      cD <- cum.cpue4(cpue[n])
      D <- cum.cpue4(cpue)
      if(smooth.cpue==TRUE) cD <- mean(cum.cpue4(cpue[n:n-n.cpue+1]))
      if(simple.empir ==TRUE){
        cD <- simple_ecdf(cpue,cpue[n])
        D <- simple_ecdf_seq(cpue)
        if(smooth.cpue==TRUE){
          tmp.recent.cpue<-c()
          for(i in 1:n.cpue){
            tmp.recent.cpue<-c(tmp.recent.cpue,simple_ecdf(cpue,cpue[n-i+1]))
          }
          cD <- mean(tmp.recent.cpue)
        }
        if(cD <= min(D)) cat("alpha = 0 because current cpue is min(cpue)\n")
      }
    }
  }else{
    if(smooth.cpue==TRUE) cD <- cum.cpue3(mean.cpue.current,cpue)
    if(empir.dist==TRUE){
      cD <- cum.cpue4(target.cpue)
      D <- cum.cpue4(cpue)
      if(smooth.cpue==TRUE) cD <- mean(cum.cpue4(ccdata$cpue[n:n-n.cpue+1]))
      if(simple.empir ==TRUE){
        D <- simple_ecdf_seq(cpue)
        if(target.cpue > min(D)) {
          cD <- simple_ecdf(cpue,target.cpue)
          if(smooth.cpue) {
            tmp.recent.cpue<-c()
            for(i in 1:n.cpue){
              tmp.recent.cpue<-c(tmp.recent.cpue,ifelse(cpue[n-i+1]>min(D),simple_ecdf(cpue.ori,cpue[n-i+1]),min(D)))
              if(cpue[n-i+1]<min(D)) cat(ccdata$year[n-i+1],"years' cpue is replaced min(cpue) because it is less than any cpue(year <=BTyear)")
            }
            cD <- mean(tmp.recent.cpue)
          }
        }else cD <- min(D)

        if(cD <= min(D)) cat("alpha <= 0 because current cpue is min(cpue) or less than any cpue (year <= BTyear) \n")
      }

    }
  }

  icum.cpue <- function(x) sD*qnorm(x,0,1)+mD   # inverse function from D to CPUE
  if(empir.dist==TRUE) {icum.cpue <- function(x) as.numeric(quantile(cpue,x))   # inverse function from empirical dist D to CPUE
  if(simple.empir==TRUE) icum.cpue <- function(x) inv_simple_ecdf(cpue,x) # inverse function from simple empirical dist D to CPUE
  }

  if (delta3 > 0){
    if(AAV=="auto"){
      AAV <- aav.f(ccdata$cpue)
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

  if(is.null(BTyear)){
    if(!(empir.dist)) alpha <- type2_func(cD,cpue[n],BT=BT,PL=PL,PB=PB,AAV=AAV,tune.par=tune.par,beta)
    else alpha <- type2_func_empir(cD,cpue,simple=simple.empir,BT=BT,PL=PL,PB=PB,AAV=AAV,tune.par=tune.par,beta)
  }else{
    if(!(empir.dist)) alpha <- type2_func(cD,target.cpue,BT=BT,PL=PL,PB=PB,AAV=AAV,tune.par=tune.par,beta)
    else alpha <- type2_func_empir(cD,cpue,simple=simple.empir,BT=BT,PL=PL,PB=PB,AAV=AAV,tune.par=tune.par,beta)
  }
  if(smooth.cpue) {
    if(!is.null(BTyear)){
      if(!(empir.dist)) alpha <- type2_func(cD,mean.cpue,BT=BT,PL=PL,PB=PB,AAV=AAV,tune.par=tune.par,beta)
    }else if(!(empir.dist)) alpha <- type2_func(cD,mean.cpue.current,BT=BT,PL=PL,PB=PB,AAV=AAV,tune.par=tune.par,beta)
    else alpha <- type2_func_empir(cD,smooth.cpue,simple=simple.empir,BT=BT,PL=PL,PB=PB,AAV=AAV,tune.par=tune.par,beta)
  }


  if(is.null(D2alpha)){
    alphafromD <- NULL
  }else{

    if(is.null(BTyear)){
      if(!(empir.dist)) alphafromD <- type2_func(D2alpha,cpue[n],BT=BT,PL=PL,PB=PB,AAV=AAV,tune.par=tune.par,beta)
      else alphafromD<- type2_func_empir(D2alpha,cpue,simple=simple.empir,BT=BT,PL=PL,PB=PB,AAV=AAV,tune.par=tune.par,beta)
    }else{
      if(!(empir.dist)) alphafromD <- type2_func(D2alpha,target.cpue,BT=BT,PL=PL,PB=PB,AAV=AAV,tune.par=tune.par,beta)
      else alphafromD<- type2_func_empir(D2alpha,cpue,simple=simple.empir,BT=BT,PL=PL,PB=PB,AAV=AAV,tune.par=tune.par,beta)
    }
  }

  if(is.null(BTyear)){
    if(!(empir.dist)) alphafromD01 <- type2_func(0.1,cpue[n],BT=BT,PL=PL,PB=PB,AAV=AAV,tune.par=tune.par,beta)
    else alphafromD01 <- type2_func_empir(0.1,cpue,simple=simple.empir,BT=BT,PL=PL,PB=PB,AAV=AAV,tune.par=tune.par,beta)
    if(!(empir.dist)) alphafromD005 <- type2_func(0.05,cpue[n],BT=BT,PL=PL,PB=PB,AAV=AAV,tune.par=tune.par,beta)
    else alphafromD005 <- type2_func_empir(0.05,cpue,simple=simple.empir,BT=BT,PL=PL,PB=PB,AAV=AAV,tune.par=tune.par,beta)
  }else{
    if(!(empir.dist)) alphafromD01 <- type2_func(0.1,target.cpue,BT=BT,PL=PL,PB=PB,AAV=AAV,tune.par=tune.par,beta)
    else alphafromD01 <- type2_func_empir(0.1,cpue,simple=simple.empir,BT=BT,PL=PL,PB=PB,AAV=AAV,tune.par=tune.par,beta)
    if(!(empir.dist)) alphafromD005 <- type2_func(0.05,target.cpue,BT=BT,PL=PL,PB=PB,AAV=AAV,tune.par=tune.par,beta)
    else alphafromD005 <- type2_func_empir(0.05,cpue,simple=simple.empir,BT=BT,PL=PL,PB=PB,AAV=AAV,tune.par=tune.par,beta)

  }
  alphafromD01 <- type2_func(0.1,cpue[n],BT=BT,PL=PL,PB=PB,AAV=AAV,tune.par=tune.par,beta)
  alphafromD005 <- type2_func(0.05,cpue[n],BT=BT,PL=PL,PB=PB,AAV=AAV,tune.par=tune.par,beta)

  # ABC
  ABC <- mean.catch * alpha
  resp_flag<-0
  if(!is.null(resp)) { # resp != NULL
    if(!is.na(catch[n])){
      if( ABC > (1+resp)*catch[n] ) {
        ABC <- catch[n]*(1+resp)
        resp_flag<-1
      }
      if( ABC < (1-resp)*catch[n] ) {
        ABC <- catch[n]*(1-resp)
        resp_flag<-2
      }
    }else{ # if latest catch == NA
      if( ABC > (1+resp)*catch[n-1] ) {
        ABC <- catch[n-1]*(1+resp)
        resp_flag<-1
      }
      if( ABC < (1-resp)*catch[n-1] ) {
        ABC <- catch[n-1]*(1-resp)
        resp_flag<-2
      }
    }
  }

  Obs_BRP <- c(icum.cpue(BT), icum.cpue(BL), icum.cpue(BB))
  Obs_percent <- icum.cpue(c(0.05,seq(from=0.1,to=0.9,by=0.1),0.95))
  Obs_percent_even <- icum.cpue(c(0.05,seq(from=0.2,to=0.8,by=0.2),0.95))

  if(is.null(BTyear)){
    Current_Status <- c(D[n],cpue[n])
  }else{ # BTyear!=NULLではHC level計算年の状態
    Current_Status <- c(cD,target.cpue)
  }
  names(Current_Status) <- c("Level","CPUE")
  Recent_Status <- c(cD,mean.cpue)


  names(Recent_Status) <- c("Level","CPUE")

  names(BRP) <- names(Obs_BRP) <- c("Target","Limit","Ban")

  if(summary_abc){ # summary_abc=Tなら以下の結果を自動で書く
    cat("---------------------\n")
    if(is.null(BTyear)){
      cat(stringr::str_c("Target CPUE value and Level: ",round(Obs_BRP[1],2)," and ", round(BRP[1],2) ,"\n",
                         "Limit CPUE value and Level: ",round(Obs_BRP[2],2)," and ", round(BRP[2],2) ,"\n",
                         "Histrical low CPUE value and Level: ",round(min(cpue),3)," and ", round(min(D),3), "  (", ccdata[ccdata$cpue==min(cpue),]$year, ")", "\n"))
    }else{
      cat(stringr::str_c("Target CPUE value and Level: ",round(Obs_BRP[1],2)," and ", round(BRP[1],2) ,"\n",
                         "Limit CPUE value and Level: ",round(Obs_BRP[2],2)," and ", round(BRP[2],2) ,"\n",
                         "Histrical low CPUE value and Level(",min(ccdata_fixedBT$year),"-",max(ccdata_fixedBT$year),"): ",round(min(cpue),3)," and ", round(min(D),3), "  (", ccdata_fixedBT[ccdata_fixedBT$cpue==min(cpue),]$year, ")", "\n"))

    }
    if(smooth.cpue==FALSE && smooth.dist==FALSE){
      if(is.null(BTyear)) cat(stringr::str_c("Last year's CPUE value and Level: ",round(cpue[n],3)," and ",round(D[n],3),"\n"))
      else cat(stringr::str_c("HC-level calculation year(",BTyear,")'s CPUE value and Level: ",round(cpue[length(cpue)],3)," and ",round(cD,3),"\n"))
    }
    else {
      if(is.null(BTyear)) cat(stringr::str_c("Recent ", n.cpue, " years' average CPUE value and Level: ",round(mean.cpue.current,3)," and ",round(cD,3),"\n"))
      else cat(stringr::str_c("Recent ", n.cpue, " years' (",BTyear-n.cpue+1,"-",BTyear,") average CPUE value and Level: ",round(mean.cpue,3)," and ",round(cD,3),"\n"))
    }
    cat(stringr::str_c("AAV of CPUE: ",round(AAV,3),"\n",
                       "alpha: ",round(alpha,3),"\n",
                       "Average catch: ",round(mean.catch,3),"\n"))
    if(!timelag0){
      cat(stringr::str_c("ABC in ",max(ccdata$year,na.rm=T)+2,": ",round(ABC,3),"\n"))
    }else{
      cat(stringr::str_c("ABC in ",max(ccdata$year,na.rm=T)+1,": ",round(ABC,3),"\n"))
    }
    if(resp_flag==1) cat(stringr::str_c("ABC was replaced by ",(1+resp)*100,"% of the Latest catch \n"))
    if(resp_flag==2) cat(stringr::str_c("ABC was replaced by ",(1-resp)*100,"% of the Latest catch \n"))
    cat(stringr::str_c("CPUE Level and alpha: 0.1  and  ",round(alphafromD01,3),"\n",
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

  if(smooth.cpue==TRUE || smooth.dist==TRUE) output <- list(BRP=BRP,Obs_BRP=Obs_BRP,Current_Status=Recent_Status,
                                                            AAV=AAV,tune.par=tune.par,ABC=ABC,arglist=arglist,
                                                            mean.catch=mean.catch,Obs_percent=Obs_percent,Obs_percent_even=Obs_percent_even,
                                                            D=D,
                                                            alpha=alpha,beta=beta,D2alpha=alphafromD)
  return(output)
}

#' 2系のCPUEデータに対してABCを返す関数
#'
#' @param cD cpueの値
#' @param BT 目標水準
#' @param PL 目標水準に対する限界水準の割合
#' @param PB 目標水準に対する禁漁水準の割合
#' @param AAV AAV
#' @param tune.par チューニングパラメタ
#' @param beta 保守的漁獲の割合
#'
#' @export
#'
#'
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
  if(cD >= BL) alpha <- exp(delta1*(cD-BT))
  assertthat::assert_that(is.numeric(alpha))
  return(alpha*beta)
  # cpue.nは必要か？
  #    k <- ifelse(cD > BB, delta1+(cD <= BL)*delta2*exp(delta3*log(AAV^2+1))*(BL-cD)/(cD-BB), Inf)    #  calculation of k
  #    ifelse(cD > BB & cpue.n > 0, exp(k*(cD-BT)), 0)    # calculation of ABC
}

#' 2系資源計算を経験分布で計算する時のCPUEデータに対してABCを返す関数
#'
#' @param cD cpueの値
#' @param cpue cpue時系列
#' @param simple 旧2系的な経験分布にする場合TRUE
#' @param BT 目標水準
#' @param PL 目標水準に対する限界水準の割合
#' @param PB 目標水準に対する禁漁水準の割合
#' @param AAV AAV
#' @param tune.par チューニングパラメタ
#' @param beta 保守的漁獲の割合
#'
#' @export
#'
#'
type2_func_empir <- function(cD,cpue,simple=FALSE,BT=0.8,PL=0.7,PB=0,AAV=0.4,tune.par=c(0.5,0.5,0.4),beta=1.0){
  delta1 <- tune.par[1]   # velocity to go to BT
  delta2 <- tune.par[2]   # correction factor when D <= BL
  delta3 <- tune.par[3]   # tuning parameter for updating BT
  BT <- BT      # Btarget
  BL <- PL*BT      # Blimit
  BB <- PB*BT      # Bban

  if(!simple){
    cum.cpue <- ecdf(cpue)
    cpue.order <- sort(unique(cpue))
    cpue.prob <-cum.cpue(cpue.order)

  }else{
    cum.cpue <- simple_ecdf_seq(cpue)
    cpue.prob <- sort(unique(cum.cpue))
  }
  trans_empir_prob<-function(prob.seq,cpue.prob){
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

  cD <- trans_empir_prob(cD,cpue.prob)

  #if(simple.empir & cD==0) cD <- (max(cpue)-min(cpue))/1000

  if(cD <= BB) alpha <- 0
  if(BB < cD & cD < BL){
    k <- delta1 + delta2* exp(delta3*log(AAV^2+1)) * (BL-cD)/(cD-BB)
    alpha <- exp(k*(cD-BT))
  }
  if(cD >= BL) alpha <- exp(delta1*(cD-BT))
  assertthat::assert_that(is.numeric(alpha))
  return(alpha*beta)
  # cpue.nは必要か？
  #    k <- ifelse(cD > BB, delta1+(cD <= BL)*delta2*exp(delta3*log(AAV^2+1))*(BL-cD)/(cD-BB), Inf)    #  calculation of k
  #    ifelse(cD > BB & cpue.n > 0, exp(k*(cD-BT)), 0)    # calculation of ABC
}

#' 2系のCPUEの確率点に対して連続的にABCを返す関数
#'
#' @export
type2_func_wrapper <- function(DL,type=NULL,...){
  if(type=="%") DL <- DL/100
  purrr::map_dbl(DL,type2_func,...)
}

#' 2系を経験分布で計算する時、CPUEの確率点に対して連続的にABCを返す関数
#'
#' @export
type2_func_empir_wrapper <- function(DL,cpue,simple,type=NULL,...){
  if(type=="%") DL <- DL/100
  purrr::map_dbl(DL,type2_func_empir,cpue=cpue,simple=simple,...)
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

simple_ecdf <- function(cpue, x){
  percent <- (x-min(cpue))/(max(cpue)-min(cpue))
  if(percent<0) percent <- 0
  return(percent)
}

simple_ecdf_seq<-function(cpue){
  cum.cpue <-c()
  i<-1
  while(i<=length(cpue)){
    cum.cpue.tmp <- simple_ecdf(cpue,cpue[i])
    cum.cpue <- c(cum.cpue,cum.cpue.tmp)
    i <- i +1
  }
  return(cum.cpue)
}

inv_simple_ecdf <- function(cpue, x){
  inv.value <- x*(max(cpue)-min(cpue))+min(cpue)
  return(inv.value)
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
#' @param x_right_space  資源量指標値のトレンドの図で、水準の%やラベルを書き込むの幅の設定（デフォルトは5年分）
#' @param abc4  北海道東部の跨り資源で使用する図を描画（TRUEなら使用、デフォルトはFALSE）
#' @param fillarea  資源量指標値の図にkobeプロットに似た色を塗る（TRUEなら塗る、デフォルトはFALSE）
#' @param cpueunit  資源量指標値の縦軸見出しに追記したい指標値の単位（例えば"（トン/網）"のように指定する）
#' @param catchunit 漁獲量トレンドの縦軸見出しに追記したい単位（例えば"（万トン）"のように指定する、デフォルトは"（トン）"）
#' @param catchdividedegit 漁獲量トレンドの縦軸を任意の桁で割る（例えば"catchdividedegit=2"で100で割るように指定する、デフォルトはNULLでそのまま）
#' @param leftalign  資源量指標値の時系列の長さが漁獲量に比べて短い時、データが無い範囲の空間を削除する（TRUEなら使用、デフォルトはFALSE）
#' @param RP  資源量指標値/年のプロットでReference Point（目標・限界管理基準線）を載せる・載せない（デフォルトはTRUE、FALSEでは直近年の資源量指標値をポイントでハイライトする）
#' @param hcrhscale HCRのプロットで縦軸の目盛幅をいくつ刻むか（sparseで0.5刻み、middleで0.25刻み、denseで0.2刻み）
#' @param hcrhline HCRのプロットで縦軸の漁獲量を増減させる係数に補助線を入れる（noneで補助線なし、oneでy=1、hscaleで縦軸目盛に補助線をひく）
#' @param plotexactframe HCRのプロットで縦軸横軸の0をプロット枠とあわせるか、余裕を持たせるか（デフォルトはFALSEで枠いっぱい）
#' @param ignore_naCatch_point ABC算出に使う最近年の漁獲量にNAが入っている場合、表示上NAとなる年のポイントと年数を引く
#'
#' @export
#'

plot_abc2 <- function(res, stock.name=NULL, fishseason=0, detABC=2, abc4=FALSE, fillarea=FALSE, cpueunit="", catchunit="(トン)", catchdividedegit=NULL, RP=TRUE, leftalign=FALSE, proposal=TRUE, hcrdist=FALSE, BThcr=FALSE,hcrhline="none",hcrhscale="middle",plotexactframe=FALSE,ignore_naCatch_point=FALSE,latest_Catch_na=FALSE){
  # abc4は北海道東部海域の「跨り資源」で資源量指標値の平均水準・過去最低値を描画する際に使用する。その際、calc_abc2の引数BTは0.5に設定すること。

  # 漁期年/年設定 ----
  ifelse(fishseason==1, year.axis.label <- "漁期年", year.axis.label <- "年")
  # Setting cpue dist
  empir.dist <- res$arglist$empir.dist
  simple.empir <- res$arglist$simple.empir
  smooth.cpue <- res$arglist$smooth.cpue
  smooth.dist <- res$arglist$smooth.dist
  # plot

  ccdata <- res$arglist$ccdata
  ccdata_fixedBT <- res$arglist$ccdata
  mean.catch <-res$mean.catch

  if(is.null(catchdividedegit)){
    ABC <- res$ABC
  }
  else{
    ccdata$catch <- res$arglist$ccdata$catch/(10^catchdividedegit)
    ccdata_fixedBT$catch <- res$arglist$ccdata$catch/(10^catchdividedegit)
    ABC <- res$ABC/(10^catchdividedegit)
    mean.catch <-res$mean.catch/(10^catchdividedegit)
  }
  BTyear <- res$arglist$BTyear
  if(is.null(BTyear) && BThcr==TRUE) stop("BThcr option works if BTyear is not NULL. See the argments in calc_abc2./n")
  if(!is.null(res$arglist$BTyear)) ccdata_fixedBT <- ccdata[which(ccdata$year <= BTyear),]
  n.catch <- res$arglist$n.catch
  #latestCatchna <- res$arglist$latestCatchna
  #if(!is.null(latestCatchna)) ignore_naCatch_point<-TRUE
  years <- ccdata$year
  last.year <- rev(years)[1]

  BT <- res$arglist$BT
  PL <- res$arglist$PL
  PB <- res$arglist$PB
  tune.par <- res$arglist$tune.par
  beta <- res$arglist$beta

  catch.abc.na<-0
  if(ignore_naCatch_point || latest_Catch_na){
    mean.catch.abc <- ccdata$catch[(length(ccdata$catch)-n.catch+1):length(ccdata$catch)]
    catch.abc.na <- sum(as.numeric(is.na(mean.catch.abc)))
    if(prod(!is.na(mean.catch.abc))) stop("ignore_naCatch_point option (latest_Catch_na option) works if catch[lastyear-n.catch+1:lastyear] contains NA.")
  }
  if(latest_Catch_na && res$arglist$timelag0==FALSE) stop("latest_Catch_na option works if timela0 option for calc_abc2 is TRUE.")

  if(BThcr==FALSE){
    if(!res$arglist$timelag0){ # 2年後ABC算出
      data_catch <- tibble(year=c((last.year-res$arglist$n.catch+1):last.year,last.year+2),
                           catch=c(rep(mean.catch,res$arglist$n.catch),ABC),
                           type=c(rep(str_c(res$arglist$n.catch-catch.abc.na,"年平均漁獲量"),n.catch),"ABC"))
    }else{ # 1年後ABC算出
      if(latest_Catch_na){ #最新年CatchがNAの時
        data_catch <- tibble(year=c((last.year-res$arglist$n.catch+1):(last.year-1),last.year+1),
                             catch=c(rep(mean.catch,res$arglist$n.catch-1),ABC),
                             type=c(rep(str_c(res$arglist$n.catch-catch.abc.na,"年平均漁獲量"),n.catch-1),"ABC"))
      }else{
        data_catch <- tibble(year=c((last.year-res$arglist$n.catch+1):last.year,last.year+1),
                             catch=c(rep(mean.catch,res$arglist$n.catch),ABC),
                             type=c(rep(str_c(res$arglist$n.catch-catch.abc.na,"年平均漁獲量"),n.catch),"ABC"))

      }
    }
  }
  else{
    res.nullBTyear <- calc_abc2(ccdata=res$arglist$ccdata,BT=BT,PL=PL,PB=PB,tune.par = tune.par, AAV=res$arglist$AAV,n.catch=res$arglist$n.catch,n.cpue=res$arglist$n.catch,smooth.cpue = res$arglist$smooth.cpue,smooth.dist = res$arglist$smooth.dist,empir.dist = res$arglist$empir.dist,simple.empir = res$arglist$simple.empir,beta = res$arglist$beta,D2alpha = res$arglist$D2alpha,BTyear = NULL,summary_abc=FALSE)
    nullBTyearABC<-res.nullBTyear$ABC
    if(!is.null(catchdividedegit))nullBTyearABC<-nullBTyearABC/(10^catchdividedegit)
    if(!res$arglist$timelag0){
      data_catch <- tibble(year=c((last.year-res$arglist$n.catch+1):last.year,last.year+2,last.year+2),catch=c(rep(mean.catch,res$arglist$n.catch),ABC,nullBTyearABC),ABC,nullBTyearABC,                              type=c(rep(str_c(res$arglist$n.catch-catch.abc.na,"年平均漁獲量"),n.catch),"ABC","入力データ最終年算出ABC"))
    }
    else{
      if(latest_Catch_na){
        data_catch <- tibble(year=c((last.year-res$arglist$n.catch+1):(last.year-1),last.year+1,last.year+1),catch=c(rep(mean.catch,res$arglist$n.catch-1),ABC,nullBTyearABC),                              type=c(rep(str_c(res$arglist$n.catch-catch.abc.na,"年平均漁獲量"),n.catch-1),"ABC","入力データ最終年算出ABC"))
      }else{
        data_catch <- tibble(year=c((last.year-res$arglist$n.catch+1):last.year,last.year+1,last.year+1),catch=c(rep(mean.catch,res$arglist$n.catch),ABC,nullBTyearABC),                              type=c(rep(str_c(res$arglist$n.catch-catch.abc.na,"年平均漁獲量"),n.catch),"ABC","入力データ最終年算出ABC"))
      }
    }
  }
  if(catch.abc.na!=0) {
    data_catch2 <- data_catch
    data_catch2$catch[which(is.na(ccdata$catch[(length(ccdata$catch)-n.catch+1):length(ccdata$catch)]))] <-NA
  }

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
  #label.y.position<-c(0.5,0.4,0.8)
  #label.y.nudge<-c(0.1,-0.1,0.1)

  linetype.set <- c("dashed","dotdash","solid")

  legend.labels2 <-c(str_c(res$arglist$n.catch-catch.abc.na,"年平均漁獲量"),"ABC")
  legend.labels2bt <-c(str_c(res$arglist$n.catch-catch.abc.na,"年平均漁獲量"),"ABC","最終年データ利用のABC")
  legend.labels2.1 <-c(str_c(res$arglist$n.catch-catch.abc.na,"年平均漁獲量"),"算定漁獲量")
  legend.labels2bt.1 <-c(str_c(res$arglist$n.catch-catch.abc.na,"年平均漁獲量"),"算定漁獲量","最終年データ利用の算定")
  if(!res$arglist$timelag0){
    legend.labels2.2 <-c(str_c(res$arglist$n.catch-catch.abc.na,"年平均漁獲量"),paste(max(years)+2,"年",gsub("年","",year.axis.label),"の予測値",sep=""))
    legend.labels2bt.2 <-c(str_c(res$arglist$n.catch-catch.abc.na,"年平均漁獲量"),paste(max(years)+2,"年",gsub("年","",year.axis.label),"の予測値",sep=""),"最終年データ利用の予測値")
  } else{
    legend.labels2.2 <-c(str_c(res$arglist$n.catch-catch.abc.na,"年平均漁獲量"),paste(max(years)+1,"年",gsub("年","",year.axis.label),"の予測値",sep=""))
    legend.labels2bt.2 <-c(str_c(res$arglist$n.catch-catch.abc.na,"年平均漁獲量"),paste(max(years)+1,"年",gsub("年","",year.axis.label),"の予測値",sep=""),"最終年データ利用の予測値")
  }
  col.BRP.hcr <- col.BRP
  data_BRP_hcr <- tibble(BRP=names(res$BRP),value_obs=res$Obs_BRP, value_ratio=res$BRP)

  # PB=0の時の禁漁水準削除設定 ----
  data_BRP2 <- data_BRP
  if(res$BRP[3] == 0) {
    if(proposal==TRUE){
      legend.labels <- c("目標管理基準値（目標水準）案","限界管理基準値（限界水準）案")
    }else{
      legend.labels <- c("目標管理基準値（目標水準）","限界管理基準値（限界水準）")
    }
    linetype.set <- c("dashed","dotdash")
    if(abc4==TRUE){
      col.BRP <- c("blue","red")
    }else{
      col.BRP <- c("#00533E","#edb918")
    }
    data_BRP <- tibble(BRP=names(res$BRP[-3]),value_obs=res$Obs_BRP[-3],value_ratio=res$BRP[-3])
  }else{
    if(abc4==TRUE){
      col.BRP <- c("blue","red","orange")
    }else{
      col.BRP <- c("#00533E","#edb918","#C73C2E")
    }

    #label.y.position<-c(0.5,0.4)
    #label.y.nudge<-c(0.1,-0.1)
  }

  # ABC決定可能/不可能設定 ----
  if(detABC==1){
    g.catch.title <- "漁獲量のトレンドと算定漁獲量"
    g.catch.abcpoint <- "算定漁獲量"
    if(BThcr==T) legend.labels2 <- legend.labels2bt.1
    else legend.labels2 <- legend.labels2.1
    #legend.labels2 <- legend.labels2.1
  }else if(detABC==2){
    g.catch.title <- "漁獲量の推移と予測値"
    g.catch.abcpoint <- "予測値"
    if(BThcr==T) legend.labels2 <- legend.labels2bt.2
    else legend.labels2 <- legend.labels2.2
    #legend.labels2 <- legend.labels2.2
  }else{
    g.catch.title <- "漁獲量のトレンドとABC"
    g.catch.abcpoint <- "ABC"
    if(BThcr==T) legend.labels2 <- legend.labels2bt
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

  if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){ ## plot 設定 for mac----
    g.cpue <- ccdata %>% ggplot() +
      geom_polygon(data=tibble(x=c(minyears,max(years)+4,max(years)+4,minyears), y=c(data_BRP2$value_obs[1],data_BRP2$value_obs[1],max(ccdata$cpue,na.rm=T)*1.05,max(ccdata$cpue,na.rm=T)*1.05)), aes(x=x,y=y), fill=colfill[1]) +
      geom_polygon(data=tibble(x=c(minyears,max(years)+4,max(years)+4,minyears), y=c(data_BRP2$value_obs[2],data_BRP2$value_obs[2],data_BRP2$value_obs[1],data_BRP2$value_obs[1])), aes(x=x,y=y), fill=colfill[2]) +
      geom_polygon(data=tibble(x=c(minyears,max(years)+4,max(years)+4,minyears), y=c(data_BRP2$value_obs[3],data_BRP2$value_obs[3],data_BRP2$value_obs[2],data_BRP2$value_obs[2])), aes(x=x,y=y), fill=colfill[3]) +
      geom_polygon(data=tibble(x=c(minyears,max(years)+4,max(years)+4,minyears), y=c(0,0,data_BRP2$value_obs[3],data_BRP2$value_obs[3])), aes(x=x,y=y), fill=colfill[4]) +
      geom_hline(yintercept=res$Obs_percent_even,color="gray",linetype=2)+
      geom_text(data=data_percent_even,aes(x=x,y=y*1.05,label=label))+
      geom_text(aes(x=max(years)+3,y=min(data_percent_even$y[data_percent_even$y >= 0])*0.75,family=font_MAC,label="(資源水準)"),size=4)

    if(RP==TRUE){
      g.cpue <- g.cpue +
        geom_hline(data=data_BRP, mapping=aes(yintercept=value_obs, color=rev(col.BRP), linetype=rev(linetype.set)), size = 0.9*1.5)+
        #ggrepel::geom_label_repel(data=data_BRP, mapping=aes(x=min(years)+0.5, y=value_obs+0.5, label=legend.labels,family = font_MAC), box.padding=0.5, nudge_x=1)+
        scale_linetype_manual(name="", values=rev(c(linetype.set)), labels=rev(c(legend.labels))) +
        scale_color_manual(name="",values=rev(c(col.BRP)),labels=rev(c(legend.labels)))
    }else{
      g.cpue <- g.cpue +
        geom_point(mapping=aes(x=rev(year)[1], y=rev(ccdata$cpue)[1], color="red"),size=4, show.legend =TRUE)+
        scale_color_manual(name="",values="red",labels="直近年の資源量指標値")
    }
    g.cpue <- g.cpue +
      geom_path(data=ccdata, aes(x=year,y=cpue),size=1)+geom_point(aes(x=year,y=cpue),size=3) +
      theme_bw()+ylab(paste("資源量指標値",cpueunit))+xlab(year.axis.label)+
      ylim(0,max(ccdata$cpue,na.rm=T)*1.05)+theme_custom()

    g.cpue <- g.cpue +
      ggtitle("")+
      theme(legend.position="top",legend.justification = c(1,0), legend.spacing=unit(0.25,'lines'), legend.key.width = unit(2.0, 'lines')) +
      theme(text = element_text(family = font_MAC))
    if(leftalign==TRUE){
      g.cpue <- g.cpue + xlim(minyears,max(ccdata[!is.na(ccdata$cpue),]$year)+4)
    }
  }else{
    g.cpue <- ccdata %>% ggplot() +

      geom_polygon(data=tibble(x=c(minyears,max(years)+4,max(years)+4,minyears), y=c(data_BRP2$value_obs[1],data_BRP2$value_obs[1],max(ccdata$cpue,na.rm=T)*1.05,max(ccdata$cpue,na.rm=T)*1.05)), aes(x=x,y=y), fill=colfill[1]) +
      geom_polygon(data=tibble(x=c(minyears,max(years)+4,max(years)+4,minyears), y=c(data_BRP2$value_obs[2],data_BRP2$value_obs[2],data_BRP2$value_obs[1],data_BRP2$value_obs[1])), aes(x=x,y=y), fill=colfill[2]) +
      geom_polygon(data=tibble(x=c(minyears,max(years)+4,max(years)+4,minyears), y=c(data_BRP2$value_obs[3],data_BRP2$value_obs[3],data_BRP2$value_obs[2],data_BRP2$value_obs[2])), aes(x=x,y=y), fill=colfill[3]) +
      geom_polygon(data=tibble(x=c(minyears,max(years)+4,max(years)+4,minyears), y=c(0,0,data_BRP2$value_obs[3],data_BRP2$value_obs[3])), aes(x=x,y=y), fill=colfill[4]) +
      geom_hline(yintercept=res$Obs_percent_even,color="gray",linetype=2)+
      geom_text(data=data_percent_even, aes(x=x,y=y*1.05,label=label))+
      geom_text(aes(x=max(years)+3,y=min(data_percent_even$y[data_percent_even$y >= 0])*0.75,label="(資源水準)"),size=4)

    if(RP==TRUE){
      g.cpue <- g.cpue +
        geom_hline(data=data_BRP, mapping=aes(yintercept=value_obs, color=rev(col.BRP), linetype=rev(linetype.set)), size = 0.9*1.5)+
        #ggrepel::geom_label_repel(data=data_BRP, mapping=aes(x=min(years)+0.5, y=value_obs+0.5, label=legend.labels), box.padding=0.5, nudge_x=1)+
        scale_linetype_manual(name="", values=rev(c(linetype.set)), labels=rev(c(legend.labels))) +
        scale_color_manual(name="",values=rev(c(col.BRP)),labels=rev(c(legend.labels)))
    }else{
      g.cpue <- g.cpue +
        geom_point(mapping=aes(x=rev(year)[1], y=rev(ccdata$cpue)[1], color="red"),size=3, show.legend =TRUE)+
        scale_color_manual(name="",values="red",labels="直近年の資源量指標値")
    }
    g.cpue <- g.cpue +
      geom_path(data=ccdata, aes(x=year,y=cpue),size=1)+geom_point(aes(x=year,y=cpue),size=4) +
      theme_bw()+ylab(paste("資源量指標値",cpueunit))+xlab(year.axis.label)+
      ylim(0,max(ccdata$cpue,na.rm=T)*1.05)+theme_custom()
    g.cpue <- g.cpue +
      ggtitle("") + theme(legend.position="top", legend.spacing=unit(0.25,'lines'), legend.key.width = unit(5, 'lines'),legend.justification=c(1,0))
    if(leftalign==TRUE){
      g.cpue <- g.cpue + xlim(minyears, max(ccdata[!is.na(ccdata$cpue),]$year)+4)
    }

    }
  g.cpue <- g.cpue %>% apply_minor_ticks_type2()

  if(isTRUE(abc4)){
    hanrei_label <- rev(c(paste(min(ccdata[!is.na(ccdata$cpue),]$year),"～",max(ccdata[!is.na(ccdata$cpue),]$year),"年", gsub("年","",year.axis.label), "の平均水準",sep=""),"過去最低値"))  ##OS200702
    g.cpue4 <- ccdata %>% ggplot() +
      geom_hline(yintercept=res$Obs_percent_even,color="gray",linetype=2)+
      geom_text(data=data_percent_even,aes(x=x+x_right_space,y=y*1.05,label=label))+
      geom_text(aes(x=max(years)+x_right_space-5,y=min(data_percent_even$y)*0.75,label="(指標値の水準)"),size=4)+
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
  

    if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){ ## plot 設定 for mac----
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
    g.cpue4 <- g.cpue4 %>% apply_minor_ticks_type2()
  }

  #漁獲管理規則案 HCR ----
  ifelse(is.null(BTyear),ccdata.plot<-ccdata,
         ccdata.plot<-ccdata_fixedBT)
  #プロットの順番；枠、資源水準vsアルファ、管理水準縦線、ラベル、軸ラベル、abc計算に用いる現状ポイント
  g.hcr <- ggplot(data=data.frame(X=c(0,100)), aes(x=X)) #プロット枠
  if(!empir.dist){
    g.hcr <- g.hcr +
      #stat_function(fun=type2_func_wrapper,
      #        args=list(BT=BT,PL=0,PB=PB,tune.par=tune.par,beta=beta,AAV=res$AAV,type="%"),
      #              color="gray")+
      stat_function(fun=type2_func_wrapper,
                    args=list(BT=BT,PL=PL,PB=PB,tune.par=tune.par,beta=beta,AAV=res$AAV,type="%"),
                    color="black",size=1)
    ## BTyear!=NULLの場合、最新年のデータまで入れたデータでの結果を同時出力
    if(BThcr) g.hcr <- g.hcr +
        stat_function(fun=type2_func_wrapper,
                      args=list(BT=BT,PL=PL,PB=PB,tune.par=tune.par,beta=beta,AAV=res.nullBTyear$AAV,type="%"), color="grey",size=0.5,linetype="dashed")

  }else{ # 経験分布利用 empir.dist=T ----
    g.hcr <- g.hcr +
      stat_function(fun=type2_func_empir_wrapper,
                    args=list(BT=BT,PL=PL,PB=PB,tune.par=tune.par,beta=beta,AAV=res$AAV,cpue=ccdata.plot$cpue,simple=simple.empir,type="%"),
                    color="black",size=1)
    ## BTyear!=NULLの場合、最新年のデータまで入れたデータでの結果を同時出力
    if(BThcr) g.hcr <- g.hcr +
        stat_function(fun=type2_func_empir_wrapper,
                      args=list(BT=BT,PL=PL,PB=PB,tune.par=tune.par,beta=beta,AAV=res.nullBTyear$AAV,cpue=ccdata$cpue,simple=simple.empir,type="%"), color="grey",size=1,linetype="dashed")
  }

  g.hcr <- g.hcr +
    geom_vline(data=data_BRP,mapping=aes(xintercept=value_ratio*100,color=BRP), size = 0.9*1.5, linetype = linetype.set)

  if(hcrhscale=="sparse") hlinebreaks <- c(0,0.5,1.0)
  if(hcrhscale=="middle") hlinebreaks <- c(0,0.25,0.5,0.75,1.0)
  if(hcrhscale=="dense") hlinebreaks <- c(0,0.2,0.4,0.6,0.8,1.0)

  if(!plotexactframe) g.hcr <- g.hcr + scale_y_continuous(breaks = hlinebreaks)
  else g.hcr <- g.hcr + scale_x_continuous(expand = c(0,0),limits = c(0,100)) + scale_y_continuous(expand = c(0,0),breaks = hlinebreaks)

  if(hcrhline=="none") hcrAuxiliaryhline <- c()
  if(hcrhline=="one") hcrAuxiliaryhline <- c(1.0)
  if(hcrhline=="sparse") hcrAuxiliaryhline <- c(0,0.5,1.0)
  if(hcrhline=="middle") hcrAuxiliaryhline <- c(0,0.25,0.5,0.75,1.0)
  if(hcrhline=="dense") hcrAuxiliaryhline <- c(0,0.2,0.4,0.6,0.8,1.0)
  if(hcrhline=="hscale") hcrAuxiliaryhline <- hlinebreaks

  g.hcr <- g.hcr +
    geom_hline(yintercept=hcrAuxiliaryhline,color="gray",linetype=2)

  if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){ ## 図中ラベルと軸ラベルの設定 mac ----
    if(res$BRP[3]==0) #禁漁水準=0の時
      g.hcr <- g.hcr +
        ggrepel::geom_label_repel(data=data_BRP,
                                  mapping=aes(x=value_ratio*100, y=c(0.55,0.35), label=legend.labels,family = font_MAC),
                                  box.padding=0.5,nudge_y=c(0.1,-0.1) )
    else  g.hcr <- g.hcr +
        ggrepel::geom_label_repel(data=data_BRP,
                                  mapping=aes(x=value_ratio*100, y=c(0.55,0.35,0.7), label=legend.labels,family = font_MAC),
                                  box.padding=0.5,nudge_y =c(0.1,-0.1,0.1) )

    g.hcr <- g.hcr+
      scale_color_manual(name="",values=rev(c(col.BRP)), guide="none")+ #,labels=rev(c(legend.labels)))+
      theme_bw()+theme_custom()+
      ggtitle("")+
      xlab("資源水準(%)")+ylab(str_c("漁獲量を増減させる係数"))+
      theme(legend.position="top",legend.justification = c(1,0)) +
      theme(text = element_text(family = font_MAC))
  }else{ ## 図中ラベルと軸ラベルの設定 mac以外 ----
    if(res$BRP[3]==0) #禁漁水準=0の時
      g.hcr <- g.hcr +
        ggrepel::geom_label_repel(data=data_BRP,
                                  mapping=aes(x=value_ratio*100, y=c(0.55,0.35), label=legend.labels),
                                  box.padding=0.5, nudge_y=c(0.1,-0.1))
    else g.hcr <- g.hcr +
        ggrepel::geom_label_repel(data=data_BRP,
                                  mapping=aes(x=value_ratio*100, y=c(0.55,0.35,0.7), label=legend.labels),
                                  box.padding=0.5, nudge_y=c(0.1,-0.1,0.1))
    g.hcr <- g.hcr +
      scale_color_manual(name="",values=rev(c(col.BRP)), guide="none")+#,labels=rev(c(legend.labels)))+
      theme_bw()+theme_custom()+
      ggtitle("")+
      xlab("資源水準(%)")+ylab(str_c("漁獲量を増減させる係数"))+
      theme(legend.position="top",legend.justification = c(1,0))
  }

  g.hcr <- g.hcr +
    geom_point(aes(x=res$Current_Status[1]*100,y=res$alpha),color=2,size=4)

  ## BTyear!=NULLの場合、BTyear==NULLの現状ポイントを出力 ----
  if(BThcr){
    if(!empir.dist) g.hcr <- g.hcr +
        geom_point(aes(x=res.nullBTyear$Current_Status[1]*100,y=res.nullBTyear$alpha),color=3,size=4)
    else g.hcr <- g.hcr +
        geom_point(aes(x=res.nullBTyear$Current_Status[1]*100,y=res.nullBTyear$alpha),color=3,size=4)
  }

  if(plotexactframe) g.hcr <- g.hcr + theme(plot.margin = margin(0,15,0,10))


  #漁獲管理規則案 HCR.Dist ----
  current_index_col <- "#1A4472"
  # cpueの桁数に応じてHCR.Distのcpue刻み幅をかえる
  model_dist<-c()
  ifelse( floor(log10(max(ccdata.plot$cpue,na.rm = T))+1) < 4, model_dist <- data.frame(cpue=seq(0, max(ccdata.plot$cpue,na.rm=T), by=0.1),  dens=NA), model_dist <- data.frame(cpue=seq(0, max(ccdata.plot$cpue,na.rm=T), by=10^floor(log10(max(ccdata.plot$cpue,na.rm = T))+1)/1000),  dens=NA) )

  if(!empir.dist) model_dist$dens <- dnorm(model_dist$cpue,mean = mean(ccdata.plot$cpue,na.rm=T),sd=sd(ccdata.plot$cpue,na.rm = T))
  else{ # empir.dist = T で累積確率から個々の確率を求めて総和(1)で割って密度にする
    if(!simple.empir){
      cum.cpue4 <- ecdf(na.omit(ccdata.plot$cpue))
      cum.probs <- cum.cpue4(model_dist$cpue)
    }
    else{
      cum.probs<-c()
      for(i in 1:length(model_dist$cpue)){
        cum.probs <- c(cum.probs,simple_ecdf(ccdata.plot$cpue,model_dist$cpue[i]))
      }
    }
    probs<-c(cum.probs[1])
    for(j in 2:length(cum.probs)){
      if(cum.probs[j-1]!=cum.probs[j]) tmp <- cum.probs[j]-cum.probs[j-1]
      else tmp <- probs[j-1]
      probs<-c(probs,tmp)
    }
    #plot(model_dist$cpue,probs)
    model_dist$dens<-probs/sum(probs)
  }

  g.hcr.dist <- ggplot(data=model_dist)+
    #stat_function(fun=dnorm,args=list(mean=mean(ccdata.plot$cpue),sd=sd(ccdata.plot$cpue)),color="black",size=1)
    geom_line(aes(x=cpue,y=dens))+
    geom_area(data=filter(model_dist, cpue < res$Current_Status[2]), aes(x=cpue, y=dens), fill="grey")

  g.hcr.dist <-  g.hcr.dist +
    geom_vline(data=data_BRP,mapping=aes(xintercept=value_obs,color=BRP), size = 0.9*1.5, linetype = linetype.set) +
    scale_linetype_manual(name="", values=rev(c(linetype.set)), labels=rev(c(legend.labels))) +
    scale_color_manual(name="",values=rev(c(col.BRP)),labels=rev(c(legend.labels)))+
    guides(colour="none")+
    coord_flip()

  if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){ ## plot 設定 for mac----
    g.hcr.dist <- g.hcr.dist +
      geom_vline(data=data_BRP,mapping=aes(xintercept=res$Current_Status[2]),color=current_index_col,size=1,linetype="dashed")+
      geom_text(aes(x=ifelse(res$Current_Status[2]<mean(ccdata.plot$cpue)/3,mean(ccdata.plot$cpue)/2,mean(ccdata.plot$cpue)/4),y=max(dens)*0.85,family=font_MAC,label="(現在の資源水準)"),color=current_index_col,size=4)
  }else{
    g.hcr.dist <- g.hcr.dist +
      geom_vline(data=data_BRP,mapping=aes(xintercept=res$Current_Status[2]),color=current_index_col,size=1,linetype="dashed")+
      geom_text(aes(x=ifelse(res$Current_Status[2]<mean(ccdata.plot$cpue)/3,mean(ccdata.plot$cpue)/2,mean(ccdata.plot$cpue)/4),y=max(dens)*0.85,label="(現在の資源水準)"),color=current_index_col,size=4)
  }

  if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){ ## plot 設定 for mac----
    g.hcr.dist <- g.hcr.dist +  ggtitle("")+
      scale_x_continuous(limits=c(0,max(ccdata$cpue,na.rm=T)*1.05)) +
      #scale_y_continuous(limits=c(0,max(ccdata$cpue,na.rm=T)*1.05)) +
      xlab("資源量指標値")+ylab("")+
      theme_bw()+theme_custom()+
      theme(legend.position="top",legend.justification = c(1,0), legend.spacing=unit(0.25,'lines'), legend.key.width = unit(2.0, 'lines'),axis.text.x = element_blank()) +
      theme(text = element_text(family = font_MAC))
  }else{
    g.hcr.dist <- g.hcr.dist +  ggtitle("")+
      scale_x_continuous(limits=c(0,max(ccdata$cpue,na.rm=T)*1.05)) +
      #scale_y_continuous(limits=c(0,max(ccdata$cpue,na.rm=T)*1.05)) +
      xlab("資源量指標値")+ylab("")+
      theme_bw()+theme_custom()+
      theme(legend.position="top",legend.justification = c(1,0), legend.spacing=unit(0.25,'lines'), legend.key.width = unit(2.0, 'lines'),axis.text.x = element_blank())
  }

  #漁獲量のトレンドとABC/算定漁獲量 ----

  if(BThcr==T) CatchABC<-c(1,2,3)
  else CatchABC<-c(1,2)

  if(!ignore_naCatch_point) { #
    g.catch <- ccdata %>% ggplot() +
      geom_path(data=data_catch,mapping=aes(x=year,y=catch,color=type),size=2)+
      geom_point(data=data_catch,mapping=aes(x=year,y=catch,color=type),size=3)
  }else{
    g.catch <- ccdata %>% ggplot() +
      geom_path(data=data_catch,mapping=aes(x=year,y=catch,color=type),size=2)+
      geom_point(data=data_catch2,mapping=aes(x=year,y=catch,color=type),size=3)
  }

  #ggrepel::geom_label_repel(data=data_catch,
  #                          mapping=aes(x=max(year)-5, y=catch, label=legend.labels2),
  #                          box.padding=0.5, nudge_y=1)+
  #scale_color_manual(name="",values=c(1,2),labels=legend.labels2)+
  g.catch <- g.catch +scale_color_manual(name="",values=CatchABC,labels=legend.labels2)+
    # geom_point(data=dplyr::filter(data_catch,type=="ABC"),
    #                    mapping=aes(x=year,y=catch),lwd=2,color=1)+
    #         geom_line(data=dplyr::filter(data_catch,type!="ABC"),
    #                    mapping=aes(x=year,y=catch),lwd=2,color="gray")+
    geom_path(aes(x=year,y=catch),size=1)+
    ylab(paste("漁獲量",catchunit))+xlab(year.axis.label)+
    ggtitle("")+
    ylim(0,NA)+ theme_custom()+
    theme(legend.position="top",legend.justification = c(1,0))

  if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){# plot 設定 for mac
    if(!ignore_naCatch_point) { #
      g.catch <- ccdata %>% ggplot() +
        geom_path(data=data_catch,mapping=aes(x=year,y=catch,color=type),size=2)+
        geom_point(data=data_catch,mapping=aes(x=year,y=catch,color=type),size=3)
    }else{
      g.catch <- ccdata %>% ggplot() +
        geom_path(data=data_catch,mapping=aes(x=year,y=catch,color=type),size=2)+
        geom_point(data=data_catch2,mapping=aes(x=year,y=catch,color=type),size=3)
    }
    #ggrepel::geom_label_repel(data=data_catch,
    #                          mapping=aes(x=max(year)-5,y=catch[1],label=legend.labels2,family=font_MAC),
    #                          box.padding=0.5, nudge_y=1)+
    #scale_color_manual(name="",values=c("black","red"),labels=legend.labels2)+
    g.catch <- g.catch + scale_color_manual(name="",values=CatchABC,labels=legend.labels2)
    # geom_point(data=dplyr::filter(data_catch,type=="ABC"),
    #                    mapping=aes(x=year,y=catch),lwd=2,color=1)+
    #         geom_line(data=dplyr::filter(data_catch,type!="ABC"),
    #                    mapping=aes(x=year,y=catch),lwd=2,color="gray")+

    g.catch <- g.catch+
      geom_path(aes(x=year,y=catch),size=1)+
      ylab(paste("漁獲量",catchunit))+xlab(year.axis.label)+
      ggtitle("")+
      ylim(0,NA)+ theme_custom()+
      theme(legend.position="top",legend.justification = c(1,0)) +
      theme(text = element_text(family = font_MAC))
  }

  g.catch <- g.catch + geom_point(aes(x = year, y = catch), shape=1, size = 2)
  g.catch <- g.catch %>% apply_minor_ticks_type2()

  # 出力設定 ----
  if(isTRUE(hcrdist)){
    if(isTRUE(abc4)){
      graph.component <- list(g.cpue4,g.cpue,g.hcr.dist,g.hcr,g.catch)
      graph.combined <- gridExtra::grid.arrange(g.cpue4,g.cpue,g.hcr.dist,g.hcr,g.catch,ncol=3,top=stock.name)
      return(list(graph.component=graph.component,graph.combined=graph.combined))
    }else{
      graph.component <- list(g.cpue,g.hcr.dist,g.hcr,g.catch)
      graph.combined <- gridExtra::grid.arrange(g.cpue,g.hcr.dist,g.hcr,g.catch,ncol=2,top=stock.name)
      return(list(graph.component=graph.component,graph.combined=graph.combined))
    }
  }else{
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
}


#' 3系のABC計算結果のをプロットするための関数
#'
#' @param res calc_abc3の返り値
#'
#' @export
#'

plot_abc3 <- function(res,stock.name=NULL,fishseason=0,detABC=0,proposal=TRUE,catchunit="(トン)"){
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
    linetype.set <- c("dashed","dotdash","solid")
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
     geom_path(data=data_catch,mapping=aes(x=year,y=catch,color=type),size=2)+
     geom_point(data=data_catch,mapping=aes(x=year,y=catch,color=type),size=3)+
     #         geom_point(data=dplyr::filter(data_catch,type=="ABC"),
     #                    mapping=aes(x=year,y=catch),lwd=2,color=1)+
     #         geom_line(data=dplyr::filter(data_catch,type!="ABC"),
     #                    mapping=aes(x=year,y=catch),lwd=2,color="gray")+
     geom_path(aes(x=year,y=catch),size=1)+
     geom_text(data=data_percent,aes(x=x,y=y,label=label))+
     geom_text(aes(x=min(ccdata$year)+2,y=min(data_percent$y)*0.75,label="(漁獲量水準)"),size=4)+
     theme_bw()+ylab(paste("漁獲量",catchunit))+xlab(year.axis.label)+theme_custom()+
     geom_hline(data=data_BRP,mapping=aes(yintercept=value_obs,color=BRP), size = 0.9*2, linetype = linetype.set)+
     scale_color_manual(name="",values=c(1,2,rev(col.BRP)),labels=legend.labels2)+
     ylim(0,NA)+xlim(min(ccdata$year)-1,NA)+
     ggtitle(g.catch.title)+
     theme(legend.position="top",legend.justification = c(1,0)))

  if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){ # plot 設定 for mac
    col.set <- c("#000000","#FF0000",rev(col.BRP))
    g.catch <- ccdata %>% ggplot() +
      theme_bw(base_family = font_MAC)+
      geom_path(data=data_catch,mapping=aes(x=year,y=catch,color=type),size=2)+
      geom_point(data=data_catch,mapping=aes(x=year,y=catch,color=type),size=3)+
      #          geom_point(data=dplyr::filter(data_catch,type=="ABC"),
      #                      mapping=aes(x=year,y=catch),lwd=2,color="red")+
      #          geom_line(data=dplyr::filter(data_catch,type!="ABC"),
      #                      mapping=aes(x=year,y=catch),lwd=3,color="black")+
      geom_path(aes(x=year,y=catch),size=1)+
      ylab(paste("漁獲量",catchunit))+xlab(year.axis.label)+theme_custom()+geom_text(data=data_percent,aes(x=x,y=y,label=label),family = font_MAC)+
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
  linetype.set <- c("dashed","dotdash","solid")
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

intersection_hcrs <- function(res.list){

  if(length(res.list)>2)stop("This function calculates the instersection between 2 hcrs, so length(res.list) must be 2.\n")
  out <- NULL

  BT1 <- res.list[[1]]$arglist$BT
  BL1 <- res.list[[1]]$arglist$PL*BT1
  BB1 <- res.list[[1]]$arglist$PB*BT1
  delta1 <- res.list[[1]]$arglist$tune.par
  AAV1<-res.list[[1]]$AAV
  beta1 <- res.list[[1]]$arglist$beta

  BT2 <- res.list[[2]]$arglist$BT
  BL2 <- res.list[[2]]$arglist$PL*BT2
  BB2 <- res.list[[2]]$arglist$PB*BT2
  delta2 <- res.list[[2]]$arglist$tune.par
  AAV2<-res.list[[2]]$AAV
  beta2 <- res.list[[2]]$arglist$beta

  # max(BL1,BL2) <= x に交点を持つケース
  D.l.BL <- alpha.l.BL <- NULL
  if(beta1==beta2){
    if((BT2*delta2[1]-BT1*delta1[1])/(delta2[1]-delta1[1])>max(BL1,BL2)){
      D.l.BL <- ((BT2*delta2[1]-BT1*delta1[1])/(delta2[1]-delta1[1]))
    }
  }else if(beta1!=1 && beta2!=1){ #beta1!=beta2
    if( (BT2*delta2[1]-BT1*delta1[1]+log(beta1/beta2))/(delta2[1]-delta1[1]) > max(BL1,BL2)){
      D.l.BL <- ((BT2*delta2[1]-BT1*delta1[1]+log(beta1/beta2))/(delta2[1]-delta1[1]))
    }
  }

  # min(BL1,BL2) < x < max(BL1,BL2) で交点を持つケース
  D1.m.BL <- D2.m.BL <- alpha1.m.BL <- alpha2.m.BL <-NULL

  if(BL1>BL2){
    Delta2.1<-(delta1[2]*exp(delta1[3]*log(AAV1^2+1)))
    Delta2.2<-0
  }
  if(BL2>=BL1){
    Delta2.1<-0
    Delta2.2<-(delta2[2]*exp(delta2[3]*log(AAV2^2+1)))
  }

  A.1<-(delta1[1]-Delta2.1)
  A.2<-(delta2[1]-Delta2.2)

  B.1<-log(beta1)+(Delta2.1*(BT1+BL1)-delta1[1]*BT1)
  B.2<-log(beta2)+(Delta2.2*(BT2+BL2)-delta2[1]*BT2)

  C.1<- (-1*Delta2.1*BT1*BL1)
  C.2<- (-1*Delta2.2*BT2*BL2)

  A<-A.1-A.2
  B<-B.1-B.2
  C<-C.1-C.2

  if(A.1!=A.2){
    if( (B^2-4*A*C)>0){
      D1.m.BL<-max((-B+sqrt(B^2-4*A*C))/(2*A),((-B-sqrt(B^2-4*A*C))/(2*A)))
      D2.m.BL<-min((-B+sqrt(B^2-4*A*C))/(2*A),((-B-sqrt(B^2-4*A*C))/(2*A)))
    }
  }else if(B.1!=B.2){
    D1.m.BL<- -C/B
  }

  # max(BB1,BB2) < x < min(BL1,BL2) に交点を持つケース
  D1.s.BL <- D2.s.BL <- D3.s.BL <-NULL
  if(BB1==0 && BB2==0){
    Delta2.1<-(delta1[2]*exp(delta1[3]*log(AAV1^2+1)))
    Delta2.2<-(delta2[2]*exp(delta2[3]*log(AAV2^2+1)))

    A.1<-(delta1[1]-Delta2.1)
    A.2<-(delta2[1]-Delta2.2)

    B.1<-log(beta1)+(Delta2.1*(BT1+BL1)-delta1[1]*BT1)
    B.2<-log(beta2)+(Delta2.2*(BT2+BL2)-delta2[1]*BT2)

    C.1<- (-1*Delta2.1*BT1*BL1)
    C.2<- (-1*Delta2.2*BT2*BL2)

    A<-A.1-A.2
    B<-B.1-B.2
    C<-C.1-C.2

    if(A.1!=A.2){
      if( (B^2-4*A*C)>0){
        D1.s.BL<-max((-B+sqrt(B^2-4*A*C))/(2*A),((-B-sqrt(B^2-4*A*C))/(2*A)))
        D2.s.BL<-min((-B+sqrt(B^2-4*A*C))/(2*A),((-B-sqrt(B^2-4*A*C))/(2*A)))
      }
    }else if(B.1!=B.2){
      D1.s.BL<- -C/B
    }

  }
  else{ # Bban!=0
    Delta2.1<-(delta1[2]*exp(delta1[3]*log(AAV1^2+1)))
    Delta2.2<-(delta2[2]*exp(delta2[3]*log(AAV2^2+1)))

    A.1<-delta1[1]-Delta2.1
    A.2<-delta2[1]-Delta2.2

    B.1<-Delta2.1*(BT1+BL1+BB2) -delta1[1]*(BT1+BB1+BB2) + log(beta1)
    B.2<-Delta2.2*(BT2+BL2+BB1) -delta2[1]*(BT2+BB2+BB1) + log(beta2)

    C.1<-delta1[1]*(BT1*BB2+BT1*BB1+BB1*BB2)-Delta2.1*(BT1*BB2+BT1*BL1+BL1*BB2)
    -(BB1+BB2)*log(beta1)
    C.2<-delta2[1]*(BT2*BB1+BT2*BB2+BB2*BB1)-Delta2.2*(BT2*BB1+BT2*BL2+BL2*BB1)-(BB2+BB1)*log(beta2)

    D.1<-Delta2.1*BT1*BL1*BB2-delta1[1]*BT1*BB1*BB2+BB1*BB2*log(beta1)
    D.2<-Delta2.2*BT2*BL2*BB1-delta2[1]*BT2*BB2*BB1+BB2*BB1*log(beta1)

    a <- A.1-A.2
    b <- B.1-B.2
    c <- C.1-C.2
    d <- D.1-D.2

    if(a!=0){
      # a x^3 + b x^2 + c x + d = 0を解く
      # 解析的に解くのは難しそうなので、数値計算
      f <- function(x) a*x^3 + b*x^2 + c*x + d
      #curve(f,c(0,1));abline(h = 0)
      search_eps<-0.000001
      search_range<-c(max(BB1,BB2),min(BL1,BL2))
      DS<-c()
      for(i in 1:100){
        delta<-(max(search_range)-min(search_range))
        if(ifelse(f(min(search_range+(delta/100)*(i-1)))>0,1,-1)*ifelse(f(min(search_range+(delta/100)*i))>0,1,-1) <0 ) DS<-c(DS,uniroot(f,c(min(search_range)+(delta/100)*(i-1),min(search_range)+(delta/100)*i))$root)
      }
      D1.s.BL<-unique(DS)[1]
      if(length(unique(DS))==2) D2.s.BL<-DS[order(unique(DS),decreasing = F)[2]]
      if(length(unique(DS))==3) {
        D2.s.BL<-DS[order(unique(DS),decreasing = F)[2]]
        D3.s.BL<-DS[order(unique(DS),decreasing = F)[3]]
      }

    }else{#a=0
      AA <- B.1-B.2
      BB <- C.1-C.2
      CC <- D.1-D.2

      if(B.1!=B.2){
        if( (BB^2-4*AA*CC)>0){
          D1.s.BL<-max((-BB+sqrt(BB^2-4*AA*CC))/(2*AA),((-BB-sqrt(BB^2-4*AA*CC))/(2*AA)))
          D2.s.BL<-min((-BB+sqrt(BB^2-4*AA*CC))/(2*AA),((-BB-sqrt(BB^2-4*AA*CC))/(2*AA)))
        }
      }else if(C.1!=C.2){
        D1.s.BL<- -CC/BB
      }
    }
  }

  Dl <- Dm <- Ds <- NULL
  if(!is.null(D.l.BL) && (D.l.BL > max(BL1,BL2)) ) Dl <- D.l.BL
  if(!is.null(D1.m.BL)  && (D1.m.BL < max(BL1,BL2)) && (D1.m.BL > min(BL1,BL2)) ) Dm <- D1.m.BL
  if(!is.null(D2.m.BL)  && (D2.m.BL < max(BL1,BL2)) && (D2.m.BL > min(BL1,BL2)) ) Dm <- c(Dm,D2.m.BL)
  if(!is.null(D1.s.BL)  && (D1.s.BL > max(BB1,BB2)) && (D1.s.BL < min(BL1,BL2)) ) Ds <- D1.s.BL
  if(!is.null(D2.s.BL)  && (D2.s.BL > max(BB1,BB2)) && (D2.s.BL < min(BL1,BL2)) ) Ds <- c(Ds,D2.s.BL)
  if(!is.null(D3.s.BL)  && (D3.s.BL > max(BB1,BB2)) && (D3.s.BL < min(BL1,BL2)) ) Ds <- c(Ds,D3.s.BL)

  D.intersect<-c(Dl,Dm,Ds)
  alpha.intersect<-c()
  if(!is.null(D.intersect)){
    for(i in 1:length(D.intersect)){
      alpha.intersect<- c(alpha.intersect,calc_abc2(ccdata = res.list[[1]]$arglist$ccdata,BT = res.list[[1]]$arglist$BT, PL = res.list[[1]]$arglist$PL, PB = res.list[[1]]$arglist$PB, tune.par = res.list[[1]]$arglist$tune.par, AAV = res.list[[1]]$arglist$AAV, n.catch = res.list[[1]]$arglist$n.catch,n.cpue = res.list[[1]]$arglist$n.cpue, smooth.cpue = res.list[[1]]$arglist$smooth.cpue, smooth.dist = res.list[[1]]$arglist$smooth.dist, empir.dist = res.list[[1]]$arglist$empir.dist, simple.empir = res.list[[1]]$arglist$simple.empir, beta = res.list[[1]]$arglist$beta, D2alpha = D.intersect[i], BTyear = res.list[[1]]$arglist$BTyear, timelag0 = res.list[[1]]$arglist$timelag0,resp = res.list[[1]]$arglist$resp, summary_abc = F)$D2alpha)
    }
  }

  out<- data.frame(D=D.intersect,alpha=alpha.intersect)
  return(out)
}

#' 2系のHCRを比較するための関数
#'
#' @param res.list calc_abc2の返り値のリスト
#' @param stock.name
#'
#' @export
#'

plot_hcr2 <- function(res.list,stock.name=NULL,proposal=TRUE, hline="none", hscale="middle",plotexactframe=FALSE, vline=TRUE, vline.listnum=1,vlineBan=TRUE,vline.text=TRUE,label.list=NULL,is_point=TRUE,select_point=NULL,change_ps=NULL,one_point=FALSE,intersection=FALSE){

  font_MAC <- "HiraginoSans-W3"#"Japan1GothicBBB"#
  if(vline.listnum>length(res.list)) stop("vline.listnum must not be larger than length(res.list).\n")
  if(!is.null(label.list) & length(res.list)!=length(label.list)) stop("length(label.list) must be identical to length(res.list).\n")

  if(!is.null(select_point) && length(res.list)!=length(select_point)) stop("length(list) must be identical to length(select_point).\n")
  if(!is.null(select_point) && one_point==T) warning("one_point and select_point options do not work simulteneously. plot_hcr2 prefers one_point option. \n")
  if(!is.null(select_point) && !is.null(change_ps)) warning("select_point option can change point size by itself. plot_hcr2 enables both select_point and change ps options. \n")

  if(proposal==TRUE){
    legend.labels.hcr <-c("目標管理基準値（目標水準）案","限界管理基準値（限界水準）案","禁漁水準案")
  }else{
    legend.labels.hcr <-c("目標管理基準値（目標水準）","限界管理基準値（限界水準）","禁漁水準")
  }
  linetype.set <- c("dashed","dotdash","solid")

  if(vlineBan==FALSE) {
    if(proposal==TRUE){
      legend.labels.hcr <- c("目標管理基準値（目標水準）案","限界管理基準値（限界水準）案")
    }else{
      legend.labels.hcr <- c("目標管理基準値（目標水準）","限界管理基準値（限界水準）")
    }
    linetype.set <- c("dashed","dotdash")
    col.BRP <- c("#00533E","#edb918")
  }

  if("arglist"%in%names(res.list)) res.list <- list(res.list)

  Currentalphas<-c()
  for(i in 1:length(res.list)){
    Currentalphas<-rbind(Currentalphas,tibble(x=res.list[[i]]$Current_Status[1]*100,y=res.list[[i]]$alpha))
  }
  col.hcr.points <- seq(2,1+length(res.list))

  # 枠設定
  if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){
    g.hcr <- ggplot(data=data.frame(X=c(0,100)), aes(x=X)) +
      theme_bw(base_family = font_MAC)+theme_custom()+
      ggtitle("")+
      xlab("資源水準(%)")+ylab(str_c("漁獲量を増減させる係数"))+
      theme(legend.position="top",legend.justification = c(1,0))+
      theme(text = element_text(family = font_MAC))
  }else{
    g.hcr <- ggplot(data=data.frame(X=c(0,100)), aes(x=X)) +
      theme_bw()+theme_custom()+
      ggtitle("")+
      xlab("資源水準(%)")+ylab(str_c("漁獲量を増減させる係数"))+
      theme(legend.position="top",legend.justification = c(1,0))
  }

  # HCR曲線
  for(i in 1:length(res.list)){
    res <- res.list[[i]]
    if(vlineBan==TRUE) data_BRP <- tibble(BRP=names(res$BRP),value_obs=res$Obs_BRP,
                                          value_ratio=res$BRP)
    else data_BRP <- tibble(BRP=names(res$BRP[-3]),value_obs=res$Obs_BRP[-3],
                            value_ratio=res$BRP[-3])
    BT <- res$arglist$BT
    PL <- res$arglist$PL
    PB <- res$arglist$PB
    tune.par <- res$arglist$tune.par
    beta <- res$arglist$beta
    empir.dist<- res$arglist$empir.dist
    simple.empir<-res$arglist$simple.empir
    if(is.null(res$arglist$BTyear)) ccdata.plot<-res$arglist$ccdata
    else ccdata.plot<-res$arglist$ccdata[which(res$arglist$ccdata$year <= res$arglist$BTyear),]

    if(!empir.dist) g.hcr <- g.hcr +
      #            stat_function(fun=type2_func_wrapper,
      #                          args=list(BT=BT,PL=0,PB=PB,tune.par=tune.par,AAV=res$AAV,type="%"),
      #                       color="gray")+
      stat_function(fun=type2_func_wrapper,
                    args=list(BT=BT,PL=PL,PB=PB,tune.par=tune.par,beta=beta,AAV=res$AAV,type="%"),
                    color="black",size=1,linetype=i)
    else g.hcr <- g.hcr +
      stat_function(fun=type2_func_empir_wrapper,
                    args=list(BT=BT,PL=PL,PB=PB,tune.par=tune.par,beta=beta,AAV=res$AAV,cpue=ccdata.plot$cpue,simple=simple.empir,type="%"),
                    color="black",size=1,linetype=i)
  }

  # 目盛設定
  if(hscale=="sparse") hlinebreaks <- c(0,0.5,1.0)
  if(hscale=="middle") hlinebreaks <- c(0,0.25,0.5,0.75,1.0)
  if(hscale=="dense") hlinebreaks <- c(0,0.2,0.4,0.6,0.8,1.0)

  if(hline=="none") hcrAuxiliaryhline <- c()
  if(hline=="one") hcrAuxiliaryhline <- c(1.0)
  if(hline=="sparse") hcrAuxiliaryhline <- c(0,0.5,1.0)
  if(hline=="middle") hcrAuxiliaryhline <- c(0,0.25,0.5,0.75,1.0)
  if(hline=="dense") hcrAuxiliaryhline <- c(0,0.2,0.4,0.6,0.8,1.0)
  if(hline=="hscale") hcrAuxiliaryhline <- hlinebreaks

  if(!plotexactframe) g.hcr <- g.hcr + scale_y_continuous(breaks = hlinebreaks)
  else g.hcr <- g.hcr + scale_x_continuous(expand = c(0,0),limits = c(0,100)) + scale_y_continuous(expand = c(0,0),breaks = hlinebreaks)

  g.hcr <-  g.hcr + geom_hline(yintercept=hcrAuxiliaryhline,color="gray",linetype=2)

  # 管理水準縦線描画
  if(vline==TRUE){
    # 一つのHCRの水準線を選択
    if(vline.listnum!=0){
      vline.num<-vline.listnum
      if(!is.null(label.list)) legend.labels.hcr <-paste0(label.list[vline.num]," ",legend.labels.hcr)
      res <- res.list[[vline.num]]
      if(vlineBan==TRUE) data_BRP <- tibble(BRP=names(res$BRP),value_obs=res$Obs_BRP,
                                            value_ratio=res$BRP)
      else data_BRP <- tibble(BRP=names(res$BRP[-3]),value_obs=res$Obs_BRP[-3],
                              value_ratio=res$BRP[-3])

      if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){
        if(vlineBan){
          g.hcr <- g.hcr + geom_vline(data=data_BRP,mapping=aes(xintercept=value_ratio*100,color=BRP), size = 0.9, linetype = linetype.set)
          if(vline.text) g.hcr <- g.hcr + ggrepel::geom_label_repel(data=data_BRP,
                                                                    mapping=aes(x=value_ratio*100, y=c(0.5,1.15,0.8), label=legend.labels.hcr,family=font_MAC),
                                                                    box.padding=0.5)
        }else{
          g.hcr <- g.hcr + geom_vline(data=data_BRP,mapping=aes(xintercept=value_ratio*100,color=BRP), size = 0.9, linetype = linetype.set)
          if(vline.text) g.hcr <- g.hcr + ggrepel::geom_label_repel(data=data_BRP,
                                                                    mapping=aes(x=value_ratio*100, y=c(0.5,1.15), label=legend.labels.hcr,family=font_MAC),
                                                                    box.padding=0.5)
        }

      }else{
        if(vlineBan){
          g.hcr <- g.hcr + geom_vline(data=data_BRP,mapping=aes(xintercept=value_ratio*100,color=BRP), size = 0.9, linetype = linetype.set)
          if(vline.text) g.hcr <- g.hcr + ggrepel::geom_label_repel(data=data_BRP,
                                                                    mapping=aes(x=value_ratio*100, y=c(0.5,1.15,0.8), label=legend.labels.hcr),
                                                                    box.padding=0.5)
        }else{
          g.hcr <- g.hcr + geom_vline(data=data_BRP,mapping=aes(xintercept=value_ratio*100,color=BRP), size = 0.9, linetype = linetype.set)
          if(vline.text) g.hcr <- g.hcr + ggrepel::geom_label_repel(data=data_BRP,
                                                                    mapping=aes(x=value_ratio*100, y=c(0.5,1.15), label=legend.labels.hcr),
                                                                    box.padding=0.5)
        }
      }
    }else{ # 複数の水準線を同時描画 vline.listnum==0

      data_BRPs <-list()
      legend.labels<-list()
      linetype.sets<-list()
      line.sizes<-list()
      col.BRPs <- list()
      label.hlevels<-list()
      data_BRP<-c()
      if(vlineBan) label.hlevel<-c(0.8,0.4,1.15) #rep
      else label.hlevel<-c(0.8,0.4) #rep
      for(i in 1:length(res.list)){
        res<-res.list[[i]]
        if(vlineBan==T) data_BRPs[[i]] <- tibble(reslist=i,BRP=names(res$BRP),value_obs=res$Obs_BRP,
                                                 value_ratio=res$BRP)
        else data_BRPs[[i]] <- tibble(reslist=i,BRP=names(res$BRP[-3]),value_obs=res$Obs_BRP[-3],
                                      value_ratio=res$BRP[-3])
        if(is.null(label.list)) legend.labels[[i]] <-paste0(legend.labels.hcr,i)
        else legend.labels[[i]] <-paste0(label.list[i]," ",legend.labels.hcr)
        linetype.sets[[i]] <- rep((i+1),nrow(data_BRPs[[i]]))
        line.sizes[[i]]<- rep(0.5*length(res.list)/i,nrow(data_BRPs[[i]]))
        col.BRPs[[i]]<-col.BRP
        label.hlevels[[i]]<-label.hlevel-0.1*(i-1)
        data_BRP<-rbind(data_BRP,data_BRPs[[i]])
      }
      if(vlineBan==T) {
        flagBan<-c()
        for(i in 1:length(res.list)){
          flagBan <- c(flagBan, data_BRPs[[i]]$value_ratio[3])
        }
        if(sum(flagBan)==0){
          legend.labels[[1]][3]<-legend.labels.hcr[3]
          linetype.sets[[1]][3] <-1
          for(i in 2:length(res.list)){
            legend.labels[[i]][3]<-""
          }
        }
      }
      legend.labels.hcr<-unlist(legend.labels)
      linetype.set<-unlist(linetype.sets)
      line.size <- unlist(line.sizes)
      label.hlevel <- unlist(label.hlevels)
      boxpaddings<-0.5 #resごとに1ずつずらす
      if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){ # for mac
        if(vlineBan){
          g.hcr <- g.hcr +
            geom_vline(data=data_BRP,mapping=aes(xintercept=value_ratio*100,color=BRP), size = line.size, linetype = linetype.set)
          if(vline.text) g.hcr <- g.hcr + ggrepel::geom_label_repel(data=data_BRP,
                                                                    mapping=aes(x=value_ratio*100, y=label.hlevel,label=legend.labels.hcr,family=font_MAC),
                                                                    box.padding=boxpaddings)
        }else{
          g.hcr <- g.hcr + geom_vline(data=data_BRP,mapping=aes(xintercept=value_ratio*100,color=BRP), size = line.size, linetype = linetype.set)
          if(vline.text) g.hcr <- g.hcr + ggrepel::geom_label_repel(data=data_BRP,
                                                                    mapping=aes(x=value_ratio*100, y=label.hlevel, label=legend.labels.hcr,family=font_MAC),
                                                                    box.padding=boxpaddings)
        }
      }else{ # for !mac
        if(vlineBan){
          g.hcr <- g.hcr +
            geom_vline(data=data_BRP,mapping=aes(xintercept=value_ratio*100,color=BRP), size = line.size, linetype = linetype.set)
          if(vline.text) g.hcr <- g.hcr + ggrepel::geom_label_repel(data=data_BRP,
                                                                    mapping=aes(x=value_ratio*100, y=label.hlevel,label=legend.labels.hcr),
                                                                    box.padding=boxpaddings)
        }else{
          g.hcr <- g.hcr + geom_vline(data=data_BRP,mapping=aes(xintercept=value_ratio*100,color=BRP), size = line.size, linetype = linetype.set)
          if(vline.text) g.hcr <- g.hcr + ggrepel::geom_label_repel(data=data_BRP,
                                                                    mapping=aes(x=value_ratio*100, y=label.hlevel, label=legend.labels.hcr),
                                                                    box.padding=boxpaddings)
        }      }
    }
  }
  if(intersection){
    list.combn<-combn(x = length(res.list),m=2)
    intersect.dat<-NULL
    for(j in 1:ncol(list.combn)){
      res.list.inters<-list(res.list[[list.combn[1,j]]],res.list[[list.combn[2,j]]])

      intersect <- intersection_hcrs(res.list.inters)
      if(!is.null(intersect)){
        if(nrow(intersect)==0) cat(stringr::str_c("no intersection betweenHCR",list.combn[1,j]," and HCR",list.combn[2,j],"\n"))
        else cat(stringr::str_c("intersection between HCR",list.combn[1,j]," and HCR",list.combn[2,j],";D=",round(intersect$D,3),", alpha=",round(intersect$alpha,3),"\n"))

        if(!nrow(intersect)==0){
          if(res.list.inters[[1]]$arglist$BT==res.list.inters[[2]]$arglist$BT) intersect<-intersect[-1,]
          intersect$D <- intersect$D*100
          if(length(which(intersect$D > 100))>0) intersect <- intersect[-which(intersect$D > 100),]
          if(length(which(intersect$D < 0))>0) intersect <- intersect[-which(intersect$D < 0),]
          intersect.dat<-rbind(intersect.dat,intersect)
        }
      }
    }
    if(!is.null(intersect.dat)){
      if(!nrow(intersect.dat)==0) g.hcr <- g.hcr + geom_vline(data=intersect.dat,mapping = aes(xintercept=D),size=0.5,linetype=3)
    }
  }#intersection
  if(one_point){
    points.size.magnify <- c(1)
    for(k in 2:nrow(Currentalphas)){
      points.size.magnify <- c(points.size.magnify,NA)
    }
    g.hcr <- g.hcr +
      geom_point(data=Currentalphas,aes(x=x,y=y),color=col.hcr.points,size=4*points.size.magnify) +
      scale_color_manual(name="",values=rev(c(col.BRP)),guide="none") #label=rev(legend.labels.hcr))
  }else if(is_point){
    if(is.null(change_ps)){
      if(is.null(select_point)) g.hcr <- g.hcr +
          geom_point(data=Currentalphas,aes(x=x,y=y),color=col.hcr.points,size=4)
      else g.hcr <- g.hcr +
          geom_point(data=Currentalphas,aes(x=x,y=y),color=col.hcr.points,size=4*select_point)
      g.hcr <- g.hcr + scale_color_manual(name="",values=rev(c(col.BRP)),guide="none") #label=rev(legend.labels.hcr))
    }else{
      points.size.magnify <- c(1)
      point.size <- 1
      for(k in 1:(nrow(Currentalphas)-1)){
        points.size.magnify <- c(points.size.magnify,point.size*change_ps)
        point.size<-point.size*change_ps
      }
      if(!is.null(select_point)) points.size.magnify<-points.size.magnify*select_point
      g.hcr <- g.hcr +
        geom_point(data=Currentalphas,aes(x=x,y=y),color=col.hcr.points,size=4*points.size.magnify) +
        scale_color_manual(name="",values=rev(c(col.BRP)),guide="none") #label=rev(legend.labels.hcr))
    }
  }else{
    g.hcr <- g.hcr + scale_color_manual(name="",values=rev(c(col.BRP)),guide="none")
  }

  return(g.hcr)
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

#' 複数の2系ABC計算結果を同時プロットするための関数
#'
#' @param res.list calc_abc2の返り値。比較結果が多すぎるとみづらくなるので一度に比較するのは５つまで。
#' @param fishseason  X軸のラベルを変更（0なら年、1なら漁期年）
#' @param abc4  北海道東部の跨り資源で使用する図を描画（TRUEなら使用、デフォルトはFALSE））
#' @param cpueunit  資源量指標値の縦軸見出しに追記したい指標値の単位（例えば"（トン/網）"のように指定する）
#' @param catchunit 漁獲量トレンドの縦軸見出しに追記したい単位（例えば"（万トン）"のように指定する、デフォルトは"（トン）"）
#' @param leftalign  資源量指標値の時系列の長さが漁獲量に比べて短い時、データが無い範囲の空間を削除する（TRUEなら使用、デフォルトはFALSE）
#' @param RP  資源量指標値/年のプロットでReference Point（目標・限界管理基準線）を載せる・載せない（デフォルトはTRUE、FALSEでは直近年の資源量指標値をポイントでハイライトする）
#' @param ignore_naCatch_point ABC算出に使う最近年の漁獲量にNAが入っている場合、表示上NAとなる年のポイントと年数を引く
#' @export
#'
plot_abc2_multires <- function(res.list, stock.name=NULL, fishseason=0, detABC=0, abc4=FALSE, cpueunit="", catchunit="（トン）",catchdividedegit =NULL,fillarea=FALSE, RP=TRUE, leftalign=FALSE, proposal=TRUE, hcrdist=FALSE,BThcr=FALSE,hcrhline="none",hcrhscale="middle",hcrvlineBan=FALSE,plotexactframe=FALSE,ignore_naCatch_point=FALSE,abclegend=NULL){
  font_MAC <- "HiraginoSans-W3"#"Japan1GothicBBB"#

  #結果比較の限界は５個まで
  if(length(res.list)>5) stop("The max number in res.list is 5.\n")

  # 漁期年/年設定 ----
  ifelse(fishseason==1, year.axis.label <- "漁期年", year.axis.label <- "年")

  catch.abc.na<-0
  if(ignore_naCatch_point){
    mean.catch.abc <- res.list[[1]]$arglist$ccdata$catch[(length(ccdata$catch)-n.catch+1):length(ccdata$catch)]
    catch.abc.na <- sum(as.numeric(is.na(mean.catch.abc)))
    if(prod(!is.na(mean.catch.abc))) stop("ignore_naCatch_point option works if catch[lastyear-n.catch+1:lastyear] contains NA.")
  }

  # 漁獲量とABC出力設定 ----
  years <- res.list[[1]]$arglist$ccdata$year
  last.year <- rev(years)[1]
  labels2<-labels2.1<-labels2.2<-c()

  if(!is.null(abclegend) & (length(res.list)!=length(abclegend))) stop("the lengh of abclegend must be identical to the length of res.list!")

  for(i in 1:length(res.list)){
    if(!res.list[[1]]$arglist$timelag0){
      if(i==1) {
        if(is.null(abclegend)){
          if(is.null(catchdividedegit)) data_catch<- tibble(year=c((last.year-res.list[[i]]$arglist$n.catch+1):last.year,last.year+2),
                                                            catch=c(rep(res.list[[i]]$mean.catch,res.list[[i]]$arglist$n.catch),res.list[[i]]$ABC),
                                                            type=c(rep(str_c("平均漁獲量(",res.list[[1]]$arglist$n.catch-catch.abc.na,"年平均)"),res.list[[i]]$arglist$n.catch),paste0(i,"番目ABC")))
          else data_catch<- tibble(year=c((last.year-res.list[[i]]$arglist$n.catch+1):last.year,last.year+2),
                                   catch=c(rep(res.list[[i]]$mean.catch/(10^catchdividedegit),res.list[[i]]$arglist$n.catch),res.list[[i]]$ABC/(10^catchdividedegit)),
                                   type=c(rep(str_c("平均漁獲量(",res.list[[1]]$arglist$n.catch-catch.abc.na,"年平均)"),res.list[[i]]$arglist$n.catch),paste0(i,"番目ABC")))
        }
        else{
          if(is.null(catchdividedegit)) data_catch<- tibble(year=c((last.year-res.list[[i]]$arglist$n.catch+1):last.year,last.year+2),
                                                            catch=c(rep(res.list[[i]]$mean.catch,res.list[[i]]$arglist$n.catch),res.list[[i]]$ABC),
                                                            type=c(rep(str_c("平均漁獲量(",res.list[[1]]$arglist$n.catch-catch.abc.na,"年平均)"),res.list[[i]]$arglist$n.catch),abclegend[i]))
          else data_catch<- tibble(year=c((last.year-res.list[[i]]$arglist$n.catch+1):last.year,last.year+2),
                                   catch=c(rep(res.list[[i]]$mean.catch/(10^catchdividedegit),res.list[[i]]$arglist$n.catch),res.list[[i]]$ABC/(10^catchdividedegit)),
                                   type=c(rep(str_c("平均漁獲量(",res.list[[1]]$arglist$n.catch-catch.abc.na,"年平均)"),res.list[[i]]$arglist$n.catch),abclegend[i]))
        }
      }
      else {
        if(is.null(abclegend)){
          if(is.null(catchdividedegit)) data_catch <- rbind(data_catch,tibble(year=last.year+2,
                                                                              catch=c(res.list[[i]]$ABC),
                                                                              type=c(paste0(i,"番目ABC"))))
          else data_catch <- rbind(data_catch,tibble(year=last.year+2,
                                                     catch=c(res.list[[i]]$ABC/(10^catchdividedegit)),
                                                     type=c(paste0(i,"番目ABC"))))
        }
        else{
          if(is.null(catchdividedegit)) data_catch <- rbind(data_catch,tibble(year=last.year+2,
                                                                              catch=c(res.list[[i]]$ABC),
                                                                              type=c(abclegend[i])))
          else data_catch <- rbind(data_catch,tibble(year=last.year+2,
                                                     catch=c(res.list[[i]]$ABC/(10^catchdividedegit)),
                                                     type=c(abclegend[i])))
        }
      }
    }else{ #timelag0=T
      if(i==1) {
        if(is.null(abclegend)){
          if(is.null(catchdividedegit)) data_catch<- tibble(year=c((last.year-res.list[[i]]$arglist$n.catch+1):last.year,last.year+1),
                                                            catch=c(rep(res.list[[i]]$mean.catch,res.list[[i]]$arglist$n.catch),res.list[[i]]$ABC),
                                                            type=c(rep(str_c("平均漁獲量(",res.list[[1]]$arglist$n.catch-catch.abc.na,"年平均)"),res.list[[i]]$arglist$n.catch),paste0(i,"番目ABC")))
          else data_catch<- tibble(year=c((last.year-res.list[[i]]$arglist$n.catch+1):last.year,last.year+1),
                                   catch=c(rep(res.list[[i]]$mean.catch/(10^catchdividedegit),res.list[[i]]$arglist$n.catch),res.list[[i]]$ABC/(10^catchdividedegit)),
                                   type=c(rep(str_c("平均漁獲量(",res.list[[1]]$arglist$n.catch-catch.abc.na,"年平均)"),res.list[[i]]$arglist$n.catch),paste0(i,"番目ABC")))
        }
        else{
          if(is.null(catchdividedegit)) data_catch<- tibble(year=c((last.year-res.list[[i]]$arglist$n.catch+1):last.year,last.year+1),
                                                            catch=c(rep(res.list[[i]]$mean.catch,res.list[[i]]$arglist$n.catch),res.list[[i]]$ABC),
                                                            type=c(rep(str_c("平均漁獲量(",res.list[[1]]$arglist$n.catch-catch.abc.na,"年平均)"),res.list[[i]]$arglist$n.catch),abclegend[i]))
          else data_catch<- tibble(year=c((last.year-res.list[[i]]$arglist$n.catch+1):last.year,last.year+1),
                                   catch=c(rep(res.list[[i]]$mean.catch/(10^catchdividedegit),res.list[[i]]$arglist$n.catch),res.list[[i]]$ABC/(10^catchdividedegit)),
                                   type=c(rep(str_c("平均漁獲量(",res.list[[1]]$arglist$n.catch-catch.abc.na,"年平均)"),res.list[[i]]$arglist$n.catch),abclegend[i]))
        }
      }
      else {
        if(is.null(abclegend)) {
          if(is.null(catchdividedegit)) data_catch <- rbind(data_catch,tibble(year=last.year+1,
                                                                              catch=c(res.list[[i]]$ABC),
                                                                              type=c(paste0(i,"番目ABC"))))
          else  data_catch <- rbind(data_catch,tibble(year=last.year+1,
                                                      catch=c(res.list[[i]]$ABC/(10^catchdividedegit)),
                                                      type=c(paste0(i,"番目ABC"))))
        }
        else{
          if(is.null(catchdividedegit)) data_catch <- rbind(data_catch,tibble(year=last.year+1,
                                                                              catch=c(res.list[[i]]$ABC),
                                                                              type=c(abclegend[i])))
          else data_catch <- rbind(data_catch,tibble(year=last.year+1,
                                                     catch=c(res.list[[i]]$ABC/(10^catchdividedegit)),
                                                     type=c(abclegend[i])))
        }
      }
    }

    labels2 <- c(labels2,paste0(i,"番目ABC"))
    labels2.1 <- c(labels2.1,paste0(i,"番目算定漁獲量"))
    labels2.2 <- c(labels2.2,paste(i,"番目",max(years)+2,"年",gsub("年","",year.axis.label),"の予測値",sep=""))
  }

  if(!is.null(abclegend)) labels2<- labels2.1<-labels2.2<-abclegend

  if(catch.abc.na!=0) {
    data_catch2 <- data_catch
    data_catch2$catch[which(is.na(ccdata$catch[(length(ccdata$catch)-n.catch+1):length(ccdata$catch)]))] <-NA
  }

  legend.labels2 <-c(str_c("平均漁獲量(",res.list[[1]]$arglist$n.catch-catch.abc.na,"年平均)"),rev(labels2))
  legend.labels2.1 <-c(str_c("平均漁獲量(",res.list[[1]]$arglist$n.catch-catch.abc.na,"年平均)"),rev(labels2.1))
  legend.labels2.2 <-c(str_c("平均漁獲量(",res.list[[1]]$arglist$n.catch-catch.abc.na,"年平均)"),rev(labels2.2))

  if(hcrvlineBan) {
    col.BRP.hcr <- col.BRP
    data_BRP_hcr <- tibble(BRP=names(res.list[[1]]$BRP),value_obs=res.list[[1]]$Obs_BRP, value_ratio=res.list[[1]]$BRP)
    legend.labels.hcr <-c("目標管理基準値（目標水準）","限界管理基準値（限界水準）","禁漁水準")
  }else {
    col.BRP.hcr <- col.BRP[-3]
    data_BRP_hcr <- tibble(BRP=names(res.list[[1]]$BRP[-3]),value_obs=res.list[[1]]$Obs_BRP[-3], value_ratio=res.list[[1]]$BRP[-3])
    legend.labels.hcr <-c("目標管理基準値（目標水準）","限界管理基準値（限界水準）")
  }

  #   # PB=0の時の禁漁水準削除設定 ----
  #   if(res.list[[1]]$BRP[3] == 0) {
  #     if(proposal==TRUE){
  #       legend.labels <- c("目標管理基準値（目標水準）案","限界管理基準値（限界水準）案")
  #     }else{
  #       legend.labels <- c("目標管理基準値（目標水準）","限界管理基準値（限界水準）")
  #     }
  #     linetype.set <- c("22","41")
  #     if(abc4==TRUE){
  #       col.BRP <- c("blue","red")
  #     }else{
  #       col.BRP <- c("#00533E","#edb918")
  #     }
  #     data_BRP2 <- data_BRP
  #     data_BRP <- tibble(BRP=names(res.list[[1]]$BRP[-3]),value_obs=res.list[[1]]$Obs_BRP[-3],value_ratio=res.list[[1]]$BRP[-3])
  #   }else{
  #     if(abc4==TRUE){
  #       col.BRP <- c("blue","red","orange")
  #     }else{
  #       col.BRP <- c("#00533E","#edb918","#C73C2E")
  #     }
  #   }

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

  linetype.set <- c("dashed","dotdash","solid")

  # 資源量指標値のトレンド ----
  if(isTRUE(abc4)){
    g.cpue4<-plot_abc2(res = res.list[[1]],stock.name=stock.name, fishseason=fishseason, abc4=abc4, fillarea=fillarea, cpueunit=cpueunit, RP=RP, leftalign=leftalign, hcrdist=hcrdist)$graph.component[[1]]
    g.cpue4<-plot_abc2(res = res.list[[1]],stock.name=stock.name, fishseason=fishseason, abc4=abc4, fillarea=fillarea, cpueunit=cpueunit, RP=RP, leftalign=leftalign, hcrdist=hcrdist)$graph.component[[2]]
  }else{
    g.cpue <- plot_abc2(res = res.list[[1]],stock.name=stock.name, fishseason=fishseason, abc4=abc4, fillarea=fillarea, cpueunit=cpueunit, RP=RP, leftalign=leftalign, hcrdist=hcrdist)$graph.component[[1]]
  }

  # 漁獲管理規則 HCR ----
  if("arglist"%in%names(res.list)) res.list <- list(res.list)
  if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){
    g.hcr <- ggplot(data=data.frame(X=c(0,100)), aes(x=X)) +
      theme_bw(base_family = font_MAC)+theme_custom()+
      ggtitle("")+
      xlab("資源水準(%)")+ylab(str_c("漁獲量を増減させる係数"))+
      theme(legend.position="top",legend.justification = c(1,0))+
      theme(text = element_text(family = font_MAC))
  }else{
    g.hcr <- ggplot(data=data.frame(X=c(0,100)), aes(x=X)) +
      theme_bw()+theme_custom()+
      ggtitle("")+
      xlab("資源水準(%)")+ylab(str_c("漁獲量を増減させる係数"))+
      theme(legend.position="top",legend.justification = c(1,0))
  }

  Currentalphas<-c()
  size.hcr.points<-c()
  for(i in 1:length(res.list)){
    Currentalphas<-rbind(Currentalphas,tibble(x=res.list[[i]]$Current_Status[1]*100,y=res.list[[i]]$alpha))
    size.hcr.points<- c(size.hcr.points, (5-(i-1)))
  }
  col.hcr.points <- seq(2,1+length(res.list))

  for(i in 1:length(res.list)){
    res <- res.list[[i]]
    if(hcrvlineBan) data_BRP <- tibble(BRP=names(res$BRP),value_obs=res$Obs_BRP,
                                       value_ratio=res$BRP)
    else data_BRP <- tibble(BRP=names(res$BRP[-3]),value_obs=res$Obs_BRP[-3],
                            value_ratio=res$BRP[-3])
    BT <- res$arglist$BT
    PL <- res$arglist$PL
    PB <- res$arglist$PB
    tune.par <- res$arglist$tune.par
    beta <- res$arglist$beta
    empir.dist<- res$arglist$empir.dist
    simple.empir<-res$arglist$simple.empir
    if(is.null(res$arglist$BTyear)) ccdata.plot<-res$arglist$ccdata
    else ccdata.plot<-res$arglist$ccdata[which(res$arglist$ccdata$year <= res$arglist$BTyear),]

    if(!empir.dist) g.hcr <- g.hcr +
      stat_function(fun=type2_func_wrapper,
                    args=list(BT=BT,PL=PL,PB=PB,tune.par=tune.par,beta=beta,AAV=res$AAV,type="%"),
                    color=rgb(0+((i-1)/5),0+((i-1)/5),0+((i-1)/5)),size=1,linetype="solid")
    else g.hcr <- g.hcr +
      stat_function(fun=type2_func_empir_wrapper,
                    args=list(BT=BT,PL=PL,PB=PB,tune.par=tune.par,beta=beta,AAV=res$AAV,cpue=ccdata.plot$cpue,simple=simple.empir,type="%"),
                    color=rgb(0+((i-1)/5),0+((i-1)/5),0+((i-1)/5)),size=1,linetype="solid")
    g.hcr <- g.hcr +
      geom_vline(data=data_BRP,mapping=aes(xintercept=value_ratio*100,color=BRP), size = 1-((i-1)/5), linetype = ifelse(i==1,"solid",i*11))

  }

  if(hcrhscale=="sparse") hlinebreaks <- c(0,0.5,1.0)
  if(hcrhscale=="middle") hlinebreaks <- c(0,0.25,0.5,0.75,1.0)
  if(hcrhscale=="dense") hlinebreaks <- c(0,0.2,0.4,0.6,0.8,1.0)

  if(!plotexactframe) g.hcr <- g.hcr + scale_y_continuous(breaks = hlinebreaks)
  else g.hcr <- g.hcr + scale_x_continuous(expand = c(0,0),limits = c(0,100)) + scale_y_continuous(expand = c(0,0),breaks = hlinebreaks)

  if(hcrhline=="none") hcrAuxiliaryhline <- c()
  if(hcrhline=="one") hcrAuxiliaryhline <- c(1.0)
  if(hcrhline=="sparse") hcrAuxiliaryhline <- c(0,0.5,1.0)
  if(hcrhline=="middle") hcrAuxiliaryhline <- c(0,0.25,0.5,0.75,1.0)
  if(hcrhline=="dense") hcrAuxiliaryhline <- c(0,0.2,0.4,0.6,0.8,1.0)
  if(hcrhline=="hscale") hcrAuxiliaryhline <- hlinebreaks
  g.hcr <- g.hcr +
    geom_hline(yintercept=hcrAuxiliaryhline,color="gray",linetype=2)

  g.hcr <- g.hcr +
    geom_point(data=Currentalphas,aes(x=x,y=y),color=col.hcr.points,size=size.hcr.points)

  if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){
    if(hcrvlineBan) g.hcr <- g.hcr+
        ggrepel::geom_label_repel(data=data_BRP,                                              mapping=aes(x=value_ratio*100, y=c(0.5,1.15,0.8), label=legend.labels.hcr,family = font_MAC),
                                  box.padding=0.5)+
        scale_color_manual(name="",values=rev(c(col.BRP.hcr)),guide="none") #label=rev(legend.labels.hcr))
    else g.hcr <- g.hcr+
        ggrepel::geom_label_repel(data=data_BRP,                                              mapping=aes(x=value_ratio*100, y=c(0.5,1.15), label=legend.labels.hcr,family = font_MAC),
                                  box.padding=0.5)+
        scale_color_manual(name="",values=rev(c(col.BRP.hcr)),guide="none") #label=rev(legend.labels.hcr))
  }else{
    if(hcrvlineBan) g.hcr <- g.hcr+
        ggrepel::geom_label_repel(data=data_BRP,                                              mapping=aes(x=value_ratio*100, y=c(0.5,1.15,0.8), label=legend.labels.hcr),
                                  box.padding=0.5)+
        scale_color_manual(name="",values=rev(c(col.BRP.hcr)),guide="none") #label=rev(legend.labels.hcr))
    else g.hcr <- g.hcr+
        ggrepel::geom_label_repel(data=data_BRP,                                              mapping=aes(x=value_ratio*100, y=c(0.5,1.15), label=legend.labels.hcr),
                                  box.padding=0.5)+
        scale_color_manual(name="",values=rev(c(col.BRP.hcr)),guide="none") #label=rev(legend.labels.hcr))
  }

  # 漁獲量のトレンドとABC ----
  CatchABC<-c(1,rev(seq(2,(length(res.list)+1))))
  if(!ignore_naCatch_point){
    g.catch <- res.list[[1]]$arglist$ccdata %>% ggplot() +
      geom_path(data=data_catch,mapping=aes(x=year,y=catch,color=type),size=2)+
      geom_point(data=data_catch,mapping=aes(x=year,y=catch,color=type),size=3)+
      scale_color_manual(name="",values=rev(CatchABC),labels=rev(legend.labels2))

  }else{
    g.catch <- res.list[[1]]$arglist$ccdata %>% ggplot() +
      geom_path(data=data_catch,mapping=aes(x=year,y=catch,color=type),size=2)+
      geom_point(data=data_catch2,mapping=aes(x=year,y=catch,color=type),size=3)+
      scale_color_manual(name="",values=rev(CatchABC),labels=rev(legend.labels2))
  }

  if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){# plot 設定 for mac
    if(is.null(catchdividedegit)) g.catch <- g.catch +
        geom_path(aes(x=year,y=catch),size=1)+
        ylab(paste("漁獲量",catchunit))+xlab(year.axis.label)+
        ggtitle("")+
        ylim(0,NA)+ theme_custom()+
        theme(legend.position="top",legend.justification = c(1,0)) +
        theme(text = element_text(family = font_MAC))
    else g.catch <- g.catch +
        geom_path(aes(x=year,y=catch/(10^catchdividedegit)),size=1)+
        ylab(paste("漁獲量",catchunit))+xlab(year.axis.label)+
        ggtitle("")+
        ylim(0,NA)+ theme_custom()+
        theme(legend.position="top",legend.justification = c(1,0)) +
        theme(text = element_text(family = font_MAC))
  }else{
    if(is.null(catchdividedegit)) g.catch <- g.catch +
        geom_path(aes(x=year,y=catch),size=1)+
        ylab(paste("漁獲量",catchunit))+xlab(year.axis.label)+
        ggtitle("")+
        ylim(0,NA)+ theme_custom()+
        theme(legend.position="top",legend.justification = c(1,0))
    else g.catch <- g.catch +
        geom_path(aes(x=year,y=catch/(10^catchdividedegit)),size=1)+
        ylab(paste("漁獲量",catchunit))+xlab(year.axis.label)+
        ggtitle("")+
        ylim(0,NA)+ theme_custom()+
        theme(legend.position="top",legend.justification = c(1,0))
  }
  g.catch <- g.catch + geom_point(aes(x = year, y = catch), shape=1, size = 2)
  g.catch <- g.catch %>% apply_minor_ticks_type2()
  # 出力設定 ----
  if(isTRUE(abc4)){
    graph.component <- list(g.cpue4,g.cpue,g.hcr.g.catch)
    graph.combined <- gridExtra::grid.arrange(g.cpue4,g.cpue,g.hcr,g.catch,ncol=2,top=stock.name)
    return(list(graph.component=graph.component,graph.combined=graph.combined))
  }  else{
    graph.component <- list(g.cpue,g.hcr,g.catch)
    graph.combined <- gridExtra::grid.arrange(g.cpue,g.hcr,g.catch,ncol=3,top=stock.name)
    return(list(graph.component=graph.component,graph.combined=graph.combined))
  }
}

#' 2系のABC計算をBTyearオプションありで計算（資源水準導出のためのCPUE時系列制御）した結果をBTyearなしの場合まで連続的に比較してプロットするための関数
#'
#' @param res calc_abc2の返り値、ただしBTyear!=NULL
#' @param fishseason  X軸のラベルを変更（0なら年、1なら漁期年)
#' @param abc4  北海道東部の跨り資源で使用する図を描画（TRUEなら使用、デフォルトはFALSE）
#' @param fillarea  資源量指標値の図にkobeプロットに似た色を塗る（TRUEなら塗る、デフォルトはFALSE）
#' @param cpueunit  資源量指標値の縦軸見出しに追記したい指標値の単位（例えば"（トン/網）"のように指定する）
#' @param catchunit 漁獲量トレンドの縦軸見出しに追記したい単位（例えば"（万トン）"のように指定する、デフォルトは"（トン）"）
#' @param catchdividedegit 漁獲量トレンドの縦軸を任意の桁で割る（例えば"catchdividedegit=2"で100で割るように指定する、デフォルトはNULLでそのまま）
#' @param leftalign  資源量指標値の時系列の長さが漁獲量に比べて短い時、データが無い範囲の空間を削除する（TRUEなら使用、デフォルトはFALSE）
#' @param RP  資源量指標値/年のプロットでReference Point（目標・限界管理基準線）を載せる・載せない（デフォルトはTRUE、FALSEでは直近年の資源量指標値をポイントでハイライトする）
#' @export
#'

plot_abc2_fixTerminalCPUE_seqOut <- function(res, stock.name=NULL, fishseason=0, abc4=FALSE, fillarea=FALSE, cpueunit="", catchunit="（トン）",RP=TRUE, leftalign=FALSE, hcrdist=FALSE){
  # abc4は北海道東部海域の「跨り資源」で資源量指標値の平均水準・過去最低値を描画する際に使用する。その際、calc_abc2の引数BTは0.5に設定すること。

  # 漁期年/年設定 ----
  ifelse(fishseason==1, year.axis.label <- "漁期年", year.axis.label <- "年")
  # Setting cpue dist
  empir.dist <- res$arglist$empir.dist
  simple.empir <- res$arglist$simple.empir
  smooth.cpue <- res$arglist$smooth.cpue
  smooth.dist <- res$arglist$smooth.dist
  # plot
  ccdata <- res$arglist$ccdata
  ccdata_fixedBT <- res$arglist$ccdata
  BTyear <- res$arglist$BTyear
  if(is.null(BTyear)) stop("This function works if BTyear was set in calc_abc2.\n")
  ccdata_fixedBT <- ccdata[which(ccdata$year <= BTyear),]
  n.catch <- res$arglist$n.catch
  years <- ccdata$year
  last.year <- rev(years)[1]

  BT <- res$arglist$BT
  PL <- res$arglist$PL
  PB <- res$arglist$PB
  tune.par <- res$arglist$tune.par
  beta <- res$arglist$beta

  res.multiBTyear<-list()
  ABCs<-c()
  ABClabels<-c()
  ccdata.plotbt<-list()
  for(i in 0:(last.year-BTyear-1)){
    if(i==0) res.multiBTyear[[i+1]]<- calc_abc2(ccdata=res$arglist$ccdata,BT=BT,PL=PL,PB=PB,tune.par = tune.par, AAV=res$arglist$AAV,n.catch=res$arglist$n.catch,n.cpue=res$arglist$n.catch,smooth.cpue = res$arglist$smooth.cpue,smooth.dist = res$arglist$smooth.dist,empir.dist = res$arglist$empir.dist,simple.empir = res$arglist$simple.empir,beta = res$arglist$beta,D2alpha = res$arglist$D2alpha,BTyear = NULL,summary_abc=FALSE)
    else res.multiBTyear[[i+1]]<- calc_abc2(ccdata=res$arglist$ccdata,BT=BT,PL=PL,PB=PB,tune.par = tune.par, AAV=res$arglist$AAV,n.catch=res$arglist$n.catch,n.cpue=res$arglist$n.catch,smooth.cpue = res$arglist$smooth.cpue,smooth.dist = res$arglist$smooth.dist,empir.dist = res$arglist$empir.dist,simple.empir = res$arglist$simple.empir,beta = res$arglist$beta,D2alpha = res$arglist$D2alpha,BTyear = last.year-i,summary_abc=FALSE)
    ABCs<-c(ABCs,res.multiBTyear[[i+1]]$ABC)
    label <-paste0(i,"年前基準ABC")
    ABClabels<-c(ABClabels,label)
    ifelse(i==0,ccdata.plotbt[[i+1]]<-res.multiBTyear[[i+1]]$arglist$ccdata,ccdata.plotbt[[i+1]]<-res.multiBTyear[[i+1]]$arglist$ccdata[which(res.multiBTyear[[i+1]]$arglist$ccdata$year <= res.multiBTyear[[i+1]]$arglist$BTyear),])
  }

  data_catch_ori <- tibble(year=c((last.year-res$arglist$n.catch+1):last.year,last.year+2),catch=c(rep(res$mean.catch,res$arglist$n.catch),res$ABC),                              type=c(rep(str_c("平均漁獲量(",res$arglist$n.catch,"年平均)"),n.catch),"ABC"))

  btlabel<-paste0(last.year-BTyear,"年前基準ABC")

  data_catch <- tibble(year=c((last.year-res$arglist$n.catch+1):last.year,rep(last.year+2,last.year-BTyear+1)),catch=c(rep(res$mean.catch,res$arglist$n.catch),res$ABC,rev(ABCs)),                              type=c(rep(str_c("平均漁獲量(",res$arglist$n.catch,"年平均)"),n.catch),btlabel,rev(ABClabels)))

  data_BRP <- tibble(BRP=names(res$BRP),value_obs=res$Obs_BRP,
                     value_ratio=res$BRP)
  data_percent <- tibble(x=rep(max(years)+2,11),
                         y=res$Obs_percent,
                         label=str_c(c(0.05,seq(from=0.1,to=0.9,by=0.1),0.95)*100,"%"))
  data_percent_even <- tibble(x=rep(max(years)+2,6),
                              y=res$Obs_percent_even,
                              label=str_c(c(0.05,seq(from=0.2,to=0.8,by=0.2),0.95)*100,"%"))
  font_MAC <- "HiraginoSans-W3"#"Japan1GothicBBB"#

  legend.labels <-c("目標管理基準値（目標水準）","限界管理基準値（限界水準）","禁漁水準")

  linetype.set <- c("dashed","dotdash","solid")

  legend.labels2 <-c(str_c("平均漁獲量(",res$arglist$n.catch,"年平均)"),"ABC")
  legend.labels2bt <-c(str_c("平均漁獲量(",res$arglist$n.catch,"年平均)"),btlabel,rev(ABClabels))

  col.BRP.hcr <- col.BRP
  data_BRP_hcr <- tibble(BRP=names(res$BRP),value_obs=res$Obs_BRP, value_ratio=res$BRP)

  # PB=0の時の禁漁水準削除設定 ----
  if(res$BRP[3] == 0) {

    legend.labels <- c("目標管理基準値（目標水準）","限界管理基準値（限界水準）")

    linetype.set <- c("dashed","dotdash")
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
  g.catch.title <- "漁獲量のトレンドとABC"
  g.catch.abcpoint <- "ABC"
  legend.labels2 <- legend.labels2bt

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

  if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){ # plot 設定 for mac----
    g.cpue <- ccdata %>% ggplot() +
      geom_polygon(data=tibble(x=c(minyears,max(years)+4,max(years)+4,minyears), y=c(data_BRP2$value_obs[1],data_BRP2$value_obs[1],max(ccdata$cpue,na.rm=T)*1.05,max(ccdata$cpue,na.rm=T)*1.05)), aes(x=x,y=y), fill=colfill[1]) +
      geom_polygon(data=tibble(x=c(minyears,max(years)+4,max(years)+4,minyears), y=c(data_BRP2$value_obs[2],data_BRP2$value_obs[2],data_BRP2$value_obs[1],data_BRP2$value_obs[1])), aes(x=x,y=y), fill=colfill[2]) +
      geom_polygon(data=tibble(x=c(minyears,max(years)+4,max(years)+4,minyears), y=c(data_BRP2$value_obs[3],data_BRP2$value_obs[3],data_BRP2$value_obs[2],data_BRP2$value_obs[2])), aes(x=x,y=y), fill=colfill[3]) +
      geom_polygon(data=tibble(x=c(minyears,max(years)+4,max(years)+4,minyears), y=c(0,0,data_BRP2$value_obs[3],data_BRP2$value_obs[3])), aes(x=x,y=y), fill=colfill[4]) +
      geom_hline(yintercept=res$Obs_percent_even,color="gray",linetype=2)+
      geom_text(data=data_percent_even,aes(x=x,y=y*1.05,label=label))+
      geom_text(aes(x=max(years)-1,y=min(data_percent_even$y)*0.75,family=font_MAC,label="(資源水準)"),size=4)
    if(RP==TRUE){
      g.cpue <- g.cpue +
        geom_hline(data=data_BRP, mapping=aes(yintercept=value_obs, color=rev(col.BRP), linetype=rev(linetype.set)), size = 0.9*1.5)+
        #ggrepel::geom_label_repel(data=data_BRP, mapping=aes(x=min(years)+0.5, y=value_obs+0.5, label=legend.labels,family = font_MAC), box.padding=0.5, nudge_x=1)+
        scale_linetype_manual(name="", values=rev(c(linetype.set)), labels=rev(c(legend.labels))) +
        scale_color_manual(name="",values=rev(c(col.BRP)),labels=rev(c(legend.labels)))
    }else{
      g.cpue <- g.cpue +
        geom_point(mapping=aes(x=rev(year)[1], y=rev(ccdata$cpue)[1], color="red"),size=4, show.legend =TRUE)+
        scale_color_manual(name="",values="red",labels="直近年の資源量指標値")
    }
    g.cpue <- g.cpue +
      geom_path(aes(x=year,y=cpue),size=1)+
      theme_bw()+ylab(paste("資源量指標値",cpueunit))+xlab(year.axis.label)+
      ylim(0,max(ccdata$cpue,na.rm=T)*1.05)+theme_custom()

    g.cpue <- g.cpue +
      ggtitle("")+
      theme(legend.position="top",legend.justification = c(1,0), legend.spacing=unit(0.25,'lines'), legend.key.width = unit(2.0, 'lines')) +
      theme(text = element_text(family = font_MAC))
    if(leftalign==TRUE){
      g.cpue <- g.cpue + xlim(minyears,max(ccdata[!is.na(ccdata$cpue),]$year)+4)
    }
  }else{
    g.cpue <- ccdata %>% ggplot() +
      geom_polygon(data=tibble(x=c(minyears,max(years)+4,max(years)+4,minyears), y=c(data_BRP2$value_obs[1],data_BRP2$value_obs[1],max(ccdata$cpue,na.rm=T)*1.05,max(ccdata$cpue,na.rm=T)*1.05)), aes(x=x,y=y), fill=colfill[1]) +
      geom_polygon(data=tibble(x=c(minyears,max(years)+4,max(years)+4,minyears), y=c(data_BRP2$value_obs[2],data_BRP2$value_obs[2],data_BRP2$value_obs[1],data_BRP2$value_obs[1])), aes(x=x,y=y), fill=colfill[2]) +
      geom_polygon(data=tibble(x=c(minyears,max(years)+4,max(years)+4,minyears), y=c(data_BRP2$value_obs[3],data_BRP2$value_obs[3],data_BRP2$value_obs[2],data_BRP2$value_obs[2])), aes(x=x,y=y), fill=colfill[3]) +
      geom_polygon(data=tibble(x=c(minyears,max(years)+4,max(years)+4,minyears), y=c(0,0,data_BRP2$value_obs[3],data_BRP2$value_obs[3])), aes(x=x,y=y), fill=colfill[4]) +
      geom_hline(yintercept=res$Obs_percent_even,color="gray",linetype=2)+
      geom_text(data=data_percent_even, aes(x=x,y=y*1.05,label=label))+
      geom_text(aes(x=max(years)-1,y=min(data_percent_even$y)*0.75,label="(資源水準)"),size=4)
    if(RP==TRUE){
      g.cpue <- g.cpue +
        geom_hline(data=data_BRP, mapping=aes(yintercept=value_obs, color=rev(col.BRP), linetype=rev(linetype.set)), size = 0.9*1.5)+
        #ggrepel::geom_label_repel(data=data_BRP, mapping=aes(x=min(years)+0.5, y=value_obs+0.5, label=legend.labels), box.padding=0.5, nudge_x=1)+
        scale_linetype_manual(name="", values=rev(c(linetype.set)), labels=rev(c(legend.labels))) +
        scale_color_manual(name="",values=rev(c(col.BRP)),labels=rev(c(legend.labels)))
    }else{
      g.cpue <- g.cpue +
        geom_point(mapping=aes(x=rev(year)[1], y=rev(ccdata$cpue)[1], color="red"),size=4, show.legend =TRUE)+
        scale_color_manual(name="",values="red",labels="直近年の資源量指標値")
    }
    g.cpue <- g.cpue +
      geom_path(data=ccdata, aes(x=year,y=cpue),size=1)+
      theme_bw()+ylab(paste("資源量指標値",cpueunit))+xlab(year.axis.label)+
      ylim(0,max(ccdata$cpue,na.rm=T)*1.05)+theme_custom()
    g.cpue <- g.cpue +
      ggtitle("") + theme(legend.position="top", legend.spacing=unit(0.25,'lines'), legend.key.width = unit(2.0, 'lines'),legend.justification=c(1,0))
    if(leftalign==TRUE){
      g.cpue <- g.cpue + xlim(minyears, max(ccdata[!is.na(ccdata$cpue),]$year)+4)
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

    if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){ ## plot 設定 for mac----
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

  #漁獲管理規則案 HCR ----
  ifelse(is.null(BTyear),ccdata.plot<-ccdata,ccdata.plot<-ccdata_fixedBT)
  if(!empir.dist){
    g.hcr <- ggplot(data=data.frame(X=c(0,100)), aes(x=X)) +
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
      xlab("資源水準(%)")+ylab(str_c("漁獲量を増減させる係数"))+
      theme(legend.position="top",legend.justification = c(1,0))

    ## plot 設定 for mac----
    if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){
      g.hcr <- ggplot(data=data.frame(X=c(0,100)), aes(x=X)) +
        #stat_function(fun=type2_func_wrapper,
        #        args=list(BT=BT,PL=0,PB=PB,tune.par=tune.par,beta=beta,AAV=res$AAV,type="%"),
        #                  color="gray")+
        stat_function(fun=type2_func_wrapper,
                      args=list(BT=BT,PL=PL,PB=PB,tune.par=tune.par,beta=beta,AAV=res$AAV,type="%"),
                      color="black",size=1)+
        geom_point(aes(x=res$Current_Status[1]*100,y=res$alpha),color="red",size=4) +

        geom_vline(data=data_BRP,mapping=aes(xintercept=value_ratio*100,color=BRP), size = 0.9*1.5, linetype = linetype.set)+
        ggrepel::geom_label_repel(data=data_BRP,
                                  mapping=aes(x=value_ratio*100, y=1.1, label=legend.labels,family = font_MAC),
                                  box.padding=0.5, nudge_y=1)+
        scale_color_manual(name="",values=rev(c(col.BRP)), guide="none")+ #,labels=rev(c(legend.labels)))+
        theme_bw()+theme_custom()+
        ggtitle("")+
        xlab("資源水準(%)")+ylab(str_c("漁獲量を増減させる係数"))+
        theme(legend.position="top",legend.justification = c(1,0)) +
        theme(text = element_text(family = font_MAC))
    }

    for(i in 1:(last.year-BTyear)){
      g.hcr <- g.hcr +
        stat_function(fun=type2_func_wrapper,
                      args=list(BT=BT,PL=PL,PB=PB,tune.par=tune.par,beta=beta,AAV=res.multiBTyear[[i]]$AAV,type="%"), color="gray",size=0.5,linetype=i+1)+
        geom_point(aes(x=res.multiBTyear[[i]]$Current_Status[1]*100,y=res.multiBTyear[[i]]$alpha),color=i+2,size=3,shape=i+1)
    }

  }else{
    # empir.dist=T ----
    g.hcr <-ggplot(data=data.frame(X=c(0,100)), aes(x=X)) +
      #stat_function(fun=type2_func_wrapper,
      #        args=list(BT=BT,PL=0,PB=PB,tune.par=tune.par,beta=beta,AAV=res$AAV,type="%"),
      #              color="gray")+
      stat_function(fun=type2_func_empir_wrapper,
                    args=list(BT=BT,PL=PL,PB=PB,tune.par=tune.par,beta=beta,AAV=res$AAV,cpue=ccdata.plot$cpue,simple=simple.empir,type="%"),
                    color="black",size=1) +
      geom_point(aes(x=res$Current_Status[1]*100,y=res$alpha),color="red",size=4)+
      geom_vline(data=data_BRP,mapping=aes(xintercept=value_ratio*100,color=BRP), size = 0.9*1.5, linetype = linetype.set)+
      ggrepel::geom_label_repel(data=data_BRP,
                                mapping=aes(x=value_ratio*100, y=1.1, label=legend.labels,family = font_MAC),
                                box.padding=0.5, nudge_y=1)+
      scale_color_manual(name="",values=rev(c(col.BRP)), guide="none" )+ #,labels=rev(c(legend.labels)))+
      theme_bw()+theme_custom()+
      ggtitle("")+
      xlab("資源水準(%)")+ylab(str_c("漁獲量を増減させる係数"))+
      theme(legend.position="top",legend.justification = c(1,0))

    ## plot 設定 for mac----
    if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){
      g.hcr <-ggplot(data=data.frame(X=c(0,100)), aes(x=X)) +
        #stat_function(fun=type2_func_wrapper,
        #        args=list(BT=BT,PL=0,PB=PB,tune.par=tune.par,beta=beta,AAV=res$AAV,type="%"),
        #              color="gray")+
        stat_function(fun=type2_func_empir_wrapper,
                      args=list(BT=BT,PL=PL,PB=PB,tune.par=tune.par,beta=beta,AAV=res$AAV,cpue=ccdata.plot$cpue,simple=simple.empir,type="%"),
                      color="black",size=1) +
        geom_point(aes(x=res$Current_Status[1]*100,y=res$alpha),color="red",size=4)+
        geom_vline(data=data_BRP,mapping=aes(xintercept=value_ratio*100,color=BRP), size = 0.9*1.5, linetype = linetype.set)+
        ggrepel::geom_label_repel(data=data_BRP,
                                  mapping=aes(x=value_ratio*100, y=1.1, label=legend.labels,family = font_MAC),
                                  box.padding=0.5, nudge_y=1)+
        scale_color_manual(name="",values=rev(c(col.BRP)), guide="none" )+ #,labels=rev(c(legend.labels)))+
        theme_bw()+theme_custom()+
        ggtitle("")+
        xlab("資源水準(%)")+ylab(str_c("漁獲量を増減させる係数"))+
        theme(legend.position="top",legend.justification = c(1,0)) +
        theme(text = element_text(family = font_MAC))

    }

    for(i in 1:(last.year-BTyear)){
      g.hcr <- g.hcr +
        stat_function(fun=type2_func_empir_wrapper,
                      args=list(BT=BT,PL=PL,PB=PB,tune.par=tune.par,beta=beta,AAV=res.multiBTyear[[i]]$AAV,cpue=ccdata.plotbt[[i]]$cpue,simple=simple.empir,type="%"), color="gray",size=0.5,linetype=i+1) +
        geom_point(aes(x=res.multiBTyear[[i]]$Current_Status[1]*100,y=res.multiBTyear[[i]]$alpha),color=i+2,size=3,shape=i+1)#color=rgb(0,0,0,alpha = (i/10)),size=0.5,linetype="dashed")+
    }
  }

  #漁獲管理規則案 HCR.Dist ----
  current_index_col <- "#1A4472"
  ccdata.plot<- ccdata_fixedBT
  model_dist<-c()
  ifelse( floor(log10(max(ccdata.plot$cpue,na.rm = T))+1) < 4, model_dist <- data.frame(cpue=seq(0, max(ccdata.plot$cpue,na.rm=T), by=0.1),  dens=NA), model_dist <- data.frame(cpue=seq(0, max(ccdata.plot$cpue,na.rm=T), by=10^floor(log10(max(ccdata.plot$cpue,na.rm = T))+1)/1000),  dens=NA) )
  if(!empir.dist) model_dist$dens <- dnorm(model_dist$cpue,mean = mean(ccdata.plot$cpue,na.rm=T),sd=sd(ccdata.plot$cpue,na.rm=T))
  else{ # empir.dist = T で累積確率から個々の確率を求めて総和(1)で割って密度にする
    if(!simple.empir){
      cum.cpue4 <- ecdf(ccdata.plot$cpue)
      cum.probs <- cum.cpue4(model_dist$cpue)
    }
    else{
      cum.probs<-c()
      for(i in 1:length(model_dist$cpue)){
        cum.probs <- c(cum.probs,simple_ecdf(ccdata.plot$cpue,model_dist$cpue[i]))
      }
    }
    probs<-c(cum.probs[1])
    for(j in 2:length(cum.probs)){
      if(cum.probs[j-1]!=cum.probs[j]) tmp <- cum.probs[j]-cum.probs[j-1]
      else tmp <- probs[j-1]
      probs<-c(probs,tmp)
    }
    #plot(model_dist$cpue,probs)
    model_dist$dens<-probs/sum(probs)
  }

  g.hcr.dist <- ggplot(data=model_dist)+
    #stat_function(fun=dnorm,args=list(mean=mean(ccdata.plot$cpue),sd=sd(ccdata.plot$cpue)),color="black",size=1)
    geom_line(aes(x=cpue,y=dens))+
    geom_area(data=filter(model_dist, cpue < res$Current_Status[2]), aes(x=cpue, y=dens), fill="grey")

  g.hcr.dist <-  g.hcr.dist +
    geom_vline(data=data_BRP,mapping=aes(xintercept=value_obs,color=BRP), size = 0.9*1.5, linetype = linetype.set) +
    scale_linetype_manual(name="", values=rev(c(linetype.set)), labels=rev(c(legend.labels))) +
    scale_color_manual(name="",values=rev(c(col.BRP)),labels=rev(c(legend.labels)))+
    guides(colour="none")+
    coord_flip()

  if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){ # plot 設定 for mac----
    g.hcr.dist <- g.hcr.dist +
      geom_vline(data=data_BRP,mapping=aes(xintercept=res$Current_Status[2]),color=current_index_col,size=1,linetype="dashed")+
      geom_text(aes(x=ifelse(res$Current_Status[2]<mean(ccdata.plot$cpue)/3,mean(ccdata.plot$cpue)/2,mean(ccdata.plot$cpue)/4),y=max(dens)*0.85,family=font_MAC,label="(現在の資源水準)"),color=current_index_col,size=4)
  }else{
    g.hcr.dist <- g.hcr.dist +
      geom_vline(data=data_BRP,mapping=aes(xintercept=res$Current_Status[2]),color=current_index_col,size=1,linetype="dashed")+
      geom_text(aes(x=ifelse(res$Current_Status[2]<mean(ccdata.plot$cpue)/3,mean(ccdata.plot$cpue)/2,mean(ccdata.plot$cpue)/4),y=max(dens)*0.85,label="(現在の資源水準)"),color=current_index_col,size=4)
  }

  if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){ # plot 設定 for mac----
    g.hcr.dist <- g.hcr.dist +  ggtitle("")+
      scale_x_continuous(limits=c(0,max(ccdata$cpue,na.rm=T)*1.05)) +
      #scale_y_continuous(limits=c(0,max(ccdata$cpue,na.rm=T)*1.05)) +
      xlab("資源量指標値")+ylab("")+
      theme_bw()+theme_custom()+
      theme(legend.position="top",legend.justification = c(1,0), legend.spacing=unit(0.25,'lines'), legend.key.width = unit(2.0, 'lines'),axis.text.x = element_blank()) +
      theme(text = element_text(family = font_MAC))
  }else{
    g.hcr.dist <- g.hcr.dist +  ggtitle("")+
      scale_x_continuous(limits=c(0,max(ccdata$cpue,na.rm=T)*1.05)) +
      #scale_y_continuous(limits=c(0,max(ccdata$cpue,na.rm=T)*1.05)) +
      xlab("資源量指標値")+ylab("")+
      theme_bw()+theme_custom()+
      theme(legend.position="top",legend.justification = c(1,0), legend.spacing=unit(0.25,'lines'), legend.key.width = unit(2.0, 'lines'),axis.text.x = element_blank())
  }

  #漁獲量のトレンドとABC/算定漁獲量 ----

  CatchABC<-seq((last.year-BTyear+2):1)
  #CatchABC<-c(1,2)
  g.catch <- ccdata %>% ggplot() +
    geom_path(data=data_catch,mapping=aes(x=year,y=catch,color=type),size=2)+
    geom_point(data=data_catch,mapping=aes(x=year,y=catch,color=type),size=3)+
    #ggrepel::geom_label_repel(data=data_catch,
    #                          mapping=aes(x=max(year)-5, y=catch, label=legend.labels2),
    #                          box.padding=0.5, nudge_y=1)+
    #scale_color_manual(name="",values=c(1,2),labels=legend.labels2)+
    scale_color_manual(name="",values=rev(CatchABC),labels=rev(legend.labels2))+
    # geom_point(data=dplyr::filter(data_catch,type=="ABC"),
    #                    mapping=aes(x=year,y=catch),lwd=2,color=1)+
    #         geom_line(data=dplyr::filter(data_catch,type!="ABC"),
    #                    mapping=aes(x=year,y=catch),lwd=2,color="gray")+
    geom_path(aes(x=year,y=catch),size=1)+
    ylab(paste("漁獲量",catchunit))+xlab(year.axis.label)+
    ggtitle("")+
    ylim(0,NA)+ theme_custom()+
    theme(legend.position="top",legend.justification = c(1,0))

  if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){# plot 設定 for mac
    g.catch <- ccdata %>% ggplot() +
      geom_path(data=data_catch,mapping=aes(x=year,y=catch,color=type),size=2)+
      geom_point(data=data_catch,mapping=aes(x=year,y=catch,color=type),size=3)+
      #ggrepel::geom_label_repel(data=data_catch,
      #                          mapping=aes(x=max(year)-5,y=catch[1],label=legend.labels2,family=font_MAC),
      #                          box.padding=0.5, nudge_y=1)+
      #scale_color_manual(name="",values=c("black","red"),labels=legend.labels2)+
      scale_color_manual(name="",values=rev(CatchABC),labels=rev(legend.labels2))+
      # geom_point(data=dplyr::filter(data_catch,type=="ABC"),
      #                    mapping=aes(x=year,y=catch),lwd=2,color=1)+
      #         geom_line(data=dplyr::filter(data_catch,type!="ABC"),
      #                    mapping=aes(x=year,y=catch),lwd=2,color="gray")+
      geom_path(aes(x=year,y=catch),size=1)+
      ylab(paste("漁獲量",catchunit))+xlab(year.axis.label)+
      ggtitle("")+
      ylim(0,NA)+ theme_custom()+
      theme(legend.position="top",legend.justification = c(1,0)) +
      theme(text = element_text(family = font_MAC))
  }
  g.catch <- g.catch + geom_point(aes(x = year, y = catch), shape=1, size = 2)
  g.catch <- g.catch %>% apply_minor_ticks_type2()

  # 出力設定 ----
  #if(outABCs) print(ABCs)
  outABC <- data.frame(label=c(ABClabels,btlabel),ABC=c(ABCs,res$ABC))
  if(isTRUE(hcrdist)){
    if(isTRUE(abc4)){
      graph.component <- list(g.cpue4,g.cpue,g.hcr.dist,g.hcr,g.catch)
      graph.combined <- gridExtra::grid.arrange(g.cpue4,g.cpue,g.hcr.dist,g.hcr,g.catch,ncol=3,top=stock.name)
      return(list(ABC=outABC,graph.component=graph.component,graph.combined=graph.combined))
    }else{
      graph.component <- list(g.cpue,g.hcr.dist,g.hcr,g.catch)
      graph.combined <- gridExtra::grid.arrange(g.cpue,g.hcr.dist,g.hcr,g.catch,ncol=2,top=stock.name)
      return(list(ABC=outABC,graph.component=graph.component,graph.combined=graph.combined))
    }
  }else{
    if(isTRUE(abc4)){
      graph.component <- list(g.cpue4,g.cpue,g.hcr,g.catch)
      graph.combined <- gridExtra::grid.arrange(g.cpue4,g.cpue,g.hcr,g.catch,ncol=2,top=stock.name)
      return(list(ABC=outABC,graph.component=graph.component,graph.combined=graph.combined))
    }else{
      graph.component <- list(g.cpue,g.hcr,g.catch)
      graph.combined <- gridExtra::grid.arrange(g.cpue,g.hcr,g.catch,ncol=3,top=stock.name)
      return(list(ABC=outABC,graph.component=graph.component,graph.combined=graph.combined))
    }
  }
}

#'
#' サブ目盛りの追加
#'
#' 2系ルールの図に最適化された設定です
#'
#' @export
#'

apply_minor_ticks_type2 <- function(plot, minor_breaks=1){
  plot +   # サブ目盛の設定
    guides(x=guide_axis(minor.ticks=TRUE), # guideは凡例を制御するための関数。目盛りのスタイルを設定するのはtheme関数だが、どんな目盛りをつけるかはguidesの範疇になる？
           y=guide_axis(minor.ticks=TRUE)) +
    # サブ目盛りをつけるので、目盛りの長さを少し長くし、線幅を狭くする
    theme(axis.ticks.length=unit(0.3,"cm"), # default=unit(0.15,"cm")
          axis.ticks=element_line(linewidth=0.4), # default=0.5
          axis.minor.ticks.length = rel(0.5)) +
    # サブ目盛りの間隔の設定
    scale_x_continuous(minor_breaks=scales::breaks_width(minor_breaks),
                       breaks      =scales::breaks_pretty())
  #    scale_x_continuous(minor_breaks=scales::breaks_width(minor_breaks))
}
