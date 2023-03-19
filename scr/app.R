#************************
#********Hydro-Model*****
#************************


library(shiny)
library(ggplot2)
library(ggthemes)
library(shinythemes)
library(RColorBrewer)
library(Cairo)   # For nicer ggplot2 output when deployed on Linux
#----------------全局函数------------

#-----三层蒸散发以及蓄满产流模型-----
Dunne <- function(P,E,DT,Imp,Kc,WM,WUM,WLM,B,C,Ex,SM,Ki,Kg,FR0){
  # DT = 24
  # Imp = 0.02
  # Kc = 0.56
  # WM = 180
  # WUM = 40
  # WLM = 40
  # B = 1.2
  # C = 0.12
  # Ex = 1.2
  # SM = 12
  # Ki = 0.4
  # Kg = 0.3
  # FR0 = 0.002
  
  #------------------------------
  # XAJ是新安江的运行程序,用于单纯形和遗传算法调用,也用于新安江模型的预报 
  # Imp   流域不透水面积比
  # P      降水，输入量
  # E      水面蒸发，输入量
  # DAREA  流域面积
  # DT     计算时段
  
  # Kc      流域蒸散发折算系数多年总径流量决定 
  # WUM     流域上层蓄水容量 
  # WLM     流域下层蓄水容量  
  # WDM     流域深层蓄水容量 
  # WM = WUM + WLM + WDM    WM:流域平均蓄水容量 
  # B       流域蓄水容量分布曲线指数 
  # C       流域深层蒸发系数 
  # Ex      流域自由水分布曲线指数 
  # SM      流域自由水平均蓄水容量 
  # Ki      自由水箱壤中流出流系数（日） 
  # Kg      自由水箱地下水出流系数（日）
  # FR0     初始时刻产流面积比例
  
  #次洪决定WMBImp 
  #由于日模型与次洪模型的计算时段长不同，参数值不能全部通用，但K、WM、WUM、WLM、B、IMP、EX、C与时段长无关可以直接引用 
  #Kc SM、Kg、KSS、Cs、Ci、Cg与时段长相关不能直接引用需要另外率定 
  
  
  D=24/DT  #可以划分的时段数，一天可以划分的时段数
  
  KSSD = (1-(1-(Ki+Kg))^(1/D))/(1+Kg/Ki) #壤中流时段出流系数；包为民主编《水文预报（第四版）》P149，式5-33
  KGD = KSSD*Kg/Ki    #地下径流时段出流系数，P149，式5-33
  
  Epp=Kc*E   #Epp为蒸发能力 
  
  EU = NULL  #上层蒸发量
  EL = NULL  #下层蒸发量
  ED = NULL  #深层蒸发量
  WL = NULL  #下层土壤含水量
  WU = NULL  #上层土壤含水量
  WD = NULL  #深层土壤含水量
  W = NULL   #流域平均土壤含水量,中间时段量 
  R = NULL  #产流量，mm
  
  WDM = WM-WUM-WLM
  
  RS = NULL  #地表径流出流径流深，mm
  RSS = NULL #壤中流出流径流深，mm
  RG = NULL  #地下径流出流径流深，mm
  PE = NULL
  WU[1] = WUM/4  #设置土壤含水的初始状态
  WL[1] = WLM/9*7
  WD[1] = WDM
  W[1] = WU[1]+WL[1]+WD[1]
  len = length(P)
  S0 = SM/4  #初始自由水箱含水量（深），mm
  a = NULL
  for (t in 1:len){         
    #t以时段为单位计算,len为数据系列长度 
    ##-------------三层蒸散发计算------------- 
    if ((WU[t]  + P[t]) >= Epp[t] ){
      EU[t] = Epp[t]  #上层蒸发                                            #Epp为EM 
      EL[t] = 0    #中层 
      ED[t] = 0    #下层 
      
    } else {
      EU[t] = WU[t]  + P[t]                      #Ww[1] + P为EU 
      EL[t] = (Epp[t] - EU[t])*WL[t]/WLM       #要求计算的下层蒸发量与剩余蒸散发能力之比不小于深层蒸散发系数c 
      ED[t] = 0 
      if (WL[t] <= C*WLM) { #第二层水量小于蒸散发能力 
        if( WL[t] >= C * (Epp[t]-EU[t])) {     #'要求计算的下层蒸发量与剩余蒸散发能力之比小于深层蒸散发系数c 
          EL[t] = C * (Epp[t]-EU[t])            
          ED[t] = 0 
        }  else {
          EL[t] = WL[t] 
          ED[t] = 0
          if ( WL[t] < C*(Epp[t]-EU[t]) ) {
            ED[t] = C*(Epp[t]-EU[t])-EL[t]
          } else {
            ED[t] = WD[t]
          }
        }
      } 
    }             
    PE[t] = P[t]-EU[t]-EL[t]-ED[t]  #净雨:要么有净雨，要么为0；
    if (PE[t] < 0) {PE[t] = 0}
    #-----------产流计算部分--------------
    
    R[t] = 0       #总产流量预设为0，mm
    Rimp = 0       #不透水面积产生的直接径流量，mm
    WMM = (1+B)*WM/(1-Imp)    # Imp不透水面积比;
    #WMM为蓄水容量极值  
    #B为蓄水容量曲线的指数
    if( PE[t]>0 ) {                
      Rimp = PE[t]*Imp          #不透水面积产生的直接径流,Imp < 1 
      PE[t] = PE[t]-Rimp        #对PE[t]作出调整，透水面积上的净雨
      if (abs(WM - W[t]) <= 0.01 ){  
        a[t] = WMM       #WMM为蓄水容量极值,单点最大值
      } else {
        a[t] = WMM*(1-(1-W[t]/WM)^(1/(1 + B)))     #a为与W[t]对应的在蓄水容量曲线的纵坐标 
      }
      
      if ((PE[t] + a[t]) <= WMM){           #流域部分产流   
        R[t] = PE[t]-WM+W[t] + WM * ((1- (PE[t] + a[t])/ WMM)^(1 + B)) 
      } else {
        R[t] = PE[t] - (WM - W[t])     #流域全面产流 
      }
      
      if (abs(R[t] - PE[t]) <= 0.01) {
        R[t] = PE[t] 
      } 
      
      WU[t+1]  = WU[t]+ P[t]-R[t]- EU[t]   #第一层蓄水变化 
      WL[t+1]  = WL[t]- EL[t]              #第二层蓄水变化 
      WD[t+1]  = WD[t]- ED[t]             #第三层蓄水变化 
    }  else  { #PE < 0时，净雨深为负时，产流为0；
      Rimp = 0
      R[t] = 0
      a[t] = WMM*(1-(1-W[t]/WM)^(1/(1 + B)))
      WU[t+1]  = WU[t]- EU[t]       #第一层蓄水变化 
      WL[t+1]  = WL[t]- EL[t]             #第二层蓄水变化 
      WD[t+1]  = WD[t]- ED[t]              #第三层蓄水变化 
    }
    #water balance correction
    if (WU[t+1] > WUM)                
    {
      WL[t+1] = WL[t]+WU[t]-WUM     
      WU[t+1] = WUM 
    } 
    if (WU[t+1] <0 ){
      W[t+1] = 0
    }
    if( WL[t+1] > WLM)  
    {
      WD[t+1]  = WD[t]+WL[t]-WLM 
      WL[t+1] = WLM 
    } else if ( WL[t+1] <0 ) {
      WL[t+1] = 0
    }
    
    if (WD[t+1] < 0) {
      WD[t+1]  = 0 
    } else if (WD[t+1] > WDM) {
      WD[t+1] = WDM
      R[t] = WD[t+1] - WDM + R[t]
    } 
    W[t+1] = WU[t+1]+WL[t+1]+WD[t+1]
    if (W[t+1] > WM) {
      R[t] = W[t+1]-WM + R[t]
      W[t+1] = WM
    } else if (W[t+1] < 0) {
      W[t+1] = 0
    }
    
#-----------------水源划分--------------- 
    X = FR0                           #FR0为时段初产流面积比例 
    if (PE[t] <= 0)                   #无净雨   
    { 
      RS[t] = 0
      RSS[t] = S0*KSSD*FR0                                                                                                                                                                                                                                                                                                                                      #依靠时段初自由水箱含水量计算出流；      
      #地下水出流出流系数 
      RG[t] = S0*KGD*FR0 
      S0 = S0-(RSS[t]+RG[t])/FR0   #s表示自由水在产流面积上的平均蓄水深 
    } else { 
      FR0 = R[t]/PE[t]              #用流量除以单位面积上的净雨可以理解为产流深即得产流面，产流面积时段末产流面积
      S0 = X*S0/FR0              #产流面积变化的影响，S0时段末的径流深，
      SS = S0
      
      Q = R[t]/FR0                  #为产流面积上的平均值 
      NN = floor(Q/5) + 1        #每次入流按5毫米分成并取整数NN为了消除前向差分误差,截尾取整
      Q = Q/NN                   # 一天分为CSng[NN]个时段 
      Kssdd = (1-(1-(KGD + KSSD))^(1/NN))/(1+KGD/KSSD) 
      Kgdd = Kssdd*KGD/KSSD 
      
      RS[t] = 0 
      RSS[t] = 0 
      RG[t] = 0 
      # SM流域的平均自由水蓄水容量 
      Smm = (1 + Ex) * SM             #Smm全流域最大单点的自由水蓄水容量 
      if (Ex < 0.001) 
      {Smmf = Smm}                     #Smmf表示*产流面积*最大一点的自由水蓄水容量 
      else {
        Smmf = Smm*(1-(1-FR0)^(1/Ex))   #Ex表示流域自有水容水容量曲线的指数 
      }  
      
      Smf = Smmf/(1+Ex)     # Smf表示*产流面积*上一点的自由水平均蓄水容量深 
      Rsd = NULL
      Rssd = NULL
      Rgd = NULL
      
      for (j in 1:NN) { 
        if (S0 >= Smf)         #s0表示自由水在产流面积上的平均蓄水深 
        {
          S0 = Smf   #Smf为流域平均自由水蓄水容量深，自由水蓄水容量-面积曲线包络下的面积，S0<=Smf;
        } 
        
        AU = Smmf*(1-(1-S0/Smf)^(1/(1 + Ex)))   #Smmf为流域单点最大的自由水蓄水容量
        
        if (Q+AU <= 0)  
        {
          Rsd[j] = 0                #当径流与此时刻的平均蓄水深之和小于0时不产流 
          Rssd[j] = 0 
          Rgd[j] = 0 
          S0 = 0 
        }
        else if  (Q+ AU>= Smmf)           #当径流与此时刻的平均蓄水深之和大于最大平均蓄水深全面产壤中流 
        {
          Rsd[j] = (Q +S0-Smf)*FR0        #Rsd中d为分段的地面流，P149，式（5-27）：FR0=FR 
          Rssd[j] = Smf * Kssdd * FR0     #Rsd中d为分段的壤中流 
          Rgd[j] = Smf * Kgdd * FR0         #Rsd中d为分段的地下径流 
          S0 = Smf - (Rssd[j] + Rgd[j]) / FR0  # s表示自由水在产流面积上的平均蓄水深 
        }else if (Q + AU < Smmf)            #当径流与此时刻的平均蓄水深之和大于最大平均蓄水深部分产壤中流 
        {
          Rsd[j] = (Q - Smf + S0 + Smf *(1-(Q + AU)/Smmf)^(1+Ex))*FR0 
          Rssd[j] = (S0 + Q - Rsd[j]/FR0)*Kssdd * FR0 
          Rgd[j] = (S0 + Q - Rsd[j]/FR0)*Kgdd * FR0 
          S0 = S0 + Q - (Rsd[j] + Rssd[j] + Rgd[j]) / FR0 
        }
        
        RS[t] = RS[t] + Rsd[j]       #累计三流 
        RSS[t] = RSS[t] + Rssd[j]    #累计 
        RG[t] = RG[t] + Rgd[j] 
      }
      
    } 
    

    RS[t] = RS[t] + Rimp
    
    
  } #循环计算结束
  out = data.frame(RS,RSS,RG)
  out = round(out,digits = 2)
  out = `colnames<-`(out,c("RS","RI","RG"))
  return(out)
}

#----超渗产流-----
Horton <- function(P,E,Kc,A,B,Wmax) {
  # ft = A+B*t(-1/2),t位计算时段
  # A      #菲利普下渗公式参数
  # B       #菲利普下渗公式参数
  # 
  # Wmax   #土壤蓄水容量，mm
  E = E*Kc
  W0 = Wmax/4    #代表上一时刻土壤含水量，mm，此处赋予初始值
  W = NULL    #单位：mm,土壤含水量，中间计算值
  Rs = NULL   #单位：mm
  ft = NULL    #菲利普下渗公式，时段下渗能力
  
  for (t in 1:length(P)) {
    
    ft[t] = B^2*(1+sqrt(1+A*W0/B^2))/W0+A
    ft[t] = round(ft[t],digits = 1)
    if (P[t] - E[t] <= 0) {
      Rs[t] = 0
      W[t] = W0+P[t]-E[t]
    }
    if (P[t] - E[t] <= ft[t] & P[t] - E[t]>0) {
      Rs[t] = 0
      W[t] = W0+P[t]-E[t]
    }
    if (P[t] - E[t] > ft[t] ) {
      Rs[t] = P[t] - E[t] - ft[t]
      W[t] = W0+ft[t]
    }
    if (W[t] < 0) {
      W[t] = 0
      Rs[t] = 0
    }
    if (W[t] > Wmax) {
      W[t] = Wmax
      Rs[t] = Rs[t] + W[t]-Wmax
    }
    W0 = W[t]
  }
  dat = data.frame(P,W[1:length(P)],ft,Rs)
  dat = round(dat,digits = 1)
  dat = `colnames<-`(dat,c("P","W","f","RS"))
  return(dat)
}

#---垂向混合产流------
BWM <- function(P,E,Kc,B,WM,KF,fc,BF,KI,KG){
  # B  ： 包气带蓄水容量分布曲线指数，人机交互率定
  # WM ： 流域平均蓄水容量,mm
  # KF ： 渗透系数，人机交互率定
  # DT：  #计算时段
  # fc ： 稳定下渗率，人机交互率定
  # BF ： 下渗分布曲线指数，人机交互率定
  # f_mean ：流域平均下渗率,mm/delta_t，人机交互率定
  # KI ： 壤中流出流系数
  # KG ： 地下径流出流系数
  E = E*Kc
  DT = 1
  #要选择降雨强度大、降雨历时短的一次降雨过程，可按照实测径流深，用平割法得到平均下渗率
  W = NULL
  #f_max = f_mean*(1+BF)  #流域内最大单点下渗能力
  WMM = WM*(1+B)   #单点最大包气带蓄水容量值,mm
  W[1] = WM/2  #初始时刻土壤含水量，mm
  #----以下定义计算的中间变量-----
  FA = NULL       #流域时段下渗水量
  f = NULL       #流域平均下渗能力;f=fc*(1+KF*(WM-W)/WM) #W为土壤含水量，WM流域平均蓄水容量
  R = NULL
  RS = NULL   #地表径流
  RR = NULL   #下渗径流总和
  RI = NULL   #壤中流，mm
  RG = NULL   #地下径流，mm
  S0 = 5     #初始敞开式自由水箱自由水量，用于划分地下径流水源
  S = NULL  #敞开式自由水箱自由水量
  
  a = NULL
  for (t in 1:length(P)) {
    f[t] = fc*(1+KF*(WM-W[t])/WM) #WW[t]下的流域平均下渗能力
    f_mean = f[t]  
    f_max = f_mean*(1+BF) #流域内最大单点下渗能力
    if ((P[t]-E[t]) > 0) {
      
      if ( (P[t]-E[t]) >= (f_max*DT) ) {
        FA[t] = f_mean*DT
      } else {
        FA[t] = DT*(f_mean-f_mean*(1-(P[t]-E[t])/f_max)^(BF+1))
      }
      if (FA[t] > (P[t]-E[t])) {
        FA[t] = P[t]-E[t]  #FA <= PE
      }
      RS[t] = P[t]-E[t]-FA[t]
      a[t] = WMM*(1-(1-W[t]/WM)^(1/(1+B)))
      if (a[t] < 0) {a[t]=0}
      if (a[t] > WMM) {
        a[t] = WMM
      }
      if ( (FA[t]+a[t]) >= WMM ) {
        RR[t] = FA[t]+W[t]-WM
      } else {
        RR[t] = FA[t]+W[t]-WM+WM*(1-(FA[t]+a[t])/WMM)^(B+1)
      }
      if(RR[t] > FA[t]){
        RR[t] = FA[t]
      }
      R[t] = RS[t]+RR[t] 
      #dW = FA[t] - RR[t] = P[t]-E[t]-RS[t]-RR[t] = P[t]-E[t]-R[t]
      W[t+1] = W[t]+P[t]-E[t]-R[t]
      if (W[t+1] > WM) {
        RS[t] = RS[t]+W[t+1] - WM
        W[t+1] = WM
      }
      #敞开式水箱划分地下径流水源---
      if (t == 1) {S_p = S0} 
      else {S_p = S[t-1]}
      S[t] = S_p+RR[t]
      RI[t] = KI*S[t]
      RG[t] = KG*S[t]
      S[t] = S[t] - RI[t] - RG[t]
      
    } else {
      
      FA[t] = 0 
      RS[t] = 0
      RR[t] = 0
      RI[t] = 0
      R[t] = RR[t]+RS[t]
      a[t] = WMM*(1-(1-W[t]/WM)^(1/(1+B)))
      W[t+1] = W[t] + P[t]-E[t]
      if (W[t+1] < 0) {W[t+1] = 0}
      
      if (t == 1) {S_p = S0} 
      else {S_p = S[t-1]}
      S[t] = S_p
      RI[t] = 0
      RG[t] = 0
    } 
  }
  
  out = data.frame(P,E,R,RR,RS,RI,RG)
  out = round(out,digits = 2)
  out = `colnames<-`(out,c("P","E","R","RR","RS","RI","RG"))
  return(out)
}

#-----三水源以线性水库进行坡面汇流------
LinearReservoir <- function(DAREA,DT,RS,RI,RG,Cs,Ci,Cg,ini){
  # Cs:地表水线性水库汇流（消退）系数 
  # Ci:壤中流线性水库汇流系数 
  # Cg:地下水线性水库汇流系数 
  # ini:用于调节计算的初始流量
  U=DAREA/(DT*3.6)  #单位转换 ,DAREA为流域面积，DT为计算时段单位，24小时或者1小时，delta t
  Qrs = NULL #地表径流出流量
  Qri = NULL
  Qrg = NULL
  Qrt = NULL  #总流量，m3/s
  Qrs0 = mean(RS)*U*ini
  Qri0 = mean(RI)*U*ini
  Qrg0 = mean(RG)*U*ini
  for (t in 1:length(RS)) {
    Qrs[t] = RS[t]*U*(1-Cs) + Qrs0 * Cs   #地面水线性水库汇流系数CS 
    Qri[t] = RI[t]*U*(1-Ci) + Qri0 * Ci         #壤中流线性水库汇流系数CI 
    Qrg[t] = RG[t]*U*(1-Cg) + Qrg0 * Cg            #地下水线性水库汇流系数Cg 
    
    Qrt[t] = Qrs[t] + Qri[t] + Qrg[t]   #total
    
    Qrs0 = Qrs[t] 
    Qri0 = Qri[t] 
    Qrg0 = Qrg[t] 
    #Rs0 = RS[t]
  }
  out = data.frame(Qrs,Qri,Qrg,Qrt)
  out = `colnames<-`(out,c("Qrs","Qri","Qrg","Qrt"))
  return(out)
}
#---------时段单位线汇流-----------
UnitHrydo = function(DAREA,DT,RS,RI,RG,unithydro,Ci,Cg,ini){
  q = unithydro[,2] #单位线
  n = length(q) #单位线长度
  m = length(RS) #净雨时段数
  Qrs = NULL #单位：m3/s
  for (t in 1:length(RS)){
    if (t < n) {
      k1 = 1
    } else if (t>=n) {
      k1 = t-n+1
    }
    if (t < m) {
      k2 = t
    } else {
      k2 = m
    }
    Qrs[t] = 0
    for (j in k1:k2) {
      Qrs[t] = Qrs[t] + RS[j]/10*q[t-j+1] #RS[j]/10,:单位净雨的倍数，此处单位净雨为10mm。
    }
   
  }
  U=DAREA/(DT*3.6)  #单位转换 ,DAREA为流域面积，DT为计算时段单位，24小时或者1小时，delta t
  Qri0 = RI[1]*U*ini #单位：m3/s,ini < 1
  Qrg0 = RG[1]*U*ini
  Qri = NULL  #单位：m3/s
  Qrg = NULL  #单位：m3/s

  for (t in 1:length(RS)) {
    Qri[t] = RI[t]*U*(1-Ci) + Qri0 * Ci         #壤中流线性水库汇流系数CI 
    Qrg[t] = RG[t]*U*(1-Cg) + Qrg0 * Cg            #地下水线性水库汇流系数Cg 
    Qri0 = Qri[t] 
    Qrg0 = Qrg[t] 
  }
  Qrt = Qrs+Qri+Qrg
  out = data.frame(Qrs,Qri,Qrg,Qrt)
  return(out)
}
#-----------马斯京根洪水演算-----
MasKing <- function(I,Ke,Xe,DT){
  C0 = (0.5*DT-Ke*Xe)/(0.5*DT+Ke-Ke*Xe)
  C1 = (0.5*DT+Ke*Xe)/(0.5*DT+Ke-Ke*Xe)
  C2 = 1-C0-C1
  Qout = NULL
  Qout[1] = I[1]
  for (i in 2:length(I)) {
    Qout[i] = C0*I[i]+C1*I[i-1]+C2*Qout[i-1]
  }
  out = data.frame(I,Qout)
  out = `colnames<-`(out,c("I","O"))
  return(out)
}
#------------滞后演算---------------
LagRounting <- function(I,CR,L) {
  Qout = NULL
  if (L != 0){
    for (i in 1:L) {
      Qout[i] = I[i]
    }
    
    for (i in (L+1):length(I)) {
      Qout[i] = CR*Qout[i-1]+(1-CR)*I[i-L]
    }
  } else {
    Qout[1] = I[1]
    for (i in 2:length(I)){
      Qout[i] = CR*Qout[i-1]+(1-CR)*I[i-L]
    }
  }
  out = data.frame(I,Qout)
  out = `colnames<-`(out,c("I","O"))
  return(out)
}


ui <- tagList(
  shinythemes::themeSelector(),
  navbarPage("HyDro-Model",
                 
                 # ---------------Introduction------------------
                 tabPanel("HyDro-Model",
                          fluidRow(column(12,
                                          wellPanel(
                                            h3("HyDro-Model软件前言",align = "center")       
                                          ))),
                          fluidRow(column(12,
                                          wellPanel(
                                            p("水文模型是指用模拟方法将复杂的水文现象和过程经概化所给出的近似的科学模型，
其在进行水文规律研究和解决生产实际问题中起着重要的作用，
随着现代科学技术的飞速发展，以计算机和通信为核心的信息技术在水文水资源及水利工程科学领域的广泛应用，
使水文模型的研究得到迅速发展，并广泛应用于水文基本规律研究、水资源评价与开发利用、气候变化及人类活动对水资源和水环境影响分析等领域。
因此，水文模型的开发研究具有重要的科学意义和应用价值。
概念性流域水文模型是以水文现象的物理概念和一些经验公式为基础构造的水文模型，
它将流域的物理基础(如下垫面等)进行概化，再结合水文经验公式(如下渗曲线、汇流单位线、蒸散发公式等)来近似地模拟流域水流过程。
HyDro-Model系统软件考虑了不同产流机制（蓄满产流、超渗产流和垂向混合产流）、
水源划分以及洪水演算方法，概化流域水文过程，模型理论以及算法实现主要参考包为民主编《水文预报（第四版）》。
使用过程中，如有问题请与软件开发者邮件联系，作者保留所有权利。
                                              "),
                                            br(),
                                            p("Programmer:Shane Guan",align = "right"),
                                            p("E-mail:xxguan@hhu.edu.cn",align = "right")  
                                            
                                          ))),
                          p(em("Address:College of hydrology and water resource, Hohai University,1,Xikang Road,NanJing210098"),align = "center")
                 ),
             #-----------预报方案------------------
                 tabPanel("预报方案",
                          sidebarLayout(
                            sidebarPanel(
                              h4("基本信息"),
                              textInput("Basinname","流域名:",
                                        value = ""),
                              numericInput("DAREA",
                                           "流域面积(km2):",
                                           step = 50,
                                           value = 4766),
                              numericInput("Imp",
                                           "流域不透水面积占比:",
                                           min = 0,max = 1,
                                           value = 0.02,step = 0.001)
                            ),
                     
                            mainPanel(
                              #h3("Comparison of simulated and recorded monthly discharge"),

                              fluidRow(
                                column(6,
                                       selectInput('chanliu', '产流方案', c("蓄满产流","超渗产流","垂向混合产流")
                                                   , selectize=TRUE)
                                ),
                                column(6,
                                       selectInput('RShuiliu', '地表坡面汇流', c("线性水库","时段单位线")
                                                     , selectize=TRUE),
                                       helpText("壤中流和地下径流默认使用线性水库汇流")
                                )
                                ),
                              hr(),
                              fluidRow(
                                column(6,
                                       selectInput('shuiyuan', '水源划分', c("三水源","二水源")
                                                   , selectize=TRUE),
                                       hr()
                                ),
                                column(6,
                                       selectInput('hewanghuiliu', '河网汇流', c("滞后演算法","马斯京根法")
                                                   , selectize=TRUE),
                                       hr()
                                )
                              )
                              )
                            )
                          
                      ),
             #--------------数据准备-----------------------
                 tabPanel("数据准备",
                          sidebarLayout(
                            sidebarPanel(
                              fileInput("file", "选择数据文件(.csv)",
                                        multiple = TRUE,
                                        accept = c("text/csv",
                                                   "text/comma-separated-values,text/plain",
                                                   ".csv")),
                              helpText("请选择本地逗号分隔符数据文件(.csv);
                                       第一行留设为字段标题：年，月，日，时，降水(mm)，蒸发(mm)，观测流量(m3/s)"),
                              numericInput("DT","计算时段(尺度)/h",value = 24,min = 1,max = 24,step =1)
                            ),
                            mainPanel(
                              tabsetPanel(
                              tabPanel("数据列表",
                                       DT::dataTableOutput("datafile")
                              ),
                              tabPanel("数据图示", 
                                       
                                       plotOutput("precipitationplot",height = 300),
                                       hr(),
                                       plotOutput("dischargeplot",height = 300)
                                       ),
                              tabPanel("数据统计", 
                                       h4("数据时限内月加总平均值统计表："),
                                       tableOutput("datastatistable")
                                       )
                              )
                            )
                          )
                        ),
                          
                          
            #----------模型率定------------------
                 tabPanel("模型率定",
                          fluidRow(
                            column(7,h3("模拟结果"),
                                   wellPanel(
                                     sliderInput("period", "Select baseline period", min = 1940, 
                                                 max = 2020, value = c(1955, 1970),step = 1),
                                     h4(em("模拟精度表:")),
                                     tableOutput("Accuracy"),
                                     hr(),
                                     h4(em("模拟与实测流量过程对比图:")),
                                     plotOutput("Simulatedplot",height = 300,
                                                brush = brushOpts(
                                                id = "Observe_brush",
                                                resetOnNew = TRUE)
                                     ),
                                     tabsetPanel(
                                       tabPanel("三水源过程",
                                                h4(em("分水源流量过程，可用于次洪率定:")),
                                                plotOutput("QTSIGplot",height = 300)      
                                       ),
                                       tabPanel("Precipitation",
                                                h4(em("逐日降水过程:")),
                                                plotOutput("Precipitation_plot",height = 300)      
                                       ),
                                       tabPanel("Monthly:Qr&Qs", 
                                                h4(em("逐月月均流量过程:")),
                                                plotOutput("Q_monthly_plot",height = 300)  
                                       ),
                                       tabPanel("For Observe", 
                                                h4(em("Brush区域显示:")),
                                                plotOutput("Observe_plot",height = 300)  
                                       )
                                     ),
                                     fluidRow(
                                       column(6,
                                         wellPanel(
                                           selectInput("dataset_download", 
                                                       "选择下载数据集:",
                                                       choices = c("Runoff generation", 
                                                                   "Simulated Process",
                                                                   "Simulated monthly Process"))
                                           
                                         )
                                       ),
                                       column(6,
                                         wellPanel(
                                           downloadButton("downloadData","点击执行下载")
                                         )
                                       )
                                     )
                                   )),
                            column(5,
                                   h3("参数模块"),
                                   wellPanel( 
                                     fluidRow("产流参数",
                                       column(12,
                                         wellPanel(
                                           conditionalPanel("input.chanliu == '蓄满产流'",
                                                            fluidRow(
                                                                     column(4,
                                                                            wellPanel(
                                                                              numericInput("SM","自由水平均蓄水容量SM:",min = 0,value = 12.3,step = 0.1),
                                                                              numericInput("Ki","自由水箱壤中流出流系数Ki:",min = 0,value = 0.399,max = 1,step = 0.001),
                                                                              numericInput("Kg","自由水箱地下水出流系数Kg:",min = 0,value = 0.495,max = 1,step = 0.001),
                                                                              helpText("Attention:Ki+Kg<=1")
                                                                            )
                                                                     ),
                                                                     column(4,
                                                                            wellPanel(
                                                                              numericInput("B1","蓄水容量分布曲线指数B:",min = 0,value = 1.2,step = 0.01),
                                                                              numericInput("C","深层蒸发系数C:",min = 0,value = 0.12,step = 0.01),
                                                                              numericInput("Ex","自由水分布曲线指数Ex:",min = 0,value = 1.58,step = 0.1),
                                                                              numericInput("FR0","初始产流面积比例:",min = 0,value = 0.002,max = 1,step = 0.001)
                                                                            )
                                                                     ),
                                                                     column(4,
                                                                            wellPanel(
                                                                              numericInput("Kc","蒸散发折算系数Kc:",min = 0,value = 0.26,step = 0.01),
                                                                              numericInput("WM1","平均蓄水容量WM/mm:",min = 0,value = 400,step = 2),
                                                                              numericInput("WUM","上层蓄水容量WUM/mm:",min = 0,value = 70,step = 2),
                                                                              numericInput("WLM","下层蓄水容量WLM/mm:",min = 0,value = 30,step = 2)
                                                                              
                                                                            )
                                                                     )
                                                            ) 
                                                         
                                                            
                                           ),
                                           conditionalPanel("input.chanliu == '超渗产流'",
                                                            fluidRow(
                                                              column(6,
                                                                wellPanel(
                                                                  numericInput("Kc2","蒸散发折算系数Kc:",min = 0,value = 0.8,step = 0.001),
                                                                  numericInput("Wmax","土壤蓄水容量/mm:",min = 0,value = 200,step = 2)
                                                                )
                                                              ),
                                                              column(6,
                                                                wellPanel(
                                                                  numericInput("A","菲利普下渗公式参数:A",value = 0.4,step = 0.01),
                                                                  numericInput("BB","菲利普下渗公式参数:B",value = 0.3,step = 0.01)
                                                                )
                                                              )
                                                            )
                                                            
                                                            
                                                            
                                           ),
                                           conditionalPanel("input.chanliu == '垂向混合产流'",
                                                            fluidRow(
                                                                     column(6,
                                                                            wellPanel(
                                                                              numericInput("Kc3","蒸散发折算系数Kc:",min = 0,value = 0.8,step = 0.001),
                                                                              numericInput("KF","渗透系数KF:",value = 0.4,step = 0.01),
                                                                              numericInput("fc","稳定下渗率fc:",value = 0.4,min = 0,max = 2,step = 0.01),
                                                                              numericInput("BF","下渗分布曲线指数BF:",min = 0,value = 1.2,step = 0.01)
                                                                              #numericInput("f_mean","流域平均下渗率:",value = 0.3,step = 0.01)
                                                                              
                                                                            )
                                                                     ),
                                                                     column(6,
                                                                            wellPanel(
                                                                              numericInput("WM2","平均蓄水容量WM/mm:",min = 0,value = 180,step = 2),
                                                                              numericInput("B2","蓄水容量分布曲线指数B:",min = 0,value = 1.2,step = 0.01),
                                                                              numericInput("KI","自由水箱壤中流出流系数Ki:",min = 0,value = 0.4,max = 1,step = 0.001),
                                                                              numericInput("KG","自由水箱地下水出流系数Kg:",min = 0,value = 0.3,max = 1,step = 0.001),
                                                                              helpText("Attention:Ki+Kg<=1")
                                                                            )
                                                                     )
                                                            )
                                                            
                                           )
                                         )
                                       )
                                     ),
                                     
                                  fluidRow("坡面汇流参数",
                                    column(12,
                                       wellPanel(
                                             conditionalPanel("input.RShuiliu == '线性水库'",
                                                              fluidRow(
                                                                column(6,
                                                                       wellPanel(
                                                                         numericInput("Cs","地表水线性水库汇流系数Cs:",min = 0,value = 0.25,max = 1,step = 0.01),
                                                                         numericInput("Ci","壤中流线性水库汇流系数Ci:",min = 0,value = 0.964,max = 1,step = 0.001)
                                                                         
                                                                       )
                                                                       ),
                                                                column(6,
                                                                       wellPanel(
                                                                         numericInput("Cg","地下水线性水库汇流系数Cg:",min = 0,value = 0.973,max = 1,step = 0.001),
                                                                         numericInput("ini","初始流量调节系数:",min = 0,value = 0.3,max = 1,step = 0.01)
                                                                       )
                                                                       )
                                                              )
                                                       ),
                                             conditionalPanel("input.RShuiliu == '时段单位线'",
                                                              fileInput("UnitHydroGraph", "选择10mm单位线文件(.csv)",
                                                                        multiple = TRUE,
                                                                        accept = c("text/csv",
                                                                                   "text/comma-separated-values,text/plain",
                                                                                   ".csv")),
                                                              helpText(em("注意时段单位线UH的时段转换！")),
                                                              numericInput("ini2","初始流量调节系数:",min = 0,value = 0.3,max = 1,step = 0.01)
                                                              )
                                          )
                                     )
                                  ),
                                  fluidRow("河网汇流参数",
                                    column(12,
                                           wellPanel(
                                             conditionalPanel("input.hewanghuiliu == '马斯京根法'",
                                                              fluidRow(
                                                                column(6,
                                                                  wellPanel(
                                                                    numericInput("Ke","河段传播时间Ke:",value= 0.5,min = 0,step = 0.001)
                                                                  )
                                                                ),
                                                                column(6,
                                                                       wellPanel(
                                                                         numericInput("Xe","流量比重系数Xe:",value = 0.45,min = 0,max=1,step = 0.001)
                                                                       )
                                                                 )
                                                               )
                                                              ),
                                             conditionalPanel("input.hewanghuiliu == '滞后演算法'",
                                                              fluidRow(
                                                                column(6,
                                                                       wellPanel(
                                                                         numericInput("CR","河网蓄水消退系数CR:",value= 0.496,min = 0,max = 1,step =0.001)
                                                                       )
                                                                ),
                                                                column(6,
                                                                       wellPanel(
                                                                         numericInput("L","滞后时段L:",value = 0,min = 0,step = 1) 
                                                                       )
                                                                       )
                                                              )                         
                                                 )
                                           )
                                      )
                                  )
                            )
                          )
                     )      
                 )
            #--------待添加功能--------
            
            # tabPanel("NULL"
            #   
            # )
                 
   )
)

#-------------------server-----------------
server <- function(input, output) {
  #-----------data input and display--------------
  UnitHydroGraph <- reactive({
    req(input$UnitHydroGraph)
    df <- read.csv(input$UnitHydroGraph$datapath,
                    header = T,
                    sep = ",",
                    quote = '"')
    `colnames<-`(df,c("ID","q"))
  })
  
  filedata <- reactive({
    req(input$file)
    df1 <- read.csv(input$file$datapath,
                    header = T,
                    sep = ",",
                    quote = '"')
    

    `colnames<-`(df1,c("year","month","day","hour","P","E","Q"))
  })
  output$datafile <- DT::renderDataTable({   #upload data
    df = filedata()
    
    DT::datatable(df)
  })
  output$precipitationplot <- renderPlot({
    dat = filedata()
    dat$x = 1:length(dat[,1])
    ggplot(dat,aes(x,y = P,fill = P)) + geom_bar(stat = "identity")+
      xlab("Time")+ylab("Precipitation/mm")+theme(legend.position = "top")+
      theme(legend.title = element_blank())#+theme(guide_legend(direction = "horizontal"))
  })
  output$dischargeplot <- renderPlot({
    dat = filedata()
    dat$x = 1:length(dat[,1])
    ggplot(dat,aes(x,y = Q))+geom_line()+geom_area(alpha = 1/4)+
      xlab("Time")+
      ylab("Discharge/(m3/s)")+
        theme(legend.title = element_blank())
  })
  output$datastatistable <- renderTable({
    dat = filedata()
    monthly <- function(year,mon,value){
      out = NULL
      for (i in 1:12){
        out[i] = sum(value[mon == i],na.rm = T)/length(names(table(year)))
      }
     return(out) 
    }
    P_mon = monthly(dat[,1],dat[,2],dat$P)
    Q_mon = monthly(dat[,1],dat[,2],dat$Q)
    mon = 1:12
    out = data.frame(mon,P_mon,Q_mon)
    `colnames<-`(out,c("月份","雨量/mm","流量/(m3/s)"))
  })
  #---------------simulation-------------------
  #------------runoffgeneration-----------
  runoffgeneration <- reactive({
    dat = filedata()
    a = input$period[1]
    b = input$period[2]
    dat = dat[(dat[,1] >= a & dat[,1]<=b),] 
    if (input$chanliu == "蓄满产流") {
     out = Dunne(dat$P,dat$E,input$DT,input$Imp,input$Kc,
            input$WM1,input$WUM,input$WLM,input$B1,input$C,
            input$Ex,input$SM,input$Ki,input$Kg,input$FR0)
    } else if (input$chanliu == "超渗产流") {
      out = Horton(dat$P,dat$E,input$Kc2,input$A,input$BB,input$Wmax)
    } else if (input$chanliu == "垂向混合产流") {
      out = BWM(dat$P,dat$E,input$Kc3,input$B2,input$WM2,input$KF,
                input$fc,input$BF,input$KI,input$KG)
    }
    out
  })
  FlowCon <- reactive({
    dat = runoffgeneration()
    if (input$chanliu == "超渗产流") {
      RS = dat$RS
      RI = rep(0,length(RS))
      RG = rep(0,length(RS))
    } else if (input$chanliu == "蓄满产流" | input$chanliu == "垂向混合产流"){
      RS = dat$RS
      RI = dat$RI
      RG = dat$RG
    }
    if (input$RShuiliu == "线性水库") {
      out = LinearReservoir(input$DAREA,input$DT,RS,RI,RG,input$Cs,input$Ci,input$Cg,input$ini)
    } else if (input$RShuiliu == "时段单位线") {
      unitHydrograph = UnitHydroGraph() #单位线文件数据框
      out = UnitHrydo(input$DAREA,input$DT,RS,RI,RG,unitHydrograph,input$Ci,input$Cg,input$ini2)
      #DAREA,DT,RS,RI,RG,unithydro,Ci,Cg,ini
    }
    out
  })
  #------------Rounting------------------
  Rounting <- reactive({
    Observed = filedata() #"year","month","day","hour","P","E","Q"

    a = input$period[1]
    b = input$period[2]
    Observed = Observed[(Observed[,1] >= a & Observed[,1]<=b),] 
    dat = FlowCon()  #Qrs,Qri,Qrg,Qrt
    if (input$hewanghuiliu == "马斯京根法") {
      Qsimulated = MasKing(dat$Qrt,input$Ke,input$Xe,input$DT)[,2]
    } else {
      Qsimulated = LagRounting(dat$Qrt,input$CR,input$L)[,2]
    }
    out = data.frame(Qsimulated,Observed$Q)
    `colnames<-`(out,c("Qsimulated","Qrecorded"))
  })
  #------------------Accuracy Table------------
  output$Accuracy <- renderTable({
    dat = Rounting()      #"Qsimulated","Qrecorded"
    rawdata = filedata()  #"year","month","day","hour","P","E","Q"
    a = input$period[1]
    b = input$period[2]
    rawdata = rawdata[(rawdata[,1] >= a & rawdata[,1]<=b),] 
    #**************日尺度精度统计*******************
    R2 = round(1-sum((dat$Qrecorded-dat$Qsimulated)^2,na.rm = T)/sum((dat$Qrecorded-mean(dat$Qrecorded))^2,na.rm = T),digits = 3)
    Re_runoff = round(sum(dat$Qsimulated-dat$Qrecorded,na.rm = T)/sum(dat$Qrecorded,na.rm = T)*100,digits = 3) #径流深相对误差（%）
    Qm_si = max(dat$Qsimulated)
    Qm_re = max(dat$Qrecorded)
    Re_discharge = round((Qm_si-Qm_re)/Qm_re*100,digits = 3) #洪峰流量相对误差（%）
    gap = which.max(dat$Qsimulated)-which.max(dat$Qrecorded) #峰现时差
    #***************月尺度精度统计*******************
    Yuedata = data.frame(rawdata$year,rawdata$month,rawdata$day,dat$Qsimulated,dat$Qrecorded)
    colnames(Yuedata) = c("year","month","day","Qsimulated","Qrecorded") #流量的单位是:m3/s
    Yuedata[,4][Yuedata[,4] < 0] = NA #流量小于0则为无效值
    Yuedata[,5][Yuedata[,5] < 0] = NA
    Y_Qi = Yuedata[1,1]
    Y_Mo = Yuedata[length(Yuedata[,1]),1]
    Q_re_mon = NULL
    Q_si_mon = NULL
    Year = Y_Qi:Y_Mo
    for (i in 1:length(Y_Qi:Y_Mo)) {
      for (j in 1:12) {
        #计算月均流量，m3/s
        Q_re_mon[(i-1)*12+j] = mean(Yuedata$Qrecorded[Yuedata$year == Year[i] & Yuedata$month == j],na.rm = T)
        Q_si_mon[(i-1)*12+j] = mean(Yuedata$Qsimulated[Yuedata$year == Year[i] & Yuedata$month == j],na.rm = T)
        
      }
    }
    R2_mon = round(1-sum((Q_re_mon - Q_si_mon)^2,na.rm = T)/sum((Q_re_mon - mean(Q_re_mon))^2,na.rm = T),digits = 3)
    Re_runoff_mon = round(sum(Q_si_mon - Q_re_mon,na.rm = T)/sum(Q_re_mon,na.rm = T)*100,digits = 3) #径流深相对误差（%）
    Qm_si_mon = max(Q_si_mon,na.rm = T)
    Qm_re_mon = max(Q_re_mon,na.rm = T)
    Re_discharge_mon = round((Qm_si_mon-Qm_re_mon)/Qm_re_mon*100,digits = 3) #洪峰流量相对误差（%）
    #****************************
    out = data.frame(c(R2,R2_mon),c(Re_runoff,Re_runoff_mon),c(Re_discharge,Re_discharge_mon),c(gap,NA))
    
    `colnames<-`(out,c("确定性系数","径流深相对误差/%","洪峰流量相对误差/%","峰现时差"))
    
  })
  #----------------Simulatedplot---------------
  ranges <- reactiveValues(x = NULL, y = NULL)
  
  output$Simulatedplot <- renderPlot({
    dat = FlowCon()       #Qrs,Qri,Qrg,Qrt
    Observed = filedata() #"year","month","day","hour","P","E","Q"
    a = input$period[1]
    b = input$period[2]
    Observed = Observed[(Observed[,1] >= a & Observed[,1]<=b),] 
    
    if (input$hewanghuiliu == "马斯京根法") {
      Qsimulated = MasKing(dat$Qrt,input$Ke,input$Xe,input$DT)[,2]
    } else if (input$hewanghuiliu == "滞后演算法") {
      Qsimulated = LagRounting(dat$Qrt,input$CR,input$L)[,2]
    }
    Q = c(Qsimulated,Observed$Q)
    x = rep(c(1:length(Qsimulated)),2)

    label = c(rep("Simulate",length(Qsimulated)),rep("Record",length(Observed$Q)))
    
    graph = data.frame(x,Q,label)
    ggplot(graph,aes(x = x,y = Q,linetype = label,colour = label))+geom_line(lwd = 0.8)+
      xlab("Time")+ylab("Discharge/(m3/s)") + theme(legend.position = "bottom")+
      guides(linetype = guide_legend(title = NULL),colour = guide_legend(title = NULL))+
      scale_color_brewer(palette = "Set1") 

  })
  output$Precipitation_plot <- renderPlot({
    dat = filedata() #"year","month","day","hour","P","E","Q"
    a = input$period[1]
    b = input$period[2]
    dat = dat[(dat[,1] >= a & dat[,1]<=b),] 
    
    x = 1:length(dat$P)
    graph = data.frame(x,y = dat$P)
    ggplot(data = graph,aes(x =x,y=y)) + geom_bar(stat="identity")+
    guides(fill = guide_legend(title = NULL))+
    xlab("Time") + ylab("Rainfall/mm") +
    coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)
  })
  
  graph_month_data <- reactive({
    #******月尺度数据********
    dat = Rounting()      #"Qsimulated","Qrecorded"
    rawdata = filedata()  #"year","month","day","hour","P","E","Q"
    a = input$period[1]
    b = input$period[2]
    rawdata = rawdata[(rawdata[,1] >= a & rawdata[,1]<=b),] 
    
    Yuedata = data.frame(rawdata$year,rawdata$month,rawdata$day,dat$Qsimulated,dat$Qrecorded)
    colnames(Yuedata) = c("year","month","day","Qsimulated","Qrecorded") #流量的单位是:m3/s
    Yuedata[,4][Yuedata[,4]<0] = NA #流量小于0则为无效值
    Yuedata[,5][Yuedata[,5]<0] = NA
    Y_Qi = Yuedata[1,1]
    Y_Mo = Yuedata[length(Yuedata[,1]),1]
    Q_re_mon = NULL
    Q_si_mon = NULL
    Year = Y_Qi:Y_Mo
    for (i in 1:length(Y_Qi:Y_Mo)) {
      for (j in 1:12) {
        #计算月均流量，m3/s
        Q_re_mon[(i-1)*12+j] = mean(Yuedata$Qrecorded[Yuedata$year == Year[i] & Yuedata$month == j],na.rm = T)
        Q_si_mon[(i-1)*12+j] = mean(Yuedata$Qsimulated[Yuedata$year == Year[i] & Yuedata$month == j],na.rm = T)
        
      }
    }
    YY = NULL
    for (i in 1:length(Year)) {
      YY = c(YY,rep(Year[i],12))
    }
    MM = rep(1:12,length(Year))
    Simulation_Monthly = data.frame(YY,MM,Q_re_mon,Q_si_mon)
    x = rep(1:length(Q_re_mon),2)
    #tick_label = rep(paste0(Simulation_Monthly$YY,"-",Simulation_Monthly$MM),2)
    Value = c(Q_si_mon,Q_re_mon)
    label = c(rep("Simulate",length(Q_si_mon)),rep("Record",length(Q_re_mon)))
    data.frame(x,Value,label)
  })
  output$Q_monthly_plot <- renderPlot({
    
    graph = graph_month_data()
    ggplot(data = graph,aes(x=x,y=Value,color = label))+geom_line()+
      xlab("Time")+ylab("Mongthly-Average Diacharge/(m3/s)")+ theme(legend.position = "bottom")+
      guides(linetype = guide_legend(title = NULL),colour = guide_legend(title = NULL))+
      scale_color_brewer(palette = "Set1")
  })
  
  observe({
    brush <- input$Observe_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })

  output$Observe_plot <- renderPlot({
    dat = FlowCon()       #Qrs,Qri,Qrg,Qrt
    Observed = filedata() #"year","month","day","hour","P","E","Q"
    a = input$period[1]
    b = input$period[2]
    Observed = Observed[(Observed[,1] >= a & Observed[,1]<=b),] 
    
    if (input$hewanghuiliu == "马斯京根法") {
      Qsimulated = MasKing(dat$Qrt,input$Ke,input$Xe,input$DT)[,2]
    } else if (input$hewanghuiliu == "滞后演算法") {
      Qsimulated = LagRounting(dat$Qrt,input$CR,input$L)[,2]
    }
    Q = c(Qsimulated,Observed$Q)

    x = rep(c(1:length(Qsimulated)),2)
    label = c(rep("Simulate",length(Qsimulated)),rep("Record",length(Observed$Q)))
    
    graph = data.frame(x,Q,label)
    ggplot(graph,aes(x = x,y = Q,linetype = label,colour = label))+geom_line(lwd = 0.8)+
      xlab("Time")+ylab("Discharge/(m3/s)") + theme(legend.position = "bottom")+
      guides(linetype = guide_legend(title = NULL),colour = guide_legend(title = NULL))+
      scale_color_brewer(palette = "Set1") +
      coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)
  })
  #-----------------QTSIGplot--------------------------
  output$QTSIGplot <- renderPlot({
    Observed = filedata()#"year","month","day","hour","P","E","Q"
    a = input$period[1]
    b = input$period[2]
    Observed = Observed[(Observed[,1] >= a & Observed[,1]<=b),] 
    
    dat = FlowCon() #Qrs,Qri,Qrg,Qrt
    if (input$hewanghuiliu == "马斯京根法"){
      Qrs = MasKing(dat$Qrs,input$Ke,input$Xe,input$DT)[,2]
      Qri = MasKing(dat$Qri,input$Ke,input$Xe,input$DT)[,2]
      Qrg = MasKing(dat$Qrg,input$Ke,input$Xe,input$DT)[,2]
      Qrt = MasKing(dat$Qrt,input$Ke,input$Xe,input$DT)[,2]
    } else {
      Qrs = LagRounting(dat$Qrs,input$CR,input$L)[,2]
      Qri = LagRounting(dat$Qri,input$CR,input$L)[,2]
      Qrg = LagRounting(dat$Qrg,input$CR,input$L)[,2]
      Qrt = LagRounting(dat$Qrt,input$CR,input$L)[,2]
    }
    Q = c(Observed$Q,Qrs,Qri,Qrg,Qrt)
    x = rep(c(1:length(Qrs)),5)
    # tick_label = paste0(Observed$year,"-",Observed$month,"-",Observed$day)
    # tick_label = as.Date(tick_label)
    
    Label = c(rep("Record",length(Observed$Q)),rep("Qs",length(Qrs)),rep("Qi",length(Qri)),rep("Qg",length(Qrg)),rep("Qt",length(Qrt)))
    graph = data.frame(x,Q,Label)
    ggplot(graph,aes(x = x,y = Q,colour = Label,linetype = Label))+geom_line(lwd = 0.8)+
      xlab("Time")+ylab("Discharge/(m3/s)") + theme(legend.position = "bottom")+
      guides(linetype = guide_legend(title = NULL),colour = guide_legend(title = NULL))+
      scale_color_brewer(palette = "Dark2") +
      coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)
  })
  #-------------data for download-------------
  datasetOutput <- reactive({
    Chanliu_Result = runoffgeneration()
    Simulation_Result = Rounting()
    
    rawdata = filedata()
    a = input$period[1]
    b = input$period[2]
    rawdata = rawdata[(rawdata[,1] >= a & rawdata[,1]<=b),] 
    
    Simulation_Result = cbind(rawdata[,c(1,2,3)],Simulation_Result)
    #******月尺度数据********
    dat = Rounting()      #"Qsimulated","Qrecorded"
    rawdata = filedata()  #"year","month","day","hour","P","E","Q"
    Yuedata = data.frame(rawdata$year,rawdata$month,rawdata$day,dat$Qsimulated,dat$Qrecorded)
    colnames(Yuedata) = c("year","month","day","Qsimulated","Qrecorded") #流量的单位是:m3/s
    Yuedata[,4][Yuedata[,4]<0] = NA #流量小于0则为无效值
    Yuedata[,5][Yuedata[,5]<0] = NA
    Y_Qi = Yuedata[1,1]
    Y_Mo = Yuedata[length(Yuedata[,1]),1]
    Q_re_mon = NULL
    Q_si_mon = NULL
    Year = Y_Qi:Y_Mo
    for (i in 1:length(Y_Qi:Y_Mo)) {
      for (j in 1:12) {
        #计算月均流量，m3/s
        Q_re_mon[(i-1)*12+j] = mean(Yuedata$Qrecorded[Yuedata$year == Year[i] & Yuedata$month == j],na.rm = T)
        Q_si_mon[(i-1)*12+j] = mean(Yuedata$Qsimulated[Yuedata$year == Year[i] & Yuedata$month == j],na.rm = T)
        
      }
    }
    YY = NULL
    for (i in 1:length(Year)) {
      YY = c(YY,rep(Year[i],12))
    }
    MM = rep(1:12,length(Year))
    Simulation_Monthly = data.frame(YY,MM,Q_re_mon,Q_si_mon)
    colnames(Simulation_Monthly) = c("年","月","实测月均流量（m3/s）","模拟月均流量（m3/s）")
    #***********************************
    switch(input$dataset_download,
           "Runoff generation" = Chanliu_Result,
           "Simulated Process" = Simulation_Result,
           "Simulated monthly Process" = Simulation_Monthly)
    
  })  #select dataset to download
  
  #----- Downloadable csv of selected dataset ----
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$dataset_download,".csv", sep = "")
      },
    content = function(file) {
      write.csv(datasetOutput(), file, row.names = FALSE,
                col.names = T,quote = F)
    }
  )
}


shinyApp(ui = ui, server = server)

