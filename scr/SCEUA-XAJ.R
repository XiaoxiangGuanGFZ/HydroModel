#--------SCEUA-XAJ---------------
# 输入资料：流域面积，流域面均降水量（mm）；
# 流域潜在蒸散发量（mm）；出口断面实测流量（m3/s）；
# 所有资料都为日尺度
# programmer:Shane Guan
# e-mail:xxguan@hhu.edu.cn

#由于日模型与次洪模型的计算时段长不同，参数值不能全部通用，但K、WM、WUM、WLM、B、IMP、EX、C与时段长无关可以直接引用 
#Kc SM、Kg、KSS、Cs、Ci、Cg与时段长相关不能直接引用需要另外率定
library(rtop)
#---------------------SCE-UA-XAJ model-------------------
Dunne <- function(P,E,DT,Imp,Kc,WDM,WUM,WLM,B,C,Ex,SM,Ki,Kg){ #,FR0
  
  WM = WUM + WLM + WDM
  D=24/DT  #可以划分的时段数，一天可以划分的时段数
  KSSD = (1-(1-(Ki+Kg))^(1/D))/(1+Kg/Ki) #壤中流时段出流系数；包为民主编《水文预报（第四版）》P149，式5-33
  KGD = KSSD*Kg/Ki    #地下径流时段出流系数，P149，式5-33
  Epp = Kc*E   #Epp为蒸发能力 
  PE = P - Epp
  Smm = (1 + Ex) * SM   #函数内的全局变量，Smm为流域单点最大的自由水蓄水容量，mm          
  EU = NULL  #上层蒸发量
  EL = NULL  #下层蒸发量
  ED = NULL  #深层蒸发量
  WL = NULL  #下层土壤含水量
  WU = NULL  #上层土壤含水量
  WD = NULL  #深层土壤含水量
  W = NULL   #流域平均土壤含水量,中间时段量 
  R = NULL  #产流量，mm
  RS = NULL  #地表径流出流径流深，mm
  RSS = NULL #壤中流出流径流深，mm
  RG = NULL  #地下径流出流径流深，mm
  
  WU = WUM/4  #设置土壤含水的初始状态
  WL = WLM/9*7
  WD = WDM
  W = WU+WL+WD
  len = length(P)
  S0 = SM/4  #初始自由水箱含水量（深），mm
  a = NULL
  
  for (t in 1:len){
    #---------*先*产流计算部分--------------
    R[t] = 0       #总产流量预设为0，mm
    Rimp = 0       #不透水面积产生的直接径流量，mm
    WMM = (1+B)*WM/(1-Imp)    # Imp不透水面积比;
    #WMM：流域最大单点土壤含水量，为蓄水容量极值  
    #B为蓄水容量曲线的指数
    if( PE[t] > 0.001 ) {                
      Rimp = PE[t]*Imp          #不透水面积产生的直接径流,Imp < 1 
      PE[t] = PE[t]-Rimp        #对PE[t]作出调整，透水面积上的净雨
      if (abs(WM - W) <= 0.001 ){  
        a = WMM       #WMM为蓄水容量极值,单点最大值
      } else {
        a = WMM*(1-(1-W/WM)^(1/(1 + B)))     #a为与W[t]对应的在蓄水容量曲线的纵坐标 
      }
      
      if ((PE[t] + a) <= WMM){           #流域部分产流   
        R[t] = PE[t]-WM+W + WM * ((1- (PE[t] + a)/ WMM)^(1 + B)) 
      } else {
        R[t] = PE[t] - (WM - W)     #流域全面产流 
      }
      
      if (abs(R[t] - PE[t]) <= 0.001) {
        R[t] = PE[t] 
      } 
      
    }  else  { #PE < 0时，净雨深为负时，产流为0；
      Rimp = 0
      R[t] = 0
    }
    
    ##-------------三层蒸散发计算------------- 
    if (PE[t] > 0) {
      if (WU + PE[t] - R[t] > WUM) {
        if (WU + WL + PE[t] - R[t] - WUM > WLM ) {
          WU = WUM
          WL = WLM
          WD = W + PE[t] - R[t] - WU - WL
        } else {
          WL = WU + WL +PE[t] - R[t] - WUM
          WU = WUM
        }
      } else {
        WU = WU + PE[t] - R[t]
      }
      EU[t] = Epp[t]
      EL[t] = 0.0
      ED[t] = 0.0
    } else {
      if (WU + PE[t] >= 0) {
        EU[t] = Epp[t]
        EL[t] = 0.0
        ED[t] = 0.0
        WU = WU + PE[t]
      } else {
        EU[t] = WU + P[t]
        WU = 0
        if (WL >= C*WLM) {
          EL[t] = (Epp[t]-EU[t]) * WL/WLM
          WL = WL - EL[t]
          ED[t] = 0
        } else {
          if (WL >= C*(Epp[t] - EU[t])) {
            EL[t] = C*(Epp[t] - EU[t])
            WL = WL - EL[t]
            ED[t] = 0
            
          } else {
            EL[t] = WL
            WL = 0
            ED[t] = C*(Epp[t] - EU[t])-EL[t]
            WD = WD - ED[t] 
          }
        }
      }
    }
    E = EU[t] + EL[t] + ED[t]

    if (WU > WUM) {WU = WUM}
    if (WU <0 ) {WU = 0}
    if (WL > WLM) {WL = WLM}
    if (WL <0 ) {WL = 0}
    if (WD > WDM) {WD = WDM}
    if (WD <0 ) {WD = 0}
    W = WU + WL + WD
    R[t] = R[t]*(1-Imp)
    if (W > WM) {W = WM}
    #-----------水源划分--------------------------- 
    
    if (t == 1) {
      S = S0
    }
    RS[t] = 0  #给t时段，地表径流设置初值，为0
    RSS[t] = 0 #给t时段，壤中流设置初值，为0
    RG[t] = 0  #给t时段，地下径流设置初值，为0
    #FR0为时段初(上一时段的)产流面积 
    if (PE[t] <= 0 | abs(PE[t]) < 1e-3 | R[t] <= 0.0005)                   #无净雨   
    { 
      #FR = 0   #此时段，产流面积为0
      RS[t] = 0
      RSS[t] = S*KSSD                                                                                                                                                                                                                                                                                                                                      #依靠时段初自由水箱含水量计算出流；      
      #地下水出流出流系数 
      RG[t] = S*KGD
      S = S-(RSS[t]+RG[t])   #s表示自由水在产流面积上的平均蓄水深
    } else { 
      #FR = R[t]/PE[t]         #用流量除以单位面积上的净雨可以理解为产流深即得产面面，产流面积时段末产流面积
      Q = R[t] 
      NN = floor(Q/5) + 1     #每次入流按5毫米分成并取整数NN为了消除前向差分误差,截尾取整
      if (NN == 1) { #Q还不足5mm，差分造成的误差影响不大
        Kssdd = KSSD
        Kgdd = KGD
        dQ = Q
      } else {
        Kssdd = (1-(1-(KGD + KSSD))^(1/NN))/(1+KGD/KSSD) 
        Kgdd = Kssdd*KGD/KSSD 
        dQ = round(Q/NN,digits = 3)                 #Q分为NN个时段,dQ为每个小时段流入自由水箱的量
      }
      #SM流域的平均自由水蓄水容量 
      #Smm全流域最大单点的自由水蓄水容量 
      Rsd = 0
      Rssd = 0
      Rgd = 0
      for (j in 1:NN) { 
        if (S >= SM)   {S = SM}
        
        #AU = Smm*(1-(1-S*FR0/(FR+0.0001)/SM + 0.0001)^(1/(1 + Ex)))   #Smm为流域单点最大的自由水蓄水容量
        AU = Smm*(1-(1-S/SM + 0.0001)^(1/(1 + Ex)))   #Smm为流域单点最大的自由水蓄水容量
        #如果FR = 0，则S*FR0/FR/SM在算术上不成立，因为被除数FR和SM不能为0；而上一时刻FR0可以为0
        
        if  (as.logical(dQ + AU >= Smm))           #当径流与此时刻的平均蓄水深之和大于最大平均蓄水深全面产壤中流 
        {
          #Rsd = (dQ +S * FR0 /(FR+0.0001) - SM)*FR
          Rsd = dQ + S - SM
          Rssd = S * Kssdd #* FR 
          Rgd = S * Kgdd #* FR 
          S = SM * (1-Kssdd-Kgdd)        #此时刻，自由水箱已经蓄满了
        }else if (as.logical(dQ + AU < Smm))            #当径流与此时刻的平均蓄水深之和大于最大平均蓄水深部分产壤中流 
        {
          #Rsd = FR*(dQ+ S*FR0/(FR+0.0001) - SM + SM *(1-(dQ + AU)/Smm)^(1+Ex))
          Rsd = dQ+ S - SM + SM *(1-(dQ + AU)/(Smm+0.001) )^(1+Ex)
          Rssd = S * Kssdd# * FR 
          Rgd = S * Kgdd #* FR 
          S = dQ + S * (1 - Kssdd - Kgdd)  #S*FR0/(FR+0.0001) + (R[t]-Rsd)/(FR+0.0001)  #本时段自由水蓄量
        }
        # Rssd = S * Kssdd * FR 
        # Rgd = S * Kgdd * FR 
        # S = dQ + S * (1 - Kssdd - Kgdd)    #下一时段自由水蓄量
        RS[t] = RS[t] + Rsd        
        RSS[t] = RSS[t] + Rssd     
        RG[t] = RG[t] + Rgd 
      }
    } 
    RS[t] = RS[t] + Rimp
    #FR0 = FR   #作为下一时段计算的初始产流面积
    
  } #循环计算结束
  RS = round(RS,digits = 3)
  RSS = round(RSS,digits = 3)
  RG = round(RS,digits = 3)
  out = data.frame(RS,RSS,RG)
  out = `colnames<-`(out,c("RS","RI","RG"))
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
#------------滞后演算---------------
LagRounting <- function(I,CR,L) {
  L = floor(L)
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
#-------XAJ模型计算函数--------------
XAJ = function(P, E, Q, DAREA, DT,
               Kc, SM, B, C,Ex,WUM, WLM, WDM, 
               Kg, Ki, 
               Cg, Cs, Ci, CR,L) {
  #FR0 = 0.01
  Imp = 0.02
  runoff <- Dunne(P,E,DT,Imp,Kc,WDM,WUM,WLM,B,C,Ex,SM,Ki,Kg) #,FR0
  RS = runoff$RS
  RI = runoff$RI
  RG = runoff$RG
  ini = 0.125
  basinrouting <- LinearReservoir(DAREA,DT,RS,RI,RG,Cs,Ci,Cg,ini) #"Qrs","Qri","Qrg","Qrt"
  I = basinrouting$Qrt    #坡面汇流之后的流量
  out <- LagRounting(I,CR,L)
  Qout = out$O
  return(Qout) #模拟的径流过程，单位m3/s
}

fit = function(x) {
  Kc = x[1]; SM = x[2];WUM = x[6];WLM = x[7];WDM = x[9];
  B = x[3];C = x[4];Ex = x[5];Kg = x[9];Ki = x[10];
  Cs = x[11];Ci = x[12];Cg = x[13];CR = x[14];L = x[15]
  out = XAJ(Validinput_P, Validinput_E, Validinput_Q, Validinput_DAREA, Validinput_DT,
            Kc, SM, B, C,Ex,WUM, WLM, WDM, 
            Kg, Ki, 
            Cg, Cs, Ci, CR,L) 
  Qsi = out
  Qre = Validinput_Q
  R2 = 1-sum((Qre-Qsi)^2,na.rm = T)/sum((Qre-mean(Qre))^2,na.rm = T) #NSE
  return(1-R2)
}
#************test****************
dat = read.table("G:/分析降水资料对模型率定的影响/闽江金溪池潭水库流域数据/Chitan_AA_stations日1982-2013_4766.csv",header = T,sep = ",")
dat = dat[dat[,1]<=2002,]
P = dat[,4]
E = dat[,5]
Q = dat[,6]
DAREA = 4766
DT = 24

# Imp = 0.019
# Kc = 0.6
# WDM = 10
# WUM = 180
# WLM = 40
# B = 0.4
# C = 0.1
# Ex = 0.5
# SM = 20
# Ki = 0.4
# Kg = 0.3
# FR0 = 0.002
# Cg = 0.9
# Ci = 0.9
# Cs = 0.4
# L = 0
# CR = 0.3

#-----定义全局变量---------
Validinput_P <<- P
Validinput_E <<- E
Validinput_Q <<- Q
Validinput_DAREA <<- DAREA
Validinput_DT <<- 24
#-------------1    2    3     4    5     6     7    8     9    10   11   12   13   14  15
#-------------Kc, SM,   B,    C   Ex   WUM   WLM  WDM    Kg    Ki   Cs   Ci   Cg   CR  L
paralower = c(0.1,5  ,  0.1, 0.05, 1.0,  40,  20,  5, 0.25, 0.25, 0.1, 0.8, 0.8, 0.1, 0)
paraupper = c(1  ,50 ,  1.5,  0.2, 1.4, 300, 150, 60, 0.49, 0.49,   6,   1,   1, 0.6, 3)
parainitial=c(0.6,20 ,  0.4,  0.1, 1.2, 180,  40, 10,  0.3,  0.4, 0.4, 0.9, 0.9, 0.3, 0)
#x = parainitial
sceua_xaj = sceua(fit, pars = parainitial, lower = paralower, 
                  upper = paraupper, maxn = 1000, pcento = 0.00001,iprint = 0)
fit(sceua_xaj$par)
#----------------------------
