
CorePath <- "C:/Users/mhort/OneDrive/Documents/Lammps/lammps-stable/lammps-2Aug2023/Testing/Verification Cases/2D Annular Thermal Induced Stress/"
library(ggplot2)
library(dplyr)
library(data.table)

dtmod = 0.1
fileincrement=10
Ti = 1
To = 0
A = 1   #inner radius
B = 11.9  #outer radius
mmax = 1000
kappa = 0.1
spacing = 0.1
rho = 3
rmin = 1
rmax = 11.9
cp = 2
r = seq(rmin,rmax,spacing)
t = c(2000) #just end time of file
E = 100
Nu = 0.3
alpha =  0.0001

xytor<-function(x,y){
  r=sqrt(x^2 + y^2)
  return(r)
}

Tget<-function(r,A,B,Ti,To){
  Tfinal=(Ti*log(B/r) + To*log(r/A))/log(B/A)

  return(Tfinal)
}
Qget<-function(r,A,B,Ti,To,kappa,rho){
  Qfinal=-kappa/rho *(To - Ti)/r*log(B/A)
  
  return(Qfinal)
}
StressRget<-function(r,A,B,Ti,To,E,Nu,alpha){
  StressRfinal=(E*alpha/((1-Nu)*log(B/A)))*((To/2)*(log(r/A)+(B^2)*(1-(A^2)/(r^2))*log(B/A)/(B^2 - A^2))-(Ti/2)*(log(B/r)+(A^2)*(1-(B^2)/(r^2))*log(B/A)/(B^2 - A^2)))
  #StressRfinal= StressRfinal - ((Ti*log(B/r) + To*log(r/A))/log(B/A))*alpha*E /((1-2*Nu))
  return(StressRfinal)
}

StressThetaget<-function(r,A,B,Ti,To,E,Nu,alpha){
  StressThetafinal=(E*alpha/((1-Nu)*log(B/A)))*((To/2)*(1-log(r/A)+(B^2)*(1+(A^2)/(r^2))*log(B/A)/(B^2 - A^2))-(Ti/2)*(-1+log(B/r)+(A^2)*(1+(B^2)/(r^2))*log(B/A)/(B^2 - A^2)))
  
  return(StressThetafinal)
}




  Tarray =c()
  Qarray =c()
  Srarray =c()
  Starray =c()
  for(j in r){
    Tarray[length(Tarray)+1]<- Tget(j,A+0.5,B-0.4,Ti,To)
    Qarray[length(Qarray)+1]<- Qget(j,A+0.5,B-0.4,Ti,To,kappa,rho)
    Srarray[length(Srarray)+1]<- StressRget(j,A,B,Ti,To,E,Nu,alpha)
    Starray[length(Starray)+1]<- StressThetaget(j,A,B,Ti,To,E,Nu,alpha)
  }
  df = data.frame("r"=r, "T" = Tarray, "Q" = Qarray, "Sr" = Srarray, "St" = Starray)

  ggplot() + geom_line(data=df,mapping = aes(x=r,y=Sr))
 




altprocessFile = function(Path,t,dtmod,increment) {
  con = file(Path, "r")
  linenum = 0
  listdata=list()
  
  #check for end of file
  line = readLines(con, n = 1)
  linenum = linenum + 1
  #read time step
  timestep = dtmod*as.numeric(readLines(con, n = 1))
  linenum = linenum + 1
  #skip atom header line
  line = readLines(con, n = 1)
  linenum = linenum + 1
  #read atoms
  atoms = as.numeric(readLines(con, n = 1))
  linenum = linenum + 1
  #skip box dimensions
  line = readLines(con, n = 4)
  linenum = linenum + 4
  #read column headers
  line = readLines(con, n = 1)
  linenum = linenum + 1
  line = gsub("ITEM: ATOMS ","",line)
  headers= strsplit(line,split=" ")
  for (i in 1:length(headers)){headers[[i]]<-gsub("\\[|\\]", "", headers[[i]])}
  start = 1
  if(0 %in% t){
    start = 2
  data <- read.table(Path, skip = linenum, sep = ' ', nrows=atoms, as.is = TRUE, col.names=headers[[1]])
  data <- data %>%
    mutate(t=as.character(timestep))
  listdata[[length(listdata)+1]]<-data
  }
  #skip through atom lines
  line = readLines(con, n = atoms)
  prelinenum = linenum
  linenum = linenum + atoms
  close(con)
  
  for (i in start:length(t)){
    print(t[i])
    lines = prelinenum + linenum*t[i]/(dtmod*increment)
    data <- fread(Path, skip = lines, sep = ' ', nrows=atoms, col.names=headers[[1]],data.table=getOption("datatable.fread.datatable", FALSE))
    data<- as.data.frame(data)
    data <- data %>%
      mutate(t=as.character(t[i]))
    listdata[[length(listdata)+1]]<-data
    rm(data)
  }
  dflammps<-bind_rows(listdata) %>%
    group_by(t)
  
  return(dflammps)
}

errorget = function(dflammps,A,B,Ti,To,kappa,rho, E, Nu, alpha){
  
  dflammps<- dflammps %>%
    mutate(TError = c_t - Tget(r,A,B,Ti,To)) %>%
    mutate(QError = Qr - Qget(r,A,B,Ti,To,kappa,rho)) %>%
    mutate(SrError = Sr - StressRget(r,A,B,Ti,To, E, Nu, alpha)) %>%
    mutate(StError = St - StressThetaget(r,A,B,Ti,To, E, Nu, alpha)) %>%
    mutate(TErrorsq = TError^2) %>%
    mutate(QErrorsq = QError^2) %>%
    mutate(SrErrorsq = SrError^2) %>%
    mutate(StErrorsq = StError^2)
  
  return(dflammps)
} 

#path loop
listrms = list()
dparray<-c(2,1.5,1,0.8,0.5,0.4,0.3,0.2)
#dparray<-c(2,1.5,1,0.8,0.5)
for(i in 1:length(dparray)){
dp = dparray[i]
#for(j in 1:2){
#if(j==1){
  #dp = dp
#}else{
  #dp = paste0(dp," Alt")
#}
dpPath<-paste0(CorePath,"dp=",dp,"/")
print(dpPath)
dflammps <- altprocessFile(paste0(dpPath,"dump.LAMMPS"),t,dtmod,fileincrement)

dflammps <- dflammps %>%
  mutate(r=xytor(x,y)) %>%
  mutate(theta =atan2(y,x)) %>%
  mutate(Qr = c_Q1*cos(theta) + c_Q2*sin(theta)) %>%
  mutate(Sr = c_S1*cos(theta)^2 + c_S2*sin(theta)^2 + 2*c_S4*sin(theta)*cos(theta)) %>%
  mutate(St = c_S1*sin(theta)^2 + c_S2*cos(theta)^2 - 2*c_S4*sin(theta)*cos(theta))

dflammps <- errorget(dflammps,A,B,Ti,To,kappa,rho, E, Nu, alpha)

dflammps <- dflammps %>%
  filter(!r<rmin) %>%
  filter(!r>rmax)

maxT<-max(dflammps$c_t)
minT<-min(dflammps$c_t)
maxQ<-max(dflammps$Qr)
minQ<-min(dflammps$Qr)
maxSr<-max(dflammps$Sr)
minSr<-min(dflammps$Sr)
maxSt<-max(dflammps$St)
minSt<-min(dflammps$St)

fig<- ggplot() +
  geom_line(data=df, aes(x=r, y=T)) +
  geom_point(data=dflammps, aes(x=r, y=c_t), colour="turquoise") +
  ylim(1.1*minT,1.1*maxT) +
  xlab(bquote('r'~(mm))) +
  ylab(bquote('T'~(K))) 

  
fig

ggsave(filename = "Temperature Profile.png",plot=fig,path=dpPath,height = 5, width = 10)

fig<- ggplot() +
  geom_line(data=df, aes(x=r, y=Q)) +
  geom_point(data=dflammps, aes(x=r, y=Qr), colour="turquoise") +
  ylim(1.1*minQ,1.1*maxQ)+
  xlab(bquote('r'~(mm))) +
  ylab(bquote(Q[r]~(g/s^3))) 


fig

ggsave(filename = "Heatflux Profile.png",plot=fig,path=dpPath,height = 5, width = 10)

fig<- ggplot() +
  geom_line(data=df, aes(x=r, y=Sr)) +
  geom_point(data=dflammps, aes(x=r, y=Sr), colour="turquoise") +
  #ylim(1.1*minSr,1.1*maxSr)+
  xlab(bquote('r'~(mm))) +
  ylab(bquote(sigma[r]~(g/mms^2))) 


fig

ggsave(filename = "Stress Radial Profile.png",plot=fig,path=dpPath,height = 5, width = 10)

fig<- ggplot() +
  geom_line(data=df, aes(x=r, y=St)) +
  geom_point(data=dflammps, aes(x=r, y=St), colour="turquoise") +
  #ylim(1.1*minSt,1.1*maxSt)+
  xlab(bquote('r'~(mm))) +
  ylab(bquote(sigma[theta]~(g/mms^2))) 


fig

ggsave(filename = "Stress Tangential Profile.png",plot=fig,path=dpPath,height = 5, width = 10)


fig<- ggplot() +
  geom_hline(yintercept = 0, colour="red", linetype="dashed") +
  geom_point(data=dflammps, aes(x=r, y=TError)) +
  xlab(bquote('r'~(mm))) +
  ylab(bquote('Error in T'~(K)))
fig

ggsave(filename = "Temperature Residual.png",plot=fig,path=dpPath,height = 5, width = 10)

fig<- ggplot() +
  geom_hline(yintercept = 0, colour="red", linetype="dashed") +
  geom_point(data=dflammps, aes(x=r, y=QError)) +
  xlab(bquote('r'~(mm))) +
  ylab(expression(paste("Error in ",Q[r]~(g/s^3))))
fig

ggsave(filename = "Heatflux Residual.png",plot=fig,path=dpPath,height = 5, width = 10)


fig<- ggplot() +
  geom_hline(yintercept = 0, colour="red", linetype="dashed") +
  geom_point(data=dflammps, aes(x=r, y=SrError)) +
  xlab(bquote('r'~(mm))) +
  ylab(expression(paste("Error in ",sigma[r]~(g/s^3))))
fig

ggsave(filename = "Stress Radial Residual.png",plot=fig,path=dpPath,height = 5, width = 10)

fig<- ggplot() +
  geom_hline(yintercept = 0, colour="red", linetype="dashed") +
  geom_point(data=dflammps, aes(x=r, y=StError)) +
  xlab(bquote('r'~(mm))) +
  ylab(expression(paste("Error in ",sigma[theta]~(g/s^3)))) 
fig

ggsave(filename = "Stress Tangential Residual.png",plot=fig,path=dpPath,height = 5, width = 10)

summarisedlammps<-dflammps%>%
  summarise(rms =sqrt(mean(TErrorsq)),rmsQ =sqrt(mean(QErrorsq)),rmsSr =sqrt(mean(SrErrorsq)),rmsSt =sqrt(mean(StErrorsq))) %>%
  mutate(dp = dp)

listrms[[length(listrms)+1]]<-summarisedlammps
}

dfrms<-bind_rows(listrms) %>%
  mutate(logrms = log10(rms)) %>%
  mutate(logrmsQ = log10(rmsQ)) %>%
  mutate(logrmsSr = log10(rmsSr)) %>%
  mutate(logrmsSt = log10(rmsSt)) %>%
  mutate(logdp = log10(dp))



lm_fit <- lm(logdp ~ logrms, data=dfrms)
lm_fitQ <- lm(logdp ~ logrmsQ, data=dfrms)
lm_fitSr <- lm(logdp ~ logrmsSr, data=dfrms)
lm_fitSt <- lm(logdp ~ logrmsSt, data=dfrms)
fitdata <-data.frame(x=min(dfrms$logdp)+0.25*(max(dfrms$logdp)-min(dfrms$logdp)),
                     y=min(dfrms$logrms)+0.75*(max(dfrms$logrms)-min(dfrms$logrms)),
                     text=paste0("y=",format(coef(lm_fit)[[2]], digits=3),"x+",format(coef(lm_fit)[[1]], digits=3)))
fitdataQ <-data.frame(x=min(dfrms$logdp)+0.25*(max(dfrms$logdp)-min(dfrms$logdp)),
                      y=min(dfrms$logrmsQ)+0.75*(max(dfrms$logrmsQ)-min(dfrms$logrmsQ)),
                      text=paste0("y=",format(coef(lm_fitQ)[[2]], digits=3),"x+",format(coef(lm_fitQ)[[1]], digits=3)))

fitdataSr <-data.frame(x=min(dfrms$logdp)+0.25*(max(dfrms$logdp)-min(dfrms$logdp)),
                       y=min(dfrms$logrmsSr)+0.75*(max(dfrms$logrmsSr)-min(dfrms$logrmsSr)),
                       text=paste0("y=",format(coef(lm_fitSr)[[2]], digits=3),"x+",format(coef(lm_fitSr)[[1]], digits=3)))

fitdataSt <-data.frame(x=min(dfrms$logdp)+0.25*(max(dfrms$logdp)-min(dfrms$logdp)),
                       y=min(dfrms$logrmsSt)+0.75*(max(dfrms$logrmsSt)-min(dfrms$logrmsSt)),
                       text=paste0("y=",format(coef(lm_fitSt)[[2]], digits=3),"x+",format(coef(lm_fitSt)[[1]], digits=3)))



fig<- ggplot() +
  geom_point(data=dfrms, aes(x=logdp, y=logrms, colour = t)) +
  geom_smooth(data=dfrms, aes(x=logdp, y=logrms),colour="black",linetype="dashed",method = "lm", se = FALSE) +
  geom_text(data=fitdata, aes(x=x,y=y,label=text)) +
  xlab(bquote('log(Particle Spacing)')) +
  ylab(bquote('log(RMS Error in T)')) 

fig

ggsave(filename = paste0("Temperature Convergence.png"),plot=fig,path=CorePath,height = 5, width = 10)

fig<- ggplot() +
  geom_point(data=dfrms, aes(x=logdp, y=logrmsQ, colour = t)) +
  geom_smooth(data=dfrms, aes(x=logdp, y=logrmsQ),colour="black",linetype="dashed",method = "lm", se = FALSE) +
  geom_text(data=fitdataQ, aes(x=x,y=y,label=text)) +
  xlab(bquote('log(Particle Spacing)')) +
  ylab(bquote('log(RMS Error in Q)')) 

fig

ggsave(filename = paste0("Heatflux Convergence.png"),plot=fig,path=CorePath,height = 5, width = 10)

fig<- ggplot() +
  geom_point(data=dfrms, aes(x=logdp, y=logrmsSr, colour = t)) +
  geom_smooth(data=dfrms, aes(x=logdp, y=logrmsSr),colour="black",linetype="dashed",method = "lm", se = FALSE) +
  geom_text(data=fitdataSr, aes(x=x,y=y,label=text)) +
  xlab(bquote('log(Particle Spacing)')) +
  ylab(expression(paste("log(RMS Error in ",sigma[r],")"))) 

fig

ggsave(filename = paste0("Stress Radial Convergence.png"),plot=fig,path=CorePath,height = 5, width = 10)

fig<- ggplot() +
  geom_point(data=dfrms, aes(x=logdp, y=logrmsSt, colour = t)) +
  geom_smooth(data=dfrms, aes(x=logdp, y=logrmsSt),colour="black",linetype="dashed",method = "lm", se = FALSE) +
  geom_text(data=fitdataSt, aes(x=x,y=y,label=text)) +
  xlab(bquote('log(Particle Spacing)')) +
  ylab(expression(paste("log(RMS Error in ",sigma[theta],")"))) 

fig

ggsave(filename = paste0("Stress Tangential Convergence.png"),plot=fig,path=CorePath,height = 5, width = 10)

totalfitdata <-data.frame(m=c(coef(lm_fit)[[2]],coef(lm_fitQ)[[2]],coef(lm_fitSr)[[2]],coef(lm_fitSt)[[2]]),
                          c=c(coef(lm_fit)[[1]],coef(lm_fitQ)[[1]],coef(lm_fitSr)[[1]],coef(lm_fitSt)[[1]]),
                          t=c(t[1],t[1],t[1],t[1]),
                          Stat = c("T","Q","Sr","St"))

write.csv(dfrms,file=paste0(CorePath,"Convergence_Data.csv"),row.names = FALSE)
write.csv(totalfitdata,file=paste0(CorePath,"Fit_Data.csv"),row.names = FALSE)

