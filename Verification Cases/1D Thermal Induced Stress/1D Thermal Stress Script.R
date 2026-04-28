
CorePath <- "C:/Users/mhort/OneDrive/Documents/Lammps/lammps-stable/lammps-2Aug2023/Testing/Verification Cases/1D Thermal Induced Stress/"
library(ggplot2)
library(dplyr)
library(data.table)

dtmod = 1
fileincrement=10
TH = 10
TC = 0
A = 2   #inner radius
B = 10  #outer radius
mmax = 1000
kappa = 0.1
spacing = 0.1
rho = 3
xmin = -10
xmax = 10
L=xmax-xmin
cp = 2
x = seq(xmin,xmax,spacing)
t = c(20000) #just end time of file
E = 100
Nu = 0.3
alpha =  0.001


Tget<-function(x,TH,TC,L,xmax){
  Tfinal=(TH-TC)*x/L + (TH+TC)/2

  return(Tfinal)
}
Qget<-function(TH,TC,L,kappa){
  Qfinal=-kappa*(TH-TC)/L
  
  return(Qfinal)
}
Sget<-function(x,E,alpha,TH,TC,xmax,xmin,L){
  S=-alpha*E*((TH-TC)*x/L + (TH+TC)/2)
  
  return(S)
}



  Tarray =c()
  Qarray =c()
  Sarray =c()
  for(j in x){
    Tarray[length(Tarray)+1]<- Tget(j,TH,TC,L,xmax)
    Qarray[length(Qarray)+1]<- Qget(TH,TC,L,kappa)
    Sarray[length(Sarray)+1]<- Sget(j,E,alpha,TH,TC,xmax,xmin,L)
  }
  df = data.frame("x"=x, "T" = Tarray, "Q" = Qarray, "S" = Sarray)







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

errorget = function(dflammps,x,L,E,alpha,TH,TC,xmax,xmin){
  
  dflammps<- dflammps %>%
    mutate(TError = c_t - Tget(x,TH,TC,L,xmax)) %>%
    mutate(QError = c_Q1 - Qget(TH,TC,L,kappa)) %>%
    mutate(SError = c_S1 - Sget(x,E,alpha,TH,TC,xmax,xmin,L)) %>%
    mutate(TErrorsq = TError^2) %>%
    mutate(QErrorsq = QError^2) %>%
    mutate(SErrorsq = SError^2) 
  
  return(dflammps)
} 

#path loop
listrms = list()
dparray<-c(2,1.5,1,0.8,0.5,0.4,0.3,0.2)
#dparray<-c(2,1.5,1,0.8,0.5)
#dparray<-c(1)
for(i in dparray){
dp = i
dpPath<-paste0(CorePath,"dp=",dp,"/")

dflammps <- altprocessFile(paste0(dpPath,"dump.LAMMPS"),t,dtmod,fileincrement)


dflammps <- errorget(dflammps,x,L,E,alpha,TH,TC,xmax,xmin)

dflammps <- dflammps %>%
  filter(!x<xmin) %>%
  filter(!x>xmax)

maxT<-max(dflammps$c_t)
minT<-min(dflammps$c_t)
maxQ<-max(dflammps$c_Q1)
minQ<-min(dflammps$c_Q1)
maxS<-max(dflammps$c_S1)
minS<-min(dflammps$c_S1)

fig<- ggplot() +
  geom_line(data=df, aes(x=x, y=T)) +
  geom_point(data=dflammps, aes(x=x, y=c_t)) +
  ylim(1.1*minT,1.1*maxT) +
  xlab(bquote('x'~(mm))) +
  ylab(bquote('T'~(K))) 

  
fig

ggsave(filename = "Temperature Profile.png",plot=fig,path=dpPath,height = 5, width = 10)

fig<- ggplot() +
  geom_line(data=df, aes(x=x, y=Q)) +
  geom_point(data=dflammps, aes(x=x, y=c_Q1)) +
  ylim(1.1*minQ,1.1*maxQ)+
  xlab(bquote('x'~(mm))) +
  ylab(bquote('Q_r'~(g/s^3))) 


fig

ggsave(filename = "Heatflux Profile.png",plot=fig,path=dpPath,height = 5, width = 10)

fig<- ggplot() +
  geom_line(data=df, aes(x=x, y=S)) +
  geom_point(data=dflammps, aes(x=x, y=c_S1)) +
  ylim(1.1*minS,1.1*maxS)+
  xlab(bquote('x'~(mm))) +
  ylab(expression(paste(sigma," (",g/s^3,")"))) 


fig

ggsave(filename = "Stress Profile.png",plot=fig,path=dpPath,height = 5, width = 10)


fig<- ggplot() +
  geom_hline(yintercept = 0, colour="red", linetype="dashed") +
  geom_point(data=dflammps, aes(x=x, y=TError)) +
  xlab(bquote('x'~(mm))) +
  ylab(bquote('Error in T'~(K)))
fig

ggsave(filename = "Temperature Residual.png",plot=fig,path=dpPath,height = 5, width = 10)

fig<- ggplot() +
  geom_hline(yintercept = 0, colour="red", linetype="dashed") +
  geom_point(data=dflammps, aes(x=x, y=QError)) +
  xlab(bquote('x'~(mm))) +
  ylab(bquote('Error in Q'~(g/s^3))) 
fig

ggsave(filename = "Heatflux Residual.png",plot=fig,path=dpPath,height = 5, width = 10)


fig<- ggplot() +
  geom_hline(yintercept = 0, colour="red", linetype="dashed") +
  geom_point(data=dflammps, aes(x=x, y=SError)) +
  xlab(bquote('x'~(mm))) +
  ylab(expression(paste('Error in ',sigma," (",g/s^3,")")))
fig

ggsave(filename = "Stress Residual.png",plot=fig,path=dpPath,height = 5, width = 10)

summarisedlammps<-dflammps%>%
  summarise(rms =sqrt(mean(TErrorsq)),rmsQ =sqrt(mean(QErrorsq)),rmsS =sqrt(mean(SErrorsq))) %>%
  mutate(dp = dp)

listrms[[length(listrms)+1]]<-summarisedlammps
}

dfrms<-bind_rows(listrms) %>%
  mutate(logrms = log10(rms)) %>%
  mutate(logrmsQ = log10(rmsQ)) %>%
  mutate(logrmsS = log10(rmsS)) %>%
  mutate(logdp = log10(dp))



lm_fit <- lm(logdp ~ logrms, data=dfrms)
lm_fitQ <- lm(logdp ~ logrmsQ, data=dfrms)
lm_fitS <- lm(logdp ~ logrmsS, data=dfrms)
fitdata <-data.frame(x=min(dfrms$logdp)+0.25*(max(dfrms$logdp)-min(dfrms$logdp)),
                     y=min(dfrms$logrms)+0.75*(max(dfrms$logrms)-min(dfrms$logrms)),
                      text=paste0("y=",format(coef(lm_fit)[[2]], digits=3),"x+",format(coef(lm_fit)[[1]], digits=3)))
fitdataQ <-data.frame(x=min(dfrms$logdp)+0.25*(max(dfrms$logdp)-min(dfrms$logdp)),
                     y=min(dfrms$logrmsQ)+0.75*(max(dfrms$logrmsQ)-min(dfrms$logrmsQ)),
                     text=paste0("y=",format(coef(lm_fitQ)[[2]], digits=3),"x+",format(coef(lm_fitQ)[[1]], digits=3)))

fitdataS <-data.frame(x=min(dfrms$logdp)+0.25*(max(dfrms$logdp)-min(dfrms$logdp)),
                      y=min(dfrms$logrmsS)+0.75*(max(dfrms$logrmsS)-min(dfrms$logrmsS)),
                      text=paste0("y=",format(coef(lm_fitS)[[2]], digits=3),"x+",format(coef(lm_fitS)[[1]], digits=3)))


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
  geom_point(data=dfrms, aes(x=logdp, y=logrmsS, colour = t)) +
  geom_smooth(data=dfrms, aes(x=logdp, y=logrmsS),colour="black",linetype="dashed",method = "lm", se = FALSE) +
  geom_text(data=fitdataS, aes(x=x,y=y,label=text)) +
  xlab(bquote('log(Particle Spacing)')) +
  ylab(expression(paste("log(RMS Error in ",sigma,")"))) 

fig

ggsave(filename = paste0("Stress Convergence.png"),plot=fig,path=CorePath,height = 5, width = 10)

totalfitdata <-data.frame(m=c(coef(lm_fit)[[2]],coef(lm_fitQ)[[2]],coef(lm_fitS)[[2]]),
                          c=c(coef(lm_fit)[[1]],coef(lm_fitQ)[[1]],coef(lm_fitS)[[1]]),
                          t=c(t[1],t[1],t[1]),
                          Stat = c("T","Q","S"))

write.csv(dfrms,file=paste0(CorePath,"Convergence_Data.csv"),row.names = FALSE)
write.csv(totalfitdata,file=paste0(CorePath,"Fit_Data.csv"),row.names = FALSE)
