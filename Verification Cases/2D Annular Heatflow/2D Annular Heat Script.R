
CorePath <- "C:/Users/mhort/OneDrive/Documents/Lammps/lammps-stable/lammps-2Aug2023/Testing/Verification Cases/2D Annular Heatflow/"
library(ggplot2)
library(dplyr)
library(data.table)

dtmod = 0.1
fileincrement=10
Ti = 2
To = 1
A = 2   #inner radius
B = 10  #outer radius
mmax = 1000
kappa = 0.1
spacing = 0.1
rho = 3
rmin = 2
rmax = 10
cp = 2
r = seq(rmin,rmax,spacing)
t = c(3000) #just end time of file


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


  Tarray =c()
  Qarray =c()
  for(j in r){
    Tarray[length(Tarray)+1]<- Tget(j,A,B,Ti,To)
    Qarray[length(Qarray)+1]<- Qget(j,A,B,Ti,To,kappa,rho)
  }
  df = data.frame("r"=r, "T" = Tarray, "Q" = Qarray)







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

errorget = function(dflammps,A,B,Ti,To,kappa,rho){
  
  dflammps<- dflammps %>%
    mutate(TError = c_t - Tget(r,A,B,Ti,To)) %>%
    mutate(QError = Qr - Qget(r,A,B,Ti,To,kappa,rho)) %>%
    mutate(TErrorsq = TError^2) %>%
    mutate(QErrorsq = QError^2)
  
  return(dflammps)
} 

#path loop
listrms = list()
dparray<-c(2,1.5,1,0.8,0.5,0.4,0.3,0.2)
for(i in dparray){
dp = i
dpPath<-paste0(CorePath,"dp=",dp,"/")

dflammps <- altprocessFile(paste0(dpPath,"dump.LAMMPS"),t,dtmod,fileincrement)

dflammps <- dflammps %>%
  mutate(r=xytor(x,y)) %>%
  mutate(theta =atan2(y,x)) %>%
  mutate(Qr = c_Q1*cos(theta) + c_Q2*sin(theta))

dflammps <- errorget(dflammps,A,B,Ti,To,kappa,rho)

dflammps <- dflammps %>%
  filter(!r<rmin) %>%
  filter(!r>rmax)

maxT<-max(dflammps$c_t)
minT<-min(dflammps$c_t)
maxQ<-max(dflammps$Qr)
minQ<-min(dflammps$Qr)

fig<- ggplot() +
  geom_line(data=df, aes(x=r, y=T)) +
  geom_point(data=dflammps, aes(x=r, y=c_t)) +
  ylim(1.1*minT,1.1*maxT) +
  xlab(bquote('r'~(mm))) +
  ylab(bquote('T'~(K))) 

  
fig

ggsave(filename = "Temperature Profile.png",plot=fig,path=dpPath,height = 5, width = 10)

fig<- ggplot() +
  geom_line(data=df, aes(x=r, y=Q)) +
  geom_point(data=dflammps, aes(x=r, y=Qr)) +
  ylim(1.1*minQ,1.1*maxQ)+
  xlab(bquote('r'~(mm))) +
  ylab(bquote('Q_r'~(g/s^3))) 


fig

ggsave(filename = "Heatflux Profile.png",plot=fig,path=dpPath,height = 5, width = 10)

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
  ylab(bquote('Error in Q'~(g/s^3))) 
fig

ggsave(filename = "Heatflux Residual.png",plot=fig,path=dpPath,height = 5, width = 10)

summarisedlammps<-dflammps%>%
  summarise(rms =sqrt(mean(TErrorsq)),rmsQ =sqrt(mean(QErrorsq))) %>%
  mutate(dp = dp)

listrms[[length(listrms)+1]]<-summarisedlammps
}

dfrms<-bind_rows(listrms) %>%
  mutate(logrms = log10(rms)) %>%
  mutate(logrmsQ = log10(rmsQ)) %>%
  mutate(logdp = log10(dp))



lm_fit <- lm(logdp ~ logrms, data=dfrms)
lm_fitQ <- lm(logdp ~ logrmsQ, data=dfrms)
fitdata <-data.frame(x=min(dfrms$logdp)+0.25*(max(dfrms$logdp)-min(dfrms$logdp)),
                     y=min(dfrms$logrms)+0.75*(max(dfrms$logrms)-min(dfrms$logrms)),
                      text=paste0("y=",format(coef(lm_fit)[[2]], digits=3),"x+",format(coef(lm_fit)[[1]], digits=3)))
fitdataQ <-data.frame(x=min(dfrms$logdp)+0.25*(max(dfrms$logdp)-min(dfrms$logdp)),
                     y=min(dfrms$logrmsQ)+0.75*(max(dfrms$logrmsQ)-min(dfrms$logrmsQ)),
                     text=paste0("y=",format(coef(lm_fitQ)[[2]], digits=3),"x+",format(coef(lm_fitQ)[[1]], digits=3)))


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

totalfitdata <-data.frame(m=c(coef(lm_fit)[[2]],coef(lm_fitQ)[[2]]),
                          c=c(coef(lm_fit)[[1]],coef(lm_fitQ)[[1]]),
                          t=c(t[1],t[1]),
                          Stat = c("T","Q"))


write.csv(dfrms,file=paste0(CorePath,"Convergence_Data.csv"),row.names = FALSE)
write.csv(totalfitdata,file=paste0(CorePath,"Fit_Data.csv"),row.names = FALSE)