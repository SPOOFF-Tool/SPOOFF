
CorePath <- "C:/Users/mhort/OneDrive/Documents/Lammps/lammps-stable/lammps-2Aug2023/Testing/Verification Cases/1D Transient Heatflow/"
library(ggplot2)
library(dplyr)
library(data.table)

dtmod = 0.1
fileincrement=10
T0 = 1
mmax = 1000
kappa = 0.1
spacing = 0.1
rho = 3
xmin = -10
xmax = 10
cp = 2
x = seq(xmin,xmax,spacing)
L = xmax-xmin
endtime = 1000
t = c(1 %o% 10^(1:log10(endtime)))

Tget<-function(x,t,T0,kappa,mmax,L,rho,cp){
  Tfinal=0
  for(m in 1:mmax){
    Tfinal= Tfinal+(((-1)^(m-1))/((2*m)-1)) * cos(((2*m)-1)*pi*x/L) * exp((-1*(((2*m)-1)^2)*(pi^2)*(kappa/(rho*cp))*t)/(L^2))
  }
  Tfinal=Tfinal*4*T0/pi
  return(Tfinal)
}
Qget<-function(x,t,T0,kappa,mmax,L,rho,cp){
  Qfinal=0
  for(m in 1:mmax){
    Qfinal= Qfinal+((-1)^(m)) * sin(((2*m)-1)*pi*x/L) * exp((-1*(((2*m)-1)^2)*(pi^2)*(kappa/(rho*cp))*t)/(L^2))
  }
  Qfinal=Qfinal*-4*kappa*T0/L
  return(Qfinal)
}
TimeArray <- list()
count=0
for (i in t){
  Tarray =c()
  Qarray =c()
  for(j in x){
    Tarray[length(Tarray)+1]<- Tget(j,i,T0,kappa,mmax,L,rho,cp)
    Qarray[length(Qarray)+1]<- Qget(j,i,T0,kappa,mmax,L,rho,cp)
  }
  d = data.frame("x"=x, "T" = Tarray, "Q" = Qarray, "t"=rep(i,length(x)))
  count=count+1
  TimeArray[[count]]<-d
}

df<-bind_rows(TimeArray) %>%
  mutate(t=as.character(t)) %>%
  group_by(t)




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

errorget = function(dflammps,T0,kappa,mmax,L,rho,cp){
  
  dflammps<- dflammps %>%
    mutate(TError = c_t - Tget(x,as.numeric(t),T0,kappa,mmax,L,rho,cp)) %>%
    mutate(QError = c_Q1 - Qget(x,as.numeric(t),T0,kappa,mmax,L,rho,cp)) %>%
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

dflammps <- errorget(dflammps,T0,kappa,mmax,L,rho,cp)

dflammps <- dflammps %>%
  filter(!x<xmin) %>%
  filter(!x>xmax)

maxT<-max(dflammps$c_t)
minT<-min(dflammps$c_t)
maxQ<-max(dflammps$c_Q1)
minQ<-min(dflammps$c_Q1)

fig<- ggplot() +
  geom_line(data=df, aes(x=x, y=T, colour = t, linetype = t)) +
  geom_point(data=dflammps, aes(x=x, y=c_t, colour = t)) +
  ylim(1.1*minT,1.1*maxT) +
  xlab(bquote('x'~(mm))) +
  ylab(bquote('T'~(K))) +
  scale_colour_discrete(name=bquote('t'~(s))) +
  scale_linetype_discrete(name=bquote('t'~(s))) 
  
fig

ggsave(filename = "Temperature Profile.png",plot=fig,path=dpPath,height = 5, width = 10)

fig<- ggplot() +
  geom_line(data=df, aes(x=x, y=Q, colour = t, linetype = t)) +
  geom_point(data=dflammps, aes(x=x, y=c_Q1, colour = t)) +
  ylim(1.1*minQ,1.1*maxQ)+
  xlab(bquote('x'~(mm))) +
  ylab(bquote('Q'~(g/s^3))) +
  scale_colour_discrete(name=bquote('t'~(s))) +
  scale_linetype_discrete(name=bquote('t'~(s))) 

fig

ggsave(filename = "Heatflux Profile.png",plot=fig,path=dpPath,height = 5, width = 10)

fig<- ggplot() +
  geom_hline(yintercept = 0, colour="red", linetype="dashed") +
  geom_point(data=dflammps, aes(x=x, y=TError, colour = t, shape=t)) +
  xlab(bquote('x'~(mm))) +
  ylab(bquote('Error in T'~(K))) +
  scale_colour_discrete(name=bquote('t'~(s)))  +
  scale_shape_discrete(name=bquote('t'~(s)))
fig

ggsave(filename = "Temperature Residual.png",plot=fig,path=dpPath,height = 5, width = 10)

fig<- ggplot() +
  geom_hline(yintercept = 0, colour="red", linetype="dashed") +
  geom_point(data=dflammps, aes(x=x, y=QError, colour = t, shape=t)) +
  xlab(bquote('x'~(mm))) +
  ylab(bquote('Error in Q'~(g/s^3))) +
  scale_colour_discrete(name=bquote('t'~(s)))  +
  scale_shape_discrete(name=bquote('t'~(s)))
fig

ggsave(filename = "Heatflux Residual.png",plot=fig,path=dpPath,height = 5, width = 10)

summarisedlammps<-dflammps%>%
  group_by(t) %>%
  summarise(rms =sqrt(mean(TErrorsq)),rmsQ =sqrt(mean(QErrorsq))) %>%
  mutate(dp = dp)

listrms[[length(listrms)+1]]<-summarisedlammps
}

dfrms<-bind_rows(listrms) %>%
  group_by(t) %>%
  mutate(logrms = log10(rms)) %>%
  mutate(logrmsQ = log10(rmsQ)) %>%
  mutate(logdp = log10(dp))
listdata<-list()
for(i in t){
dfrmst<-dfrms%>%
  filter(t==as.character(i))

lm_fit <- lm(logdp ~ logrms, data=dfrmst)
lm_fitQ <- lm(logdp ~ logrmsQ, data=dfrmst)
fitdata <-data.frame(x=min(dfrmst$logdp)+0.25*(max(dfrmst$logdp)-min(dfrmst$logdp)),
                     y=min(dfrmst$logrms)+0.75*(max(dfrmst$logrms)-min(dfrmst$logrms)),
                      text=paste0("y=",format(coef(lm_fit)[[2]], digits=3),"x+",format(coef(lm_fit)[[1]], digits=3)))
fitdataQ <-data.frame(x=min(dfrmst$logdp)+0.25*(max(dfrmst$logdp)-min(dfrmst$logdp)),
                     y=min(dfrmst$logrmsQ)+0.75*(max(dfrmst$logrmsQ)-min(dfrmst$logrmsQ)),
                     text=paste0("y=",format(coef(lm_fitQ)[[2]], digits=3),"x+",format(coef(lm_fitQ)[[1]], digits=3)))


fig<- ggplot() +
  geom_point(data=dfrmst, aes(x=logdp, y=logrms, colour = t)) +
  geom_smooth(data=dfrmst, aes(x=logdp, y=logrms),colour="black",linetype="dashed",method = "lm", se = FALSE) +
  geom_text(data=fitdata, aes(x=x,y=y,label=text)) +
  xlab(bquote('log(Particle Spacing)')) +
  ylab(bquote('log(RMS Error in T)')) +
  scale_colour_discrete(name=bquote('t'~(s))) +
  scale_linetype_discrete(name=bquote('t'~(s))) 

fig

ggsave(filename = paste0("Temperature Convergence t=",i,".png"),plot=fig,path=CorePath,height = 5, width = 10)

fig<- ggplot() +
  geom_point(data=dfrmst, aes(x=logdp, y=logrmsQ, colour = t)) +
  geom_smooth(data=dfrmst, aes(x=logdp, y=logrmsQ),colour="black",linetype="dashed",method = "lm", se = FALSE) +
  geom_text(data=fitdataQ, aes(x=x,y=y,label=text)) +
  xlab(bquote('log(Particle Spacing)')) +
  ylab(bquote('log(RMS Error in Q)')) +
  scale_colour_discrete(name=bquote('t'~(s))) +
  scale_linetype_discrete(name=bquote('t'~(s))) 

fig

ggsave(filename = paste0("Heatflux Convergence t=",i,".png"),plot=fig,path=CorePath,height = 5, width = 10)

totalfitdata <-data.frame(m=c(coef(lm_fit)[[2]],coef(lm_fitQ)[[2]]),
                          c=c(coef(lm_fit)[[1]],coef(lm_fitQ)[[1]]),
                          t=c(t,t),
                          Stat = c("T","Q"))
listdata[[length(listdata)+1]]<-totalfitdata
}
totalfitdata<-bind_rows(listdata)
write.csv(dfrms,file=paste0(CorePath,"Convergence_Data.csv"),row.names = FALSE)
write.csv(totalfitdata,file=paste0(CorePath,"Fit_Data.csv"),row.names = FALSE)