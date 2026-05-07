
CorePath <- "C:/Users/mhort/OneDrive/Documents/Lammps/lammps-stable/lammps-2Aug2023/Testing/Verification Cases/Elastic Tensile Failure/"
library(ggplot2)
library(dplyr)
library(data.table)

dtmod = 0.01
fileincrement=200
ymin = -8
ymax = 8
L0=ymax-ymin
Boundary = 1 #thickness of one boundary region
Area = 4 #assume unit thickness
cp = 1
x = seq(xmin,xmax,spacing)
t = c(0,1000*0.01,2000*0.01,3000*0.01,4000*0.01,5000*0.01,6000*0.01,7000*0.01,8000*0.01,9000*0.01,10000*0.01,12000*0.01,15000*0.01,17000*0.01,20000*0.01,22000*0.01,25000*0.01,27000*0.01,30000*0.01,40000*0.01) #just end time of file
#t = c(5000*0.01) 
E=2.3E11
Yeild = 1E-1






  Earray =seq(0,Yeild,0.01)
  Sarray =Earray*E







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
    print(timestep)
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

errorget = function(dflammps){
  
  dflammps<- dflammps %>%
    mutate(SError = c_S2 - E*c_E2) %>%
    mutate(SErrorsq = SError^2) 
  
  return(dflammps)
} 

#path loop
listrms = list()
dparray<-c(2,1.5,1,0.8,0.5,0.4,0.3,0.2)
#dparray<-c(2,1.5,1,0.8,0.5,0.4,0.2)
massarray=(dparray^2)/1
numberarray = dparray
totalmassarray  = dparray
correctedheatarray  = dparray
correctedLfarray  = dparray
#dparray<-c(2,1.5,1,0.8,0.5)
#dparray<-c(1,0.8,0.2)
#dparray<-c(0.5)
icount=0
for(i in dparray){
icount=icount+1
dp = i
dpPath<-paste0(CorePath,"dp=",dp,"/")

dflammps <- altprocessFile(paste0(dpPath,"dump.LAMMPS"),t,dtmod,fileincrement)


df = data.frame("e"=Earray, "S" = Sarray)

#dflammps <- dflammps %>%
#  mutate(c_E2=abs(c_E2)) %>%
#  mutate(c_S2=abs(c_S2))

dflammps <- errorget(dflammps)

StressStrainData =data.frame("S"=c(0),"E"=c(0),"t"=c(0))
for (k in t){
  temp<-dflammps%>%
    filter(t==k)
  if(k==0){
    L0 = max(temp$y) - min(temp$y)
    Top = max(temp$y)
    Bottom = min(temp$y)
    temp<-temp%>%
      filter(y<7)%>%
      filter(y>-7)
    Boundary = (L0 - (max(temp$y) - min(temp$y)))/2
    TopBoundary = Top - max(temp$y)
    BottomBoundary = min(temp$y) - Bottom
    StressStrainData =data.frame("S"=c(0),"E"=c(0),"t"=c(0))
  }else{ 
  extension = ((max(temp$y) - min(temp$y) - (2*Boundary)) - (L0 - (2*Boundary)))/(L0 - (2*Boundary))
  temp<-temp%>%
    filter(y<max(temp$y)-TopBoundary) %>%
    filter(y>max(temp$y)-(TopBoundary+dp))
  print(paste0("Time: ",k," dp: ",i," Particles in Stress Region: ",length(temp$c_S2)))
  Force=0
  for(m in temp$c_S2){
    Force = Force + m
  }
  Force = Force/length(temp$c_S2)
  StressStrainData = rbind(StressStrainData,data.frame("S"=c(Force),"E"=c(extension),"t"=c(k)))
  }
}

dflammps <- dflammps %>%
  filter(c_E2<Yeild) %>%
  filter(c_E2>0)

maxS<-max(dflammps$c_S2)
minS<-min(dflammps$c_S2)
maxE<-max(dflammps$c_E2)
minE<-min(dflammps$c_E2)
#print(df)

fig<- ggplot() +
  geom_line(data=df, aes(x=e, y=S)) +
  geom_point(data=StressStrainData, aes(x=E, y=S, colour = t)) +
  #ylim(0.9*minT,1.1*maxT) +
  xlab(expression(paste(epsilon," (",mm,")"))) +
  ylab(expression(paste(sigma," (",g/s^3,")"))) 


fig

ggsave(filename = "Stress Strain Curve.png",plot=fig,path=dpPath,height = 5, width = 10)

fig<- ggplot() +
  geom_line(data=df, aes(x=e, y=S)) +
  geom_point(data=dflammps, aes(x=c_E2, y=c_S2, colour = t)) +
  #ylim(0.9*minT,1.1*maxT) +
  xlab(expression(paste(epsilon," (",mm,")"))) +
  ylab(expression(paste(sigma," (",g/s^3,")"))) 

  
fig

ggsave(filename = "Particle Stress Strain Curve.png",plot=fig,path=dpPath,height = 5, width = 10)

fig<- ggplot() +
  geom_hline(yintercept = 0, colour="red", linetype="dashed") +
  geom_point(data=dflammps, aes(x=c_E2, y=SError)) +
  xlab(expression(paste(epsilon," (",mm,")"))) +
  ylab(expression(paste("Error in ",sigma," (",g/s^3,")"))) 
fig

ggsave(filename = "Particle Stress Strain Residual.png",plot=fig,path=dpPath,height = 5, width = 10)


summarisedlammps<-dflammps%>%
  summarise(rms =sqrt(mean(SErrorsq))) %>%
  mutate(dp = dp)

listrms[[length(listrms)+1]]<-summarisedlammps
}

dfrms<-bind_rows(listrms) %>%
  mutate(logrms = log10(rms)) %>%
  mutate(logdp = log10(dp))



lm_fit <- lm(logdp ~ logrms, data=dfrms)
fitdata <-data.frame(x=min(dfrms$logdp)+0.25*(max(dfrms$logdp)-min(dfrms$logdp)),
                     y=min(dfrms$logrms)+0.75*(max(dfrms$logrms)-min(dfrms$logrms)),
                      text=paste0("y=",format(coef(lm_fit)[[2]], digits=3),"x+",format(coef(lm_fit)[[1]], digits=3)))


fig<- ggplot() +
  geom_point(data=dfrms, aes(x=logdp, y=logrms, colour = t)) +
  geom_smooth(data=dfrms, aes(x=logdp, y=logrms),colour="black",linetype="dashed",method = "lm", se = FALSE) +
  geom_text(data=fitdata, aes(x=x,y=y,label=text)) +
  xlab(bquote('log(Particle Spacing)')) +
  ylab(expression(paste("log(RMS Error in ",sigma,")"))) 

fig

ggsave(filename = paste0("Stress Strain Convergence.png"),plot=fig,path=CorePath,height = 5, width = 10)



totalfitdata <-data.frame(m=c(coef(lm_fit)[[2]]),
                          c=c(coef(lm_fit)[[1]]),
                          t=c(t[1]),
                          Stat = c("S"))

write.csv(dfrms,file=paste0(CorePath,"Convergence_Data.csv"),row.names = FALSE)
write.csv(totalfitdata,file=paste0(CorePath,"Fit_Data.csv"),row.names = FALSE)
