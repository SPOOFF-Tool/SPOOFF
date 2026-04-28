
CorePath <- "C:/Users/mhort/OneDrive/Documents/Lammps/lammps-stable/lammps-2Aug2023/Testing/Verification Cases/Centreline Convection/"
library(ggplot2)
library(dplyr)
library(data.table)

dtmod = 0.01
fileincrement=100
fuelradius = 4.1 #mm (This is to the edge of the boundary, not system)
rho = 10970 #kg/m3
endtime = 500
t = c(endtime)

enigmadata<-as.data.frame(na.omit(read.csv(paste0(CorePath,"data CSV.csv"))))


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



#path loop
listrms = list()
heatarray<-c(1,5,10,20,30,40,50,60,70,80,90,100,110,120)
for(i in heatarray){
heat = i #in W/g
heatPath<-paste0(CorePath,"Heat=",heat,"Wperg/")

dflammps <- altprocessFile(paste0(heatPath,"dump.LAMMPS"),t,dtmod,fileincrement)

if(heat==50){
  dflammps50 <- dflammps
}

dflammps <- dflammps %>%
  filter(c_r1<fuelradius)


summarisedlammps<-dflammps%>%
  summarise(Centre =max(c_t),Outer =min(c_t)) %>%
  mutate(LinearHeat = (heat*rho*pi*(fuelradius/1000)^2)) #Kw/m

listrms[[length(listrms)+1]]<-summarisedlammps
}

df<-bind_rows(listrms) 

Key<-"SPOOFF Centre Line Temperature"
name2<-"SPOOFF Outer Pellet Temperature"
name3<-"ENIGMA Centre Line Temperature"
name4<-"ENIGMA Outer Pellet Temperature"

fig<- ggplot() +
  geom_point(data=df, aes(x=LinearHeat, y=Centre, colour = Key, shape = Key)) +
  geom_point(data=df, aes(x=LinearHeat, y=Outer, colour = name2, shape = name2)) +
  geom_line(data=df, aes(x=LinearHeat, y=Centre, colour = Key, linetype = Key)) +
  geom_line(data=df, aes(x=LinearHeat, y=Outer, colour = name2, linetype = name2)) +
  geom_point(data=enigmadata, aes(x=Linear_Heat, y=Centre_Line, colour = name3, shape = name3)) +
  geom_point(data=enigmadata, aes(x=Linear_Heat, y=Outer, colour = name4, shape = name4)) +
  geom_line(data=enigmadata, aes(x=Linear_Heat, y=Centre_Line, colour = name3, linetype = name3)) +
  geom_line(data=enigmadata, aes(x=Linear_Heat, y=Outer, colour = name4, linetype = name4)) +
  xlab('Linear Heat Rate kW/m') +
  ylab('Temperature K') + theme(legend.position="bottom")

fig

ggsave(filename = paste0("Centre Line Temperature Data.png"),plot=fig,path=CorePath,height = 5, width = 10)


write.csv(df,file=paste0(CorePath,"Centre Line Temperature Data.csv"),row.names = FALSE)

contourdata=data.frame()
for(i in seq(min(dflammps50$x)-0.2,max(dflammps50$x)+0.2,0.4)){
  for(j in seq(min(dflammps50$y)-0.2,max(dflammps50$y)+0.2,0.4)){
    tempdata=data.frame("x"=i,"y"=j,"T"=0)
    for(k in 1:length(dflammps50$y)){
      if(abs(dflammps50$y[k]-j)<0.001&abs(dflammps50$x[k]-i)<0.001){
        tempdata$T[1]=dflammps50$c_t[k]
      }
    }
    contourdata=rbind(contourdata,tempdata)
  }
}

dflammps50<-dflammps50 %>%
  mutate(Type =ifelse(c_r1<4.1,"Fuel",ifelse(c_r1<5.0,"Cladding","Coolant")))

fig<- ggplot() +
  geom_point(data=dflammps50, mapping =aes(x=x,y=y,colour = Type)) 

fig

ggsave(filename = paste0("Centre Line Mesh Data.png"),plot=fig,path=CorePath,height = 5, width = 6)

