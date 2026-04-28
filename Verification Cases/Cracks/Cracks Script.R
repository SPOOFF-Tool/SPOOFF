
CorePath <- "C:/Users/mhort/Documents/lammps-2Aug2023/Verification Cases/Cracks/"
library(ggplot2)
library(dplyr)
library(data.table)

dtmod = 0.001
fileincrement=1000
fuelradius = 4.1 #mm (This is to the edge of the boundary, not system)
rho = 10970 #kg/m3
endtime = (50000+400000+50000)*0.001
t = c(endtime)




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
    print(lines)

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
heatarray<-c(10,20,30,40,50,60,70,80,90,100,110,120)
for(i in heatarray){
heat = i #in W/g
heatPath<-paste0(CorePath,"Heat=",heat,"Wperg/")

if(heat==50 ){
  dflammps <- altprocessFile(paste0(heatPath,"dump.LAMMPS"),t-(10000*dtmod),dtmod,fileincrement)
}else{
dflammps <- altprocessFile(paste0(heatPath,"dump.LAMMPS"),t,dtmod,fileincrement)
}


summarisedlammps<-dflammps%>%
  summarise(Damage =mean(c_dam1)*100) %>%
  mutate(LinearHeat = heat*rho*pi*(fuelradius/1000)^2) #Kw/m

listrms[[length(listrms)+1]]<-summarisedlammps
}

df<-bind_rows(listrms) 

fig<- ggplot() +
  geom_point(data=df, aes(x=LinearHeat, y=Damage)) +
  geom_line(data=df, aes(x=LinearHeat, y=Damage)) +
  xlab('Linear Heat Rate kW/m') +
  ylab('Damage %') 

fig

ggsave(filename = paste0("Average Damage Data.png"),plot=fig,path=CorePath,height = 5, width = 10)


write.csv(df,file=paste0(CorePath,"Centre Line Temperature Data.csv"),row.names = FALSE)