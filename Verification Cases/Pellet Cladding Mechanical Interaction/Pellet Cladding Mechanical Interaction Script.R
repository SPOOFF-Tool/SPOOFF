
CorePath <- "C:/Users/mhort/Documents/lammps-2Aug2023/Verification Cases/Pellet Cladding Mechanical Interaction/"
library(ggplot2)
library(dplyr)
library(data.table)

dtmod = 0.001
fileincrement=1000
fuelradius = 4.1 #mm (This is to the edge of the boundary, not system)
rho = 10970 #kg/m3
endtime = 200
t = seq(10,endtime,10)




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
heatarray<-c(120)
for(i in heatarray){
heat = i #in W/g
heatPath<-paste0(CorePath,"Heat=",heat,"Wperg/")

if(heat==50 ){
  dflammps <- altprocessFile(paste0(heatPath,"dump.LAMMPS"),t-(10000*dtmod),dtmod,fileincrement)
}else{
dflammps <- altprocessFile(paste0(heatPath,"dump.LAMMPS"),t,dtmod,fileincrement)
}

dflammps<- dflammps%>%
  mutate(mindistance=10000)

for(i in 1:length(dflammps$c_t)){
  t_get=dflammps$t[i]
  r_get = c(dflammps$x[i],dflammps$y[i])
  mindistance = 100000
  if(dflammps$type[i]==1){
  lammpstemp<-dflammps %>%
    filter(t==t_get)%>%
    filter(type==2)
  }else{
    lammpstemp<-dflammps %>%
      filter(t==t_get)%>%
      filter(type==1)
  }
  for(j in 1:length(lammpstemp$c_t)){
    r_new = c(lammpstemp$x[j],lammpstemp$y[j])
    mindistance = min(mindistance, sqrt((r_get[1]-r_new[1])^2 +(r_get[2]-r_new[2])^2))
  }
  dflammps$mindistance[i] = mindistance

}
summarisedlammps<-dflammps%>%
  group_by(t) %>%
  summarise(CladdingGap =min(mindistance)*1000) %>%
  mutate(LinearHeat = heat*rho*pi*(fuelradius/1000)^2) #Kw/m

listrms[[length(listrms)+1]]<-summarisedlammps
}

df<-bind_rows(listrms) 
df<- df %>%
  mutate(t=as.numeric(t))
fig<- ggplot() +
  geom_point(data=df, aes(x=t, y=CladdingGap)) +
  geom_line(data=df, aes(x=t, y=CladdingGap)) +
  xlab('Time Step') +
  ylab('Cladding Gap mm') 

fig

ggsave(filename = paste0("Cladding Gap Data.png"),plot=fig,path=CorePath,height = 5, width = 10)


write.csv(df,file=paste0(CorePath,"Cladding Gap Data.csv"),row.names = FALSE)