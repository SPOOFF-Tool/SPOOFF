
CorePath <- "C:/Users/mhort/Documents/lammps-2Aug2023/Verification Cases/Pellet Cladding Mechanical Interaction/"
library(ggplot2)
library(dplyr)
library(data.table)
library(viridis)

dtmod = 0.001
fileincrement=1000
fuelradius = 4.1 #mm (This is to the edge of the boundary, not system)
rho = 10970 #kg/m3
endtime = 650
#t = seq(10,endtime,10)
t= c(endtime)
atom_num = 6995
spacing = 0.5
Area = spacing^2 





processbondFile = function(Path,t,dtmod,increment) {
  Header_len= 9
  lines=0
  step = 0
  stepprotect=1000
  listdata=list()
  while(TRUE){
    step = step + 1
    if(step>stepprotect){break}
    TimeStep = as.numeric(fread(Path, skip = 1+lines, sep = ' ', nrows=1))
    Time= TimeStep*dtmod
    print(Time)
    atoms = as.numeric(fread(Path, skip = 3+lines, sep = ' ', nrows=1))
    if(Time %in% t){
    headers= fread(Path, skip = 8+lines, sep = ' ', nrows=1, header=FALSE)
    headers = headers[,-1]
    headers = headers[,-1]
    data <- fread(Path, skip = Header_len+lines, sep = ' ', nrows=atoms, header=FALSE,data.table=getOption("datatable.fread.datatable", FALSE))
    data<- as.data.frame(data)
    colnames(data)<-headers
    data <- data %>%
      mutate(t=as.character(Time))
    listdata[[length(listdata)+1]]<-data
    }
    lines=Header_len+lines+atoms
    if(Time==endtime){break}
    
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

dflammps <- processbondFile(paste0(heatPath,"dump2.LAMMPS"),t,dtmod,fileincrement)
dfpairs<-dflammps

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


dtmod = 0.001
fileincrement=1000
fuelradius = 4.1 #mm (This is to the edge of the boundary, not system)
rho = 10970 #kg/m3
endtime = 200
t = c(endtime)

dflammps <- altprocessFile(paste0(heatPath,"dump.LAMMPS"),t,dtmod,fileincrement)
dflammps<- dflammps%>%
  mutate(Fragment=0)

Fragment_list = c()
Atom_list =c(rep(0,atom_num))
Atom_Fuel=c(rep(0,atom_num))
Atom_Count = 0

for(k in 1:length(dfpairs$patom1)){
  i = dfpairs$patom1[k]
  j = dfpairs$patom2[k]
  if(dflammps$c_r1[i]<4.25&dflammps$c_r1[j]<4.25){
    Atom_Fuel[i]=1
    Atom_Fuel[j]=1
  if(dfpairs$degredation[k]<0.1){
  if(Atom_list[i] == 0 & Atom_list[j] == 0){
    if(length(Fragment_list)==0){
      fragment_num=1
      }else{
        adjust=0
        fragment_num=0
        while(fragment_num==0){
          fragment_num = Fragment_list[length(Fragment_list)-adjust]
          adjust = adjust + 1
        }
    fragment_num = fragment_num+1
    }
    Fragment_list[fragment_num] = fragment_num
    Atom_list[i] = Atom_list[j] = fragment_num
  }else if(Atom_list[i] == 0){
    fragment_num = Atom_list[j]
    Atom_list[i] = fragment_num
  }else if(Atom_list[j] == 0){
    fragment_num = Atom_list[i]
    Atom_list[i] = fragment_num
  }else if(Atom_list[j]!=Atom_list[i]){
    print(paste0("Found link Between Fragments: ",Atom_list[i],Atom_list[j],"... Ajusting Fragment List... Removing:",max(Atom_list[i],Atom_list[j])))
    fragment_num = min(Atom_list[i],Atom_list[j])
    fragment_remove= max(Atom_list[i],Atom_list[j])
    Fragment_list[fragment_remove] = 0
    for(m in 1:length(Atom_list)){
      if(Atom_list[m]==fragment_remove){Atom_list[m]=fragment_num}
    }
  }
  }
  }
}

for (i in 1: length(Atom_list)){
  if(Atom_Fuel[i]==1){
    Atom_Count=Atom_Count+1
  if(Atom_list[i]==0){
    if(length(Fragment_list)==0){
      fragment_num=1
    }else{
      adjust=0
      fragment_num=0
      while(fragment_num==0){
        fragment_num = Fragment_list[length(Fragment_list)-adjust]
        adjust = adjust + 1
      }
      fragment_num = fragment_num+1
    }
    Fragment_list[fragment_num] = fragment_num
    Atom_list[i]  = fragment_num
  }
  }
}

Atom_list_old<- Atom_list

fragments_no_zero = if(length(which(Fragment_list==0)!=0)) Fragment_list[-which(Fragment_list==0)]
num_fragments = length(fragments_no_zero)
Fragment_mass = c(rep(0,num_fragments))
Fragment_Size = c(rep(0,num_fragments))
for(m in 1:length(Fragment_mass)){
  fragment_check = fragments_no_zero[m]
for(i in 1:length(Atom_list_old)){
  if(Atom_list_old[i]==fragment_check){
    Atom_list[i] = m
    Fragment_mass[m] = Fragment_mass[m]+1
  }
}
  #assume Circle
  Fragment_Size[m]=2*sqrt((Fragment_mass[m]*Area)/pi)
}

Fragment_Ratio = (Fragment_mass/Atom_Count)*100

summarisedlammps<-data.frame("Fragment_ID"=c(seq(1,num_fragments,1)),"Fragment_mass"=Fragment_mass, "Fragment_Ratio"=Fragment_Ratio, "Fragment_Size" = Fragment_Size )

listrms[[length(listrms)+1]]<-summarisedlammps
}

df<-bind_rows(listrms) 

df_new <- data.frame("Size"=c("<0.125mm","0.125mm-0.25mm","0.25mm-0.5mm","0.5mm-1mm","1mm-2mm","2mm-4mm",">4mm"),"Num_Fragments"=c(rep(0,7)),"Percentage_Mass"=c(rep(0,7)),"Test"=c(rep("SPOOFF",7)),"Order"=c(seq(1,7,1)))

for(i in 1:length(df$Fragment_ID)){
  if(df$Fragment_Size[i]<0.125){
    j=1
  }else if(df$Fragment_Size[i]<0.25){
    j=2
  }else if(df$Fragment_Size[i]<0.5){
    j=3
  }else if(df$Fragment_Size[i]<1){
    j=4
  }else if(df$Fragment_Size[i]<2){
    j=5
  }else if(df$Fragment_Size[i]<4){
    j=6
  }else{
    j=7
  }
  
  df_new$Num_Fragments[j]=df_new$Num_Fragments[j]+1
  df_new$Percentage_Mass[j] = df_new$Percentage_Mass[j]+df$Fragment_Ratio[i]
}

Exp_Data<-read.csv(paste0(CorePath,"Fragment_Data.csv"))
Exp_Data<-as.data.frame(Exp_Data)
Size_names = c("<0.125mm","0.125mm-0.25mm","0.25mm-0.5mm","0.5mm-1mm","1mm-2mm","2mm-4mm",">4mm")
Exp_Data<-Exp_Data%>%
  mutate(Order = 0)
for (i in 1:length(Exp_Data$Size)){
  Exp_Data$Order[i]=match(Exp_Data$Size[i],Size_names)
}


df_all<-rbind(Exp_Data,subset(df_new,select=-c(Num_Fragments)))
Test_names = c("Studsvik Rod 191","Studsvik Rod 192","Studsvik Rod 193","Studsvik Rod 196","Studsvik Rod 198","SPOOF")
df_all<-df_all%>%
  mutate(Order2 = 0)
for (i in 1:length(df_all$Size)){
  df_all$Order2[i]=match(df_all$Test[i],Test_names)
}

fig<- ggplot() +
  geom_bar(data=df_all, mapping=aes(x=reorder(Size,Order), y=Percentage_Mass, fill=reorder(Test,Order2)),position="dodge",stat="identity") +
  xlab('Approxiamte Fragment Diameter mm') +
  ylab('Percentage of Total Mass %') 

fig

ggsave(filename = paste0("Fragment Size Data.png"),plot=fig,path=CorePath,height = 5, width = 10)


write.csv(df_all,file=paste0(CorePath,"Fragment Size Data.csv"),row.names = FALSE)







for(i in 1:length(dflammps$x)){
  dflammps$Fragment[i]=Atom_list[dflammps$id[i]]
}

dflammps2<-dflammps%>%
  filter(dflammps$Fragment!=0)

fig<- ggplot() +
  geom_point(data=dflammps, mapping=aes(x=x, y=y, colour=Fragment),size=2.5) +
  xlab('x mm') +
  ylab('y') +
  theme(legend.position = "none") +
  scale_color_stepsn(colours=colorRampPalette(c("blue", "green", "purple", "orange","yellow","black", "pink","lightblue", "lightgreen", "mediumpurple","gold","tan","grey", "maroon", "hotpink","cyan", "darkgreen", "orchid","goldenrod","orangered","grey45", "maroon3","blue", "green", "purple", "orange","yellow","black", "pink","lightblue", "lightgreen", "mediumpurple","gold","tan","grey", "maroon", "hotpink","cyan", "darkgreen", "orchid","goldenrod","orangered","grey45", "maroon3"))(num_fragments),n.breaks=num_fragments)

fig

ggsave(filename = paste0("Fragment Colour Map.png"),plot=fig,path=CorePath,height = 10, width = 10)

