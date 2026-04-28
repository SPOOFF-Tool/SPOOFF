Path <- "C:/Users/mhort/OneDrive/Documents/Lammps/lammps-stable/lammps-2Aug2023/Testing/Heat Source Test/dump.LAMMPS"
library(ggplot2)
library(dplyr)

dtmod = 0.1
T0 = 1
mmax = 1000
kappa = 0.1
spacing = 0.1
rho = 2
xmin = -48
xmax = 48
cp = 2
x = seq(xmin,xmax,spacing)
L = xmax-xmin
endtime = 10000
t = c(0,1 %o% 10^(1:log10(endtime)))
HP=100
qdot=HP*rho

Tget<-function(x,t,rho,kappa,mmax,L,cp,qdot){
  Tfinal=0
  for(m in 1:mmax){
    Tfinal= Tfinal +((4*L*L*((-1)^(m)))/((pi^3)*(((2*m)-1)^3))) * cos(((2*m)-1)*pi*x/L) * exp((-1*(((2*m)-1)^2)*(pi^2)*(kappa/(rho*cp))*t)/(L^2))
  }
  Tfinal= Tfinal+ L*L/8 - x*x/2
  Tfinal=Tfinal*qdot/kappa
  return(Tfinal)
}
TimeArray <- list()
count=0
for (i in t){
  Tarray =c()
  for(j in x){
    Tarray[length(Tarray)+1]<- Tget(j,i,rho,kappa,mmax,L,cp,qdot)
  }
  d = data.frame("x"=x, "T" = Tarray, "t"=rep(i,length(x)))
  count=count+1
  TimeArray[[count]]<-d
}

df<-bind_rows(TimeArray) %>%
  mutate(t=as.character(t)) %>%
  group_by(t)




processFile = function(Path,t,dtmod) {
  con = file(Path, "r")
  linenum = 0
  listdata=list()
  while ( TRUE ) {
    #check for end of file
    line = readLines(con, n = 1)
    linenum = linenum + 1
    if ( length(line) == 0 ) {
      break
    }
    #read time step
    timestep = dtmod*as.numeric(readLines(con, n = 1))
    print(timestep)
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
    if(timestep %in% t){
      data <- read.table(Path, skip = linenum, sep = ' ', nrows=atoms, as.is = TRUE, col.names=headers[[1]])
      data <- data %>%
        mutate(t=as.character(timestep))
      listdata[[length(listdata)+1]]<-data
    }
    #skip through atom lines
    readLines(con, n = atoms)
    linenum = linenum + atoms
  }
  
  close(con)
  
  dflammps<-bind_rows(listdata) %>%
    group_by(t)
  
  return(dflammps)
}

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
  data <- read.table(Path, skip = linenum, sep = ' ', nrows=atoms, as.is = TRUE, col.names=headers[[1]])
  data <- data %>%
    mutate(t=as.character(timestep))
  listdata[[length(listdata)+1]]<-data
  #skip through atom lines
  line = readLines(con, n = atoms)
  prelinenum = linenum
  linenum = linenum + atoms
  close(con)
  
  for (i in 2:length(t)){
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

#dflammps <- altprocessFile(Path,t,dtmod,10)

fig<- ggplot() +
  geom_line(data=df, aes(x=x, y=T, colour = t, linetype = t)) #+
  #geom_point(data=dflammps, aes(x=x, y=c_t, colour = t))

fig

#fig2<- ggplot() +
#  geom_point(data=dflammps, aes(x=x, y="vx", colour = t, linetype = t))

#fig2