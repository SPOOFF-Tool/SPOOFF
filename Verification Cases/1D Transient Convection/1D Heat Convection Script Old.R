Path <- "C:/Users/mhort/OneDrive/Documents/Lammps/lammps-stable/lammps-2Aug2023/Testing/Heat Source Test/dump.LAMMPS"
library(ggplot2)
library(dplyr)
library(data.table)

dtmod = 0.1
mmax = 100
kappa = 1
spacing = 1
rho = 1
xmin = -0
xmax = 100
cp = 1
x = seq(xmin,xmax,spacing)
L = xmax-xmin
endtime = 100
t = c(0,1,1 %o% 10^(1:log10(endtime)))
h = 1


Tget<-function(x,t,kappa,mmax,L,rho,cp,h){
  Tfinal=0
  DT = kappa/(rho*cp)
  for(m in 1:mmax){
    lambda=rootget(m,L,DT,h)
    an=aget(lambda,L)
    #print(paste0("n; ",m," lambda: ",lambda," a: ",an))
    Tfinal= Tfinal+an*exp(-1*(lambda^2)*DT*t)*sin(lambda*x)
  }
  return(Tfinal)
}
Qget<-function(x,t,kappa,mmax,L,rho,cp,h){
  Qfinal=0
  DT = kappa/(rho*cp)
  for(m in 1:mmax){
    lambda=rootget(m,L,DT,h)
    an=aget(lambda,L)
    #print(paste0("n; ",m," lambda: ",lambda," a: ",an))
    Qfinal= Qfinal-kappa*lambda*an*exp(-1*(lambda^2)*DT*t)*cos(lambda*x)
  }
  return(Qfinal)
}


Itterate<-function(xm,L,kappa,h){
x = xm-((tan(xm*L)+kappa*xm/h)/(L/(cos(xm*L)^2) + kappa/h))  
}
rootget<-function(n,L,kappa,h){
  epsilon = 1E-8
  tolerance = 1E-15
  start = (((2*n)-1)*pi/(2*L)) + epsilon
  end = (((2*n)+1)*pi/(2*L)) - epsilon
  xm = start
  x=Itterate(xm,L,kappa,h)
  count = 1
  error = abs(x-xm)
  while(error>tolerance){
    if(x<start|x>end){
      xm = start + (end-start)/count
      count = count +1
      x=Itterate(xm,L,kappa,h)
      error = abs(x-xm)
      xm=x
    }else{
      x=Itterate(xm,L,kappa,h)
      error = abs(x-xm)
      xm=x
    }
  }
  return(x)
}
aget<-function(lambda,L){
  an = (2*(sin(lambda*L)-lambda*L*cos(lambda*L)))/(lambda*(lambda*L-sin(lambda*L)*cos(lambda*L)))
  return(an)
}



TimeArray <- list()
count=0
for (i in t){
  Tarray =c()
  Qarray =c()
  for(j in x){
    Tarray[length(Tarray)+1]<- Tget(j,i,kappa,mmax,L,rho,cp,h)
    Qarray[length(Qarray)+1]<- Qget(j,i,kappa,mmax,L,rho,cp,h)
  }
  d = data.frame("x"=x, "T" = Tarray, "Q" = Qarray, "t"=rep(i,length(x)))
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
    line = readLines(con, n = atoms)
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

dflammps <- altprocessFile(Path,t,dtmod,10)

fig<- ggplot() +
  geom_line(data=df, aes(x=x, y=T, colour = t, linetype = t)) +
  geom_point(data=dflammps, aes(x=x, y=c_t, colour = t))

fig

fig<- ggplot() +
  geom_line(data=df, aes(x=x, y=Q,colour = t, linetype = t)) +
  geom_point(data=dflammps, aes(x=x, y=c_Q1, colour = t)) +
  ylim(0,18)

fig

#fig2<- ggplot() +
#  geom_point(data=dflammps, aes(x=x, y="vx", colour = t, linetype = t))

#fig2