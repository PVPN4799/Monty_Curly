library(compiler)
library(schoolmath)

################################################################

#forever constants....
print("Initialising constants....")
kb = 1.38064852 * (10^-23)
T = 300
H = 5*(10^(-19))
rad=5*(10^(-9))

#...not this one
max = rep_len(2,length.out = 10000)

################################################################

#initialisation of particle position in simulation box
#Size of particle: 10 nm

print("Creating simulation box...")
l = list(c())
count = 1
for (z in c(1,202)){
  for (i in seq(1,202,by=3)){
    for (j in seq(1,202,by=3)){
      l[[count]] = c(i,j,z)
      count = count+1
    }
  }}
for (i in seq(1,33,by=3)){
  for (j in seq(1,202,by=3)){
    l[[count]] = c(i,j,101)
    count = count+1
  }
}
for (j in seq(1,12,by=3)){
  l[[count]] = c(34,j,101)
  count = count+1
}

################################################################
print("Building functions...")

#acceptance probability
P <- function(dE) {
  #print("Calculating probabilities...")
  prob = min(1,(exp(-dE/(kb*T))))
  return(prob)
}

#individual particle position update
update_pos <- function(oldr,maxr) {
  newr = c()
  oldr = unlist(oldr)
  ran = c(runif(1),runif(1),runif(1))
  newr = oldr + (((2*ran)-1)*maxr)
  return(list(newr))
}

#distance from origin
dist <- function(p){
  p = unlist(p)
  d = ((p[1]^2)+(p[2]^2)+(p[3]^2))^(0.5)
  return(d*5*(10^(-9)))
}

#relative displacement of individual particle rest of system particles
disp <- function(j,st){
  gg = lapply(st,function(k,j) j - unlist(k),unlist(j))
  d = sapply(gg, dist) - (2*rad)
  d[which(d<0)]=0
  return(d)
}
# detour_disp <- function(j,k) {
#   m = max(abs((j-unlist(k))-2*rad),0)
#   return(m)
# }

#pairwise Van der Waal energy
part_en <- function(d) {
  e = -(H*rad)/(12*d)
  if(e == -Inf) return(0)
  else
  return(e)
}

#max allowed displacement update
update_max <- function(max,p){
  if (p - 0.5 > 0.01)
    return(max*1.05)
  else if (p - 0.5 < (-0.01))
    return(max*0.95)
  else
    return(max)
}

#update system configuration for 1000 time steps
time_step <- function(st,l) {
  st[[1]] = l[[1]]
  print("Initial config...")
  print("Calculating Euclidean distance...")
  dist_l = lapply(l[[1]], dist)
  print("Calculating displacements...")
  part_disp = lapply(l[[1]],disp,l[[1]])
  print("Calculating VDW energies...")
  part_en_list = c()
  step_part_en = list()
  for (k in 1:10000){
    part_en_list=c(part_en_list,sum(sapply(part_disp[[k]],part_en)))
  }
  step_part_en[[1]] = part_en_list
  
  print("100 time steps begin...")
  for (i in 1:100){
    print(i)
    #o  = st[[length(st)]]
    print("Updating positions...")
    nl = Map(update_pos,st[[length(st)]],max)
    st[[length(st)+1]] = nl
    print("Distances...")
    dist_l = lapply(st[[length(st)]], dist)
    print("Displacements...")
    part_disp = lapply(nl,disp,nl)
    part_en_list = c()
    print("Energies...")
    for (k in 1:10000){
      part_en_list=c(part_en_list,sum(sapply(part_disp[[k]],part_en)))
    }
    step_part_en[[2]] = part_en_list
    print("Comparing energies...")
    dE = step_part_en[[2]]-step_part_en[[1]]
    print("Probabilities...")
    part_P = lapply(dE,P)
    
    #update max allowed disp
    print("Constraints...")
    mx = c()
    for (k in 1:10000)
      mx[k] = update_max(max[k],part_P[[k]])
    max <<- mx
    step_part_en[[1]] = part_en_list
    step_part_en[[2]] = c()
  }
  return (st)
}

#make time step function go fast
g = cmpfun(time_step)

##################################################################

#10^3 time step updates
st = list()
o = l
for (iter in 1:10){
sprintf("%d timesteps begin",i*100)
big_iter = 1
st = g(st,o)
#saving mechanism
lo = paste0("C:/Users/buser/Desktop/temp/BigStep_",big_iter)
if(!dir.exists(lo)){
  dir.create(lo)
}
flo = paste0(lo,"/SmallStep_",iter,".RData")
st_lite = st[1:100]
o = st[101]
save(st_lite,file=flo)
steps = 100*iter
sprintf("%d timesteps end",steps)
rm(st)
gc()
st = list()
}

tot_en = sum(part_en_list)