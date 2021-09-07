library(tidyverse)

### ========= Box decomposition for strain calculation ==========
L=12  # number of boxes in each direction
a=4.0626

box_atoms <- list() # keeps atom id.s included to the box[[i]]
xbox <- 0
ybox <- 0
zbox <- 0

x <- 0
y <- 0
z <- 0
length(seq(L^3))
for (i in seq(L^3)) { # number of total boxes = L^3
  points_in_L=L+1 # number of atoms for a given 1D length along the NP box
  which_box_layer=i%/%(L+1e-5)
  # print(which_box_layer)
  which_base_block=i%/%(L^2+1e-5)
  
  base_nodes <- c(i+which_box_layer+which_base_block*points_in_L,
                  i+which_box_layer+which_base_block*points_in_L+1,
                  i+which_box_layer+which_base_block*points_in_L+points_in_L,
                  i+which_box_layer+which_base_block*points_in_L+points_in_L+1)
  
  top_nodes <- c(i+which_box_layer+which_base_block*points_in_L+points_in_L^2,
                 i+which_box_layer+which_base_block*points_in_L+points_in_L^2+1,
                 i+which_box_layer+which_base_block*points_in_L+points_in_L+points_in_L^2,
                 i+which_box_layer+which_base_block*points_in_L+points_in_L+points_in_L^2+1)
  
  box_atoms[[i]] <- c(base_nodes,top_nodes)
  
  if (i<=(points_in_L)^2) {
    which_box_layer_node=i%/%(points_in_L+1e-5)
    x[i] <- -L/2+(i-1-which_box_layer_node*points_in_L)
    y[i] <- -L/2+which_box_layer_node
    z[i] <- -L/2
  }
  
  if (i<=(L)^2) {
    which_box_layer_box=i%/%(L+1e-5)
    xbox[i] <- -L/2+(i-1-which_box_layer_box*L)+0.5
    ybox[i] <- -L/2+which_box_layer_box+0.5
    zbox[i] <- -L/2+0.5
  }
}
length(box_atoms)

x_all <- as.vector(replicate(L+1,x))
y_all <- as.vector(replicate(L+1,y))
z_all <- z


for (i in seq(L)) {
  z_all <- c(z_all,z+i)
}
x_all <- x_all*a
y_all <- y_all*a
z_all <- z_all*a

x_box_all <- as.vector(replicate(L,xbox))
y_box_all <- as.vector(replicate(L,ybox))
z_box_all <- zbox

for (i in 1:(L-1)) {
  z_box_all <- c(z_box_all,zbox+i)
}
x_box_all <- x_box_all*a
y_box_all <- y_box_all*a
z_box_all <- z_box_all*a

# library(rgl)
rgl::plot3d(x_box_all,y_box_all,z_box_all)
rgl::plot3d(x_all,y_all,z_all)

library(tidyverse)
read_delim("../Downloads/5.xyz"," ",col_names = F,skip = 2) %>% 
  select(2:4)%>% rename(X1=1,X2=2,X3=3)->dt
read_delim("../Downloads/min5.xyz"," ",col_names = F,skip = 2) %>% 
  select(2:4)%>% rename(X1=1,X2=2,X3=3)->dt_min

x_np <- dt$X1
y_np <- dt$X2
z_np <- dt$X3

x_np_min <- dt_min$X1
y_np_min <- dt_min$X2
z_np_min <- dt_min$X3

count=0
x_all_min <- 0
y_all_min <- 0
z_all_min <- 0

for (i in seq(length(x_all))) {
  does_match=FALSE
  any(x_np%in%x_all[i] & y_np%in%y_all[i] & z_np%in%z_all[i]) -> does_match
  
  if(!does_match){
    x_all[i] <- NA
    y_all[i] <- NA
    z_all[i] <- NA
    
    x_all_min[i] <- NA
    y_all_min[i] <- NA
    z_all_min[i] <- NA
  
    } else {
    which_match <- which(x_np%in%x_all[i] & y_np%in%y_all[i] & z_np%in%z_all[i])
    x_all[i] <-x_all[i]
    y_all[i] <-y_all[i]
    z_all[i] <-z_all[i]
    
    x_all_min[i] <- x_np_min[which_match]
    y_all_min[i] <- y_np_min[which_match]
    z_all_min[i] <- z_np_min[which_match]
  }
}

 # tibble(x=x_all_min,y=y_all_min,z=z_all_min) %>% drop_na()->np_match
 # plot3d(np_match$x,np_match$y,np_match$z,size = 6)

displ_gradient=list()

for (i in seq(length(box_atoms))) {
  # x=u,y=v,z=w
  node1=box_atoms[[i]][1]
  node2=box_atoms[[i]][2]
  node3=box_atoms[[i]][3]
  node4=box_atoms[[i]][4]
  node5=box_atoms[[i]][5]
  node6=box_atoms[[i]][6]
  node7=box_atoms[[i]][7]
  node8=box_atoms[[i]][8]
  
  exx1=abs(x_all_min[node1]-x_all_min[node2])
  exx2=abs(x_all_min[node3]-x_all_min[node4])
  exx3=abs(x_all_min[node5]-x_all_min[node6])
  exx4=abs(x_all_min[node7]-x_all_min[node8])
  exx_av=mean(c(exx1,exx2,exx3,exx4),na.rm = F)
  du=(exx_av-a)

  eyy1=abs(y_all_min[node1]-y_all_min[node3])
  eyy2=abs(y_all_min[node2]-y_all_min[node4])
  eyy3=abs(y_all_min[node5]-y_all_min[node7])
  eyy4=abs(y_all_min[node6]-y_all_min[node8])
  eyy_av=mean(c(eyy1,eyy2,eyy3,eyy4),na.rm = F)
  dv=(eyy_av-a)
  
  ezz1=abs(z_all_min[node1]-z_all_min[node5])
  ezz2=abs(z_all_min[node2]-z_all_min[node6])
  ezz3=abs(z_all_min[node3]-z_all_min[node7])
  ezz4=abs(z_all_min[node4]-z_all_min[node8])
  ezz_av=mean(c(ezz1,ezz2,ezz3,ezz4),na.rm = F)
  dw=(ezz_av-a)
  
  displ_gradient[[i]] <- c(du/a,du/a,du/a,dv/a,dv/a,dv/a,dw/a,dw/a,dw/a)
}


exx <- NULL
eyy <- NULL
ezz <- NULL
x_box_clean <- NULL
y_box_clean <- NULL
z_box_clean <- NULL

# for (i in 1:length(displ_gradient)) {
#     c(exx,displ_gradient[[i]][1])->exx
#     c(eyy,displ_gradient[[i]][5])->eyy
#     c(ezz,displ_gradient[[i]][9])->ezz
# }
# 
# c(exx,eyy,ezz)->af
# as.vector(na.omit(af))->af2

for (i in 1:length(displ_gradient)) {
  any(is.na(displ_gradient[[i]]))->yes_na
  if (!yes_na) {
    # print("dd")
    c(exx,displ_gradient[[i]][1])->exx
    c(eyy,displ_gradient[[i]][5])->eyy
    c(ezz,displ_gradient[[i]][9])->ezz
    c(x_box_clean,x_box_all[i])->x_box_clean
    c(y_box_clean,y_box_all[i])->y_box_clean
    c(z_box_clean,z_box_all[i])->z_box_clean
  } 
}

c(exx,eyy,ezz)->m1

# library(rgl)
tibble(x_box_clean,y_box_clean,z_box_clean,exx,eyy,ezz) %>% drop_na()->aaa
aaa  %>% mutate(el="Au") %>% select(el,1:6) %>%  write_delim("../Desktop/ssss_corner.xyz")

rgl::plot3d(x_box_clean,y_box_clean,z_box_clean)

nrow(aaa)


mean(c(exx,eyy,ezz),na.rm = T)

hist(exx)

### Finalize strain from all atoms (included 1st neighbors also)
read_delim("~/../Desktop/ssss_corner.xyz"," ",col_names = F,skip = 2)->sa

sa$X5->exs
sa$X6->eys
sa$X7->ezs

mean(c(exs,eys,ezs))

############################### New Approach #################################
library(tidyverse)
{
read_delim("D:/research/2232/np rapor/lammps data/5.xyz"," ",col_names = F,skip = 2) %>% select(-1)->ref
read_delim("D:/research/2232/np rapor/lammps data/min5.xyz"," ",col_names = F,skip = 2) %>% select(-1)->myminexact

myminexact %>% 
  rename(x=1,y=2,z=3) %>% 
  mutate(CN=0,n_list=list(0),id=1:nrow(.),def_grad=0,strain=0,lat_strain=0,lat_strain_magnitude=0,strain_inf=0,
         rotation=0,displ_gr=0) -> merdanmin
criteria_neighbor <- 3.2 #enough distance to capture all neighbors
xxd <- merdanmin[[1]]
yyd <- merdanmin[[2]]
zzd <- merdanmin[[3]]
merdann <-as.list(merdanmin)

system.time(
for (i in 1:nrow(merdanmin)) {
  center_atom_x <- xxd[i] #merdanmin$x[i]
  center_atom_y <- yyd[i] #merdanmin$y[i]
  center_atom_z <- zzd[i] #merdanmin$z[i]
  merdanmin[merdann[[1]]>(center_atom_x-criteria_neighbor) & merdann[[1]]<(center_atom_x+criteria_neighbor) &
              merdann[[2]]>(center_atom_y-criteria_neighbor) & merdann[[2]]<(center_atom_y+criteria_neighbor) &
              merdann[[3]]>(center_atom_z-criteria_neighbor) & merdann[[3]]<(center_atom_z+criteria_neighbor) &
              merdann[[6]]!=i,]->filtered
  nrow(filtered) -> merdann[[4]][i] #do -1 finally, because of self counting
  list(filtered$id)->merdann[[5]][i] #$n_list[i]

}
)
merdanmin[[4]] <- merdann[[4]]
merdanmin[[5]] <- merdann[[5]]

merdan2 <- as.list(merdanmin)
ref2 <- as.list(ref)
xxr <- ref2[[1]]
yyr <- ref2[[2]]
zzr <- ref2[[3]]


system.time(
  for (i in seq(nrow(merdanmin))) {
    n0 <- NULL # ref neighbor vectors
    n <- NULL # deformed neighbor vectors
    n_list_i <- merdan2[[5]][[i]]
    for (j in seq(length(n_list_i))) {
      neighbor_id <- n_list_i[j]
      x <- xxd[neighbor_id]-xxd[i] 
      y <- yyd[neighbor_id]-yyd[i] 
      z <- zzd[neighbor_id]-zzd[i]  
      
      x0 <- xxr[neighbor_id]-xxr[i] 
      y0 <- yyr[neighbor_id]-yyr[i] 
      z0 <- zzr[neighbor_id]-zzr[i] 
      
      n <- c(n,c(x,y,z))
      n0 <- c(n0,c(x0,y0,z0))
    }
    x_m <- matrix(n0,nrow = length(n_list_i),byrow = T)
    y_m <- matrix(n,nrow = length(n_list_i),byrow = T)
    grad_i <- lsfit(x_m,y_m,intercept = F)$coefficients
    merdan2[[7]][i] <- list(grad_i)
    merdan2[[8]][i] <-  list((grad_i%*%t(grad_i)-diag(3))/2) 
    merdan2[[9]][i] <-  list(eigen(merdan2[[8]][[i]])$values) 
    merdan2[[10]][i] <- sqrt(mean((eigen(merdan2[[8]][[i]])$values)^2))
    merdan2[[13]][i] <-  list(grad_i-diag(3)) 
    merdan2[[11]][i] <-  list((merdan2[[13]][[i]]+t(merdan2[[13]][[i]]))/2) 
    merdan2[[12]][i] <-  list((merdan2[[13]][[i]]-t(merdan2[[13]][[i]]))/2)
  }
)
merdanmin[[7]] <- merdan2[[7]]
merdanmin[[8]] <- merdan2[[8]]
merdanmin[[9]] <- merdan2[[9]]
merdanmin[[10]] <- merdan2[[10]]
merdanmin[[11]] <- merdan2[[11]]
merdanmin[[12]] <- merdan2[[12]]
merdanmin[[13]] <- merdan2[[13]]


###########################
surf <- merdanmin[[6]][merdanmin[[4]]<12]
core <- merdanmin[[6]][merdanmin[[4]]==12]
distance <- 0
xx <- merdanmin[[1]]
yy <- merdanmin[[2]]
zz <- merdanmin[[3]]
surfxx <- xx[surf]
surfyy <- yy[surf]
surfzz <- zz[surf]

system.time(
for (i in seq_along(core)) {
  xi <- xx[core[i]]
  yi <- yy[core[i]]
  zi <- zz[core[i]]
  dist_i <- sqrt((surfxx-xi)^2+(surfyy-yi)^2+(surfzz-zi)^2)
  distance[i] <- min(dist_i)
}
)
}
tibble(core,distance) %>% arrange(distance)->dt
##########
merdanmin %>% filter(CN==12)->merdanmin2

mean(unlist(merdanmin2[[8]]))*1e2

unlist(merdanmin2[[9]])*100->st
st[st>-2 & st<2]->stf

hist(stf,breaks = 1000)


mean(merdanmin2[[10]]*1e6)

sqrt(mean(unlist(merdanmin2[[9]])^2))*1e6

tibble(unlist(merdanmin[[9]])) #%>% 
  # write_delim("../Desktop/1212.dat")

# 
dt[["strain"]] <- merdanmin[dt[[1]],][[10]]
# dt[1:10,] %>% ggplot(aes(distance,strain))+
#   geom_col(fill="red")

dt %>% select(2,3) %>% group_by(distance) %>% summarize(strain=mean(strain)) %>% 
  write_delim("../Desktop/strainmin5.dat",col_names = F)

dt %>% mutate(strain=strain*1e6) %>% 
  filter(strain<10) %>% arrange(distance)

hist(distance,plot = F,breaks = 500)->mids
dt[["mids"]] <- rep(mids$mids,mids$counts)
mids$mids->midsall
dt
dt %>% count(mids) %>% View(.)
lat_mids <- 0

system.time(
for (i in seq_along(mids$mids)) {
  dt %>% mutate(mids=as.double(mids)) %>%  filter(mids==midsall[i])->dti
  merdanmin[dti[[1]],]->merdan_i
  lat_strain <- unlist(merdan_i$lat_strain)
  lat_i <- sqrt(mean(lat_strain^2))*1e6
  lat_mids[i] <- lat_i
}
)

tibble(midsall,lat_mids) %>% ggplot(aes(midsall,lat_mids))+geom_col(fill="skyblue",color="red")+
  # geom_text(aes(label = round(lat_mids,0)), vjust = -0.5,check_overlap = T)+
  theme_bw()+
  theme(axis.text = element_text(color="black",size=12))+
  labs(x="distance (A)",y="rms strain")+
  geom_hline(yintercept = 300)
  # ggsave("../Desktop/10.png")




###########################

merdanmin %>% filter(CN==12)->merdanmin2
nrow(merdanmin2)

lat_strain=unlist(merdanmin2$lat_strain)
hist(lat_strain,breaks = nrow(merdanmin2))

mean(lat_strain)*4.0626+4.0626
sqrt(mean(lat_strain^2))*1e6

merdanmin %>% mutate(type="Au",filt=if_else(CN>=9,1,0)) %>% 
  select(14,1:3,15)->tst

library(plotly)

plot_ly(tst %>% filter(z<=0),x=~x,y=~y,z=~z,size = 100,color = ~filt) %>% add_markers()

tst %>% mutate(z=round(z,3))%>% filter(z<=0)->tst2

pth="../Desktop/c_section.xyz"
write_lines(nrow(tst2),pth)
write_lines(" ",pth,append = T)
write_delim(tst2,pth,append = T,col_names=F)

aaa$mids->br
aaa$counts->cn
tibble(br,cn) %>% ggplot(aes(br,cn))+geom_col(fill="gold",color="blue")
+scale_y_log10()+scale_x_continuous(n.breaks = 10)+theme_bw()+
  labs(x="strain",y="counts",title="method 2: continuum mechanics",subtitle = "r.m.s: 0.00663;  mean=-0.0021")+
  theme(axis.text = element_text(color="black"))+
  geom_vline(xintercept = -0.002077,linetype="dashed")+ggsave("~/../Desktop/met2.jpg")

mean(m2)
sqrt(mean(m2^2))

#### m1
hist(m1)->aaa

aaa$mids->br
aaa$counts->cn
tibble(br,cn) %>% ggplot(aes(br,cn))+geom_col(fill="gold",color="blue")+scale_y_log10()+scale_x_continuous(n.breaks = 10)+theme_bw()+
  labs(x="strain",y="counts",title="method 1: direct box decomposition",subtitle = "r.m.s: 0.00693;  mean=-0.0025")+
  theme(axis.text = element_text(color="black"))+geom_vline(xintercept =  -0.002459735,linetype="dashed")+ggsave("~/../Desktop/met1.jpg")
mean(m1)
sqrt(mean(m1^2))

################ method 2 #############
library(tidyverse)
pthxx="D:/research/2232/strain/strain data/energy minimized/shengmin/30/epsi_xx.dat"
pthyy="D:/research/2232/strain/strain data/energy minimized/shengmin/30/epsi_yy.dat"
pthzz="D:/research/2232/strain/strain data/energy minimized/shengmin/30/epsi_zz.dat"
pthxy="D:/research/2232/strain/strain data/energy minimized/shengmin/30/epsi_xy.dat"
pthyz="D:/research/2232/strain/strain data/energy minimized/shengmin/30/epsi_yz.dat"
pthzx="D:/research/2232/strain/strain data/energy minimized/shengmin/30/epsi_zx.dat"

read_delim(pthxx," ",col_names = F)->exx
read_delim(pthyy," ",col_names = F)->eyy
read_delim(pthzz," ",col_names = F)->ezz
read_delim(pthxy," ",col_names = F)->exy
read_delim(pthyz," ",col_names = F)->eyz
read_delim(pthzx," ",col_names = F)->ezx

as.vector(as_vector(exx))->exx
as.vector(as_vector(eyy))->eyy
as.vector(as_vector(ezz))->ezz
as.vector(as_vector(exy))->exy
as.vector(as_vector(eyz))->eyz
as.vector(as_vector(ezx))->ezx



mean(c(exx,eyy,ezz))*4.0626+4.0626
data.frame(exx,exy,ezx,eyy,eyz,ezz)->dt
lall=as.list(dt)
len=length(exx)
pr_strains=list(rep(0,len))
tracemem(lall)
tracemem(pr_strains)
for (i in seq(len)) {
  strn<-c(lall[[1]][i],lall[[2]][i],lall[[3]][i],lall[[2]][i],lall[[4]][i],lall[[5]][i],lall[[3]][i],lall[[5]][i],lall[[6]][i])
  matrix(strn,ncol=3,byrow = T)->strn
  pr_strains[[i]] <-eigen(strn)$values
}

strs <- unlist(pr_strains)
len
mean(strs)*4.0626+4.0626
sqrt(mean(strs^2))

sqrt(mean(c(exx,eyy,ezz)^2))

size=c(5,10,15,20,30)
method1<-c(4.0542,
4.0588,
4.0608,
4.061817,
4.0624)
tibble(method="method1-CN=12",lat=method1,size)->met1

sizex=c(5,10,15,20,30)
method1_2<-c(4.0503,
           4.0568,
           4.0595,
           4.0608,
           4.0617)
tibble(method="method1-CN>=9",lat=method1_2,size=sizex)->met1_2


method2<-c(
  4.0496,
  4.0567,
  4.0594,
  4.0608,
  4.0617
)

tibble(method="method2",lat=method2,size)->met2
x_ray <- c(
  4.0533,
  4.0586,
  4.0601
)
tibble(method="rietveld",lat=x_ray,size=size[1:3])->xray

strain_method1 <- c(
  0.01096502,
  0.007548377,
  0.006041681,
  0.00509711,
  0.003817262
)
tibble(method="method 1",strain=strain_method1*1e6,size)->sm1

strain_method2 <- c(
  0.01085007,
  0.00749452,
  0.005925029,
  0.004983491,
  0.003846752
)

tibble(method="method 2",strain=strain_method2*1e6,size)->sm2

# strain_method21 <- c(
#   0.009062467,
#   0.006139527,
#   0.004762605
# )
# tibble(method="method 2 - normal",strain=strain_method21*1e6,size)->sm3

rietveld_strain <- c(
  5094.55042,
  1308.96773,
  170.32151
)
tibble(method="Rietveld",strain=rietveld_strain,size=size[1:3])->sm4

rbind(met1,met1_2,met2,xray) %>% ggplot(aes(size,lat,color=method,shape=method))+
  geom_point(size=3)+geom_line(size=1)+
  theme_bw()+
  scale_y_continuous(limits = c(4.047,4.063),n.breaks = 8)+
  labs(x="particle size (nm)",y="a (A)")+
  theme(axis.text = element_text(color="black",size=12),axis.title = element_text(color="black",size=12),
        legend.position = c(0.7,0.4))+
  geom_hline(yintercept = 4.0626,linetype="dashed",size=1)+
  scale_x_continuous(n.breaks = 6)+
  ggsave("../Desktop/lat2.png")
  
rbind(sm1,sm2,sm4) %>% ggplot(aes(size,strain,color=method,shape=method))+
  geom_point(size=4)+
  geom_line(size=1.4)+
  theme_bw()+
  theme(axis.text = element_text(color="black",size=12),legend.position = c(0.7,0.7),axis.title = element_text(color="black",size=12))+
  labs(x="Particle size (nm)", y="r.m.s strain (microstrain)")+
  scale_shape_manual(values = c(18:15))+
  scale_y_continuous(limits = c(0,13e3),n.breaks = 8)+
  scale_x_continuous(n.breaks = 6)+
  ggsave("../Desktop/str.png")


pth="../Pictures/tst/grap.data"
library(tidyverse)
read_delim(pth, " ",col_names=F,skip=16) %>% select(-6)->dt

dt %>% summarise_all(max)



############## tensile simulation---> neat graphene data for lammps ###################
size=c("1_5","2_5","5","5_10","10","15","20")
pth="../Pictures/tst/data/"

for (i in seq_along(size)) {
  pth_in=paste0(pth,size[i],".xyz")
  pth_out=paste0(pth,size[i],"_clean",".data")
  
  read_delim(pth_in," ",col_names=F,skip=2,trim_ws = T) %>% 
    mutate(id=1:nrow(.),mol=1) %>% 
    select(5,6,2:4) ->dt
  
  dt %>% summarize_all(min) %>% .[[3]] - 0.7 -> xlo
  dt %>% summarize_all(min) %>% .[[4]] - 0.7 -> ylo
  dt %>% summarize_all(min) %>% .[[5]] - 10 -> zlo
  
  dt %>% summarize_all(max) %>% .[[3]] + 0.7 -> xhi
  dt %>% summarize_all(max) %>% .[[4]] + 0.7 -> yhi
  dt %>% summarize_all(max) %>% .[[5]] + 10 -> zhi
  
  write_lines("comment line....",pth_out)
  write_lines(paste(nrow(dt),"atoms",sep = " "),pth_out,append = T)
  write_lines("1 atom types",pth_out,append = T)
  write_lines(paste(xlo,xhi,"xlo xhi",sep = " "),pth_out,append = T)
  write_lines(paste(ylo,yhi,"ylo yhi",sep = " "),pth_out,append = T)
  write_lines(paste(zlo,zhi,"zlo zhi",sep = " "),pth_out,append = T)
  write_lines("       ",pth_out,append = T)
  write_lines("Masses",pth_out,append = T)
  write_lines("       ",pth_out,append = T)
  write_lines("1 12.0107",pth_out,append = T)
  write_lines("       ",pth_out,append = T)
  write_lines("Atoms",pth_out,append = T)
  write_lines("       ",pth_out,append = T)
  dt %>% write_delim(pth_out,append = T,col_names = F)
  
}







