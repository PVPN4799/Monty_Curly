#plotting mechanism
library(rgl)
i = 1
j = 100
o = unlist(st_lite[j])
par3d(windowRect=c(1000,1000,60000,6000))
rgl.viewpoint(fov = 30, zoom = 0.75)
print(plot3d(o[1],o[2],o[3], radius = 1))
for (i in seq(4,29998, by=3)){
  print(spheres3d(o[i],o[i+2],o[i+1], color = "blue"))
}
bbox3d(color = "green")
aspect3d(1,1,1)
filename = "300_steps.png"
rgl.snapshot(filename)