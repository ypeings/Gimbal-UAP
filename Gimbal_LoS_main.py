# This code examine potential flight paths for the Gimbal UAP filmed in January 2015 by Navy aviators 
# It plots the F-18 flight path, in function of its True Air Speed (TAS), banking angle (extracted from the video every 10 frames), and estimated wind speed/direction. 
# Different trajectory configurations can be tested for the Gimbal object, and the script outputs graphs about the speed, altitude, and other parameters of a given flight path.
# Y. Peings, 04/2023, prepared for AIAA 2023

import numpy 
import numpy as np
import math 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import interp1d

######## User Parameters #############

# Wind 
Wind_speed=120 # Knots
Wind_direction=-35 # Relative to F-18 initial heading : 0=headwind ; 180=tailwind 

# UAP Path
path_type="const_head" # const_head/const_alt

# if const_head
offset=15  # Offset angle relative to wind direction, in degrees
dist_init=8. # Initial distance from the F-18 (Nautical miles)

# if constant_alt
Alt_object=18700 # Feet

# Field of view ATTFLIR pod (to refine lines of sight from background cloud motion)
FOV=0.35

# Smoothing of the LOS using background cloud motion (True/False)
smth_LOS=True

##########################################

# Open input file/Read variables
data = numpy.loadtxt(open("Data_Banking-Az-nFov", "rb"), skiprows=1)
frame=data[:,0] # Frame number
bank=data[:,1]  # Plane bank in deg
Az=data[:,2] # Extracted Azimuth angle
nFOV=data[:,3]  # number of scanned FOV at different points of the video
IAS=data[:,4]  # Indicated Air Speed (Kts)

# Interpolate nFOV
not_nan = np.logical_not(np.isnan(nFOV))
indices = np.arange(len(nFOV))
interp = interp1d(indices[not_nan],nFOV[not_nan])
nFOVi=interp(indices)

################# F-18 Flight Path #############################

# F-18 True Air Speed (Knots)
TAS = []
for i in range(104):
    TAS.append(0)
for i in range(104):
    if IAS[i] == 238 :
       TAS[i] = 363
    if IAS[i] == 239 :
       TAS[i] = 364
    if IAS[i] == 240 :
       TAS[i] = 366
    if IAS[i] == 241 :
       TAS[i] = 367
    if IAS[i] == 242 :
       TAS[i] = 369

# Frame rate (s)
dt=(1/30)

# Turn Rate
turn_rate = []
for i in range(104):
    turn_rate.append(0)
for i in range(0,104):
  turn_rate[i]=math.tan(bank[i]*3.14159/180)*1091/TAS[i]

# Initialization of variables
# Coordinates of flight path, with wind
x = []
y = []
z = []

# Coordinates of flight path, without wind
x2 = []
y2 = []
z2 = []

for i in range(104):
    x.append(0)
    y.append(0)
    z.append(0)
    x2.append(0)
    y2.append(0)
    z2.append(0)
  
heading_angle = []
for i in range(104):
    heading_angle.append(0)
heading_angle[0] = 0.

# Heading angle (does not depend on wind, only ground track does)
for i in range(1,104):
  heading_angle[i]=heading_angle[i-1]+(turn_rate[i-1]*10*dt)

# X,Y coordinates at every time step (every 10 frames)
for i in range(1,104):
  x[i]=x[i-1]+((TAS[i]/3600)*10*dt*math.sin(heading_angle[i]*3.14159/180))-(Wind_speed/3600)*10*dt*math.sin(Wind_direction*3.14159/180)
  x2[i]=x2[i-1]+((TAS[i]/3600)*10*dt*math.sin(heading_angle[i]*3.14159/180))
  y[i]=y[i-1]+((TAS[i]/3600)*10*dt*math.cos(heading_angle[i]*3.14159/180))-(Wind_speed/3600)*10*dt*math.cos(Wind_direction*3.14159/180)
  y2[i]=y2[i-1]+((TAS[i]/3600)*10*dt*math.cos(heading_angle[i]*3.14159/180))

#  heading_angle[i]=(np.arctan((x2[i]-x2[i-1])/(y2[i]-y2[i-1])))*180/3.14159

# 2D Plot
fig = plt.figure(figsize = (6,6))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # main axes
ax.plot(x,y,color='black',linewidth=2, label="-35\N{DEGREE SIGN} wind")
ax.plot(x2,y2,color='black',linewidth=2,linestyle='dashed', label="No wind")
ax.set_xlabel('X (Nm)', fontsize = 12)
ax.set_ylabel('Y (Nm)', fontsize = 12)

ax.set_title('F-18 Flight Path - Top view - Effect of wind')
ax.set_xticks([-2,-1.5,-1,-0.5,0,0.5,1])
ax.set_xticklabels(['-2','-1.5','-1','-0.5','0','0.5','1'])
ax.set_yticks([0,0.5,1,1.5,2,2.5,3])
ax.set_yticklabels(['0','0.5','1','1.5','2','2.5','3'])

plt.legend(loc="upper right",fontsize=14)
plt.savefig('Figures/Fig_F18_path.pdf', format="pdf", bbox_inches="tight")  
plt.close()

################# Lines of Sight ###############################

# Smoothed Azimuths
Az_smth = []
for i in range(104):
    Az_smth.append(0)
window=5
for i in range(0,104):
     Az_smth[i]=np.mean(-Az[i:i+window])
Az_smth[103]=6.4

# Smooth lines of sight to match the clouds -> Correction term
LoS_smth = []
sumA=0.
sumB=0.
for i in range(104):
    LoS_smth.append(0)

if smth_LOS == True :
  for i in range(1,104):
      LoS_smth[i]=LoS_smth[i-1]+(-(nFOVi[i]-nFOVi[i-1])*FOV)-(heading_angle[i]-heading_angle[i-1]+Az_smth[i]-Az_smth[i-1])
      sumA=sumA+heading_angle[i]-heading_angle[i-1]+Az_smth[i]-Az_smth[i-1]
      sumB=sumB+heading_angle[i]-heading_angle[i-1]+Az_smth[i]-Az_smth[i-1]+LoS_smth[i]-LoS_smth[i-1]
  print("Scanned FOV without correction:",round(sumA/FOV,2))
  print("Scanned FOV with correction:",round(sumB/FOV,2))

# At each time step, define a reference point for the LoS
El=np.arange(-2.,-1.95,0.05/104)

# Print position of F-18 at the end  
print("Last F-18 coordinate point: (",round(x[103],2),"/",round(y[103],2),")")

# Create points to plot the lines of sight
Lx = []
Ly = []
Lz = []
Gx = []
Gy = []
Gz = []
dist = []
distG = []
distG2 = []
Px = []
Py = []

for i in range(104):
    Lx.append(0)
    Ly.append(0)
    Lz.append(0)
    Gx.append(0)
    Gy.append(0)
    Gz.append(0)
    dist.append(dist_init+(dist_init*0.25))
    distG.append(0)
    distG2.append(0)
    Px.append(0)
    Py.append(0)

# Lines of sight
for i in range(1,104):
  Lx[i]=x[i]+(dist[i]*math.cos(El[i]*3.14159/180)*math.sin((heading_angle[i]+Az_smth[i]+LoS_smth[i])*3.14159/180))
  Ly[i]=y[i]+(dist[i]*math.cos(El[i]*3.14159/180)*math.cos((heading_angle[i]+Az_smth[i]+LoS_smth[i])*3.14159/180))
  Lz[i]=z[i]+(dist[i]*math.sin(El[i]*3.14159/180))

################ Flight Path #######################################

# Constant altitude
if path_type == 'const_alt':
  Alt_conv=(Alt_object-25000)/6076.12 # Nm above the F-18
  dist_init=Alt_conv/math.sin(El[0]*3.14159/180)
  for i in range(0,104):
    tmp_dist=0.
    delta=10.
    for j in np.arange(dist_init,dist_init+10,0.001):  
      tmp_alt=abs(z[i]+j*math.sin(El[i]*3.14159/180)-Alt_conv)
      if tmp_alt < delta:
        tmp_dist=j
        delta=tmp_alt
    distG2[i]=tmp_dist  

# Constant heading
if path_type == 'const_head':
  head_object=Wind_direction+offset # Heading relative to wind
  distb=dist_init
  for i in range(1,104):
    tmp_angle=0.
    delta=360.
    for j in np.arange(distb-0.05,distb+0.05,0.0001):  
      Px=x[i]+(j*math.cos(El[i]*3.14159/180)*math.sin((heading_angle[i]+Az_smth[i]+LoS_smth[i])*3.14159/180))
      Py=y[i]+(j*math.cos(El[i]*3.14159/180)*math.cos((heading_angle[i]+Az_smth[i]+LoS_smth[i])*3.14159/180))

      Pxb=x[i-1]+(distb*math.cos(El[i-1]*3.14159/180)*math.sin((heading_angle[i-1]+Az_smth[i-1]+LoS_smth[i-1])*3.14159/180))
      Pyb=y[i-1]+(distb*math.cos(El[i-1]*3.14159/180)*math.cos((heading_angle[i-1]+Az_smth[i-1]+LoS_smth[i-1])*3.14159/180))

      tmp_angle=abs(((np.arctan((Px-Pxb)/(Py-Pyb)))*180/3.14159)-head_object)

      if tmp_angle < delta:
        tmp_dist=j
        delta=tmp_angle
    distG2[i]=tmp_dist  
    distb=tmp_dist

print("Start Distance:",distG2[1],"Nm    End Distance:",round(distG2[103],2),"Nm    ",int(((distG2[103]/distG2[1])*100)-100),"% change")

# Flight Path    
for i in range(1,104):
  Gx[i]=x[i]+(distG2[i]*math.cos(El[i]*3.14159/180)*math.sin((heading_angle[i]+Az_smth[i]+LoS_smth[i])*3.14159/180))
  Gy[i]=y[i]+(distG2[i]*math.cos(El[i]*3.14159/180)*math.cos((heading_angle[i]+Az_smth[i]+LoS_smth[i])*3.14159/180))
  Gz[i]=z[i]+(distG2[i]*math.sin(El[i]*3.14159/180))

# Plot
import matplotlib as mpl
c = np.arange(1,104)
norm = mpl.colors.Normalize(vmin=c.min(), vmax=c.max())
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.jet)
cmap.set_array([])

fig2 = plt.figure()
ax = plt.axes(projection='3d')
p=ax.plot(x, y, z, c='black',linewidth=5.5,zorder=10,linestyle='solid')
p=ax.scatter(x, y, z, c=frame, cmap='jet', linewidth=0.1);
#ax.scatter(Gx, Gy, Gz, c=frame, cmap='jet', linewidth=1.);
for i in range(0,104):
  ax.plot([x[i],Lx[i]], [y[i],Ly[i]], [z[i],Lz[i]],c=cmap.to_rgba(i))
  
# Flight Path curve
for i in range(2,104):
    ax.plot([Gx[i-1],Gx[i]], [Gy[i-1],Gy[i]], [Gz[i-1],Gz[i]],c='hotpink',linewidth=4.5)

# Wind Vector
ax.quiver(0, 0, 0, -(Wind_speed/100)*math.sin(Wind_direction*3.14159/180),-(Wind_speed/100)*math.cos(Wind_direction*3.14159/180),0)

if dist_init <= 10 :
   nx=1
   ny=1
else :
   nx=5
   ny=5
nz=0.5
xticks = np.arange(int(min(Lx)-1),int(max(Lx)+1),nx)
yticks = np.arange(0,int(max(Ly)+1),ny)
zticks = np.arange(int(min(Lz)-1),int(max(Lz)+1),nz)

ax.set_xticks(xticks)
ax.set_yticks(yticks)
ax.set_zticks(zticks)

ax.set_box_aspect(aspect = (1.,(max(Ly)-min(Ly))/(max(Lx)-min(Lx)),2*(max(Lz)-min(Lz))/(max(Lx)-min(Lx))))

# Hide grid lines
#ax.grid(False)
#plt.axis('off')

ax.set_xticklabels(xticks, fontsize = 12)
ax.set_yticklabels(yticks, fontsize = 12)

ax.set_xlabel('X (Nm)', fontsize = 16)
ax.set_ylabel('Y (Nm)', fontsize = 16)

cbar=fig2.colorbar(p, ax=ax,shrink=0.3,pad = 0.005)
cbar.ax.tick_params(labelsize=14)
cbar.set_label(label='FRAME',size=18)

plt.savefig('Figures/Fig_overview.pdf', format="pdf", bbox_inches="tight")  
plt.show()

########################### Path characteristics ##########################

ax.tick_params(labelsize=12)

# Ground Speed
vx=[]
vy=[]
vz=[]
v2d=[]
v3d=[]

for i in range(104):
    vx.append(0)
    vy.append(0)
    vz.append(0)
    v2d.append(0)
    v3d.append(0)

key=0
for i in range(2,104):
    dx=Gx[i]-Gx[i-1]
    dy=Gy[i]-Gy[i-1]
    dz=Gz[i]-Gz[i-1]
    if ((Wind_direction<=0 and dx>0) or (Wind_direction>0 and dx<0)) and key==0:
       key=1 
       revf=frame[i]
       print("Reverse direction at frame:",int(revf)) 

    vx[i]=((dx/(dt*10))*3600)
    vy[i]=((dy/(dt*10))*3600)
    vz[i]=(dz/(dt*10))*3600

    v2d[i]=math.sqrt((vx[i]*vx[i])+(vy[i]*vy[i]))
    v3d[i]=math.sqrt((vx[i]*vx[i])+(vy[i]*vy[i])+(vz[i]*vz[i]))

v2d[0]=v2d[2]
v2d[1]=v2d[2]
v3d[0]=v3d[2]
v3d[1]=v3d[2]

# Air speed
Vx=[]
Vy=[]
Vz=[]
V2d=[]
V3d=[]

for i in range(104):
    Vx.append(0)
    Vy.append(0)
    Vz.append(0)
    V2d.append(0)
    V3d.append(0)

key=0
for i in range(2,104):
    dx=Gx[i]-Gx[i-1]
    dy=Gy[i]-Gy[i-1]
    dz=Gz[i]-Gz[i-1]

    Vx[i]=((dx/(dt*10))*3600)+Wind_speed*math.sin(Wind_direction*3.14159/180)
    Vy[i]=((dy/(dt*10))*3600)+Wind_speed*math.cos(Wind_direction*3.14159/180)
    Vz[i]=(dz/(dt*10))*3600

    V2d[i]=math.sqrt((Vx[i]*Vx[i])+(Vy[i]*Vy[i]))
    V3d[i]=math.sqrt((Vx[i]*Vx[i])+(Vy[i]*Vy[i])+(Vz[i]*Vz[i]))

V2d[0]=V2d[2]
V2d[1]=V2d[2]
V3d[0]=V3d[2]
V3d[1]=V3d[2]

# Smooth Speed
window=10
for i in range(0,104):
     v2d[i]=np.mean(v2d[i:i+window])
     V2d[i]=np.mean(V2d[i:i+window])

# Plot Graph
from matplotlib.ticker import AutoMinorLocator

fig = plt.figure(figsize=(12,4))
ax = fig.add_subplot(111)

plt.title("Ground and Air Speed")
ax.xaxis.set_minor_locator(AutoMinorLocator())

ax.plot(frame,v2d,linewidth=2.0,color="darkgreen")
ax.plot(frame,V2d,linewidth=2.0,color="darkblue")
#ax.plot(frame,Vz,linewidth=2.0,color="red")

ax.axvspan(861,920, color='orange', alpha=0.5, lw=0)
ax.axvline(720, color='orange', lw=2, alpha=0.7,linestyle='dotted')
ax.axvline(820, color='orange', lw=2, alpha=0.7,linestyle='dotted')
ax.axvline(980, color='orange', lw=2, alpha=0.7,linestyle='dotted')
#ax.axvline(int(revf), color='green', lw=3, alpha=0.7,linestyle='dashed')

ax.set_xlabel("Frame", fontsize=14)
ax.set_ylabel("Speed (Kts)", fontsize=14)
plt.legend(["Ground Speed","Air Speed"])
plt.savefig("Figures/Fig_speed.pdf")
plt.close()


# Altitude
altitude=[]
for i in range(104):
    altitude.append(0)

for i in range(0,104):
    altitude[i]=25000+(Gz[i]*6076.12)
altitude[0]=altitude[1]

from matplotlib.ticker import AutoMinorLocator

fig = plt.figure(figsize=(12,4))
ax = fig.add_subplot(111)

#plt.xlim([0,104])
#plt.ylim([0,24000])

plt.title("Altitude")
ax.xaxis.set_minor_locator(AutoMinorLocator())

ax.plot(frame,altitude,linewidth=2.0,color="darkred")

ax.set_xlabel("Frame", fontsize=14)
ax.set_ylabel("Altitude (feet)", fontsize=14)
plt.legend(["Altitude"])
plt.savefig("Figures/Fig_alt.pdf")
plt.close()

# Heading
heading=[]
for i in range(104):
    heading.append(0)

for i in range(1,104):
    heading[i]=(np.arctan((Gx[i]-Gx[i-1])/(Gy[i]-Gy[i-1])))*180/3.14159
heading[0]=heading[2]
heading[1]=heading[2]

window=10
for i in range(0,104):
     heading[i]=np.mean(heading[i:i+window])

from matplotlib.ticker import AutoMinorLocator

fig = plt.figure(figsize=(12,4))
ax = fig.add_subplot(111)

plt.title("UAP Heading")
ax.xaxis.set_minor_locator(AutoMinorLocator())

ax.plot(frame,heading,linewidth=2.0,color="darkred")

ax.set_xlabel("Frame", fontsize=14)
ax.set_ylabel("Heading (deg)", fontsize=14)
plt.legend(["Heading"])
plt.savefig("Figures/Fig_UAP_heading.pdf")
plt.close()

# Aspect angle
head_obj_angle=[]
head_f18_angle=[]
asp_angle=[]
for i in range(104):
    head_obj_angle.append(0)
    head_f18_angle.append(0)
    asp_angle.append(0)

for i in range(1,104):
    head_obj_angle[i]=np.arctan((Gx[i]-Gx[i-1])/(Gy[i]-Gy[i-1]))*180/3.14159
    head_f18_angle[i]=np.arctan((x2[i]-x2[i-1])/(y2[i]-y2[i-1]))*180/3.14159
    asp_angle[i]=head_obj_angle[i]-head_f18_angle[i]-Az_smth[i]

asp_angle[0]=asp_angle[4]
asp_angle[1]=asp_angle[4]
asp_angle[2]=asp_angle[4]
asp_angle[3]=asp_angle[4]

# Smooth
window=5
for i in range(0,104):
     asp_angle[i]=np.mean(asp_angle[i:i+window])

from matplotlib.ticker import AutoMinorLocator

fig = plt.figure(figsize=(12,4))
ax = fig.add_subplot(111)

plt.title("Aspect Angle")
ax.xaxis.set_minor_locator(AutoMinorLocator())

ax.plot(frame,asp_angle,linewidth=2.0,color="darkred")

ax.set_xlabel("Frame", fontsize=14)
ax.set_ylabel("Aspect Angle (deg)", fontsize=14)
ax.tick_params(bottom=True, top=False, left=True, right=True)

plt.legend(["Aspect Angle"])
plt.savefig("Figures/Fig_asp_angle.pdf")
plt.close()

# Distance
from matplotlib.ticker import AutoMinorLocator

fig = plt.figure(figsize=(12,4))
ax = fig.add_subplot(111)

plt.title("Distance")
ax.xaxis.set_minor_locator(AutoMinorLocator())

distG2[0]=distG2[1]
ax.plot(frame,distG2,linewidth=2.0,color="darkblue")

ax.set_xlabel("Frame", fontsize=14)
ax.set_ylabel("Distance (Nm)", fontsize=14)
ax.tick_params(bottom=True, top=False, left=True, right=True)

plt.legend(["Distance"])
plt.savefig("Figures/Fig_distance.pdf")
plt.close()

# F-18 heading and wind angle
heading2=[]
wind_ang=[]
for i in range(104):
    heading2.append(0)
    wind_ang.append(0)

for i in range(1,104):
    heading2[i]=(np.arctan((x[i]-x[i-1])/(y[i]-y[i-1])))*180/3.14159
    wind_ang[i]=Wind_direction-heading2[i]
heading2[0]=heading2[1]
wind_ang[0]=wind_ang[1]

window=10
for i in range(0,104):
     heading2[i]=np.mean(heading2[i:i+window])
     wind_ang[i]=np.mean(wind_ang[i:i+window])

from matplotlib.ticker import AutoMinorLocator

fig = plt.figure(figsize=(12,4))
ax = fig.add_subplot(111)

plt.title("F-18 Heading and wind angle")
ax.xaxis.set_minor_locator(AutoMinorLocator())

ax.plot(frame,heading2,linewidth=2.0,color="darkred")
ax.plot(frame,wind_ang,linewidth=2.0,color="darkblue")
#ax.axvline(600, color='orange', lw=2, alpha=0.7,linestyle='dotted')
ax.axvline(720, color='orange', lw=2, alpha=0.7,linestyle='dotted')

ax.axhline(0, color='black', lw=2, alpha=0.7,linestyle='dotted')

ax.set_xlabel("Frame", fontsize=14)
ax.set_ylabel("Heading (deg)", fontsize=14)
plt.legend(["Heading","Wind Angle to the F-18"])
plt.savefig("Figures/Fig_F18_heading_wind-angle.pdf")
plt.close()

# Angle between pod LOS and direction of travel
ang_LOS_grndtrack=[]
for i in range(104):
    ang_LOS_grndtrack.append(0)

for i in range(0,104):
    ang_LOS_grndtrack[i]=heading2[i]-head_f18_angle[i]-Az_smth[i]+LoS_smth[i]

from matplotlib.ticker import AutoMinorLocator

fig = plt.figure(figsize=(12,4))
ax = fig.add_subplot(111)

plt.title("Angle between pod LOS and direction of travel")
ax.xaxis.set_minor_locator(AutoMinorLocator())

ax.plot(frame,ang_LOS_grndtrack,linewidth=2.0,color="darkred")
ax.axvline(720, color='orange', lw=2, alpha=0.7,linestyle='dotted')

ax.axhline(0, color='black', lw=2, alpha=0.7,linestyle='dotted')

ax.set_xlabel("Frame", fontsize=14)
ax.set_ylabel("Angle (deg)", fontsize=14)
str1 = 'Initial Wind direction = ' + str(Wind_direction)
print(str1)
plt.legend([str1,"14L"])
plt.savefig("Figures/Fig_ang_LOS_grndtrack.pdf")
plt.close()
