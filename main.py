
import math
import random

#Constant forces from FBD
Fx = 10 #N
Fy = 20 #N
Fz = 30 #N
Mx = 0  #Nm
My = 0 #Nm
Mz = 0 #Nm

#Class with the attributes and methods for each fastener
class Fastener:

    def __init__(self, Diameter, x_coord, z_coord):
        #sets coordinates and dimensions for each fastener in our design
        self.Diameter=float(Diameter)
        self.x_coord=float(x_coord)
        self.z_coord=float(z_coord)
        self.area=(math.pi)*(self.Diameter*0.5)**2
        self.force_vectors=((0.0,0.0,0.0),(0.0,0.0,0.0),(0.0,0.0,0.0))
    # give the coordinates weighted and areas of fastener of cg calculation
    def provide_x_weighted_average(self):
        return (self.area*self.x_coord), (self.area)
    def provide_z_weighted_average(self):
        return (self.area*self.z_coord), (self.area)

#Calculate the cg location of the fastener pattern  
def cg_location():
    x_num_sum=0
    x_den_sum=0
    z_num_sum=0
    z_den_sum=0
    for item in Fasteners:
        x_num_sum+=(item.provide_x_weighted_average()[0])
        x_den_sum+=(item.provide_x_weighted_average()[1])
        z_num_sum+=(item.provide_z_weighted_average()[0])
        z_den_sum+=(item.provide_z_weighted_average()[1])
    return (x_num_sum/x_den_sum, z_num_sum/z_den_sum)
    
Fasteners=[] #create list for all fastener instances
#generate random fasteners, TO BE REPLACED WITH REAL ONES! When ready with fastener dimensions
# and coordinates, for each one append manually:
#Fasteners.append(Fastener(diameter,x coodinate, zcoordinate))
#This will create a Fastener instance for each fastener, which will be used in the cg calculation
for i in range(4):
    Fasteners.append(Fastener(0.01,random.randint(0,5),random.randint(0,5)))

#Translate froces into the cg of the fastener pattern
Fcgx = Fx
Fcgz = Fz
Mcgy = My + (Fz*cg_location()[0]) - (Fx*cg_location()[1])

#Assign forces to each fastener based on the formulas provided in 4.5
def assign_fastener_forces():
    cg_x, cg_z = cg_location()
    nf=len(Fasteners)
    area_r2_sum=sum(
        fastener.area*((fastener.x_coord-cg_x)**2+(fastener.z_coord-cg_z)**2)
        for fastener in Fasteners
    )
    for fastener in Fasteners:
        dx=fastener.x_coord-cg_x
        dz=fastener.z_coord-cg_z
        r=math.hypot(dx,dz)
        F_inplanex=(Fcgx/nf if nf else 0.0,0.0,0.0)
        F_inplanez=(0.0,0.0,Fcgz/nf if nf else 0.0)
        if area_r2_sum>0 and r>0:
            magnitude=Mcgy*fastener.area*r/area_r2_sum
            tangential=(-dz/r,dx/r) #the tangential thingy from the figure 4.5
            moment_force=(magnitude*tangential[0],0.0,magnitude*tangential[1])
        else:
            moment_force=(0.0,0.0,0.0)
        fastener.force_vectors=(F_inplanex,F_inplanez,moment_force)

assign_fastener_forces()