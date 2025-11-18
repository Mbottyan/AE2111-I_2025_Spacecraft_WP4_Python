
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
        self.A=(math.pi)*(self.Diameter*0.5)**2
        return (self.A*self.z_coord), (self.A)
        

    
Fasteners=[]
#create list for all fastener instances
#generate random fasteners, TO BE REPLACED WITH REAL ONES!
for i in range(4):
    Fasteners.append(Fastener(0.01,random.randint(0,5),random.randint(0,5)))

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
print(cg_location())