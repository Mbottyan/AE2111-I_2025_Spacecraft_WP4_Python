
import math
import random







#Constant forces from FBD
Fx = 10 #N
Fy = 20 #N
Fz = 30 #N
Mx = 5  #Nm
My = 15 #Nm
Mz = 25 #Nm




class Fastener:
    def __init__(self, Diameter, x_coord, z_coord):
        #sets coordinates and dimensions for each fastener in our design
        self.Diameter=float(Diameter)
        self.x_coord=float(x_coord)
        self.z_coord=float(z_coord)
    # give the coordinates weighted and areas of fastener of cg calculation
    def provide_x_weighted_average(self):
        self.A=(math.pi)*(self.Diameter*0.5)**2
        return (self.A*self.x_coord), (self.A)
    def provide_z_weighted_average(self):
        self.A=(math.pi)*(self.Diameter*0.5)**2
        return (self.A*self.z_coord), (self.A)
        

    
Fasteners=[]
#create list for all fastener instances
#generate random fasteners, TO BE REPLACED WITH REAL ONES! When ready with fastener dimensions
# and coordinates, for each one append manually:
#Fasteners.append(Fastener(diameter,x coodinate, zcoordinate))
#This will create a Fastener instance for each fastener, which will be used in the cg calculation
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

(cg_location())
#This produces the coordinates of the centre of gravity in (xcg, zcg) format