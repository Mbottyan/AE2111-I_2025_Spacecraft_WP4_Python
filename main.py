
import math
import random







#Constant forces from FBD
Fx = 10 #N
Fy = 20 #N
Fz = 30 #N
Mx = 0  #Nm
My = 0 #Nm
Mz = 0 #Nm
t2=0.005 #m



def Number_Of_Fasteners(w, D_2, material, N_min):
    # x-z plane defined same as in 4.1

    # defining constraints

    # edge constraints are given as a range in 4.4
    if material == " 1 ":
        edge_center_min = [2, 3]
    elif material == " 2 ":
        edge_center_min = [4, 5]
    
    # fastener spacing always the same
    center_center_min = 1.5

    N_max_check = [] # Number of fasteners check
    for e in range(0, len(edge_center_min)):

        something = 1 + (w-2*edge_center_min[e]*D_2)/(center_center_min*D_2) # unrounded N_max
        N_m = 1 + (w-2*edge_center_min[e]*D_2)//(center_center_min*D_2)
        N_max_check.append(N_m)
    
        e_test = ( N_m - 1 - w/(center_center_min*D_2) ) * -center_center_min/2 # chooses highest N and its associated edge spacing
        if e_test > edge_center_min[e]:
            N_max = max(N_max_check)
            edge_spacing = e_test
            #print(edge_spacing, edge_center_min[e], N_max, "1") # for testing

        w_new = (N_min-1) * center_center_min + 2 * min(edge_center_min)
        #print(e_test, N_max, "2") # for testing

    print("Max N, given w:",N_max,",", edge_spacing,",", "New w, given minimum N:", w_new, N_min)
    
    return N_max, edge_spacing, center_center_min, w_new, N_min
    #returns max number of fasteners, edge spacing, center spacing, new width for minimum fasteners, minimum fasteners

class Fastener:
    def __init__(self, Diameter, x_coord, z_coord):
        #sets coordinates and dimensions for each fastener in our design
        self.Diameter=float(Diameter)
        self.x_coord=float(x_coord)
        self.z_coord=float(z_coord)
        self.force_vectors=((0,0,0),(0,0,0),(0,0,0)) #will hold the force vectors assigned to each fastener
    # give the coordinates weighted and areas of fastener of cg calculation
    def provide_x_weighted_average(self):
        self.area=(math.pi)*(self.Diameter*0.5)**2
        return (self.area*self.x_coord), (self.area)
    def provide_z_weighted_average(self):
        return (self.area*self.z_coord), (self.area)
    def find_bearinhg_stresses (self):
        #calculate magnitude of z and x component forces, calculate the stress.
        x_forces=(self.fcgx[0]+self.fmom[0])
        z_forces=(self.fcgz[1]+self.fmom[1])
        self.Pi=(x_forces,z_forces)
        self.Pi_magnitude=math.sqrt(x_forces**2+z_forces**2)
        self.bearing_stress=self.Pi_magnitude/(self.Diameter*t2)
        return (self.Pi, self.Pi_magnitude, self.bearing_stress)
        #produces tuple with the (force-vector, magnitude, bearing stress) for comparison with maximum
        print('Local wall thickness should be'+str(self.Pi_magnitude/(stress_max*self.Diameter)))
        #Prints the local wall thickness required

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

# test case for Number_Of_Fastners function
# print(Number_Of_Fasteners(15, 2, " 1 ", 2))

