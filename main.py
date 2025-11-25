
import math
import random

#Functions list (able to be called whenever you need these values in later code):

    # - cg_location() - Will give you a tuple with the (x, z) coordinates of the centre of gravity of the fasteners, you can index these too!
    # - Number_Of_Fasteners((w, D_2, material, N_min) - #returns max number of fasteners, edge spacing, center spacing, new width for minimum fasteners, minimum fasteners in a tuple




#Constant forces from FBD
Fx = 97.119 #N      #plus or minus
Fy = 97.119 #N      #plus or minus
Fz = -425 #N
Mx = -386.0275  #Nm #88.2132 to -386.0275
My = 88.2132 #Nm   #plus or minus
Mz = 1.8166 #Nm     #plus or minus
w = 0 #m  (Put in the real value here)
h = 0 #m  (Put in the real value here)
t1 = 0 #m  (Put in the real value here)
t2=0.005 #m (Put in the real value here)
t3 = 0 #m  (Put in the real value here)
D_1 = 0 #m  (Put in the real value here)
D_2 = 0 #m  (Put in the real value here)
P=0 #N Make a function to find P below and use it to give this variable the correct value

Materials = {'Aluminium': {'type (metal or composite)': 1, 'Modulus': 69, 'Thermal Coefficient': 23*10^{-6}}, 'Carbon Composite': {'category (metal or composite)': 2, 'Modulus': 200}, 'Titanium': {'type (metal or composite)': 1, 'Modulus': 124},  'Thermal Coefficient': 8.6*10^{-6}}


material_used = 'Aluminium'



Fasteners=[] #create list for all fastener instances
# When ready with fastener dimensions and coordinates, for each one append manually:
#Fasteners.append(Fastener(diameter,x coodinate, zcoordinate))
#This will create a Fastener instance for each fastener, which will be used in the cg calculation
#hello


class Fastener:
    def __init__(self, Diameter, x_coord, z_coord):
        #sets coordinates and dimensions for each fastener in our design
        self.Diameter=float(Diameter)
        self.x_coord=float(x_coord)
        self.z_coord=float(z_coord)
        self.force_vectors=((0,0,0),(0,0,0),(0,0,0)) #will hold the force vectors assigned to each fastener (xforces, zforces, momentforces)
    # give the coordinates weighted and areas of fastener of cg calculation
    def provide_x_weighted_average(self):
        self.area=(math.pi)*(self.Diameter*0.5)**2
        return (self.area*self.x_coord), (self.area)
    def provide_z_weighted_average(self):
        return (self.area*self.z_coord), (self.area)
    def find_bearing_stresses (self):
        #calculate magnitude of z and x component forces, calculate the stress.
        x_forces=(self.force_vectors[0][0]+self.force_vectors[1][0]+self.force_vectors[2][0])
        z_forces=(self.force_vectors[0][2]+self.force_vectors[1][2]+self.force_vectors[2][2])
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

# calculates max ammount of fasteners that fit on the lug, in case a higher ammount of fasteners (N_min) is needed, it checks and if necessary recalculates the lug height to fit (N_min) ammount
def Number_Of_Fasteners(w, D_2, material, N_min): 
    # x-z plane defined same as in 4.1

    # defining constraints

    # edge constraints are given as a range in 4.4
    if(Materials[material_used]['type (metal or composite)']) == 1 :
        edge_center_min = [2, 3]*D_2
    elif (Materials[material_used]['type (metal or composite)']) == 2 :
        edge_center_min = [4, 5]*D_2
    
    # fastener spacing always the same
    center_center_min = 1.5*D_2

    N_max_check = [] # Number of fasteners check
    for e in range(0, len(edge_center_min)):

        N_m = 1 + (w-2*edge_center_min[e])//(center_center_min)
        N_max_check.append(N_m)
        
        e_test = ( w - ( N_m - 1 ) * center_center_min) /2 # chooses highest N and its associated edge spacing
        if e_test > float(edge_center_min[e]):
            N_max = int(max(N_max_check))
            edge_spacing = e_test
            #print(e_test, edge_center_min, N_max, "1") # for testing
    
    if N_min > N_max: # if more fasteners are needed, this will calculate the spacing
        w = (N_min-1) * center_center_min + 2 * min(edge_center_min)

        # rawdogging it to only give relevant information as output, this line doesnt make sense but trust
        N_max = N_min
        edge_spacing = min(edge_center_min)*D_2
        #print(edge_spacing, N_max, "2") # for testing

    return N_max, edge_spacing, center_center_min, w, D_2

    #returns max number of fasteners, edge spacing, center spacing, new width for minimum fasteners, minimum fasteners

# test case for Number_Of_Fastners function
#print(Number_Of_Fasteners(15, 2, " 1 ", 3))

# calculates locations of fastener centroids
def Fasteners_location(N_max: int, edge_spacing, center_center_min, w_new, h, t1, D_2): # w, h and t1 are
    # at the time of writing this code atleast the fasteners are symmetrically split across the z axis so this calculates positive x-axis fasteners and negative x-axis ones separately
    for f in range(N_max):
        
        diameter = D_2 # as of writing this, every fastener has same diameter, if this changes, this has to be rewritten

        # first half of fasteners on positive x-axis
        x = h/2 + t1 + edge_spacing  # only 1 column of fasteners per side so far, if later that seems to not be enough the calculation of x has to change
        z = f * center_center_min + edge_spacing - w_new/2
        Fasteners.append(Fastener(diameter, x , z))
        #print(f, x, z)
        # second half of fasteners on negative x-axis
        x = -h/2 - t1 - edge_spacing
        Fasteners.append(Fastener(diameter, x , z))
        #print(f, x, z)
        
    return Fasteners
#print(Fasteners_location(4, 3.0, 3.0, 15, 2, 1, 2))


#generate random fasteners, TO BE REPLACED WITH REAL ONES! When ready with fastener dimensions
# and coordinates, for each one append manually:
#Fasteners.append(Fastener(diameter,x coodinate, zcoordinate))
#This will create a Fastener instance for each fastener, which will be used in the cg calculation

#for i in range(4):
    #Fasteners.append(Fastener(0.01,random.randint(0,5),random.randint(0,5)))


#Translate foces into the cg of the fastener pattern
Fcgx = Fx
Fcgz = Fz
Mcgy = My #+ (Fz*cg_location()[0]) - (Fx*cg_location()[1])

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


#material_type = Materials[material_used]['type (metal or composite)']

def Compliance_parts(Modulus,D_outer,D_inner,thickness):
    Compliance = 4*thickness/(Modulus*math.pi*(D_outer**2 - D_inner**2))
    return Compliance

def Compliance_fastener(Modulus,Cross_area,length):
    Compliance = length/(Modulus*Cross_area)
    return Compliance

Compliance_a = Compliance_parts(Materials['Aluminium']['Modulus'],D_fo,D_fi,t2)
Compliance_b = Compliance_fastener(Materials['Titanium']['Modulus'],(0.5*D_fi)**2*math.pi,t2)

#Calculate force ratio
def force_ratio(Compliance_a, Compliance_b):
    force_ratio = Compliance_a/(Compliance_a + Compliance_b)
    return force_ratio



