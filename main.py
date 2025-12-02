import math
import random
import numpy as np

#Functions list (able to be called whenever you need these values in later code):

    # - cg_location() - Will give you a tuple with the (x, z) coordinates of the centre of gravity of the fasteners, you can index these too!
    # - Number_Of_Fasteners((w, D_2, material, N_min) - #returns max number of fasteners, edge spacing, center spacing, new width for minimum fasteners, minimum fasteners in a tuple


MS_main = []

#Constant forces from FBD
Fx = 97.119 #N      #plus or minus
Fy = 97.119 #N      #plus or minus
Fz = -425 #N
Mx = -386.0275  #Nm #88.2132 to -386.0275
My = 88.2132 #Nm   #plus or minus
Mz = 1.8166 #Nm     #plus or minus
w = 0.1 #m  (Put in the real value here)
h = 0.05 #m  (Put in the real value here)
t1 = 0.02 #m  (Put in the real value here)
t2=0.005 #m (Put in the real value here)
t3 = 0.004 #m  (Put in the real value here)
t3_list=[]
t3_2_list=[]
D_1 = 0 #m  (Put in the real value here)
D_2 = 0.01 #m  (Put in the real value here)
D_in = 0.001 #m  (Put in the real value here)
D_fo = 0.002 #m  (Put in the real value here)
P=0 #N Make a function to find P below and use it to give this variable the correct value
safety_factor=1

Materials = {'Aluminium': {'type (metal or composite)': 1, 'Modulus': 73500000000, 'Thermal Coefficient': 23*10**(-6), 'Yield Stress':345000000 }, 'Carbon Composite': {'type (metal or composite)': 2, 'Modulus': 230000000000, 'Yield Stress': 4400000000}, 'Titanium': {'type (metal or composite)': 1, 'Modulus': 124000000000, 'Yield Stress': 170000000, 'Thermal Coefficient': 8.6*10**(-6)}}


material_used = 'Aluminium'
material_used_fastener = 'Titanium'
material_used_body = 'Aluminium'

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
        self.force_vectors_inplane=((0,0,0),(0,0,0),(0,0,0)) #will hold the force vectors assigned to each fastener (xforces, zforces, momentforces)
        self.force_vectors_outofplane=((0,0,0),(0,0,0)) #will hold the force vectors assigned to each fastener (yforces, shearforces, outofplanemomentforces)
        self.Pi_magnitude=0
        self.passes_bearing=False
        self.passes_pullthrough=False

        #safey factors MS
        self.MS_t2_bearing=0
        self.MS_t2_bearing_cold=0
        self.MS_t2_bearing_hot=0
        self.MS_pullthrough_t2=0
        self.MS_t3_bearing_thermal=0
        self.MS_pullthrough_t3=0
    # give the coordinates weighted and areas of fastener of cg calculation
    def provide_x_weighted_average(self):
        self.area=(math.pi)*(self.Diameter*0.5)**2
        return (self.area*self.x_coord), (self.area)
    def provide_z_weighted_average(self):
        return (self.area*self.z_coord), (self.area)
    def find_bearing_stresses(self,bearing_allowable_stress,thermal=0):
        #calculate magnitude of z and x component forces, calculate the stress.
        x_forces=(self.force_vectors_inplane[0][0]+self.force_vectors_inplane[1][0]+self.force_vectors_inplane[2][0])
        z_forces=(self.force_vectors_inplane[0][2]+self.force_vectors_inplane[1][2]+self.force_vectors_inplane[2][2])
        self.Pi=(x_forces,z_forces)
        if self.Pi_magnitude==0:
            self.Pi_magnitude=math.sqrt(x_forces**2+z_forces**2)
        if thermal==0:
            self.bearing_stress=self.Pi_magnitude/(self.Diameter*t2)
            bearingstressfortest=self.bearing_stress
            self.MS_t2_bearing=bearing_allowable_stress/bearingstressfortest - 1
        elif thermal==1:
            self.bearing_stress_cold=self.Pi_magnitude/(self.Diameter*t2)
            bearingstressfortest=self.bearing_stress_cold
            self.MS_t2_bearing_cold=bearing_allowable_stress/bearingstressfortest - 1
        elif thermal==2:
            self.bearing_stress_hot=self.Pi_magnitude/(self.Diameter*t2)
            bearingstressfortest=self.bearing_stress_hot
            self.MS_t2_bearing_hot=bearing_allowable_stress/bearingstressfortest - 1

        
        if safety_factor*bearingstressfortest<bearing_allowable_stress:
            self.passes_bearing=True
        else:
            print('Attention: The fastener located at ('+str(self.x_coord)+', 0, '+str(self.z_coord)+') is expected to fail in bearing by a factor of '+str(safety_factor*bearingstressfortest/bearing_allowable_stress)+'!!!!')
            self.passes_bearing=False
        self.local_wall_thickness=safety_factor*self.Pi_magnitude/(bearing_allowable_stress*self.Diameter)
        #print('Local wall thickness for the fastener at ('+str(self.x_coord)+', 0, '+str(self.z_coord)+') should be '+str(self.local_wall_thickness*1000)+'mm.')

        return (self.Pi, self.Pi_magnitude, self.bearing_stress)
        #produces tuple with the (force-vector, magnitude, bearing stress) for comparison with maximum

    def check_pull_through_failure(self, yield_stress_t2, yield_stress_t3):
        # Pull-through load (magnitude)
        self.p_pull = abs(self.force_vectors_outofplane[0][1]+self.force_vectors_outofplane[1][1])  #The y component forces
        
        # Shear Area = pi * (D_fi) * t
        
        #t2 lug
        area_shear_t2 = math.pi * D_in * t2
        self.shear_stress_t2 = self.p_pull / area_shear_t2 if area_shear_t2 > 0 else 0
        
        #t3 wall
        area_shear_t3 = math.pi * D_in * t3
        self.shear_stress_t3 = self.p_pull / area_shear_t3 if area_shear_t3 > 0 else 0
        
        #Von Mises Stress (Eq 4.8)
        #Sigma_vm = sqrt(3 * tau^2)
        sigma_vm_t2 = math.sqrt(3 * self.shear_stress_t2**2)
        sigma_vm_t3 = math.sqrt(3 * self.shear_stress_t3**2)

        self.MS_pullthrough_t2 = yield_stress_t2 / sigma_vm_t2 - 1
        self.MS_pullthrough_t3 = yield_stress_t3 / sigma_vm_t3 - 1
        if (safety_factor * sigma_vm_t2 < yield_stress_t2) and (safety_factor * sigma_vm_t3 < yield_stress_t3):
            self.passes_pullthrough = True
        else:
            print('Attention: The fastener located at ('+str(self.x_coord)+', 0, '+str(self.z_coord)+') is expected to fail in pull through!!!!')




    




# calculates max ammount of fasteners that fit on the lug, in case a higher ammount of fasteners (N_min) is needed, it checks and if necessary recalculates the lug height to fit (N_min) ammount
def Number_Of_Fasteners(w, D_2, N_min): 
    # x-z plane defined same as in 4.1

    # defining constraints

    # edge constraints are given as a range in 4.4
    if(Materials[material_used]['type (metal or composite)']) == 1 :
        edge_center_min = np.array([2, 3])*D_2
    elif (Materials[material_used]['type (metal or composite)']) == 2 :
        edge_center_min = np.array([4, 5])*D_2
    
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
#works
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
#works

#generate random fasteners, TO BE REPLACED WITH REAL ONES! When ready with fastener dimensions
# and coordinates, for each one append manually:
#Fasteners.append(Fastener(diameter,x coodinate, zcoordinate))
#This will create a Fastener instance for each fastener, which will be used in the cg calculation

#for i in range(4):
    #Fasteners.append(Fastener(0.01,random.randint(0,5),random.randint(0,5)))

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
    return (x_num_sum/x_den_sum,0, z_num_sum/z_den_sum)


#Assign forces to each fastener based on the formulas provided in 4.5
def assign_fastener_forces():
    cg_x, cg_y, cg_z = cg_location()
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
        F_pi=(0.0,Fy/nf if nf else 0.0,0.0)
        if area_r2_sum>0 and r>0:
            magnitude=Mcgy*fastener.area*r/area_r2_sum
            magnitude_outofplane=Mz*fastener.area*r/area_r2_sum
            tangential=(-dz/r,dx/r) #the tangential thingy from the figure 4.5
            moment_force=(magnitude*tangential[0],0.0,magnitude*tangential[1])
            moment_outofplane_force=(0.0,magnitude_outofplane,0.0)
        else:
            moment_force=(0.0,0.0,0.0)
        fastener.force_vectors_inplane=(F_inplanex,F_inplanez,moment_force)
        fastener.force_vectors_outofplane=(F_pi, moment_outofplane_force)

#assign_fastener_forces()
#def bearing_passes_func():
    #bearing_passes=0
    #for fastn in Fasteners:
        #fastn.find_bearing_stresses(Materials[material_used]['Yield Stress'])
        #if fastn.passes_bearing==True:
            #bearing_passes+=1
    #t3_list.append(fastn.local_wall_thickness)
    #if bearing_passes==len(Fasteners):
        #print('All fasteners pass the bearing check.')


#material_type = Materials[material_used]['type (metal or composite)']

def Compliance_parts(Modulus,D_outer,D_inner,thickness):
    Compliance = 4*thickness/(Modulus*math.pi*(D_outer**2 - D_inner**2))
    return Compliance

def Compliance_fastener(Modulus,Cross_area,length):
    Compliance = length/(Modulus*Cross_area)
    return Compliance

#Compliance_a = Compliance_parts(Materials['Aluminium']['Modulus'],D_fo,D_fi,t2)
#Compliance_b = Compliance_fastener(Materials['Titanium']['Modulus'],(0.5*D_fi)**2*math.pi,t2)

#Calculate force ratio
def force_ratio(Compliance_a, Compliance_b):
    force_ratio = Compliance_a/(Compliance_a + Compliance_b)
    return force_ratio


#Thermal expansion created stress function
def thermal1():
    a_c = (Materials[material_used]['Thermal Coefficient'])
    a_f = (Materials[material_used_fastener]['Thermal Coefficient'])
    T_ref = 15
    T_operate = [-100, 130]
    psi = force_ratio(delta_a, delta_b)
    lst = np.zeros((len(Fasteners), len(T_operate)))
    thermal_failure = False

    bearing_passes=0
    for i, item in enumerate(Fasteners):
        A_sw = item.provide_x_weighted_average()[1]
        E_b = (Materials[material_used_fastener]['Modulus'])
        for j, T in enumerate(T_operate):
            delta_T = T - T_ref
            F_t = (a_c - a_f) * delta_T * E_b * A_sw * (1 - psi)
            lst[i, j] = F_t
        
            # bearing check
            #print(item.Pi_magnitude)
            #print('ft'+str(F_t))
            placeholder=item.Pi_magnitude
            item.Pi_magnitude=(item.Pi_magnitude+F_t)

            #print(item.Pi_magnitude)
            item.find_bearing_stresses(Materials[material_used]['Yield Stress'],j+1)
            if item.passes_bearing==True:
                bearing_passes+=1
            t3_2_list.append(item.local_wall_thickness)
            

            item.Pi_magnitude=placeholder

            #Stress_max = (Materials[material_used]['Yield Stress'])

            #if Stress > Stress_max:
                #thermal_failure = True
    if bearing_passes==len(Fasteners):
        print('All fasteners pass the thermal bearing check.')
    print('Minimum required wall thicknesses after thermal (in m):', max(t3_2_list))

    return thermal_failure


#                                                             #
##                                                           ##
###                                                         ###
####                                                       ####
#####                                                     #####
######                                                   ######
#######                                                 #######
######                                                   ######
#####                                                     #####
####                                                       ####
###                                                         ###
##                                                           ##
#                                                             #

         #Generate fastener coordinates
NOF=Number_Of_Fasteners(w,D_2,3)
Fasteners_location(NOF[0],NOF[1],NOF[2],NOF[3],h,t1,NOF[4])
fastener_zcoords=[]
for fastn in Fasteners:
    fastener_zcoords.append(fastn.z_coord)

plate_centre=(0,0,((min(fastener_zcoords)+max(fastener_zcoords))*0.5))
Fastener_Quantity=len(Fasteners)

#Translate foces into the cg of the fastener pattern
Fcgx = Fx
Fcgz = Fz
Mcgy = My + (Fz*cg_location()[0]-plate_centre[0]) - (Fx*cg_location()[2]-plate_centre[2])
#print(cg_location())
#print(Fcgx, Fcgz, Mcgy)


#Now we run the assign fastener forces function, which directly affects the Fastener instances
assign_fastener_forces()

#Time to do the stresss checks
#I'm using yield stress as the comparison figure, we need to double check if this is right!
#-----------------------------------------------------------------------------------------
bearing_passes=0
pull_through_passes=0
for fastn in Fasteners:
    #bearing
    fastn.find_bearing_stresses(Materials[material_used]['Yield Stress'])
    if fastn.passes_bearing==True:
        bearing_passes+=1
    #print((fastn.bearing_stress))
    t3_list.append(fastn.local_wall_thickness)

    #pull through
    fastn.check_pull_through_failure(Materials[material_used]['Yield Stress'], Materials[material_used_body]['Yield Stress'])
    if fastn.passes_pullthrough==True:
        pull_through_passes+=1

if bearing_passes==len(Fasteners):
    print('All fasteners pass the bearing check.')
if pull_through_passes==len(Fasteners):
    print('All fasteners pass the pull through check.')
print('Minimum required wall thicknesses(in m):', max(t3_list))


#Compliances
delta_a=Compliance_parts(Materials[material_used]['Modulus'],D_fo, D_in, t3)
delta_b=Compliance_fastener(Materials[material_used]['Modulus'],(math.pi*(D_in/2)**2), 0.03)

#thermal stress check
thermal1()

for fastn in Fasteners:
    MS_main.append((fastn.MS_t2_bearing, min(fastn.MS_t2_bearing_cold, fastn.MS_t2_bearing_hot), fastn.MS_pullthrough_t2, fastn.MS_t3_bearing_thermal, fastn.MS_pullthrough_t3))

print(MS_main)