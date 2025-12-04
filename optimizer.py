import math
import numpy as np
import json
# from scipy.optimize import minimize # Removed dependency

# Load parameters from JSON file
with open('parameters.json', 'r') as f:
    params = json.load(f)

# Constants
Fx = params['forces']['Fx']
Fy = params['forces']['Fy']
Fz = params['forces']['Fz']
Mx = params['forces']['Mx']
My = params['forces']['My']
Mz = params['forces']['Mz']
w_original = params['geometry']['w']
h = params['geometry']['h']
t1 = params['geometry']['t1']
D_in = params['geometry']['D_in']
D_fo = params['geometry']['D_fo']
safety_factor = params['safety_factor']

Materials = params['materials']
material_used = params['material_selection']['material_used']
material_used_fastener = params['material_selection']['material_used_fastener']
material_used_body = params['material_selection']['material_used_body']

T_ref = params['thermal']['T_ref']
T_operate = params['thermal']['T_operate']

# Material Properties
Yield_Stress_Material = Materials[material_used]['Yield Stress']
Yield_Stress_Body = Materials[material_used_body]['Yield Stress']
Modulus_Material = Materials[material_used]['Modulus']
Modulus_Fastener = Materials[material_used_fastener]['Modulus']
Thermal_Coeff_Material = Materials[material_used]['Thermal Coefficient']
Thermal_Coeff_Fastener = Materials[material_used_fastener]['Thermal Coefficient']

class Fastener:
    def __init__(self, Diameter, x_coord, z_coord):
        self.Diameter = float(Diameter)
        self.x_coord = float(x_coord)
        self.z_coord = float(z_coord)
        self.force_vectors_inplane = ((0,0,0),(0,0,0),(0,0,0))
        self.force_vectors_outofplane = ((0,0,0),(0,0,0))
        self.Pi_magnitude = 0
        self.area = (math.pi)*(self.Diameter*0.5)**2
        
        # Results
        self.MS_t2_bearing = 0
        self.MS_t2_bearing_cold = 0
        self.MS_t2_bearing_hot = 0
        self.MS_pullthrough_t2 = 0
        self.MS_t3_bearing_thermal = 0
        self.MS_pullthrough_t3 = 0

    def provide_x_weighted_average(self):
        return (self.area*self.x_coord), (self.area)
    
    def provide_z_weighted_average(self):
        return (self.area*self.z_coord), (self.area)

def Number_Of_Fasteners(w, D_2, N_min): 
    if(Materials[material_used]['type (metal or composite)']) == 1 :
        edge_center_min = np.array([2, 3])*D_2
    elif (Materials[material_used]['type (metal or composite)']) == 2 :
        edge_center_min = np.array([4, 5])*D_2
    
    center_center_min = 1.5*D_2

    N_max_check = []
    for e in range(0, len(edge_center_min)):
        N_m = 1 + (w-2*edge_center_min[e])//(center_center_min)
        N_max_check.append(N_m)
        e_test = ( w - ( N_m - 1 ) * center_center_min) /2
        if e_test > float(edge_center_min[e]):
            N_max = int(max(N_max_check))
            edge_spacing = e_test
    
    if N_min > N_max:
        w = (N_min-1) * center_center_min + 2 * min(edge_center_min)
        N_max = N_min
        edge_spacing = min(edge_center_min)*D_2

    return N_max, edge_spacing, center_center_min, w, D_2

def Fasteners_location(N_max, edge_spacing, center_center_min, w_new, h, t1, D_2):
    fasteners = []
    for f in range(int(N_max)):
        diameter = D_2
        # Positive x
        x = h/2 + t1 + edge_spacing
        z = f * center_center_min + edge_spacing - w_new/2
        fasteners.append(Fastener(diameter, x , z))
        # Negative x
        x = -h/2 - t1 - edge_spacing
        fasteners.append(Fastener(diameter, x , z))
    return fasteners

def cg_location(fasteners):
    x_num_sum=0
    x_den_sum=0
    z_num_sum=0
    z_den_sum=0
    for item in fasteners:
        x_num_sum+=(item.provide_x_weighted_average()[0])
        x_den_sum+=(item.provide_x_weighted_average()[1])
        z_num_sum+=(item.provide_z_weighted_average()[0])
        z_den_sum+=(item.provide_z_weighted_average()[1])
    return (x_num_sum/x_den_sum, 0, z_num_sum/z_den_sum)

def assign_fastener_forces(fasteners):
    cg = cg_location(fasteners)
    cg_x, cg_z = cg[0], cg[2]
    nf = len(fasteners)
    
    # Calculate plate centre for moment transfer
    z_coords = [f.z_coord for f in fasteners]
    plate_centre = (0, 0, (min(z_coords)+max(z_coords))*0.5)
    
    # Transfer moments to CG
    Fcgx = Fx
    Fcgz = Fz
    Mcgy = My + (Fz*cg_x - plate_centre[0]) - (Fx*cg_z - plate_centre[2]) # Simplified from original code logic
    # Note: Original code: Mcgy = My + (Fz*cg_location()[0]-plate_centre[0]) - (Fx*cg_location()[2]-plate_centre[2])
    # Assuming plate_centre[0] is 0 (x=0) and plate_centre[2] is center of pattern z.
    # cg_location()[0] is x_cg.
    
    area_r2_sum = sum(
        f.area*((f.x_coord-cg_x)**2+(f.z_coord-cg_z)**2)
        for f in fasteners
    )
    
    for f in fasteners:
        dx = f.x_coord - cg_x
        dz = f.z_coord - cg_z
        r = math.hypot(dx, dz)
        
        F_inplanex = (Fcgx/nf if nf else 0.0, 0.0, 0.0)
        F_inplanez = (0.0, 0.0, Fcgz/nf if nf else 0.0)
        F_pi = (0.0, Fy/nf if nf else 0.0, 0.0)
        
        if area_r2_sum > 0 and r > 0:
            magnitude = Mcgy * f.area * r / area_r2_sum
            magnitude_outofplane = Mz * f.area * r / area_r2_sum
            tangential = (-dz/r, dx/r)
            moment_force = (magnitude*tangential[0], 0.0, magnitude*tangential[1])
            moment_outofplane_force = (0.0, magnitude_outofplane, 0.0)
        else:
            moment_force = (0.0, 0.0, 0.0)
            moment_outofplane_force = (0.0, 0.0, 0.0)
            
        f.force_vectors_inplane = (F_inplanex, F_inplanez, moment_force)
        f.force_vectors_outofplane = (F_pi, moment_outofplane_force)
        
        # Calculate Pi magnitude for bearing
        x_forces = (f.force_vectors_inplane[0][0] + f.force_vectors_inplane[1][0] + f.force_vectors_inplane[2][0])
        z_forces = (f.force_vectors_inplane[0][2] + f.force_vectors_inplane[1][2] + f.force_vectors_inplane[2][2])
        f.Pi_magnitude = math.sqrt(x_forces**2 + z_forces**2)
        
        # Calculate Pull magnitude
        f.p_pull = abs(f.force_vectors_outofplane[0][1] + f.force_vectors_outofplane[1][1])

def Compliance_parts(Modulus, D_outer, D_inner, thickness):
    if thickness <= 0: return 1e9 # Penalty
    return 4*thickness/(Modulus*math.pi*(D_outer**2 - D_inner**2))

def Compliance_fastener(Modulus, Cross_area, length):
    return length/(Modulus*Cross_area)

def force_ratio(Compliance_a, Compliance_b):
    return Compliance_a/(Compliance_a + Compliance_b)

def solve_thickness(N_rows, D_2):
    # Iteratively solve for t2 and t3
    
    # Initial guess
    t2 = 0.001
    t3 = 0.001
    
    # Constants for this configuration
    NOF = Number_Of_Fasteners(w_original, D_2, N_rows)
    fasteners = Fasteners_location(NOF[0], NOF[1], NOF[2], NOF[3], h, t1, NOF[4])
    assign_fastener_forces(fasteners)
    
    # Max loads across all fasteners (Mechanical only first)
    max_Pi = 0
    max_Pull = 0
    for f in fasteners:
        if f.Pi_magnitude > max_Pi: max_Pi = f.Pi_magnitude
        if f.p_pull > max_Pull: max_Pull = f.p_pull
        
    # 1. Minimum t2 for Mechanical Bearing
    # sigma = Pi / (D * t2) <= Yield / SF
    # t2 >= Pi * SF / (D * Yield)
    t2_mech_bearing = max_Pi * safety_factor / (D_2 * Yield_Stress_Material)
    
    # 2. Minimum t2 for Pull-through
    # sigma_vm = sqrt(3) * tau = sqrt(3) * Pull / (pi * D_in * t2) <= Yield / SF
    # t2 >= sqrt(3) * Pull * SF / (pi * D_in * Yield)
    t2_mech_pull = math.sqrt(3) * max_Pull * safety_factor / (math.pi * D_in * Yield_Stress_Material)
    
    # 3. Minimum t3 for Pull-through
    t3_mech_pull = math.sqrt(3) * max_Pull * safety_factor / (math.pi * D_in * Yield_Stress_Body)
    
    # Start with mechanical minimums
    t2 = max(t2_mech_bearing, t2_mech_pull, 0.0001)
    t3 = max(t3_mech_pull, 0.0001)
    
    # Iteratively adjust for Thermal
    for _ in range(20): # Max iterations
        # Calculate Thermal Loads with current t3
        delta_a = Compliance_parts(Modulus_Material, D_fo, D_in, t3)
        delta_b = Compliance_fastener(Modulus_Material, (math.pi*(D_in/2)**2), 0.03)
        psi = force_ratio(delta_a, delta_b)
        
        a_c = Thermal_Coeff_Material
        a_f = Thermal_Coeff_Fastener
        E_b = Modulus_Fastener
        
        max_t2_req_thermal = 0
        
        for f in fasteners:
            A_sw = f.area
            # Check both temperature extremes
            for T in T_operate:
                delta_T = T - T_ref
                F_t = (a_c - a_f) * delta_T * E_b * A_sw * (1 - psi)
                
                # Total load for bearing
                # Note: F_t is added to Pi_magnitude in original code?
                # Original: item.Pi_magnitude=(item.Pi_magnitude+F_t)
                # But F_t can be negative. Pi_magnitude is positive.
                # We should take the worst case magnitude.
                # If F_t is aligned or opposed? Original code just adds them.
                # Assuming worst case alignment or scalar addition as per original code logic.
                
                total_load = f.Pi_magnitude + F_t
                # Bearing check
                # t2 >= Total_Load * SF / (D * Yield)
                req = abs(total_load) * safety_factor / (D_2 * Yield_Stress_Material)
                if req > max_t2_req_thermal:
                    max_t2_req_thermal = req
        
        # Update t2
        new_t2 = max(t2_mech_bearing, t2_mech_pull, max_t2_req_thermal)
        
        # t3 is not affected by t2 in the formulas provided (only t3 affects compliance)
        # So t3 stays at t3_mech_pull (unless there's a t3 thermal check?)
        # Original code: fastn.MS_t3_bearing_thermal=t3/max(t3_2_list...) - 1
        # Wait, t3_2_list contains local_wall_thickness from thermal check.
        # local_wall_thickness = safety_factor*Pi_magnitude/(bearing_allowable_stress*Diameter)
        # This seems to be calculating t2 thickness?
        # "MS_t3_bearing_thermal" name suggests t3 check, but the formula uses t3 / max(t2_reqs).
        # If it compares t3 to t2_reqs, that's weird.
        # Let's look at original code:
        # t3_2_list.append(item.local_wall_thickness) -> This is t2 thickness required.
        # fastn.MS_t3_bearing_thermal=t3/max(t3_2_list[2*i],t3_2_list[2*i+1],t3_list[i]) - 1
        # It seems it's checking if t3 is greater than the required t2 thickness?
        # That implies t3 must be at least as thick as the required bearing thickness?
        # Or maybe it's a typo in the original code and it should be t2?
        # But MS_t2_bearing is calculated separately.
        # Let's assume the requirement is t3 >= required_bearing_thickness (maybe for the wall?).
        # If so, we must update t3 too.
        
        new_t3 = max(t3_mech_pull, max_t2_req_thermal) # Assuming t3 must also satisfy bearing thickness?
        # Actually, let's stick to the explicit checks.
        # If the user code checks t3 against "max(t3_2_list...)", and t3_2_list are bearing thicknesses...
        # Then t3 must be >= max bearing thickness.
        
        if abs(new_t2 - t2) < 1e-6 and abs(new_t3 - t3) < 1e-6:
            t2 = new_t2
            t3 = new_t3
            break
        t2 = new_t2
        t3 = new_t3
        
    return t2, t3

def optimize():
    # Discrete variables
    N_rows_options = [3, 4, 5, 6, 7, 8] # Number of rows (N_min)
    D_2_options = [0.003, 0.004, 0.005, 0.006, 0.008, 0.010] # Standard diameters
    
    best_result = None
    best_config = None
    
    print(f"{'N_rows':<10} {'D_2 (mm)':<10} {'t2 (mm)':<10} {'t3 (mm)':<10} {'Score':<10}")
    print("-" * 60)

    for N in N_rows_options:
        for D in D_2_options:
            try:
                t2, t3 = solve_thickness(N, D)
                
                # Calculate score (weight proxy)
                score = t2 + t3
                
                print(f"{N:<10} {D*1000:<10.1f} {t2*1000:<10.4f} {t3*1000:<10.4f} {score:<10.6f}")
                
                if best_result is None or score < best_result:
                    best_result = score
                    best_config = {
                        'N_min': N,
                        'D_2': D,
                        't2': t2,
                        't3': t3
                    }
            except Exception as e:
                # print(f"Failed for {N}, {D}: {e}")
                pass

    print("\nOptimization Complete.")
    if best_config:
        print("Best Configuration Found:")
        print(json.dumps(best_config, indent=4))
        
        # Calculate final margins for verification
        print("\nVerifying Margins for Best Config:")
        ms = calculate_margins([best_config['t2'], best_config['t3']], best_config['N_min'], best_config['D_2'])
        print(f"Minimum Margin of Safety: {ms}")
        
    else:
        print("No valid configuration found.")

if __name__ == "__main__":
    optimize()
