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

# Densities (kg/m^3) - Added for Mass Optimization
Densities = {
    "Aluminium": 2780,
    "Titanium": 4430,
    "Carbon Composite": 1600, # Approximate
    "Steel": 7850
}
Rho_Plate = Densities.get(material_used, 2780)
Rho_Fastener = Densities.get(material_used_fastener, 4430)
Rho_Body = Densities.get(material_used_body, 2780)

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
    # Determine minimum edge spacing based on material type
    if(Materials[material_used]['type (metal or composite)']) == 1 :
        edge_center_min_factor = np.array([2, 3]) # e1, e2
    elif (Materials[material_used]['type (metal or composite)']) == 2 :
        edge_center_min_factor = np.array([4, 5])
    
    edge_center_min = edge_center_min_factor * D_2
    min_edge = min(edge_center_min)
    
    center_center_min = 1.5 * D_2

    # Calculate how much space N_min fasteners need
    # Width needed = (N-1)*spacing + 2*edge_margin
    required_width = (N_min - 1) * center_center_min + 2 * min_edge
    
    if required_width <= w:
        # It fits in the original width
        # Center the pattern by increasing edge spacing
        N_out = N_min
        w_out = w
        # Calculate new edge spacing to center the pattern
        # Total span of fasteners = (N-1)*s
        # Remaining space = w - span
        # Edge spacing = Remaining space / 2
        edge_spacing = (w - (N_min - 1) * center_center_min) / 2
    else:
        # It does not fit, need to expand width
        N_out = N_min
        w_out = required_width
        edge_spacing = min_edge

    return N_out, edge_spacing, center_center_min, w_out, D_2

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

def calculate_margins(x, N_rows, D_2, D_in, D_fo):
    # x = [t2, t3]
    t2 = x[0]
    t3 = x[1]
    
    if t2 <= 1e-5 or t3 <= 1e-5:
        return -1e9 # Invalid
        
    # 1. Geometry & Fasteners
    NOF = Number_Of_Fasteners(w_original, D_2, N_rows)
    fasteners = Fasteners_location(NOF[0], NOF[1], NOF[2], NOF[3], h, t1, NOF[4])
    assign_fastener_forces(fasteners)
    
    margins = []
    
    # 2. Mechanical Checks
    for f in fasteners:
        # Bearing t2
        # Use D_in for bearing check as requested
        bearing_stress = f.Pi_magnitude / (D_in * t2)
        ms_bearing = Yield_Stress_Material / (safety_factor * bearing_stress) - 1
        margins.append(ms_bearing)
        
        # Pull-through t2
        area_shear_t2 = math.pi * D_in * t2
        shear_stress_t2 = f.p_pull / area_shear_t2
        sigma_vm_t2 = math.sqrt(3 * shear_stress_t2**2)
        ms_pull_t2 = Yield_Stress_Material / (safety_factor * sigma_vm_t2) - 1
        margins.append(ms_pull_t2)
        
        # Pull-through t3
        area_shear_t3 = math.pi * D_in * t3
        shear_stress_t3 = f.p_pull / area_shear_t3
        sigma_vm_t3 = math.sqrt(3 * shear_stress_t3**2)
        ms_pull_t3 = Yield_Stress_Body / (safety_factor * sigma_vm_t3) - 1
        margins.append(ms_pull_t3)

    # 3. Thermal Checks
    delta_a = Compliance_parts(Modulus_Material, D_fo, D_in, t3)
    delta_b = Compliance_fastener(Modulus_Material, (math.pi*(D_in/2)**2), 0.03)
    
    psi = force_ratio(delta_a, delta_b)
    
    a_c = Thermal_Coeff_Material
    a_f = Thermal_Coeff_Fastener
    E_b = Modulus_Fastener
    
    for f in fasteners:
        A_sw = f.area
        for T in T_operate:
            delta_T = T - T_ref
            F_t = (a_c - a_f) * delta_T * E_b * A_sw * (1 - psi)
            
            # Thermal Bearing
            total_load = f.Pi_magnitude + F_t
            bearing_stress_thermal = total_load / (D_in * t2)
            ms_bearing_thermal = Yield_Stress_Material / (safety_factor * abs(bearing_stress_thermal)) - 1
            margins.append(ms_bearing_thermal)
            
    return min(margins)

def solve_thickness(N_rows, D_2, D_in, D_fo):
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
    t2_mech_bearing = max_Pi * safety_factor / (D_in * Yield_Stress_Material)
    
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
                total_load = f.Pi_magnitude + F_t
                # Bearing check
                # t2 >= Total_Load * SF / (D * Yield)
                req = abs(total_load) * safety_factor / (D_in * Yield_Stress_Material)
                if req > max_t2_req_thermal:
                    max_t2_req_thermal = req
        
        # Update t2
        # Add 0.1% buffer to ensure MS >= 0 despite float errors
        buffer = 1.001
        new_t2 = max(t2_mech_bearing, t2_mech_pull, max_t2_req_thermal) * buffer
        
        new_t3 = max(t3_mech_pull, max_t2_req_thermal) * buffer
        
        if abs(new_t2 - t2) < 1e-6 and abs(new_t3 - t3) < 1e-6:
            t2 = new_t2
            t3 = new_t3
            break
        t2 = new_t2
        t3 = new_t3
        
    return t2, t3
     {right_mark}      (z={z:.4f})")

def plot_top_configurations(top_results):
    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import Rectangle, Circle
    except ImportError:
        print("Matplotlib not installed.")
        return

    n_res = min(len(top_results), 3)
    # Create figure with 2 rows (Front View, Side View)
    fig, axes = plt.subplots(2, n_res, figsize=(6 * n_res, 10))
    
    # Normalize axes to 2D array [row, col]
    if n_res == 1:
        axes = axes.reshape(2, 1)
    elif len(axes.shape) == 1: # Should not happen with 2 rows, but for safety
        axes = axes.reshape(2, -1)

    for i, res in enumerate(top_results):
        if i >= n_res: break
        
        cfg = res['config']
        # Extract geometry
        N_min = cfg['N_min']
        D_in_val = cfg['D_in']
        D_2_val = cfg['D_2']
        t2_val = cfg['t2']
        t3_val = cfg['t3']
        
        # Re-calculate layout
        NOF = Number_Of_Fasteners(w_original, D_2_val, N_min)
        # NOF returns: N_out, edge_spacing, center_center_min, w_out, D_2
        N_act, edge_spacing, s, w_new, _ = NOF
        
        fasteners = Fasteners_location(N_act, edge_spacing, s, w_new, h, t1, D_2_val)
        
        # --- Front View (X-Z) ---
        ax_front = axes[0, i]
        
        # Plate Dimensions
        # Fasteners are at +/- x_f.
        # Plate width = 2 * (abs(x_f) + edge_spacing)
        if fasteners:
            x_f = abs(fasteners[0].x_coord)
            plate_width = 2 * (x_f + edge_spacing)
            # Force plate height to be w_original as requested
            plate_height = w_original 
        else:
            plate_width = w_original
            plate_height = w_original

        # Draw Plate
        rect_plate = Rectangle((-plate_width/2, -plate_height/2), plate_width, plate_height, 
                               linewidth=2, edgecolor='black', facecolor='#e0e0e0', alpha=0.5, label='Base Plate')
        ax_front.add_patch(rect_plate)
        
        # Draw Flanges (Front View footprint)
        # Left Flange
        rect_fl_l = Rectangle((-h/2 - t1, -plate_height/2), t1, plate_height, 
                              linewidth=1, edgecolor='black', facecolor='#a0a0a0', linestyle='--', label='Flange')
        ax_front.add_patch(rect_fl_l)
        # Right Flange
        rect_fl_r = Rectangle((h/2, -plate_height/2), t1, plate_height, 
                              linewidth=1, edgecolor='black', facecolor='#a0a0a0', linestyle='--')
        ax_front.add_patch(rect_fl_r)

        # Draw Fasteners
        for f in fasteners:
            circle = Circle((f.x_coord, f.z_coord), f.Diameter/2, color='blue', alpha=0.8)
            ax_front.add_patch(circle)
            ax_front.plot(f.x_coord, f.z_coord, 'k+')

        ax_front.set_title(f"Rank {i+1}: Front View (X-Z)\nN={N_min}, D={D_in_val*1000:.1f}mm")
        ax_front.set_xlabel("X (m)")
        ax_front.set_ylabel("Z (m)")
        ax_front.axis('equal')
        ax_front.grid(True, alpha=0.3)
        # Legend only on first plot to save space
        if i == 0: ax_front.legend(loc='upper right', fontsize='x-small')

        # --- Side View (Y-Z) ---
        ax_side = axes[1, i]
        
        # Coordinate system for Side View: Y is horizontal (Thickness), Z is vertical.
        # Y=0 is interface between Wall and Plate.
        # Wall is Y < 0. Plate is Y > 0.
        
        # Draw Wall
        wall_vis_width = max(t3_val, 0.005)
        rect_wall = Rectangle((-wall_vis_width, -plate_height/2 - 0.002), wall_vis_width, plate_height + 0.004,
                              hatch='///', edgecolor='black', facecolor='none', label='Wall')
        ax_side.add_patch(rect_wall)
        
        # Draw Base Plate
        rect_plate_side = Rectangle((0, -plate_height/2), t2_val, plate_height,
                                    edgecolor='black', facecolor='#e0e0e0', label='Base Plate')
        ax_side.add_patch(rect_plate_side)
        
        # Draw Lug/Flange
        # Extends from Y = t2_val.
        # Length? If D_1 is 0, assume length = plate_height (w_new) for visualization.
        D_1_param = params['geometry']['D_1']
        lug_length = D_1_param if D_1_param > 0 else w_new
        lug_y_start = t2_val
        
        # Draw simple rectangle for lug profile
        rect_lug = Rectangle((lug_y_start, -plate_height/2), lug_length, plate_height,
                             edgecolor='black', facecolor='#a0a0a0', alpha=0.7, label='Lug')
        ax_side.add_patch(rect_lug)
        
        # Draw Pin Hole
        # Center at (lug_y_start + lug_length/2, 0) if D_1 not specified?
        pin_x = lug_y_start + (lug_length * 0.75) # Heuristic
        pin_z = 0
        pin_radius = w_new * 0.2 # Heuristic
        circle_pin = Circle((pin_x, pin_z), pin_radius, color='white', ec='red', linewidth=2)
        ax_side.add_patch(circle_pin)
        
        ax_side.set_title(f"Side View (Y-Z)\nt2={t2_val*1000:.1f}mm, t3={t3_val*1000:.1f}mm")
        ax_side.set_xlabel("Y (Thickness) (m)")
        ax_side.set_ylabel("Z (m)")
        ax_side.axis('equal')
        ax_side.grid(True, alpha=0.3)
        if i == 0: ax_side.legend(loc='upper right', fontsize='x-small')

    plt.suptitle("Top 3 Optimized Fastener Configurations - Detailed Views")
    plt.tight_layout()
    plt.show()

def optimize():
    # Discrete variables
    N_rows_options = [1, 2, 3, 4, 5, 6] # Number of rows (N_min) - Added 6
    # Assumed D options in meters
    D_options = [0.0001, 0.0002, 0.0003, 0.0004, 0.0005, 0.0006, 0.0008, 0.001, 0.0015, 0.002, 0.0025, 0.003, 0.004, 0.005, 0.006, 0.008, 0.010] 
    
    results = []
    
    print(f"{'N_rows':<8} {'D(mm)':<8} {'t2(mm)':<10} {'t3(mm)':<10} {'Mass(kg)':<10} {'Score(Mass)':<12}")
    print("-" * 70)

    # Check if w_original is very small, which might cause all configs to fail
    print(f"Constraint: Plate Height (w) fixed to {w_original*1000:.2f} mm from parameters.json")

    for N in N_rows_options:
        for D in D_options:
            try:
                # Assumption: D_2 (spacing) scales with D
                D_2 = D 
                # Assumption: D_in (fastener) is D
                D_in = D
                # Assumption: D_fo (outer) is 2.0 * D
                D_fo = 2.0 * D
                
                # Check geometric constraint FIRST
                # Re-run geometry to check if it fits in w_original
                NOF = Number_Of_Fasteners(w_original, D_2, N)
                # NOF returns: N_out, edge_spacing, center_center_min, w_out, D_2
                w_required = NOF[3]
                
                if w_required > w_original * 1.0001: # Tolerance for float comparison
                    # Skip if it requires more width than available in JSON
                    # print(f"Skipping N={N}, D={D}: Required w={w_required*1000:.1f}mm > {w_original*1000:.1f}mm")
                    continue

                t2, t3 = solve_thickness(N, D_2, D_in, D_fo)
                
                # --- Mass Calculation ---
                # 1. Calculate Plate Dimensions
                # Use w_original as the fixed height
                
                # Fastener X-position (distance from center)
                # x = h/2 + t1 + edge_spacing
                # Note: edge_spacing from NOF might be calculated based on w_required if we didn't enforce it.
                # But since we enforced w_required <= w_original, NOF returned w_out = w_original (or smaller? No, NOF logic sets w_out=w if fits)
                
                N_act, edge_spacing, s, w_out, _ = NOF
                
                x_pos = h/2 + t1 + edge_spacing
                
                # Plate Width (X-direction)
                plate_width_x = 2 * (x_pos + max(edge_spacing, D_fo/2))
                
                # Plate Height (Z-direction) is fixed to w_original
                plate_height_z = w_original
                
                area_plate = plate_width_x * plate_height_z
                
                # 2. Volumes
                vol_plate_t2 = area_plate * t2
                vol_plate_t3 = area_plate * t3 
                
                # Fastener Length: t2 + t3 + extra for nut/head (approx 1.5*D)
                bolt_length = t2 + t3 + 1.5 * D_in
                vol_bolt_total = N_act * 2 * (math.pi * (D_in/2)**2 * bolt_length) 
                
                vol_holes_t2 = N_act * 2 * (math.pi * (D_in/2)**2 * t2)
                vol_holes_t3 = N_act * 2 * (math.pi * (D_in/2)**2 * t3)
                
                # 3. Mass
                mass_plate_t2 = (vol_plate_t2 - vol_holes_t2) * Rho_Plate
                mass_plate_t3 = (vol_plate_t3 - vol_holes_t3) * Rho_Body
                mass_bolts = vol_bolt_total * Rho_Fastener
                
                total_mass = mass_plate_t2 + mass_plate_t3 + mass_bolts
                
                # Score is now Mass
                score = total_mass
                
                print(f"{N:<8} {D*1000:<8.1f} {t2*1000:<10.4f} {t3*1000:<10.4f} {total_mass:<10.6f} {score:<12.6f}")
                
                results.append({
                    'score': score,
                    'mass': total_mass,
                    'config': {
                        'N_min': N,
                        'D_2': D_2,
                        'D_in': D_in,
                        'D_fo': D_fo,
                        't2': t2,
                        't3': t3
                    }
                })
            except Exception as e:
                # print(f"Failed for {N}, {D}: {e}")
                pass

    print("\nOptimization Complete.")
    
    if results:
        # Sort by score (Mass)
        results.sort(key=lambda x: x['score'])
        
        best_res = results[0]
        best_mass = best_res['mass']
        
        print("\nComparison to Optimal Configuration (Rank 1 - Lightest):")
        print(f"{'Rank':<5} {'N':<5} {'D(mm)':<8} {'t2(mm)':<10} {'t3(mm)':<10} {'Mass(kg)':<10} {'Diff(kg)':<12}")
        print("-" * 80)
        
        for i, res in enumerate(results):
            cfg = res['config']
            diff_mass = res['mass'] - best_mass
            
            print(f"{i+1:<5} {cfg['N_min']:<5} {cfg['D_in']*1000:<8.1f} {cfg['t2']*1000:<10.4f} {cfg['t3']*1000:<10.4f} {res['mass']:<10.6f} {diff_mass:<+12.6f}")

        print("\nTop 3 Configurations:")
        for i, res in enumerate(results[:3]):
            print(f"\nRank {i+1} (Score: {res['score']:.6f}):")
            print(json.dumps(res['config'], indent=4))
            
            # Calculate final margins for verification
            cfg = res['config']
            NOF = Number_Of_Fasteners(w_original, cfg['D_2'], cfg['N_min'])
            fasteners = Fasteners_location(NOF[0], NOF[1], NOF[2], NOF[3], h, t1, NOF[4])
            
            ms = calculate_margins([cfg['t2'], cfg['t3']], cfg['N_min'], cfg['D_2'], cfg['D_in'], cfg['D_fo'])
            print(f"Minimum Margin of Safety: {ms}")
            
        plot_top_configurations(results[:3])
    else:
        print("No valid configuration found.")

if __name__ == "__main__":
    optimize()
