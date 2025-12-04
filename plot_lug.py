import matplotlib.pyplot as plt
import matplotlib.patches as patches
import json
import numpy as np
import math

# Load parameters
param_file = 'parameters 6 bolts.json'
try:
    with open(param_file, 'r') as f:
        params = json.load(f)
except FileNotFoundError:
    # Fallback to 6 bolts if 8 not found
    param_file = 'parameters 6 bolts.json'
    with open(param_file, 'r') as f:
        params = json.load(f)

# Extract parameters
w = params['geometry']['w']
h_param = params['geometry']['h']
t1 = params['geometry']['t1']
t2 = params['geometry']['t2']
D_1 = params['geometry']['D_1']
D_2 = params['geometry']['D_2']
D_in = params['geometry']['D_in']
D_fo = params['geometry']['D_fo']
materials = params['materials']
material_used = params['material_selection']['material_used']

# Override h=w as in main.py
h = w 
# Note: In main.py, h is set to w. 
# However, looking at the image, h is the gap between ears. w is the height.
# If h=w, the gap is equal to the height.

# Fastener class and functions from main.py
class Fastener:
    def __init__(self, Diameter, x_coord, z_coord):
        self.Diameter = float(Diameter)
        self.x_coord = float(x_coord)
        self.z_coord = float(z_coord)

Fasteners = []

def Number_Of_Fasteners(w, D_in, N_min): 
    if(materials[material_used]['type (metal or composite)']) == 1 :
        edge_center_min = np.array([2, 3])*D_in
    elif (materials[material_used]['type (metal or composite)']) == 2 :
        edge_center_min = np.array([4, 5])*D_in
        
    edge_spacing_2 = max(edge_center_min)
    center_center_min = 1.5*D_in

    N_max_check = []
    for e in range(0, len(edge_center_min)):
        N_m = 1 + (w-2*edge_center_min[e])//(center_center_min)
        N_max_check.append(N_m)
        
        e_test = ( w - ( N_m - 1 ) * center_center_min) /2 
        if e_test > float(edge_center_min[e]):
            N_max = int(max(N_max_check))
            edge_spacing_1 = e_test
            
    if N_min != 0:
        if N_min > N_max:
            w = (N_min-1) * center_center_min + 2 * edge_spacing_1
        else:
            # Distribute fasteners evenly in the available range (keeping edge spacing fixed)
            edge_spacing_1 = max(edge_center_min)
            if N_min > 1:
                center_center_min = (w - 2 * edge_spacing_1) / (N_min - 1)
            else:
                # If only 1 fastener, center it (spacing doesn't matter, but edge spacing does)
                edge_spacing_1 = w / 2
        N_max = N_min
        edge_spacing_2 = max(edge_center_min)
    
    return N_max, edge_spacing_1, center_center_min, w, D_in, edge_spacing_2

def Fasteners_location(N_max, edge_spacing_1, center_center_min, w_new, h, t1, D_in, edge_spacing_2):
    Fasteners.clear()
    for f in range(N_max):
        diameter = D_in 
        # Positive x side
        x = h/2 + t1 + edge_spacing_2
        z = f * center_center_min + edge_spacing_1 - w_new/2
        Fasteners.append(Fastener(diameter, x , z))
        
        # Negative x side
        x = -h/2 - t1 - edge_spacing_2
        Fasteners.append(Fastener(diameter, x , z))
    return Fasteners

# Calculate Fasteners
N_min = 3 # Default from main.py
NOF = Number_Of_Fasteners(w, D_in, N_min)
# NOF returns: N_max, edge_spacing_1, center_center_min, w_new, D_2, edge_spacing_2
w_new = NOF[3]
edge_spacing_2 = NOF[5]
Fasteners_location(NOF[0], NOF[1], NOF[2], NOF[3], h, t1, NOF[4], NOF[5])

# Plotting
fig, ax = plt.subplots(figsize=(8, 6))

# Base Plate View (x-z plane)
ax.set_title("Base Plate View (Wall Attachment)")
ax.set_xlabel("X (m)")
ax.set_ylabel("Z (m)")
ax.axis('equal')

# Draw Fasteners
for fastn in Fasteners:
    circle = patches.Circle((fastn.x_coord, fastn.z_coord), fastn.Diameter/2, edgecolor='black', facecolor='none')
    ax.add_patch(circle)
    # Draw crosshair
    ax.plot([fastn.x_coord], [fastn.z_coord], 'k+')

# Draw Flanges
# Flange width in x: from h/2 + t1 to h/2 + t1 + 2*edge_spacing_2
# Flange height in z: w_new (from -w_new/2 to w_new/2)
flange_width = 2 * edge_spacing_2
flange_height = w_new

# Right Flange
right_flange_x = h/2 + t1
rect_right = patches.Rectangle((right_flange_x, -flange_height/2), flange_width, flange_height, linewidth=1, edgecolor='blue', facecolor='none', label='Flange')
ax.add_patch(rect_right)

# Left Flange
left_flange_x = -h/2 - t1 - flange_width
rect_left = patches.Rectangle((left_flange_x, -flange_height/2), flange_width, flange_height, linewidth=1, edgecolor='blue', facecolor='none')
ax.add_patch(rect_left)

# Draw Central Lug Block (Gap + Ears)
# From -h/2 - t1 to h/2 + t1
# Height w_new? Or w? The ears might have height w.
# Assuming ears have height w_new for consistency with flange.
central_width = h + 2*t1
rect_center = patches.Rectangle((-central_width/2, -flange_height/2), central_width, flange_height, linewidth=1, edgecolor='red', facecolor='none', linestyle='--', label='Lug Central Body')
ax.add_patch(rect_center)

# Annotations
# D_in label
# Move down below the fastener to avoid overlap
ax.text(Fasteners[0].x_coord, Fasteners[0].z_coord - 0.01, f'D_in={D_in*1000:.3f}mm', ha='center')

# h label (Gap)
# Draw dimension line for h
h_dim_y = -flange_height/2 - 0.005
# Extension lines for h
ax.plot([-h/2, -h/2], [-flange_height/2, h_dim_y], 'k-', linewidth=0.5)
ax.plot([h/2, h/2], [-flange_height/2, h_dim_y], 'k-', linewidth=0.5)
# Dimension line
ax.annotate('', xy=(-h/2, h_dim_y), xytext=(h/2, h_dim_y), arrowprops=dict(arrowstyle='<->', lw=0.5))
ax.text(0, h_dim_y - 0.002, f'h={h*1000:.1f}mm', ha='center', va='top', fontsize=9)

# w label (Height)
# Move to the left side to avoid overlap
w_dim_x = left_flange_x - 0.005 # Further reduced offset
# Extension lines for w
ax.plot([left_flange_x, w_dim_x], [-flange_height/2, -flange_height/2], 'k-', linewidth=0.5)
ax.plot([left_flange_x, w_dim_x], [flange_height/2, flange_height/2], 'k-', linewidth=0.5)
# Dimension line for w
ax.annotate('', xy=(w_dim_x, -flange_height/2), xytext=(w_dim_x, flange_height/2), arrowprops=dict(arrowstyle='<->', lw=0.5))
ax.text(w_dim_x - 0.002, 0, f'w={w_new*1000:.1f}mm', va='center', ha='right', rotation=90)

# Total width annotation
total_width = h + 2*t1 + 2*flange_width
dim_y = -flange_height/2 - 0.02 # Reduced offset

# Extension lines
ax.plot([left_flange_x, left_flange_x], [-flange_height/2, dim_y], 'k-', linewidth=0.5)
ax.plot([right_flange_x + flange_width, right_flange_x + flange_width], [-flange_height/2, dim_y], 'k-', linewidth=0.5)

# Dimension line
ax.annotate('', xy=(left_flange_x, dim_y), xytext=(right_flange_x + flange_width, dim_y), arrowprops=dict(arrowstyle='<->'))
ax.text(0, dim_y - 0.002, f'Total Width={total_width*1000:.001f}mm', ha='center', va='top')

# Adjust limits
ax.autoscale_view()
# Add some margin to the bottom for annotations
ylim = ax.get_ylim()
ax.set_ylim(dim_y - 0.015, ylim[1] + 0.01)
# Add margin to the left for w annotation
xlim = ax.get_xlim()
ax.set_xlim(w_dim_x - 0.01, xlim[1])

ax.grid(True, linestyle=':', alpha=0.6)

plt.tight_layout()
plt.savefig('lug_plot.png')
print("Plot saved to lug_plot.png")

