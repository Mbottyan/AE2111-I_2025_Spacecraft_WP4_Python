import numpy as np
import json




# are knowns
net_force = 436 # newtons (aka P)
F_y = 97.119 # newtons
safety_margin = 3





# pin
pin_materials = [
    'titanium',
]
tau_yield_pin = [
    550e6,  # titanium
]
density_pin = [
    4420,  # titanium
]



pin_material_count = len(pin_materials)




# bearing
bearing_materials = [
    "aluminum_alloy_transverse_grain",
]
bearing_yield = [
    276e6,  # aluminum_alloy_transverse_grain
]
bearing_materials_count = len(bearing_materials)




# lug
lug_materials = [
    "aluminum_7075",
    "steel",
    "titanium",
]
lug_yield = [
    503e6,  # aluminum_7075
    250e6,  # steel
    880e6,  # titanium
]
density_lug = [
    2810,  # aluminum_7075
    7850,  # steel
    4420,  # titanium
]
lug_materials_count = len(lug_materials)









# iterate over all pin materials to get D for each material
w_over_d = [ # w/d is just 2* e/d
    2.0,
    2.4,
    2.8,
    3.2,
    3.6,
    4.0
]


# lines of data as a fx of w/d from fig D1.13
lug_design_data = np.array([
    [0.90, 1.12, 1.27, 1.40, 1.47, 1.60], # curve 1 from fig D1.13
    [0.83, 0.92, 1.07, 1.09, 1.11, 1.13], # curve 2 from fig D1.13
])
number_of_curves = lug_design_data.shape[0]






def do_analysis(bearing_material_id, pin_material_id, lug_material_id, w_d_id, curve_id):
    # print the name of the lug material
    
    print("")
    print("")
    print("--------------------------------")
    print("Lug Material:", lug_materials[lug_material_id])
    print("Pin Material:", pin_materials[pin_material_id])
    print("Bearing Material:", bearing_materials[bearing_material_id])
    print("w/d: ", w_over_d[w_d_id])
    print("Chosen curve:", curve_id + 1)
    print("--------------------------------")

    tau_pin = tau_yield_pin[pin_material_id]
    print("tau pin yield: ", tau_pin)
    yield_lug = lug_yield[lug_material_id]
    print("yield lug: ", yield_lug)
    yield_bearing = bearing_yield[bearing_material_id]
    print("yield bearing: ", yield_bearing)
    w_over_d_value = w_over_d[w_d_id]
    print("w/d value: ", w_over_d_value)
    k_t_specific = lug_design_data[curve_id, w_d_id]
    print("k_t at given w/d: ", k_t_specific)
    
    # calculate D based on yield strength
    D = np.sqrt (
        (4 * net_force) /
        (np.pi * tau_yield_pin[pin_material_id]/safety_margin)
    )
    print("D: ", D)
    
    W = w_over_d_value * D
    print("W: ", W)
    
    A_br = net_force / (k_t_specific * yield_bearing/safety_margin)
    print("A_br: ", A_br)
    
    A_lug = F_y / (yield_lug/safety_margin)
    
    
    t = A_lug / W
    print("t: ", t)
    
    
    volume_lug = 2 * (( np.pi * (0.5*W)**2 ) - ( np.pi * (0.5*D)**2 )) * t
    volume_pin = 2 * np.pi * (0.5*D)**2 * t
    
    mass_lug = volume_lug * density_lug[lug_material_id]    
    mass_pin = volume_pin * density_pin[pin_material_id]
    
    mass_total = mass_lug + mass_pin
    mass_total_milligram = mass_total * 1e6
    
    
    data = {
        "Lug Material": lug_materials[lug_material_id],
        "Pin Material": pin_materials[pin_material_id],
        "Bearing Material": bearing_materials[bearing_material_id],
        "w/d": w_over_d[w_d_id],
        "Chosen curve": curve_id + 1,
        "tau pin yield": tau_pin,
        "yield lug": yield_lug,
        "yield bearing": yield_bearing,
        "w/d value": w_over_d_value,
        "k_t at given w/d": k_t_specific,
        "D": D,
        "W": W,
        "A_br": A_br,
        "A_lug": A_lug,
        "t": t,
        "volume_lug": volume_lug,
        "volume_pin": volume_pin,
        "mass_lug": mass_lug,
        "mass_pin": mass_pin,
        "mass_total": mass_total,
        "mass_total_milligram": mass_total_milligram
    }
    
    # add to array
    data_list.append(data)
    
    return

data_list = []
# get t as a function of D for each bearing material
for bearing_material_id in range(bearing_materials_count):
    for pin_material_id in range(pin_material_count):
        for lug_material_id in range(lug_materials_count):
            for w_d_id in range(len(w_over_d)):
                for curve_id in range(number_of_curves):
                    do_analysis(bearing_material_id, pin_material_id, lug_material_id, w_d_id, curve_id)





with open("output.json", "w") as f:
    json.dump(data_list, f, indent=4)




# calculate t as it is constraint by D and the bearing specs
