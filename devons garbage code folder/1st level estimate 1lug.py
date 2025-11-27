
Fx = 97.119 #N      #plus or minus
Fy = 97.119 #N      #plus or minus
Fz = -425 #N
Mx = -386.0275  #Nm #88.2132 to -386.0275
My = 88.2132 #Nm   #plus or minus
Mz = 1.8166 #Nm     #plus or minus


safety_factor = 2.5

# we assume
h = 0.04
t_1 = 0.0075
w = 0.04
d_pin = 0.005
area_pin = 3.1416 * (d_pin/2)**2
area_flange = 0.5*(w-d_pin) * t_1 # weakest point.
stress_concentration_factor = 3  # assumed for rectangular cutout with rounded edges

# material properties, assume aluminium 2219
yield_strength = 324e6 # Pa
shear_strength = 250e6 # Pa



# shear stress applied on pin
shear_stress_pin = safety_factor * 0.5*(Fx + Fz)/area_pin


# forces applied on each flange
sum_of_tension = 0.5 * Fy + Mz * 2 / (h+t_1)
sum_of_shear = 0.5*(Fx + Fz) + My*2/(h+t_1)

normal_stress_flange = safety_factor * stress_concentration_factor * sum_of_tension/area_flange
shear_stress_flange = safety_factor * stress_concentration_factor * sum_of_shear/area_flange

print("Sum of tension per flange:", sum_of_tension, "N")
print("Sum of shear per flange:", sum_of_shear, "N")

print("stress in flange:", normal_stress_flange , "Pa")
print("max allowed yield", yield_strength, "Pa")

print("shear stress in flange:", shear_stress_flange, "Pa")
print("max allowed shear", shear_strength, "Pa")



checkpin = shear_strength - shear_stress_pin
if checkpin > 0:
    print("Pin in shear: OK")
else:
    print("Pin in shear: NOT OK")


check1 = yield_strength - normal_stress_flange
if check1 > 0:
    print("Flange in tension: OK")
else:
    print("Flange in tension: NOT OK")

check2 = shear_strength - shear_stress_flange
if check2 > 0:
    print("Flange in shear: OK")
else:
    print("Flange in shear: NOT OK")