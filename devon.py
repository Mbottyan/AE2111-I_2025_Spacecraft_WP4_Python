import math


d1=0.03 # meter
yield_bolt = 240e6 # Pascal
static_friction_coefficient = 1.05


likely_arm_to_satelite_com = 1 # meter
area_solar_panel = 2 # meter^2
solar_sail = 4e-9 # pascal (its nanopascal)




pi = math.pi




mx_maybe = solar_sail * area_solar_panel * likely_arm_to_satelite_com



area_bolt_pin = (d1/2)*(d1/2)*pi
max_tension_bolt = area_bolt_pin * yield_bolt
friction_from_that = static_friction_coefficient * max_tension_bolt

print(friction_from_that)
print(mx_maybe)