#Bluebird Aerospace
#Flow calculations for feed system, Charley Z
#This section uses the Darcy-Weisbach equation for estimation.
#Should add in two-phase dynamics later on

#ALL UNITS METRIC
import math
import CoolProp.CoolProp as CP

g = 9.81 #m/s^2, Earth

pc = 15e5 #Pa (15 bars)

epsilon = 1.2e-6 #1.2 micrometers

rho_water = CP.PropsSI('D', 'T', 290, 'P', 101325, 'Water')


print("-- DARCY-WEISBACH CALCULATION TOOL FOR SYSTEM PRESSURE DROP -- \n")

# OXIDIZER
oxidizer = 'NitrousOxide'
T_init_ox = 273.15 #kelvin, initial temperature 
mdot_ox = 1.22055 #kg/s
rho_ox = CP.PropsSI('D', 'Q', 0, 'T', T_init_ox, oxidizer) #saturated liquid nitrous oxide density
mu_ox = 0.003237 #saturated liquid nitrous oxide dynamic viscosity, CoolProp doesn't work here for some reason
cv_ox = 32.3 #cv of valve: McMaster 47865K23
d_pipe_ox = 0.0127 #inner diameter of main ox line. This case, it is inner diameter of hose 4468K814, which is 1/2 inches
r_pipe_ox = d_pipe_ox / 2
L_ox = 0.8 #meters, feed line length of oxidizer
s_inj_ox = 0.2 #injector stiffness, oxidizer

print(f"OXIDIZER {oxidizer} \n")

dP_inj_ox = pc * s_inj_ox #∆p over injector, oxidizer
print(f"Pressure drop in injector, {oxidizer} side: {dP_inj_ox} Pa \n")

#standard gravity, oxidizer
sg_ox = rho_ox / rho_water
print(f"Standard gravity of {oxidizer}: {sg_ox} \n")

#volumetric flow rate, oxidizer
q_ox = mdot_ox / rho_ox
print(f"Volumetric flow rate of {oxidizer}: {q_ox} m^3/s \n")

#average velocity, oxidizer
v_ox = q_ox / (math.pi * (r_pipe_ox ** 2))
print(f"Average flow velocity of {oxidizer}: {v_ox} m/s \n")

#kinematic viscosity, oxidizer
nu_ox = mu_ox / rho_ox
print(f"Kinematic viscosity of {oxidizer}: {nu_ox} m^2/s \n")

#Reynolds number, oxidizer
Re_ox = (v_ox * d_pipe_ox) / (nu_ox)
print(f"Reynolds number of {oxidizer} flow: {Re_ox} \n")

# Pipe flow calcs
#Friction factor, oxidizer
f_ox = ( -1.8 * math.log( (6.9/Re_ox) +  (( epsilon / d_pipe_ox ) / 3.7) ** 1.11 ,10) ) ** -2
print(f"Friction factor of {oxidizer} feed lines: {f_ox} \n")

#head loss, oxidizer
#Darcy-Weisbach
hf_ox = f_ox * (L_ox / d_pipe_ox) * (v_ox ** 2) / (2 * g)
print(f"Head loss in {oxidizer} feed lines: {hf_ox} m \n")

#line pressure drop
dP_line_ox = hf_ox * rho_ox * g
print(f"Pressure drop in {oxidizer} line: {dP_line_ox} Pa \n")

# Valve flow calcs
gpm_ox = q_ox * 15850.323140625 #convert to gallons per minute to use the formula in NASA SP-125
dP_valve_ox = (sg_ox * (gpm_ox **2)) / (cv_ox ** 2)
print(f"Pressure drop over {oxidizer} main valve: {dP_valve_ox} Pa \n")

#Total oxidizer pressure drop
dP_total_ox = dP_inj_ox + dP_line_ox + dP_valve_ox
print(f"TOTAL PRESSURE DROP IN {oxidizer}: {dP_total_ox} Pa \n")


# FUEL
fuel = 'Ethanol'
T_init_fuel = 273.15 #kelvin, initial temperature 
mdot_fuel = 0.43591 #kg/s
rho_fuel = CP.PropsSI('D', 'Q', 0, 'T', T_init_fuel, fuel) #ethanol density
mu_fuel = CP.PropsSI('V', 'T', T_init_fuel, 'P', 101325, fuel) #ethanol dynamic viscosity
cv_fuel = 32.3 #cv of valve: McMaster 47865K23
d_pipe_fuel = 0.0127 #inner diameter of main fuel line. This case, it is inner diameter of hose 4468K814, which is 1/2 inches
r_pipe_fuel = d_pipe_fuel / 2
L_fuel = 0.8 #meters, feed line length of fuel
s_inj_fuel = 0.2 #injector stiffness, fuel

print(f"FUEL {fuel} \n")

dP_inj_fuel = pc * s_inj_fuel #∆p over injector, fuel
print(f"Pressure drop in injector, {fuel} side: {dP_inj_fuel} Pa \n")

#standard gravity, fuel
sg_fuel = rho_fuel / rho_water
print(f"Standard gravity of {fuel}: {sg_fuel} \n")

#volumetric flow rate, fuel
q_fuel = mdot_fuel / rho_fuel
print(f"Volumetric flow rate of {fuel}: {q_fuel} m^3/s \n")

#average velocity, fuel
v_fuel = q_fuel / (math.pi * (r_pipe_fuel ** 2))
print(f"Average flow velocity of {fuel}: {v_fuel} m/s \n")

#kinematic viscosity, fuel
nu_fuel = mu_fuel / rho_fuel
print(f"Kinematic viscosity of {fuel}: {nu_fuel} m^2/s \n")

#Reynolds number, fuel
Re_fuel = (v_fuel * d_pipe_fuel) / (nu_fuel)
print(f"Reynolds number of {fuel} flow: {Re_fuel} \n")

# Pipe flow calcs
#Friction factor, fuel
f_fuel = ( -1.8 * math.log( (6.9/Re_fuel) +  (( epsilon / d_pipe_fuel ) / 3.7) ** 1.11 ,10) ) ** -2
print(f"Friction factor of {fuel} feed lines: {f_fuel} \n")

#head loss, fuel
#Darcy-Weisbach
hf_fuel = f_fuel * (L_fuel / d_pipe_fuel) * (v_fuel ** 2) / (2 * g)
print(f"Head loss in {fuel} feed lines: {hf_fuel} m \n")

#line pressure drop
dP_line_fuel = hf_fuel * rho_fuel * g
print(f"Pressure drop in {fuel} line: {dP_line_fuel} Pa \n")

# Valve flow calcs
gpm_fuel = q_fuel * 15850.323140625 #convert to gallons per minute to use the formula in NASA SP-125
dP_valve_fuel = (sg_fuel * (gpm_fuel **2)) / (cv_fuel ** 2)
print(f"Pressure drop over {fuel} main valve: {dP_valve_fuel} Pa \n")

#Total fuel pressure drop
dP_total_fuel = dP_inj_fuel + dP_line_fuel + dP_valve_fuel
print(f"TOTAL PRESSURE DROP IN {fuel}: {dP_total_fuel} Pa")


