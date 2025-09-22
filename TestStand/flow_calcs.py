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

print(f"OXIDIZER {oxidizer}")

dP_inj_ox = pc * s_inj_ox #∆p over injector, oxidizer
print(f"Pressure drop in injector, {oxidizer} side: {dP_inj_ox} Pa")

#standard gravity, oxidizer
sg_ox = rho_ox / rho_water
print(f"Standard gravity of {oxidizer}: {sg_ox}")

#volumetric flow rate, oxidizer
q_ox = mdot_ox / rho_ox
print(f"Volumetric flow rate of {oxidizer}: {q_ox} m^3/s")

#average velocity, oxidizer
v_ox = q_ox / (math.pi * (r_pipe_ox ** 2))
print(f"Average flow velocity of {oxidizer}: {v_ox} m/s")

#kinematic viscosity, oxidizer
nu_ox = mu_ox / rho_ox
print(f"Kinematic viscosity of {oxidizer}: {nu_ox} m^2/s")

#Reynolds number, oxidizer
Re_ox = (v_ox * d_pipe_ox) / (nu_ox)
print(f"Reynolds number of {oxidizer} flow: {Re_ox}")

# Pipe flow calcs
#Friction factor, oxidizer
f_ox = ( -1.8 * math.log( (6.9/Re_ox) +  (( epsilon / d_pipe_ox ) / 3.7) ** 1.11 ,10) ) ** -2
print(f"Friction factor of {oxidizer} feed lines: {f_ox}")

#head loss, oxidizer
#Darcy-Weisbach
hf_ox = f_ox * (L_ox / d_pipe_ox) * (v_ox ** 2) / (2 * g)
print(f"Head loss in {oxidizer} feed lines: {hf_ox} m")

#line pressure drop
dP_line_ox = hf_ox * rho_ox * g
print(f"Pressure drop in {oxidizer} line: {dP_line_ox} Pa")

# Valve flow calcs
gpm_ox = q_ox * 15850.323140625 #convert to gallons per minute to use the NASA formula
dP_valve_ox = (sg_ox * (gpm_ox **2)) / (cv_ox ** 2)

print(f"Pressure drop over {oxidizer} main valve: {dP_line_ox} Pa")

#Total oxidizer pressure drop
dP_total_ox = dP_inj_ox + dP_line_ox + dP_valve_ox
print(f"TOTAL PRESSURE DROP IN {oxidizer}: {dP_total_ox} Pa \n")
