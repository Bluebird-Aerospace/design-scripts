# Pintle optimization calculator, modified from the work from PSAS

import math
import numpy as np
from matplotlib import pyplot as plt


# Desired flow parameters

of_ratio = 1.21134 / 0.43262

pintle_mdot = 1.21134  # kg/s  , desired Pintle mass flow rate
annulus_mdot = pintle_mdot / of_ratio  # kg/s  , desired Annulus mass flow rate

num_pintle_entrance = 1
num_annulus_entrance = 1
# Fluid properties

fluid_pintle = 'N2O'
fluid_annulus = 'Eth'

if fluid_pintle == 'water':
    pintle_rho = 999.7  # kg/m^3
    pintle_dynvisc = 0.001307  # N*s/m^2
    pintle_surface_tension = 0.0720  # N/m

elif fluid_pintle == 'Eth':
    pintle_rho = 797  # kg/m^3
    pintle_dynvisc = 0.0012318  # N*s/m^2
    pintle_surface_tension = 0.02189  # N/m

elif fluid_pintle == 'N2O':
    pintle_rho = 1002  # kg/m^3
    pintle_dynvisc = 0.003237  # N*s/m^2
    pintle_surface_tension = 0.02373  # N/m

if fluid_annulus == 'water':
    annulus_rho = 999.7  # kg/m^3
    annulus_dynvisc = 0.001307  # N*s/m^2
    annulus_surface_tension = 0.0720  # N/m

elif fluid_annulus == 'IPA':
    annulus_rho = 877.8  # kg/m^3
    annulus_dynvisc = 0.002  # N*s/m^2
    annulus_surface_tension = 0.0720 * 0.3 + 0.0209 * 0.7  # N/m 70% IPA estimate

elif fluid_annulus == 'Eth':
    annulus_rho = 797  # kg/m^3
    annulus_dynvisc = 0.0012318  # N*s/m^2
    annulus_surface_tension = 0.02189  # N/m

elif fluid_annulus == 'N2O':
    annulus_rho = 1002  # kg/m^3
    annulus_dynvisc = 0.003237  # N*s/m^2
    annulus_surface_tension = 0.02373  # N/m


# V_3.0 geometry

def inch_to_meters(inch):
    meter = inch / 39.36996
    return (meter)


pintle_entrance_dia = 0.0157988 #0.009906 for tubing
pintle_inner_dia = 0.01270669 #contracted from 1/2 inch tubing
pintle_exit_dia = 0.00149351
pintle_outer_dia = 0.02399905
pintle_number_outlets = 40
pintle_angle_deg = 45

#display(pintle_outer_dia)

annulus_entrance_dia = 0.0157988
annulus_chamber_dia = 0.04249502
annulus_channel_depth = 0.02456740
annulus_channel_width = 0.04249502
annulus_outer_dia = 0.02473185
annulus_hydrolic_dia = annulus_outer_dia - pintle_outer_dia
annular_gap = (annulus_outer_dia - pintle_outer_dia) / 2

# Area calculations

pintle_entrance_area = (math.pi * (pintle_entrance_dia / 2) ** 2) * num_pintle_entrance  # m^2 , flow area 1
pintle_inner_area = math.pi * (pintle_inner_dia / 2) ** 2  # m^2 , flow area 2
pintle_exit_area = pintle_number_outlets * math.pi * (
            pintle_exit_dia / 2) ** 2  # m^2 , flow area 3, 16 small diameter holes

annulus_entrance_area = (math.pi * (annulus_entrance_dia / 2) ** 2) # m^2 , flow area 1
annulus_manifold_area = annulus_entrance_area  # m^2 , flow area 2
annulus_channel_area = (annulus_channel_depth * annulus_channel_width) * num_annulus_entrance  # m^2 , flow area 3
annulus_chamber_area = math.pi * (annulus_chamber_dia / 2) ** 2  # m^2 , flow area 4
annulus_inner_area = math.pi * (pintle_outer_dia / 2) ** 2  # m^2
annulus_outer_area = math.pi * (annulus_outer_dia / 2) ** 2  # m^2
annulus_exit_area = annulus_outer_area - annulus_inner_area  # m^2 , flow area 5

#display(pintle_exit_area)
#display(annulus_exit_area)

# velocity calculations

pintle_entrance_velocity = pintle_mdot / (pintle_rho * pintle_entrance_area)
pintle_inner_velocity = pintle_mdot / (pintle_rho * pintle_inner_area)
pintle_exit_velocity = pintle_mdot / (pintle_rho * pintle_exit_area)

annulus_entrance_velocity = annulus_mdot / (annulus_rho * annulus_entrance_area)
annulus_manifold_velocity = annulus_mdot / (annulus_rho * annulus_manifold_area)
annulus_channel_velocity = annulus_mdot / (annulus_rho * annulus_channel_area)
annulus_chamber_velocity = annulus_mdot / (annulus_rho * annulus_chamber_area)
annulus_exit_velocity = annulus_mdot / (annulus_rho * annulus_exit_area)

# Entrance Dynamic Pressure

pintle_dynamic_pressure_in = (pintle_rho * pintle_entrance_velocity ** 2 / 2) * 0.1450377 / 1000 #Pa to psi
pintle_dynamic_pressure_out = (pintle_rho * pintle_exit_velocity ** 2 / 2) * 0.1450377 / 1000
annulus_dynamic_pressure_in = (annulus_rho * annulus_entrance_velocity ** 2 / 2) * 0.1450377 / 1000

# Reynolds numbers

pintle_re_1 = (pintle_rho * pintle_entrance_velocity * pintle_entrance_dia) / pintle_dynvisc
pintle_re_2 = (pintle_rho * pintle_inner_velocity * pintle_inner_dia) / pintle_dynvisc
pintle_re_3 = (pintle_rho * pintle_exit_velocity * pintle_exit_dia) / pintle_dynvisc

annulus_re_1 = (annulus_rho * annulus_entrance_velocity * annulus_entrance_dia) / annulus_dynvisc
annulus_re_3 = (annulus_rho * annulus_channel_velocity *
                ((4 * annulus_channel_area) / (2 * annulus_channel_depth + 2 * annulus_channel_width))
                ) / annulus_dynvisc
annulus_re_4 = (annulus_rho * annulus_chamber_velocity * annulus_chamber_dia) / annulus_dynvisc
annulus_re_5 = (annulus_rho * annulus_exit_velocity * annulus_hydrolic_dia) / annulus_dynvisc

# Loss coefficents; taken from Fluids textbook

pintle_k_1 = 0.07  # loss coefficient for first entrance region (from Munson)
pintle_k_2 = 2  # inner region for transition between centeral flow path and radial holes (guess based on Munson)
pintle_k_3 = 1  # Exit region to atmosphere (from Munson)

annulus_k_1 = 0.07  # Transition from entrance to manifold  (2.25 guess based on Munson)
annulus_k_2 = 4  # Transition between manifold and inner channels (2.8 guess based on Munson)
annulus_k_3 = 1.2  # Transition from channels into annular chamber (from Munson)
annulus_k_4 = 2  # Transition from upper fuel chamber diameter to annular gap (0.07 from Munson)
annulus_k_5 = 1.75  # exit region to atmosphere (1 from Munson)

# calculated head losses in kPa

pintle_hL_1 = (pintle_k_1 * 0.5 * pintle_rho * pintle_entrance_velocity ** 2) / 1000
pintle_hL_2 = (pintle_k_2 * 0.5 * pintle_rho * pintle_inner_velocity ** 2) / 1000
pintle_hL_3 = (pintle_k_3 * 0.5 * pintle_rho * pintle_exit_velocity ** 2) / 1000
pintle_hl_tot = pintle_hL_1 + pintle_hL_2 + pintle_hL_3

annulus_hL_1 = (annulus_k_1 * 0.5 * annulus_rho * annulus_entrance_velocity ** 2) / 1000
annulus_hL_2 = (annulus_k_2 * 0.5 * annulus_rho * annulus_manifold_velocity ** 2) / 1000
annulus_hL_3 = (annulus_k_3 * 0.5 * annulus_rho * annulus_channel_velocity ** 2) / 1000
annulus_hL_4 = (annulus_k_4 * 0.5 * annulus_rho * annulus_chamber_velocity ** 2) / 1000
annulus_hL_5 = (annulus_k_5 * 0.5 * annulus_rho * annulus_exit_velocity ** 2) / 1000
annulus_hl_tot = annulus_hL_1 + annulus_hL_2 + annulus_hL_3 + annulus_hL_4 + annulus_hL_5

# convert to psi from kPa for head loss

pintle_hL_1_psi = pintle_hL_1 * 0.1450377
pintle_hL_2_psi = pintle_hL_2 * 0.1450377
pintle_hL_3_psi = pintle_hL_3 * 0.1450377
pintle_hL_total_psi = (pintle_hL_1_psi + pintle_hL_2_psi + pintle_hL_3_psi)

annulus_hL_1_psi = annulus_hL_1 * 0.1450377
annulus_hL_2_psi = annulus_hL_2 * 0.1450377
annulus_hL_3_psi = annulus_hL_3 * 0.1450377
annulus_hL_4_psi = annulus_hL_4 * 0.1450377
annulus_hL_5_psi = annulus_hL_5 * 0.1450377
annulus_hL_total_psi = (annulus_hL_1_psi + annulus_hL_2_psi + annulus_hL_3_psi +
                        annulus_hL_4_psi + annulus_hL_5_psi)

# calculate overall loss coefficient
pintle_K_L = (2 * pintle_hl_tot * 1000) / (pintle_rho * pintle_exit_velocity ** 2)
pintle_K_L_alt = pintle_mdot / (pintle_exit_area * math.sqrt(2 * pintle_rho * pintle_hl_tot * 1000))
annulus_K_L = (2 * annulus_hl_tot * 1000) / (annulus_rho * annulus_exit_velocity ** 2)

# calulcate Velocity scale
pintle_velocity_scale = pintle_entrance_velocity * 0.1
annulus_velocity_scale = annulus_entrance_velocity * 0.1

# calculate Turbulence intensity

pintle_turbulence_intensity = math.sqrt((3 * pintle_velocity_scale ** 2) / (2 * pintle_entrance_velocity ** 2))
annulus_turbulence_intensity = math.sqrt((3 * annulus_velocity_scale ** 2) / (2 * annulus_entrance_velocity ** 2))

# calculate Turbulence length scale

pintle_turbulence_length = 0.38 * pintle_entrance_dia * 1000
annulus_turbulence_length = 0.38 * annulus_entrance_dia * 1000

# calculate exit weber numbers

pintle_exit_weber = (pintle_rho * pintle_exit_velocity ** 2 * pintle_exit_dia) / (pintle_surface_tension)
annulus_exit_weber = (annulus_rho * annulus_exit_velocity ** 2 * annular_gap) / (annulus_surface_tension)

# blockage factor

pintle_blockage_factor = (pintle_number_outlets * pintle_exit_dia) / (math.pi * pintle_outer_dia)

# Momentum Ratio

#pintle_force = pintle_rho * pintle_exit_velocity ** 2 * math.pi * (pintle_exit_dia / 2) ** 2
#annulus_force = annulus_rho * annulus_exit_velocity ** 2 * annular_gap * pintle_exit_dia

pintle_force2 = pintle_rho * pintle_exit_velocity ** 2 * math.pi * (pintle_exit_dia / 2) ** 2 * pintle_number_outlets
annulus_force2 = annulus_rho * annulus_exit_velocity ** 2 * annulus_exit_area

pintle_angle = pintle_angle_deg * (math.pi/180)

pintle_momentum = pintle_mdot * pintle_exit_velocity * math.cos(pintle_angle)
annulus_momentum = annulus_mdot * annulus_exit_velocity + pintle_mdot * pintle_exit_velocity * math.sin(pintle_angle)


force_ratio = pintle_force2 / annulus_force2


momentum_ratio = pintle_momentum / annulus_momentum

spray_angle = math.acos(1/(1+force_ratio)) * (180/math.pi)

print("Pintle V3: Calculated flow values\n")
print('Pintle calculations if using %s\n' % fluid_pintle)

print('Pintle entrance velocity        = %8.3f m/s' % pintle_entrance_velocity)
# print('Pintle inner velocity       = %8.3f m/s' % pintle_inner_velocity)
print('Pintle exit velocity            = %8.3f m/s \n' % pintle_exit_velocity)

print('Pintle entrance Reynolds #  = %8.0f' % pintle_re_1)
# print('Pintle Reynolds # 2         = %10.0f' % pintle_re_2)
print('Pintle exit Reynolds #      = %8.0f \n' % pintle_re_3)
print('Pintle headloss 1           = %8.3f psi' % pintle_hL_1_psi)
print('Pintle headloss 2           = %8.3f psi' % pintle_hL_2_psi)
print('Pintle headloss 3           = %8.3f psi\n' % pintle_hL_3_psi)

print('Pintle Total head loss           = %8.2f psi' % pintle_hL_total_psi)
print('Pintle  Total head loss          = %8.3f kPa' % pintle_hl_tot)
print('Pintle  KL total                 = %8.3f \n ' % pintle_K_L)
print('Pintle turbulence intensity      = %8.4f ' % pintle_turbulence_intensity)
print('Pintle turbulence length scale   = %8.3f mm ' % pintle_turbulence_length)
print('Pintle exit weber number         = %8.3f \n' % pintle_exit_weber)

print('Pintle exit diameter             = %8.7f m \n' % pintle_exit_dia)
print('Pintle exit area             = %8.7f m^2 \n' % pintle_exit_area)

print('Pintle Entrance Dynamic Pressure = %8.3f psi \n' % pintle_dynamic_pressure_in)
print('Pintle Exit Dynamic Pressure     = %8.3f psi \n' % pintle_dynamic_pressure_out)
print('Pintle blockage factor           = %8.3f \n' % pintle_blockage_factor)

print('Annulus calculations if using %s\n' % fluid_annulus)

print('annulus entrance velocity        = %8.3f m/s' % annulus_entrance_velocity)
# print('annulus inner velocity      = %8.3f m/s' % annulus_chamber_velocity)
print('annulus exit velocity            = %8.3f m/s \n' % annulus_exit_velocity)

print('annulus entrance Reynolds # = %8.0f' % annulus_re_1)
# print('annulus Reynolds # 2        = %10.3f' % annulus_re_2)
print('annulus exit Reynolds #     = %8.0f \n' % annulus_re_5)

print('annulus headloss 1          = %8.3f psi' % annulus_hL_1_psi)
print('annulus headloss 2          = %8.3f psi' % annulus_hL_2_psi)
print('annulus headloss 3          = %8.3f psi' % annulus_hL_3_psi)
print('annulus headloss 4          = %8.3f psi' % annulus_hL_4_psi)
print('annulus headloss 5          = %8.3f psi\n' % annulus_hL_5_psi)

print('annulus Total head loss           = %8.2f psi' % annulus_hL_total_psi)
print('annulus Total head loss           = %8.3f kPa' % annulus_hl_tot)
print('annulus KL total                  = %8.3f \n ' % annulus_K_L)
print('annulus turbulence intensity      = %8.4f ' % annulus_turbulence_intensity)
print('annulus turbulence length scale   = %8.3f mm' % annulus_turbulence_length)
print('annulus exit weber number         = %8.3f \n' % annulus_exit_weber)

print('annular gap                       = %8.7f m \n' % annular_gap)
print('annulus exit area                      = %8.7f m^2 \n' % annulus_exit_area)
print('annulus Entrance Dynamic Pressure = %8.3f \n' % annulus_dynamic_pressure_in)

print('Force Ratio (Fo/Pf)           = %8.3f \n' % force_ratio)
print('Momentum Ratio (Po/Pf)           = %8.3f \n' % momentum_ratio)
print('Pintle spray theoretical half angle = %8.3f deg' % spray_angle)

print(pintle_momentum)
print(annulus_momentum)
