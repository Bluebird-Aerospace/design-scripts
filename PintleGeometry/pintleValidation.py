#PINTLE GEOMETRY VALIDATION ALGORITHM
#Adapted from the Fenton project, United Kingdom

# Pintle injection geometry calculator using SPI, HEM & NHNE models
# Due to high vapour pressure of N20, models must be used to account accurately for pressure drop over the injector as a result of 'flashing' / caviation.
import math
import CoolProp.CoolProp as CP

print("IF NUMBER DISPLAYED BELOW IS 373.1242958476844, COOLPROP DATABASE IS CONNECTED:")
print(CP.PropsSI('T', 'P', 101325, 'Q', 0, 'Water'))  # Test to make sure CoolProp is working

# For N20 pintle:
p_inj = 2.079e+6
p_cham = 1.5e+6
mdot_N2O = 1.21134  # Desired mass flow rate for combustion chamber
Re_guess = 1e6  # Reynolds estimate
Cf = 0.0791 * math.pow(Re_guess, -0.25)
# Cd_N2O = (1 / math.sqrt((4 * Cf * 10) + 2.28)) - 0.1 #-0.1 to account for cold flow testing tolerances. Estimate of Cd from UCL ----- Make sure to update with actual orifice geometry to get more accurate result. Might not be L/D of 10.
Cd_N2O = 0.65
T = -15
# s_in = CP.PropsSI('S', 'P', p_inj, 'T', T + 273.15, 'NitrousOxide') # Entropy of saturated N20 at inlet (assumed 50 bar self pressurised) injection pressure and ambient temperature
s_in = CP.PropsSI('S', 'T', T + 273.15, 'Q', 0,
                  'NitrousOxide')  # Entropy of saturated N20 at outlet (assumed 20 bar chamber pressure)
rho_inl = CP.PropsSI('D', 'Q', 0, 'T', T + 273.15, 'NitrousOxide')  # Density of liquid N20 at the inlet

# SPI Model:
deltaP = p_inj - p_cham  # Pressure drop across injector
G_SPI = Cd_N2O * math.sqrt(2 * rho_inl * deltaP)  # Mass flux with SPI model
print("Incompressible model mass flux:", G_SPI, "kgm^-2s")

# HEM Model:
h_in = CP.PropsSI('H', 'T', T + 273.15, 'Q', 0, 'NitrousOxide')  # Enthalpy of saturated liquid N20 at inlet
h_out = CP.PropsSI('H', 'P', p_cham, 'S', s_in,
                   'NitrousOxide')  # Given that entropy equivalence is assumed in the HEM model, enthalpy at the outlet can be calculated by chamber pressure and entropy at the inlet
rho_out = CP.PropsSI('D', 'P', p_cham, 'S', s_in, 'NitrousOxide')  # Density of saturated N20 at outlet
G_HEM = Cd_N2O * rho_out * math.sqrt(2 * (
            h_in - h_out))  # Mass flux with HEM model. Simplified equation can be used as upstream vapour quality can be assumed to be liquid.
A_N20 = mdot_N2O / G_HEM
A_N20 *= 10 ** 6
print("Area of injector using HEM model:", A_N20, "mm^2 at", T, "°C")

# NHNE Model:
L = 0.03241
p_sat = CP.PropsSI("P", "T", T + 273.15, "Q", 0,
                   'NitrousOxide')  # Saturated vapour pressure at specified temperature. Should be around 50 bar but this is more accurate and can be set on the day.
tau_b = math.sqrt(3 / 2 * rho_inl / (
            p_inj - p_cham))  # Characteristic Bubble Growth Rate (Make sure that p_inj reflects the actual saturated condition and not the estimate (vapour pressure))
tau_r = L * math.sqrt(rho_inl / 2 * deltaP)  # Fluid residence time
kappa = 1 / (1 + tau_b / tau_r)  # non-equilibrium factor kappa
print("Kappa: ",kappa)

G_NHNEest = (G_SPI + G_HEM) / 2  # Under-estimate of NHNE mass flux
inj_totalarea = mdot_N2O / G_NHNEest  # Therefore over-estimate of corresponding injector area

mdot_NHNE = 0
tolerance = 1e-6
step_size = inj_totalarea * 0.00001  # step by 0.001% of injector area estimate

while abs(mdot_N2O - mdot_NHNE) >= tolerance:  # Taking under-estimate of area and constantly increasing it by a small amount until it reaches a value that corresponds to the correct mdot. Code is written to prevent runaway bug that occurs when mdot_NHNE continues to increase
    inj_totalarea += step_size
    G_NHNE = (((1-(1 / (1 + kappa))) * G_SPI) + ((1 / (1 + kappa)) * G_HEM))
    mdot_NHNE = G_NHNE * inj_totalarea
    print(mdot_NHNE)
    if mdot_NHNE >= mdot_N2O:
        break

inj_totalarea = (mdot_N2O / G_NHNE) * 10 ** 6  # Final injector area
print("FINAL injector area:", inj_totalarea, "mm^2")
num_hole = 7
hole_area = inj_totalarea / num_hole
hole_radius = round(math.sqrt(hole_area / math.pi), 2)
print("Radius of pintle N20 holes:", hole_radius, "mm,", "so use", round(hole_radius, 1), "mm")

# For IPA annulus
# Using Ethanol as approximation for IPA in CoolProp database
mdot_IPA = 0.43262  # Mass flow rate of IPA from chamber code / NASA CEA
p_in_IPA = 1.9e+6
s_in_IPA = CP.PropsSI('S', 'P', p_in_IPA, 'T', T + 20 + 273.15, 'Ethanol')  # Entropy of Ethanol at inlet
rho_in_IPA = CP.PropsSI('D', 'P', p_in_IPA, 'S', s_in_IPA, 'Ethanol')  # Density of liquid N20 at the inlet
s_out_IPA = CP.PropsSI('S', 'P', 2e+6, 'T', T + 20 + 273.15,
                       'Ethanol')  # Entropy of Ethanol at outlet assuming temperature hasn't increased yet
rho_out_IPA = CP.PropsSI('D', 'P', 2e+6, 'S', s_out_IPA,
                         'Ethanol')  # Density of liquid N20 at the outlet for calculating spray angle
# mu_IPA = 2.1e-3 #Lookup table states 2.1mPA at 25°C
# D_IPA =
# V_IPA = mdot_IPA / rho_in_IPA * IPA_area
# L_IPA =
# Re_IPA = (rho_in_IPA * V_IPA* D_IPA) / mu_IPA
# Cf_IPA = 0.0791 * math.pow(Re_IPA, -0.25)
# Cd_IPA = 1 / math.sqrt((4 * Cf_IPA * L_IPA / D_IPA) + 2.28)
Cd_IPA = 0.75 - 0.1  # Crude estimate - 0.1 to give space for adjustment

# SPI Model:
# IPA doesn't have a high vapour pressure like N20 so incompressible model can be assumed to be accurate
Annulus_deltaP = p_in_IPA - p_cham  # Pressure drop across injector
G_SPI_IPA = Cd_IPA * math.sqrt(2 * rho_in_IPA * Annulus_deltaP)  # Mass flux with SPI model
A_IPA = mdot_IPA / G_SPI_IPA
print("Area of annulus:", A_IPA * 10 ** 6, "mm^2")
circumference_N2O = ((
                                 hole_radius / 10 ** 3) * 2 * 30)  # Circumference is set to account for size N20 holes. Thickness of pintle will account for any overlap in holes (in layers)
radius_N2O = circumference_N2O / (2 * math.pi)
annular_radius = math.sqrt((A_IPA + (math.pi * radius_N2O ** 2)) / math.pi)
annular_width = annular_radius - radius_N2O
print("Annular radius (From centre):", annular_radius * 10 ** 3, "mm")
print("Annular width:", annular_width * 10 ** 3, "mm")
print("Radius of pintle:", radius_N2O * 10 ** 3, "mm")

# Spray angle
v_IPA = Cd_IPA * math.sqrt((2 * Annulus_deltaP) / rho_out_IPA)  # Velocity of IPA exiting annulus using Bernoulli's
v_N2O = Cd_N2O * math.sqrt((2 * deltaP) / rho_out)  # Velocity of N2O exiting pintle using Bernoulli's
TMR = (mdot_N2O * v_N2O) / (mdot_IPA * v_IPA)
rad = math.acos(1 / (1 + TMR))
θ = rad * 180 / math.pi
print(mdot_N2O * v_N2O)
print(mdot_IPA * v_IPA)
print(v_IPA,"m/s")
print(v_N2O,"m/s")
print("Spray angle of atomised mixture:", θ, "°")