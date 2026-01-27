# Calculate expected pressure and temperature of nitrous tank when giving a delta mass
# Isentropic cooling assumption

import CoolProp.CoolProp as CP

def print_full_properties(dm):
    fluid = 'NitrousOxide'
    V = 0.0058  # m^3
    m1 = 5.0    # kg
    T1 = 262.0  # K
    
    rho1 = m1 / V
    p1 = CP.PropsSI('P', 'T', T1, 'D', rho1, fluid)
    s1 = CP.PropsSI('S', 'T', T1, 'D', rho1, fluid)
    x1 = CP.PropsSI('Q', 'T', T1, 'D', rho1, fluid)
    
    # isentropic expansion
    m2 = m1 - dm
    rho2 = m2 / V
    t2 = CP.PropsSI('T', 'S', s1, 'D', rho2, fluid)
    p2 = CP.PropsSI('P', 'S', s1, 'D', rho2, fluid)
    x2 = CP.PropsSI('Q', 'S', s1, 'D', rho2, fluid)

    print("INITIAL:")
    print(f"Mass: {m1}kg")
    print(f"Pressure: {p1/1e5:.2f} bar")
    print(f"Quality: {x1:.4f}")
    print(f"Entropy: {s1/1000:.4f} kJ/kg K")
    
    print("FINAL:")
    print(f"Mass: {m1-dm}kg")
    print(f"Final temperature: {t2:.2f} K (dropped {T1-t2:.2f} K)")
    print(f"Final pressure: {p2/1e5:.2f} bar")
    print(f"Final quality: {x2:.4f}")

print_full_properties(1.15*2.5)