#flow rate converter
L_per_s = 0.548  # L/s flow rate actual
P_atm_psi = 14.696 #atmosphere
P_upstream_psi = 2000
P_downstream_psi = 400 
T_actual_K = 298  # operating temperature, Kelvin
T_std_K = 273.15

P_actual_psi = (P_upstream_psi + P_downstream_psi) / 2

L_per_min = L_per_s * 60

stdL_min = L_per_min * (P_atm_psi / P_actual_psi) * (T_actual_K / T_std_K)

print(stdL_min)
