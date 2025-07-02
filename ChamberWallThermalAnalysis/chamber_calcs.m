clear; clc;
%% ALL UNITS SI %%
%% === INPUTS === %%
% DE-1 (Cardinal) %
mdot = 1.61;    %mdot (kg/s)
%Chamber Dimensions%
r_i = 0.0664866;    %inner radius
t_w = 0.01;    %wall thickness
r_o = r_i + t_w;    %outer radius
A_loc = pi * (r_i)^2;
r_t = 0.0221622;    %throat radius
D_t = 2 * r_t;
A_t = pi * (r_t)^2;
L = 0.15;    %chamber length (use RPA)

% Chamber Material%
rho = 2700; %density
c_p = 896;  %specific heat (J/kg K)
k = 170;    %Thermal conductivity (W/m K)
sigma_yield = 276e6;    %yield strength, room temperature
sigma_ultimate = 310e6; %ultimate strength

P_chamber = 1.5e6;

T_gas = 1820;   %gas temperature (K)

c_star = P_chamber * A_t / mdot;

%Combustion Gas Properties (Take From RPA)%
mu = 7.359e-5;  %viscosity (kg/(m.s) or Pa.s)
c_p_gas = 1850; %specific heat, constant pressure
Pr = 0.5905;    %Prandtl number, effective
u = 700;    %gas velocity

%% === BARTZ HEAT TRANSFER COEFFICIENT === %%
L_c = t_w; %characteristic length
D = 2 * r_i;

h_bartz = 0.026 * (mu^0.2 * c_p_gas) / (Pr^0.6 * D_t^0.2) * (P_chamber / c_star)^0.8 * (D_t / (D_t/2))^0.1 * (A_t / (A_loc)) ^ 0.9;

h = h_bartz;

%% === BIOT NUMBER CHECK === %%
Bi = h * L_c / k;

fprintf('Biot Number Check:\n');
fprintf('Biot Number: %.3f\n',Bi);

if Bi < 0.1
    fprintf('Lumped mass model is valid.\n\n');
else
    fprintf('Lumped mass model is NOT valid.\n\n');
end

%% === LUMPED MASS HEAT TRANSFER MODEL === %%
A_inner = 2 * pi * r_i;
V_wall = pi * (r_o^2 - r_i^2) * L;
m_wall = V_wall*rho;
disp(m_wall);
C_th = m_wall*c_p;

dTwdt = @(t,Tw) h * A_inner * (T_gas-Tw) / (m_wall*c_p);

t_end = 2.6;    %total burn time (s)
Tw0 = 300;
[t,Tw]=ode45(dTwdt, [0,t_end],Tw0);

figure;
plot(t,Tw,'-o');
xlabel("Time(s)");
ylabel("Wall Temperature(K)");
title('Bartz Output');
grid on;

%% === HOOP STRESS === %%
ECF = 0.52;
sigma = (P_chamber * (r_o^2 + r_i^2) / (r_o^2 - r_i^2));
sf = (sigma_yield * ECF) / sigma;
disp(['Safety Factor:',num2str(sf)]);
