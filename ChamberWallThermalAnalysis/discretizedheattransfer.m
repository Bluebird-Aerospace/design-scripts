clc
clearvars


%%%%%%%%%%%%%%%%%% SIM CONTROL VALUES %%%%%%%%%%%%%%

numpts = 200;

numpts_radial = 120; % number of radial annular sections

numpts_time = 30;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu = 0.0001063; % (kg/(m*s)) viscosity

mdot_total = 2.21; % (Kg/s)

%%%%%%%%%%%%%%%%%%%%%%% Structural Parameters

k_gr = 40; % W/(M*K)

% Aluminum 6061-T6 properties
k_al = 167;      % W/(m*K)
rho_al = 2700;   % kg/m^3
cp_al = 896;     % J/(kg*K)

% Copper properties
k_cu = 391; % W/(m*K)
rho_cu = 8960; % kg/m^3
cp_cu = 385; % J/(kg*K)

Tinit = 297; % (K)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Le = 6; % length of engine in inches
Spacing = 0.3; %space between each evaluation point in inches


x = 0:Spacing:Le;

interpts = linspace(0,Le,numpts); %

Diameter_Array = [0.721;0.634;0.546;0.457;
    0.372;0.46;0.633;0.807;0.98;1.134;1.256;
    1.351;1.42;1.465;1.486;1.488;1.488;1.488;
    1.488;1.488;1.488]*1.8; %include top of chamber as final measurement

outerdia = 4; % (in)
lower_outerdia = 4; % (in)
thick_length = 0.5; % optional length of a the thicker base section (as in a flange)

Spline_Dia = spline(x,Diameter_Array,interpts); % interpolate chamber points with spline

figure(1)
plot(x,Diameter_Array/2,'k-',LineWidth=1.5) % top side of chamber
hold on
plot(x,-Diameter_Array/2,'k-',LineWidth=1.5) % bottom side
plot(interpts,Spline_Dia/2,'b-',LineWidth=1.5) % interpolated top
plot(interpts,-Spline_Dia/2,'b-',LineWidth=1.5) % interpolated bottom
title('Thrust Chamber Cross Section (in)')
grid on
axis equal
hold off

% convert to metric
Metric_Dia = Spline_Dia/(39.3700787); % interpolated inner diameter (m)
x_metric = interpts/(39.3700787); % not interpolated (m)
Le_Metric = Le/(39.3700787); % chamber length, (m)

% Gas parameter interpolation

x_hga = [0,1.5,4,6]/39.3700787; % x values for hga
hga = [1520.68,2520.07,2809.6,2809.6];  % heat gas array (K)
Prga = [0.6277,0.5063,0.4486,0.4486]; % Prandtl number array
kga = [0.142,0.2955,0.401,0.401]; % thermal conductivity of the hot gas (W/(M*K))


% interpolation with pchip
hotgasarray = pchip(x_hga,hga,x_metric);
pr_gas_array = pchip(x_hga,Prga,x_metric);
k_array = pchip(x_hga,kga,x_metric);


%Initialize and solve for Nusselt # based HT coefficient approximation
hg_array_Nusselt = zeros(size(Metric_Dia));

% Dittus Boelter equation
for i = 1:numpts
    Prantl = pr_gas_array(i);
    k = k_array(i);
    hg_array_Nusselt(i) = (((0.023*k)/(Metric_Dia(i)^1.8))*(((4*mdot_total)/(pi*mu))^0.8))*(Prantl^0.34);
end

figure;
plot(x_metric,k_array)
drawnow;
xlabel('Axial position(m)', 'FontSize',16);
ylabel('Thermal Conductivity of the gas (W/(M*K)', 'FontSize',16);
title('Axial position vs thermal conductivity of the hot gas','FontSize',16)


figure(2)
plot(x_metric,hg_array_Nusselt,'r-',Linewidth=2) 
drawnow;
xlabel('Axial position (m)','FontSize',16)
ylabel('Convective Coefficient (W/(m^2K)','FontSize',16)        
grid on
title('Reaction Heat Transfer Coefficient vs Axial Position','FontSize',16)

% this is for a flange. sets the first 0.5 inches to lower_outerdia
Outer_Dia = ones(size(x_metric))*outerdia; % (in)
thick_param = (thick_length/6)*numpts;
for i = 1:thick_param
    Outer_Dia(i) = lower_outerdia;
end

Outer_Dia = Outer_Dia/39.3700787; % convert to meters

Tg_array = ones(size(x_metric))*400;


Tf = 5; % burn duration, seconds
timepoints = linspace(0,Tf,numpts_time);

Tdiskarray = zeros([numpts,numpts_time]);

rho_gr = 2260; % (kg/m3)
Cp_gr = 720; % specific heat (J/(Kg*K))
invtau_array = ones(size(x_metric));
biot_array = ones(size(x_metric));
for i = 1:numpts
    hotgas_temp = hotgasarray(i);
    hg = hg_array_Nusselt(i);
    india = Metric_Dia(i); % inner diam
    pakistan = Outer_Dia(i); % outer diam
    V = (Le_Metric/numpts)*(pi/4)*(pakistan^2 - india^2); % volume of the ring
    invtau = (hg*(pi*india*(Le_Metric/numpts)))/(rho_gr*Cp_gr*V); % 1 / tau, the time constant in heat transfer. jank stuff is area of cylinder wall.
    biot_array(i) = (hg*V)/(k_gr*(pi*india*(Le_Metric/numpts))); % biot number
    invtau_array(i) = invtau;
    for j = 1:numpts_time
        Tdiskarray(i,j) = hotgas_temp + (297 - hotgas_temp)*exp(-invtau*timepoints(j)); %solved newtons law of cooling

    end
end


%%%%%%%%%%%%%%%%% Lumped Parameter 'Disc'-retization %%%%%%%


VolRad = ones([numpts_radial 1]); % volume column vector for each radial segment (m^2)

InRadofslice = ones([numpts_radial numpts]); % inner radius for each shell (m)

biotslice = ones([numpts numpts_radial]); % biot number for each shell 
Aradinner_array = ones(size(biotslice)); % inner surface area of each shell (m^2)
Aradouter_array = ones(size(biotslice)); % outer surface area of each shell (m^2)
darr = ones([numpts 1]); % delta radius array, stores dR (m)
disp('Discretizing Nozzle...')
for k = 1:length(x_metric)
    dR = (Outer_Dia(k) - Metric_Dia(k))/(2*numpts_radial); % thickness of one radial slice (wall thickness / numpts_radial)
    darr(k) = dR;
    n = 0;
    for i = 1:numpts_radial

        %inner rad calcs
        inner_rad_slice = Metric_Dia(k)/2 + n*dR;
        InRadofslice(i,k) = inner_rad_slice;

        %volume calc
        Vrad = (pi)*((inner_rad_slice + dR)^2 - inner_rad_slice^2)*(Le_Metric/numpts);
        VolRad(i) = Vrad;

        %surface area calcs
        Aradinner = 2*pi*(inner_rad_slice)*(Le_Metric/numpts);
        Aradinner_array(k,i) = Aradinner;
        Aradouter = 2*pi*(inner_rad_slice + dR)*(Le_Metric/numpts);
        Aradouter_array(k,i) = Aradouter;

        %biot calc
        biotslice(k,i) = (hg_array_Nusselt(k)*Vrad)/(k_gr*Aradinner);

        % repeat for next layer
        n=n+1;
    end
    
end

%{ 
    AN EXPLANATION:
    First law of thermodynamics: change in internal energy (U) is equal to
    heat added minus work done:

    dU = Q - W

    Our wall is constant and incompressible (Volume doesn't change), so
    work done is zero. Therefore:

    dU = m * c * dT, where c is specific heat. Take the derivative
    of both sides with respect to time to get: 

    dU/dt = m * c * dT/dt
    
    m * c is also equal to the thermal capacitance C: density * specific heat * volume
    Therefore, dU/dt = C * dT/dt

    dU/dt = Heat flow in - heat flow out.
    
    This can be expressed by Fourier's law of conduction, which tells us that heat flow Q_dot = G * dT,
    where G is the thermal conductance, -k * A / dx (dR in our case), and
    dT is the temperature difference between the hotter inner shell T_(i-1)
    and the current shell T_i.

    For the first layer (called Gasvec in code), heat flow is convective. By Newton's law of cooling it is 
    equal to Q_gas = h * A * (T_gas - T_wall), where h is the convective heat 
    transfer coefficient(hg in our case) and A is area. Written exactly this is:

    C * dT/dt = G_(i−1,i) * ​(Ti−1​−Ti​) − G_(i,i+1​)(Ti​−Ti+1​) + Q_Gas

    If we isolate the dT/dt (and drop the subscript) we get:

    dT/dt = C^-1 (G * T + Q_gas) = C \ G * T + C \ Q_gas.

    This looks awfully similar to the single variable ODE dT/dt = α * T + f, where f is
    a constant forcing function. The solution to this is:

    T(t) = T_0 * e^αt + 1/α * (e^αt - 1) * f

    For us, in matrix exponential form that looks like:
    
    T(t) = T_init * e^(C \ G * t) + (C \ G) \ (e^(C \ G * t) - I) * C \ Q_gas

    %}


Teatime = ones(numpts,numpts_radial,numpts_time); % T over time for each shell (K)
disp('Iterating Heat Transfer...')
for k = 1:length(x_metric)
    hg = hg_array_Nusselt(k); % convective heat transfer coefficient (W/m^2 * K)
    Nmat = zeros(numpts_radial); % Inner wall conductance (W / K)
    Omat = zeros(numpts_radial); % Outer wall conductance (W / K)
    Cmat = zeros(numpts_radial); % Thermal capacitance (J/K)

    % Calculating the conductances using G = k*A/L

    for i = 2:numpts_radial % 1st layer is Gasvec
        Nmat(i,i-1) = ((Aradinner_array(k,i))*k_gr)/(darr(k)); % hotter inner shell conductance
        Nmat(i,i) = -((Aradinner_array(k,i))*k_gr)/darr(k); % current shell
    end

    for i = 1:(numpts_radial-1) %last layer not needed
        Omat(i,i+1) = (-Aradouter_array(k,i)*k_gr)/darr(k); % cooler outer shell conductance
        Omat(i,i) = (Aradouter_array(k,i)*k_gr)/darr(k); % current shell)
    end
    
    for i = 1:numpts_radial
        Cmat(i,i) = rho_gr*Cp_gr*VolRad(i); % thermal capacitance
    end
    
    Gasvec = zeros(numpts_radial,1); % heat flow into 1st layer (W)
    Gasvec(1) = (hg*Aradinner_array(k,1))*hotgasarray(k); % heat flow Q_dot = h*A*dT

    Nmat(1,1) = -hg*Aradinner_array(k,1);
    Gmat = Nmat - Omat; % conductance total equals in - out
    Gmat = Cmat \ Gmat; % inverse C * G  (1/s)
    Gasvec = Cmat \ Gasvec; % inverse C * Gasvec (K/s)

    for t = 1:numpts_time
        Gmatdt = (timepoints(t))*Gmat; % the exponentiated term
        Teatime(k,:,t) = (expm(Gmatdt))*(Tinit*ones([numpts_radial 1])) + Gmat \ (expm(Gmatdt) - eye(numpts_radial))*Gasvec; % see explanation
    end
end


Biot_max = max(biotslice,[],"all");
disp('Maximum Biot Number:')
disp(Biot_max)

%%%%%%%%%%%%%%%%%%%%

tevalindex = numpts_time; 

Tcut = squeeze(Teatime(:,:,tevalindex))';   

% not sure how this plotting method actually works its a template but it
% looks right
Rmat = InRadofslice(:,1:numpts)/Le_Metric;

[Xmat,~] = meshgrid(x_metric/Le_Metric, 1:numpts_radial);

figure(7)
surf(Xmat, Rmat, Tcut)
hold on
surf(Xmat,-Rmat,Tcut)
shading interp
view(2)
axis equal
colormap jet
xlabel('Axial position (L_e)')
ylabel('Radius (L_e)')
title(['Temperature distribution at t = ' num2str(timepoints(tevalindex)) ' s'])

%%%%%%%%%%

figure(8)

Tcut2 = squeeze(Teatime(:,1,:))';

mesh(x_metric/Le_Metric,timepoints,Tcut2,'linewidth',1.5);
shading interp;
colormap jet;
xlabel('Axial position (L_e)');
ylabel('Time (s)');
zlabel('Temperature (K)');
title('Transient Nozzle Wall Temperature (K)')