clc;
clear;

% Constants
L = 0.0084; % Define L as a constant
n_max = 10; % Maximum terms in the summation
accuracy = 100
x = linspace(0, L, accuracy); % Discretize x
alpha = 0.000111;
% Time vector for t
t = linspace(0, 10, accuracy ...
    );
sxx = -0.0565694276772

% Predefine 3D matrix for u(x, t)
u = zeros(length(x), length(t));

% Function definitions
phi = @(x) 293 - (290 + (x / L) * 2510);
A_t = @(t) 2800 - 2510 * exp(sxx * t);
B_t = @(t) 2800;
A_t_prime = @(t) -2510 * sxx * exp(sxx * t);
% Nested loops to calculate u(x, t) for all x and t
for i = 1:length(x)
    for j = 1:length(t)
        % Variables for current x and t
        x_curr = x(i);
        t_curr = t(j);
        % Summation components
        
        sum1 = 0;
        sum2 = 0;
        for n = 1:n_max
            
            % Compute Fourier coefficients and terms
            coeff_phi = (2 / L) * integral(@(x) phi(x) .* sin(n * pi * x / L), 0, L);
            coeff_A = (2 / L) * integral(@(x) (A_t_prime(t_curr) + (x / L) * (-A_t_prime(t_curr))) .* sin(n * pi * x / L), 0, L);

            % First summation term
            sum1 = sum1 + coeff_phi * sin(n * pi * x_curr / L) * exp(-alpha * n^2 * pi^2 * t_curr / L^2);

            % Second summation term
            sum2 = sum2 + (1/(alpha * n^2 * pi^2 / L^2)) * -1 * coeff_A * sin(n * pi * x_curr / L) * exp(-alpha * n^2 * pi^2 * (t_curr) / L^2) + coeff_A * sin(n * pi * x_curr / L) * (1/(alpha * n^2 * pi^2 / L^2));
        fprintf('Current x: %.2f, Current t: %.2f, Sum1: %.4f, Sum2: %.4f, coeff_phi: %.4f, coeff_A: %.4f\n', x_curr, t_curr, sum1, sum2, coeff_phi, coeff_A);
        
        end

        % Combine terms
        u(i, j) = sum1 + sum2 + (x_curr / L) * (B_t(t_curr) - A_t(t_curr)) + A_t(t_curr);
        if isinf(u(i, j)) || isnan(u(i, j))
            fprintf('u(%d, %d) = %f at t = %.2f, x = %.2f\n', i, j, u(i, j), t_curr, x_curr);
        end
    end
end

% Plotting the 3D surface
[X, T] = meshgrid(x, t);
figure;
surf(X, T, u');
axis([0 L 0 10 ])
shading interp;
colorbar;
title('3D Surface Plot of u(x, t)');
xlabel('x');
ylabel('t');
zlabel('u(x, t)');

