clear
close all
clc

% Dimensions
n = 10;
% Iterations
N = 500;

rho_range = 0.5:0.1:1.3;
chi_range = 1.1:0.1:1.9;
gamma_range = 0.1:0.1:0.9;
sigma_range = 0.1:0.1:0.9;

steps = zeros(1, N);
res = zeros(1, N);
avg_steps = zeros(1, length(rho_range));
avg_res = zeros(1, length(rho_range));

% Select function
[f,choice]= menu();

% Parameters to check
rho = 1;
chi = 1.1; 
gamma = 0.7;
sigma = 0.6;

i = 1;

for rho = rho_range
    for j = 1:N
        [xk, fk, execution_time, k] = ...
            nelder(f, [], n, rho, chi, gamma, sigma);
        steps(j) = k;
        res(j) = fk;
    end
    avg_steps(i) = mean(steps);
    avg_res(i) = mean(res);
    i = i+1;
    disp(i);
end

% Display the average number of iteration with respect to parameter
figure;
plot(rho_range, avg_steps, '-o', 'LineWidth', 1.2);
title("Number of iterations");
xlabel('rho');
ylabel('avg steps');
grid on;

% Display the average value of f with respect to parameter
figure;
plot(rho_range, avg_res, '-o', 'LineWidth', 1.2);
title("Convergence");
xlabel('rho');
ylabel('avg f(xk)');
set(gca, 'YScale', 'log')
grid on;