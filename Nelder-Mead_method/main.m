clc
close all
clear

% Set random seed to enhance reproducibility
rng(337517)

% Open menu configuration: choice of function
[f, choice]= menu();

dimensions=[10, 25, 50];
title_plot = ["Extended Freudenstein and Roth function", "Problem 75", "Problem 76"];

% Parameters of Nelder-Mead method
rho = 1;
chi = 1.1;
gamma = 0.7;
sigma = 0.6;

k = 1;

% Store execution times for time plot
execution_times = zeros(10, length(dimensions));

% Compute Nelder-Mead method for rosenbrock function with different
% starting point
rosenbrock([1.2;1.2], rho, chi, gamma, sigma)
rosenbrock([-1.2;1], rho, chi, gamma, sigma)

for n=dimensions
    % Starting point with respect to the choosen function
    if choice==1
        x0 = 60*ones(n, 1);
        x0(1:2:end) = 90;
    elseif choice==2
        x0 = -1.2*ones(n, 1);
        x0(end) = -1;
    else
        x0 = 2* ones(n, 1);
    end

    % Creation of the hypercube
    hypercube = x0 + (rand(n, 10) * 2 - 1);

    disp("-----------------")
    fprintf("Dimension n=%g\n",n)
    disp("-----------------")
    
    % Nelder-Mead method computed with starting point x0
    [xk, fk, execution_time, iter] = nelder(f,x0, n, rho, chi, gamma, sigma);
    print_results(xk, fk, execution_time, iter, title_plot(choice));
    
    % Nelder-Mead method computed for each column vector of the hypercube
    for j= 1:size(hypercube,2)
        point=hypercube(:,j);
        [xk, fk, execution_time, iter] = nelder(f, point, n, rho, chi, gamma, sigma);
        execution_times(j, k) = execution_time;
        print_results(xk, fk, execution_time, iter, title_plot(choice));
    end

    k=k+1;
end

% Plot execution time with respect to the dimensions
figure;
hold on;
grid on;
plot(dimensions, mean(execution_times, 1), '-o', 'LineWidth', 1.5);
title('Execution time');
xlabel('Dimension n');
ylabel('Time(s)');
set(gca, 'YScale', 'log');
grid on;
hold off;
