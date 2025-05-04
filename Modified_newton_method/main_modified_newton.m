clear
close all
clc

load("global_parameters.mat")
warning('off', 'all')

% Parameters
h=1e-2;
type = "fw";
rng(random_seed)

% Menu configuration
[f, x0, grad_f, Hess_f, c1, mode, opts, title_function, choice, prec_choice]= menu_modified(h, type);

% Modified newton method on rosenbrock function in two different starting point
rosenbrock([1.2; 1.2], 100, tolgrad, 1e-4, rho, btmax, delta, "smallestreal", struct(), prec_choice);
rosenbrock([-1.2; 1], 100, tolgrad, 1e-4, rho, btmax, delta, "smallestreal", struct(), prec_choice);

% n = 1000, n = 10000, n = 100000
for n = logspace(3, 5, 3) 
    disp('------------------')
    disp(['DIMENSION:',mat2str(n)])
    disp('------------------')

    % Starting point with respect to dimension n
    x0_val = x0(n);
    
    % Execution starting time
    tic
    % Modified newton method in starting point x0
    [xk, fk, gradfk_norm, k, failure] = ...
        modified_newton_bcktrck(n, x0_val, f, grad_f, Hess_f, kmax, tolgrad, c1, rho, btmax, delta, mode, opts,prec_choice);
    % Execution ending time
    toc
    print_results(xk, fk, gradfk_norm, k, kmax, title_function, failure);

    % Generate hypercube
    random_points = x0_val + (rand(n, 10) * 2 - 1);

    % Modified newton method for each column vector of the hypercube
    for i = 1:10
        tic
        [xk, fk, gradfk_norm, k, failure] = ...
            modified_newton_bcktrck(n, random_points(:,i), f, grad_f, Hess_f, kmax, tolgrad, c1, rho, btmax, delta, mode, opts,prec_choice);
        toc
        print_results(xk, fk, gradfk_norm, k, kmax, title_function, failure);
    end
end

if n == 2
    if choice == 1
        [X, Y] = meshgrid(linspace(-6200, 10, 1000), linspace(-60, 60, 1000));
    else
        [X, Y] = meshgrid(linspace(-6, 6, 100), linspace(-6, 6, 100));
    end
    Z = arrayfun(@(x, y) f([x; y]), X, Y);

    % Contour Plot
    figure;
    contour(X, Y, Z, 30); hold on;
    plot(xseq(1, :), xseq(2, :), '--*r'); 
    title('Contour Plot with trajectory');
    xlabel('x'); ylabel('y');
    hold off;

    % Surface Plot
    figure;
    surf(X, Y, Z, 'EdgeColor', 'none'); hold on;
    plot3(xseq(1, :), xseq(2, :), arrayfun(@(i) f(xseq(:, i)), 1:size(xseq, 2)), '--*r');
    title('Function with descent direction');
    xlabel('x'); ylabel('y'); zlabel('f(x)');
    hold off;
end