function rosenbrock(x0, rho, chi, gamma, sigma)
    % Starting point and relative f value
    f = @(x) rosenbrock_function(x);
    fprintf('x0: %f;\n', x0);
    fprintf('f(x0): %f;\n', f(x0));
    
    [xk, fk, execution_time, k] = nelder(f, x0, 2, rho, chi, gamma, sigma);
    print_results(xk, fk, execution_time, k, "Rosenbrock function 2D");
end