function F = problem75_newton(x)
    n = length(x);
    
    f_k = zeros(n, 1);
    f_k(1) = x(1) - 1;
    
    for i = 2:n
        f_k(i) = 10 * (i - 1) * (x(i) - x(i - 1))^2; % Caso 1 < k ≤ n
    end
  
    F = 0.5 * sum(f_k.^2);
end
