function F_val = problem76_newton(x)
    n = length(x); 
    fk = zeros(n, 1); 
    
    for k = 1:n-1
        fk(k) = x(k) - (x(k+1)^2) / 10;
    end

    fk(n) = x(n) - (x(1)^2) / 10;
    F_val = 0.5 * sum(fk.^2);
end
