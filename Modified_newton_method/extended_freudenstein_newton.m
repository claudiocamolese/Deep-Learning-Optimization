function F = extended_freudenstein_newton(x)
    [m, ~] = size(x);

    F = 0;
    for k = 1:m
        if mod(k, 2) == 1  
            if k == m
                fk = x(k) - 13;
            else
                fk = x(k) + ((5 - x(k+1)) * x(k+1) - 2) * x(k+1) - 13;
            end
        else  
            fk = x(k-1) + ((x(k) + 1) * x(k) - 14) * x(k) - 29;
        end
         F =F+ fk^2;
    end
    F=0.5*F;
end