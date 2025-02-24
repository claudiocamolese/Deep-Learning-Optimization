function F = problem76(X)   
    [n, m] = size(X); 

    % F is a column vector whose values are f evaluations in each column of
    % the simplex
    F = zeros(m, 1); 
    
    % For each point of the simplex
    for i = 1:m
        x = X(:, i);
        fk = zeros(n, 1); 
        
        % Compute F for i-th point of the simplex
        for k = 1:n-1
            fk(k) = x(k) - (x(k+1)^2) / 10;
        end
        
        fk(n) = x(n) - (x(1)^2) / 10;
        F(i) = 0.5 * sum(fk.^2);
    end
end
