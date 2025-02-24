function simplex = generate_simplex(x, n)
    simplex = rand(n, n+1);
    simplex(:, 1) = x;
    
    % Check if points are affinely independent, if the rank of the 
    % differential matrix is not n, it regenerate the simplex
    while rank(simplex(:, 2:end) - simplex(:, 1)) < n
        simplex = rand(n, n+1);
        simplex(:, 1) = x;
    end
end