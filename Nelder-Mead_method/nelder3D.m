function [xk, fk, execution_time, k] = nelder3D(f, x, n, rho, chi, gamma, sigma)
    % Parameters
    max_iter = 500000; 
    tol = 1e-6;     
    k = 1;

    % Simplex generator
    simplex = generate_simplex(x, n);

    % check of the dimensions of input array x and value n
    if length(x) ~= n
        error("The dimension of the array and the number of dimensions n must coincide")
    end

    % Store the initial point as one of the points of the simplex (used in tuning)
    % if ~isempty(x)
    %     simplex(:, 1) = x;
    % end

    % Only for Rosenbrock function 3D plot
    if n == 2
        plot_nelder_mead_3d(f, [-2, 2], [-1, 3], [0, 2500]);
        hold on;
    end

    % Store all the f values here
    f_values_history = zeros(1,max_iter);   

    tic 

    % Stopping criteria
    while k <= max_iter && max(vecnorm(simplex - mean(simplex, 2), 2, 1)) >= tol
        f_values = f(simplex);
        [~, idx] = sort(f_values);
        % Sort simplex points based on f_values
        simplex = simplex(:, idx);
        
        % Store best f_value of the simplex (with lowest f_value) (with highest f_value)
        f_values_history(k)=f_values(1);

        % Compute the centroid of all the simplex points minus the worst
        centroid = mean(simplex(:, 1:end-1), 2);

        % Reflection
        x_r = centroid + rho * (centroid - simplex(:, end));
        f_r = f(x_r);

        % Update worst point with reflected point
        if f_r < f_values(end-1) && f_r >= f_values(1)
            simplex(:, end) = x_r;

        elseif f_r < f_values(1)
            % Expansion
            x_e = centroid + chi * (x_r - centroid);

            % Update worst point with expanded point
            if f(x_e) < f_r
                simplex(:, end) = x_e;
            else
                simplex(:, end) = x_r;
            end

        else
            if f_r < f_values(end)
                % Contraction
                x_c = centroid + gamma * (x_r - centroid); 
            else
                x_c = centroid + gamma * (simplex(:, end) - centroid);
            end

            % Update worst point with contracted point
            if f(x_c) < f_values(end)
                simplex(:, end) = x_c;
            else
                % Shrinkage
                simplex(:, 2:end) = simplex(:, 1) + sigma * (simplex(:, 2:end) - simplex(:, 1));
            end
        end

        % Only for rosenbrock function: plot simplex points
        if n == 2
            simplex_z = arrayfun(@(i) f(simplex(:, i)), 1:size(simplex, 2));
            plot3([simplex(1, :), simplex(1, 1)], ...
                  [simplex(2, :), simplex(2, 1)], ...
                  [simplex_z, simplex_z(1)], '-o', 'MarkerSize', 6);
        end

        k = k + 1;
    end

    % Final values
    xk = simplex(:, 1);
    fk = f(xk);
    execution_time=toc;
    f_values_history= f_values_history(:,1:k);

    % Objective function with respect to the k-th iteration
    figure;
    semilogy(1:k, f_values_history, '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
    title(sprintf("Objective function value with n=%g", n));
    xlabel('Iterazioni');
    ylabel('f(x)');
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    grid on;
end
