function [xk, fk, gradfk_norm, k, failure] = modified_newton_bcktrck(n, x0, f, gradf, Hessf, kmax, tolgrad, c1, rho, btmax, delta, mode, opts,prec_choice)

farmijo = @(fk, alpha, c1_gradfk_pk) fk + alpha * c1_gradfk_pk;

% Store all the values of backtracking iterations
btseq = zeros(1, kmax);
% Store all the values for f history plot
f_values_history = zeros(1,kmax);
% Sequence of xk (used for plots, to try it, add xseq to output parameters)
% xseq = zeros(n, kmax);

% Flags
failure = false;
bcktrck = false;
eigFound = false;

% Starting values
xk = x0;
fk = f(xk);
f_values_history(1) = fk;
gradfk = gradf(xk);
k = 0;
gradfk_norm = norm(gradfk);

% Stopping criteria
while k < kmax && gradfk_norm >= tolgrad
    % Compute hessian matrix
    Hk = Hessf(xk);
    Bk = Hk;
    
    % Try to compute the incomplete cholesky factorization
    try 
        L = ichol(Bk, struct('type', 'nofill'));
    % If it fails because the matrix is not positive definite, make it
    % positive definite
    catch ME
        if contains(ME.message, 'nonpositive pivot') || contains(ME.message, 'nonsingular')
            % If the eigs function do not converge to a minimum try to
            % increase the number of max iteration until an upper bound
            while opts.maxit <= 1e+5 && eigFound == false
                eigFound = true;
                % Compute smallest eigenvalue
                lambda_min = eigs(Hk,1,mode,opts); 
                % If eigs fails, it returns NaN so try it again increasing
                % maxIterations
                if isnan(lambda_min)
                    opts.maxit = opts.maxit * 10;
                    eigFound = false;
                end
            end
            % If eigs completely fails exit from the method
            if eigFound == false
                failure = true;
                disp("Error during eigs computation.")
                break;
            end

            if lambda_min <= 0
                % Correction of the Hessian matrix
                tau_k = max(0, delta - lambda_min);
                Bk = Hk + tau_k * speye(n);                
            end
            % Recompute incomplete cholesky factorization
            if prec_choice==1
                try
                    L = ichol(Bk, struct('type', 'nofill'));
                catch ME
                    failure = true;         
                    fprintf("Error during ichol computation: %s\n", ME.message);
                    break;
                end
            end
        else
            % If the error of ichol is not "the input matrix is not
            % positive definite" throw an error
            failure = true;         
            fprintf("Error during ichol computation: %s\n", ME.message);
            break;
        end
    end
    
    % If the user selects the preconditioning option
    if prec_choice==1 
        % compute precodnitioning
        tol = 1e-4;
        maxIter = 100;
        [pk,~] = pcg(Bk, -gradfk, tol, maxIter, L, L');
    else
        % compute pk without preconditioning
        pk = -Bk \ gradfk;   
    end
    
    % Compute parameters for Armijo condition for backtracking strategy
    alpha = 1;
    
    xnew = xk + alpha * pk;
    fnew = f(xnew);
    
    c1_gradfk_pk = c1 * gradfk' * pk;
    bt = 0;
   
    % Backtracking strategy
    while bt < btmax && fnew > farmijo(fk, alpha, c1_gradfk_pk)
        alpha = rho * alpha;  
        xnew = xk + alpha * pk;
        fnew = f(xnew);
        bt = bt + 1;
        % flag used to know if plot backtrack iterations
        bcktrck = true;
    end

    % If Armijo condition is not satisfied and bt index reaches the maximum
    % number of iterations exit from the method
    if bt == btmax && fnew > farmijo(fk, alpha, c1_gradfk_pk)
        failure = true;
        disp("Max backtrack iterations reached")
        break;
    end
    
    % Update all the parameters for the next iteration
    k = k + 1;
    xk = xnew;
    fk = fnew;
    f_values_history(k) = fk;
    gradfk = gradf(xk);
    gradfk_norm = norm(gradfk);
    btseq(k) = bt;
    % xseq(:, k) = xk;
end

% Reshape of the backtrack vector
btseq = btseq(1:k);
f_values_history = f_values_history(1:k);
% xseq = xseq(:, 1:k);

% If there are not errors and the backtracking strategy has been done al
% least one time plot backtrack values
if failure == false && bcktrck == true
    figure;
    bar(btseq);
    title('Backtracking iterations');
    xlabel('Iterations of modified method');
    ylabel('Iterations of backtracking strategy');
    grid on;
end

% Objective function with respect to the k-th iteration
% figure;
% semilogy(1:k, f_values_history, '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
% title(sprintf("Objective function value with n=%g", n));
% xlabel('Iterazioni');
% ylabel('f(x)');
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
% grid on;

end