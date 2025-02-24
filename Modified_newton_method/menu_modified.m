function [f, x0, grad, H, c1, mode, opts, title_function, choice, prec_choice] = menu_modified(h, type)

    % Menu choice of the function
    datasets = {'Extended Freudenstein and Roth function', 'Problem 75', 'Problem 76'};
    [choice, ok] = listdlg('PromptString', 'Select a function:', ...
                           'SelectionMode', 'single', ...
                           'ListString', datasets, ...
                           'ListSize', [250, 250], ...
                           'Name', 'Function selection', ...
                           'OKString', 'Select', ...
                           'CancelString', 'Cancel');
    
    if ~ok
        error('Abort operation');
    end
    
    % Menu choice of derivative method
    deriv_methods = {'Analytical derivatives', 'Finite differences'};
    [deriv_choice, flag] = listdlg('PromptString', 'Select a differentiation method:', ...
                                   'SelectionMode', 'single', ...
                                   'ListString', deriv_methods, ...
                                   'ListSize', [250, 250], ...
                                   'Name', 'Differentiation method', ...
                                   'OKString', 'Select', ...
                                   'CancelString', 'Cancel');
    
    if ~flag
        error('Abort operation');
    end

    % Menu choice of preconditioning method
    prec_methods = {'Use preconditioning', 'Do not use preconditioning'};
    [prec_choice, prec] = listdlg('PromptString', 'Select a preconditioning method:', ...
                                   'SelectionMode', 'single', ...
                                   'ListString', prec_methods, ...
                                   'ListSize', [250, 250], ...
                                   'Name', 'Preconditioning method', ...
                                   'OKString', 'Select', ...
                                   'CancelString', 'Cancel');

    if ~prec
        error('Abort operation');
    end
    
    % Choice of the function and its parameters
    switch choice
        case 1
            f = @(x) extended_freudenstein_newton(x);
            x0 = @(n) (60 * ones(n, 1)) .* (mod((1:n)', 2) == 0) + (90 * ones(n, 1)) .* (mod((1:n)', 2) ~= 0);
            c1=1e-6; 
            mode='smallestreal';
            opts.maxit=1000;
            opts.tol=1e-6;
            title_function = 'Extended Freudenstein and Roth function';

            if deriv_choice == 1
                grad= @(x) extended_freudenstein_grad(x);
                H = @(x) extended_freudenstein_hessian(x);
               
            else               
                grad=@(x) findiff_gradf(f, x, h, type);
                H = @(x) extended_freudenstein_findiff_hessian(f, x, sqrt(h));
            end

        case 2
            f = @(x) problem75_newton(x);
            x0 = @(n) [-1.2 * ones(n-1, 1); -1];
            c1=1e-4;
            mode="smallestabs";
            opts=struct();
            title_function = 'Problem 75';

            if deriv_choice == 1
                grad= @(x)problem75_grad(x);
                H = @(x) problem75_hessian(x);
            else
                mode = "smallestreal";
                opts.maxit=1000;
                opts.tol=1e-6;
                grad=@(x) findiff_gradf(f, x, h, type);
                H = @(x) problem75_findiff_hessian(f, x, sqrt(h));
            end

        case 3
            f = @(x) problem76_newton(x);
            x0 =@(n) 2* ones(n, 1);
            c1=1e-4;
            mode='smallestreal';
            opts.maxit=1000;
            opts.tol=1e-6;
            title_function = 'Problem 76';

            if deriv_choice == 1
                grad=@(x) problem76_grad(x);
                H = @(x) problem76_hessian(x);
            else
                grad=@(x)findiff_gradf(f, x, h, type);
                H = @(x) problem76_findiff_hessian(f, x, sqrt(h));
            end
    end
    
    disp('Function and differentiation method loaded!');
end
