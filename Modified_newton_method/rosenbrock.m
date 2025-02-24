function rosenbrock(x0, kmax, tolgrad, c1, rho, btmax, delta, mode, opts,prec_choice)
    f = @(x) 100*(x(2)-x(1).^2).^2 + (1-x(1)).^2;
    grad_f = @(x) [-400*x(1)*(x(2) - x(1)^2) + 2*(x(1) - 1);
                   200*(x(2) - x(1)^2)];
    Hess_f = @(x) sparse([1200*x(1)^2 - 400*x(2) + 2, -400*x(1);
                   -400*x(1), 200]);

    [xk, fk, gradfk_norm, k, failure] = ...
    modified_newton_bcktrck(2, x0, f, grad_f, Hess_f, kmax, tolgrad, c1, rho, btmax, delta, mode, opts,prec_choice);
    print_results(xk, fk, gradfk_norm, k, kmax, "Rosenbrock function 2D", failure);
end