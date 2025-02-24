function print_results(xk, fk, gradfk_norm, k, kmax, title, failure)
    fprintf('**** MODIFIED NEWTON METHOD: FINISHED POINT: %s\n', title);
    fprintf('**** MODIFIED NEWTON METHOD: RESULTS *****\n');

    % Print results only if there are no errors in modified newton method
    if failure == false
        % fprintf('xk: %s;\n', mat2str(xk));
        fprintf('f(xk): %g;\n', fk);
        fprintf('Gradient norm: %g;\n', gradfk_norm);
        fprintf('N. of Iterations: %d/%d;\n\n', k, kmax);
    end
end