function print_results(xk, fk, execution_time, k, title)
    fprintf('**** NELDER MEAD METHOD: FINISHED POINT: %s\n', title);
    fprintf('**** NELDER MEAD METHOD: RESULTS *****\n');
    % fprintf('xk: %f;\n', xk);
    fprintf('f(xk): %g;\n', fk);
    fprintf('N. of Iterations: %d;\n', k);
    fprintf('Elapsed time is %g seconds.\n\n', execution_time);
end