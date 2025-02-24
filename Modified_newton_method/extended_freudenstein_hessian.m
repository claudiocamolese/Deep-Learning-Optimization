function H = extended_freudenstein_hessian(x)
    n = size(x, 1); 

    % vector of principal diagonal and first upper diagonal
    vectDiag0 = zeros(1, n);
    vectDiag1 = zeros(1, n-1);

    for k = 1:n
        if mod(k, 2) == 1
                vectDiag0(k) = 2;
            if k < n
                vectDiag1(k) = 12*x(k+1)-16;
            end
        else 
            vectDiag0(k) = 30*x(k)^4 -80*x(k)^3 +12*x(k)^2 -240*x(k) +12*x(k-1) +12;
        end
    end
    H = sparse(1:n,1:n,vectDiag0,n,n) + sparse(2:n,1:n-1,vectDiag1,n,n) + sparse(1:n-1,2:n,vectDiag1,n,n);
    
end


