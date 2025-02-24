function [H] = problem76_findiff_hessian(f, x, h)
n = length(x);

% vector of principal diagonal and first upper diagonal
Hessf0 = zeros(1,n);
Hessf1= zeros(1,n-1);

for j=1:n
    % Elements on the principal diagonal
    xh_plus = x;
    xh_minus = x;
    xh_plus(j) = xh_plus(j) + h;
    xh_minus(j) = xh_minus(j) - h;
    Hessf0(j) = (f(xh_plus) - 2*f(x) + f(xh_minus))/(h^2);
   
    if j < n
        % Elements on the first upper diagonal
        i=j+1;
        xh_plus_ij = x;
        xh_plus_ij([i, j]) = xh_plus_ij([i, j]) + h;
        xh_plus_i = x;
        xh_plus_i(i) = xh_plus_i(i) + h;
        xh_plus_j = x;
        xh_plus_j(j) = xh_plus_j(j) + h;
        Hessf1(j) = (f(xh_plus_ij) - ...
            f(xh_plus_i) - f(xh_plus_j) + f(x))/(h^2);
    end
   
end

% Creation of the sparse hessian
H = sparse(1:n,1:n,Hessf0,n,n) + sparse(2:n,1:n-1,Hessf1,n,n) + sparse(1:n-1,2:n,Hessf1,n,n);

% Add elements in position (1,n) and (n,1) based on the structure of the hessian
xh_plus_ij = x;
xh_plus_ij([1, n]) = xh_plus_ij([1, n]) + h;
xh_plus_i = x;
xh_plus_i(1) = xh_plus_i(1) + h;
xh_plus_j = x;
xh_plus_j(n) = xh_plus_j(n) + h;
H(1,n) = (f(xh_plus_ij) - ...
    f(xh_plus_i) - f(xh_plus_j) + f(x))/(h^2);

H(n,1) = H(1,n);
    
end
