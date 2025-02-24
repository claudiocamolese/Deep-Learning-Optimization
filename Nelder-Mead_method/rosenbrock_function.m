function f=rosenbrock_function(x)    
    if size(x,2)==1
        % If x is a column vector
        f= 100*(x(2)-x(1).^2).^2 + (1-x(1)).^2;
    else
        % If x is a simplex
        f=100*(x(2,:)-x(1,:).^2).^2+(1-x(1,:)).^2;
        f=f';
    end
end