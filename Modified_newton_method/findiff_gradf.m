function [gradfx] = findiff_gradf(f, x, h, type)

gradfx = zeros(size(x));

switch type
    % forward finite differences
    case 'fw'
        for i=1:length(x)
            xh = x;
            xh(i) = xh(i) + h;
            gradfx(i) = (f(xh) - f(x))/ h;
            
        end
    % centered finite differences
    case 'c'
        for i=1:length(x)
            xh_plus = x;
            xh_minus = x;
            xh_plus(i) = xh_plus(i) + h;
            xh_minus(i) = xh_minus(i) - h;
            gradfx(i) = (f(xh_plus) - f(xh_minus))/(2 * h);
        end
    % backward finite differences
    otherwise 
        for i=1:length(x)
            xh = x;
            xh(i) = xh(i) + h;
            gradfx(i) = (f(xh) - f(x))/h;
        end
end
end
