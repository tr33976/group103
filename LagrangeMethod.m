function [lagPoly, x1, x2] = LagrangeMethod(h, cam, timeSpan, position)
% LagrangeMethod Constructs an interpolating polynomial using the Lagrange
% Method and returns an anonymous function for evaluating the polynomial at
% values of x
%
% Inputs:
% h: Global step size
% cam: Camera height from the river
% timeSpan: Discrete time range to integrate over
% position: Discrete range of position values of jumper
%
% Outputs:
% lagPoly: Lagrange polynomial encoded as anon function
% x1: A bracket for root finding range
% x2: B bracket for root finding range

    for i = 2:length(timeSpan)
    if i > 3 && position(i-2) < cam && position(i-1) > cam
        y0 = position(i-3) ; % y(i)
        x0 = (i-4) * h ;     % t(i)
        y1 = position(i-2) ; % y(i+1)
        x1 = (i-3) * h ;     % t(i+1)
        y2 = position(i-1) ; % y(i+2)
        x2 = (i-2) * h ;     % t(i+2)
        y3 = position(i) ;   % y(i+3)
        x3 = (i-1) * h ;     % t(i+3)
       break 
    end 
    end
    
    % Points to be used for polynomial:
    % (y(i), t(i)), (y(i+1), t(i+1)), (y(i+2), t(i+2)), (y(i+3), t(i+3))
    % (42.9978, 3.3334), (42.9993, 3.3335), (43.0008, 3.3336), (43.0022, 3.3337) 
    % 
    % The Lagrange method will be used:
    L0 = @(x) (((x) - x1)/(x0 - x1)).*(((x) - x2)/(x0 - x2)).*(((x) - x3)/(x0 - x3)) ;
    L1 = @(x) (((x) - x0)/(x1 - x0)).*(((x) - x2)/(x1 - x2)).*(((x) - x3)/(x1 - x3)) ;
    L2 = @(x) (((x) - x0)/(x2 - x0)).*(((x) - x1)/(x2 - x1)).*(((x) - x3)/(x2 - x3)) ;
    L3 = @(x) (((x) - x0)/(x3 - x0)).*(((x) - x1)/(x3 - x1)).*(((x) - x2)/(x3 - x2)) ;
    
    lagPoly = @(x) y0 .* L0(x) + y1 .* L1(x) + y2 .* L2(x) + y3 .* L3(x) ;
end