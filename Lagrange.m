function [P] = Lagrange(point_coords)
% Lagrange2 Generates a Lagrange interpolated polynomial using four points.
% 
% Requires points to be in a 2xn matrix with x-coordinates in first row and
% corresponding y-coordinates in second row. Returns the resulting 
% interpolated polynomial in the form of an anonymous function 'P.'
%
% Inputs:
% point_coords: A 2-by-n matrix of the points to be included in polynomial
%
% Outputs:
% P: Polynomial that passes through each point in 'point_coords'

x0 = point_coords(1, 1) ;
x1 = point_coords(1, 2) ;
x2 = point_coords(1, 3) ;
x3 = point_coords(1, 4) ;

y0 = point_coords(2, 1) ;
y1 = point_coords(2, 2) ;
y2 = point_coords(2, 3) ;
y3 = point_coords(2, 4) ;

L0 = @(x) (((x) - x1)/(x0 - x1)).*(((x) - x2)/(x0 - x2)).*(((x) - x3)/(x0 - x3)) ;
L1 = @(x) (((x) - x0)/(x1 - x0)).*(((x) - x2)/(x1 - x2)).*(((x) - x3)/(x1 - x3)) ;
L2 = @(x) (((x) - x0)/(x2 - x0)).*(((x) - x1)/(x2 - x1)).*(((x) - x3)/(x2 - x3)) ;
L3 = @(x) (((x) - x0)/(x3 - x0)).*(((x) - x1)/(x3 - x1)).*(((x) - x2)/(x3 - x2)) ;

P = @(x) y0 .* L0(x) + y1 .* L1(x) + y2 .* L2(x) + y3 .* L3(x) ;
end