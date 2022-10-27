function [x] = Bisection(a, b, f, error)
% Bisection Approximates the root of a function using the bisection method.
% Given some function  of x 'f', the value of x when f = 0 is
% approximated, provided this value of x lies between 'a' and 'b.'
%
% Inputs:
% a0:    Initial lower boundary (Must be less than actual value of root)
% b0:    Initial upper boundary (Must be greater than actual value of root)
% f:     Function to find roots of 
% error: Number of decimal points wanted for answer
%
% Outputs:
% x:     Final estimate of root of function

c = (a + b) / 2 ; % Midpoint between lower and upper bounds

while abs(f(c)) > error
    if f(c) < 0 && f(a) < 0
        a = c ;       % Updating a
    else
        b = c ;       % Updating b
    end
    c = (a + b) / 2 ; % Updating c
end

x = c ;
end