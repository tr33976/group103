function c = BisectionMethod(a, b, func, error)
% BisectionRoot Performs iterative root finding in a bracket of A of B for
% input function to degree of input error

% Inputs:
% a: upper bound for bisection
% b: lower bound for bisection
% func: function to root find
% error: error target to iterate until
%
% Outputs: 
% c: Mid point of bracket estimating root value to within error bound. 

    c = (a + b) / 2; % Midpoint between lower and upper bounds
    
    while abs(func(c)) > error
        if func(c) < 0 && func(a) < 0
            a = c ;       % Updating a
        else
            b = c ;       % Updating b
        end
        c = (a + b) / 2 ; % Updating c
    end
end