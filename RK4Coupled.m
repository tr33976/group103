function [position, velocity] = RK4Coupled(in_g, timeSpan, h, inital_y, initial_v)
% RK4COUPLED Solves position and velocity coupled ODE.
%
% Inputs:
% in_g: Anon function representing second ODE
% timeSpan: Time range to operate over
% h: global step size for calculations
% inital_y: starting value for position
% initual_v: starting value for velocity
% 
% Outputs:
% position: Vector of position values indexed to timespan
% velocity: Vector of velocity values indexed to timespan

    position =  zeros(1,length(timeSpan));
    velocity =  zeros(1,length(timeSpan),1);
    
    y = inital_y;
    v = initial_v;
    
    position(1) = y;
    velocity(1) = v;
    
    for iterator = 2:length(timeSpan)      
       k1 = h * v;
       r1 = h * in_g(y, v);
       
       k2 = h * (v + r1 / 2);
       r2 = h * in_g(y + k1 / 2 , v + r1 / 2);
       
       k3 = h * (v + r2 / 2);
       r3 = h * in_g(y + k2 / 2 , v + r2 / 2);
       
       k4 = h * (v + r3);
       r4 = h * in_g(y + k3, v + r3);
       
       y = y + (1/6) * (k1 + 2*k2 + 2*k3 + k4);
       v = v + (1/6) * (r1 + 2*r2 + 2*r3 + r4);

       position(iterator) = y;
       velocity(iterator) = v;

    end
 
end
