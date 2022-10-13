function dist = DefIntegral(timeSpan,velocity,h,startTime,endTime)
% DefIntegral calculates the distance travelled over a given time
% using Simpson's Rule
%
% Inputs:
% timeSpan: Time range to operate over
% velocity: Vector of Velocity values
% h: global step size for calculations
% startTime: start time in seconds
% endTime: end time in seconds
% 
% Outputs:
% dist: Distance travelled over the period of time
dist = 0;
f = @(t) velocities(t);
end