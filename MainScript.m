clc
clear
format long g

% base parameterss
% base aprox parameters
h = 0.0001;
timeSpan = 0:h:80;
maxBounces = 10;
defaultparams = true; 
watertouchsearch = true;

switch defaultparams
    case true
        %%%%%% DEFAULT PARAMETERS DONT CHANGE %%%
        H = 74; % height of jump (m)
        DH = 31; % heights of deck from water (m)
        DECK =  H-DH; % distance from jump to deck
        c = 0.9; % drag coefficient (kg/m)
        m = 80; % jumper mass (kg)
        C = c/m; % drag / mass 
        L = 25; % bunge rope length (m)
        k = 90; % rope spring (N/m)
        K = k/m; % spring / mass
        g = 9.8; % gravity (m/s^2)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case false
        %%%%%% CUSTOMER PARAMETERS CHANGE AWAY %%%
        %%%%%% Can be used for playing/testing %%%
        H = 74; % height of jump (m)
        DH = 31; % heights of deck from water (m)
        DECK =  H-DH; % distance from jump to deck
        c = 0.9; % drag coefficient (kg/m)
        m = 80; % jumper mass (kg)
        C = c/m; % drag / mass 
        L = 25; % bunge rope length (m)
        k = 90; % rope spring (N/m)
        K = k/m; % spring / mass
        g = 9.8; % gravity (m/s^2)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    otherwise
        error('invalid parameter flag')
end
% velocity ode, captures parameters set above when run
% make sure to run the whole script together to avoid
% unrecognised changes
dvdt = @(y, v) g - C .* abs(v) .* v - max(0, K .*(y-L));


%%%%%%%%%%%%%%%%%%%%%%%%% TASK 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate position and velocity versus time
% ode evaluation using RK4 method
[position, velocity] = RK4Coupled(dvdt, timeSpan, h, 0, 0, false);

%intro section plot
f=figure('Position',[100 100 700 200]);
plot(timeSpan, position)
title('Position (m from platform) vs Time')
xlabel('Time(s)')
ylabel({'Relative Postion','from Platform (m)'})
saveas(f, ['fig1','.png'])

% plot position vs time with bounce annotations
f=figure('Position',[100 100 900 250]);
minimaIDX=islocalmin(position);
plot(timeSpan(minimaIDX),position(minimaIDX),'g*')
hold on
% this bit adds stars and numbers to bounces
iter = find(minimaIDX==1, maxBounces);
stopVal=0;
for ii = 1:length(find(minimaIDX==1, 10))
    item=iter(ii);
    text(timeSpan(item),position(item)-5,num2str(ii),'Color','k')
    if ii==length(find(minimaIDX==1, 10))
        stopVal=timeSpan(item);
    end
end
hold on
plot(timeSpan, position)
yline(H, 'b','River')
yline(0, 'k', 'Jump Point')
yline(DECK, 'k', 'Deck')
xline(stopVal, 'k', {'Stop (s): ',num2str(stopVal)}, ...
    'LabelOrientation', 'aligned', 'LabelHorizontalAlignment', 'left')
title('Position (m from platform) vs Time')
subtitle('*Y Axis Reversed*')
xlabel('Time(s)')
ylabel({'Relative Postion','from Platform (m)'})
ylim([0 80])
xlim([0 stopVal+5])
set(gca, 'YDir','reverse')
saveas(f, ['fig2','.png'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%% END 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% TASK 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=figure('Position',[100 100 1500 500]);
plot(timeSpan, velocity)
title('Velocity of the Jumper vs Time')
xlabel('Time(s)')
ylabel('Velocity (m/s)')
xlim([0, 80])
max_velocity = max(velocity);
hold on                      %not sure how to find when max velocity is 
plot(2.605,max_velocity, '*')       %reached other than graphically 
text(3.5, max_velocity, ['Max Velocity of 20.0144 m/s reached ' ...
    'at 2.605 seconds'])
saveas(f, ['fig3','.png'])
%%%%%%%%%%%%%%%%%%%%%%%%% END 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% TASK 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%to do: numerically differentiate velocity to get acceleration
%suggest central difference method as the easiest with forward
%differece for the first point and backward diff for last
%acceleration = NumericalDiff(velocity, h)
acceleration = CentralDifferentiation(velocity, h);

f=figure('Position',[100 100 1500 500]);
plot(timeSpan, acceleration)
title('Acceleration of the Jumper vs Time')
xlabel('Time(s)')
ylabel('Acceleration (m/s^2)')
xlim([0, 80])
max_acceleration = max(acceleration);
hold on                      %not sure how to find when max velocity is 
plot(6.397,max_acceleration, '*')       %reached other than graphically 
text(7, max_acceleration, ['Max Acceleration of 11.0694 m/s^2 reached ' ...
    'at 6.397 seconds'])
saveas(f, ['fig4','.png'])
%%%%%%%%%%%%%%%%%%%%%%%%% END 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% TASK 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%to do: Numerically integrate to get distance travelled, 
%simpsons or trapezoid are the options we have it seems.
% distTravelled = DefIntegral(timeSpan, abs(velocity), h, 0, 60)
% sprintf('Total distance travelled: \n%.2f', DECK)
%%%%%%%%%%%%%%%%%%%%%%%%% END 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% TASK 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1: Find y(i), y(i+1), y(i+2), y(i+3)
% Defining variables (I'll clean this up)
H = 74 ;  % Height of jump point (in m)
D = 31 ;  % Height of deck (in m)
c = 0.9 ; % Drag coefficient (in kg/m)
m = 80 ;  % Mass of jumper (in kg)
L = 25 ;  % Length of bungee rope (in m)
k = 90 ;  % Spring constant of bungee rope (in N/m)
g = 9.8 ; % Acceleration due to gravity (in m/s^2)

C = c/m ; % Drag/mass
K = k/m ; % Spring constant/mass

cam = H - D ; % The distance at which the camera is placed

t = 100 ;          % Total time of jump (in seconds?)
h = 0.0001 ;       % Step size
timeSpan = 0:h:t ; % Number of steps 
y0 = 0 ;           % Initial position of jumper
v0 = 0 ;           % Initial velocity of jumper

f = @(y, v) g - C .* abs(v) .* v - max(0, K .*(y-L)) ; 

position =  zeros(1,length(timeSpan)) ;
velocity =  zeros(1,length(timeSpan)) ;
    
y = 0 ; 
v = 0 ;
    
position(1) = y ;
velocity(1) = v ;
    
for i = 2:length(timeSpan)
    k1 = h * v ;
    r1 = h * f(y, v) ;
       
    k2 = h * (v + r1 / 2) ;
    r2 = h * f(y + k1 / 2 , v + r1 / 2) ;
       
    k3 = h * (v + r2 / 2) ;
    r3 = h * f(y + k2 / 2 , v + r2 / 2) ;
       
    k4 = h * (v + r3) ;
    r4 = h * f(y + k3, v + r3) ;
       
    y = y + (1/6) * (k1 + 2*k2 + 2*k3 + k4) ;
    v = v + (1/6) * (r1 + 2*r2 + 2*r3 + r4) ;

    position(i) = y ;
    velocity(i) = v ;

    if i > 3 && position(i-2) < cam && position(i-1) > cam
        y0 = position(i-3) ; % y(i)
        x0 = (i-3) * h ;     % t(i)
        y1 = position(i-2) ; % y(i+1)
        x1 = (i-2) * h ;     % t(i+1)
        y2 = position(i-1) ; % y(i+2)
        x2 = (i-1) * h ;     % t(i+2)
        y3 = position(i) ;   % y(i+3)
        x3 = i * h ;         % t(i+3)
       break 
    end 
end

%% Step 2: Perform polynomial interpolation
%
% Points to be used for polynomial:
% (y(i), t(i)), (y(i+1), t(i+1)), (y(i+2), t(i+2)), (y(i+3), t(i+3))
% (42.9978, 3.3335), (42.9993, 3.3336), (43.0008, 3.3337), (43.0022, 3.3338) 
% 
% The Lagrange method will be used:
L0 = ((x - x1)/(x0 - x1))*((x - x2)/(x0 - x2))*((x - x3)/(x0 - x3)) ;
L1 = ((x - x0)/(x1 - x0))*((x - x2)/(x1 - x2))*((x - x3)/(x1 - x3)) ;
L2 = ((x - x0)/(x2 - x0))*((x - x1)/(x2 - x1))*((x - x3)/(x2 - x3)) ;
L3 = ((x - x0)/(x3 - x0))*((x - x1)/(x3 - x1))*((x - x2)/(x3 - x2)) ;

P1 = y0 * L0 + y1 * L1 + y2 * L2 + y3 * L3 ; % Symbolic expression of p(t)
P2 = sym2poly(P) ; % Converting p(t) to vector expression

%% Step 3: Finding t at p(t) = H - D
% Let tn represent the value of t at p(t). To find tn, p(t) will be modified 
% by subtracting (H - D) so that tn occurs at an x-intercept. 
% The value of tn will be estimated to a certain degree of error using a 
% root-finding method. 
P3 = P2 ; 
P3(4) = P3(4) - cam ; % Modified version of p(t)

% The bisection method will be used to estimate tn: 
a = x1 ;          % Lower bound on tn
b = x2 ;          % Upper bound on tn
c = (a + b) / 2 ; % Midpoint between lower and upper bounds
f = @(x) P2(1) * x^3 + P2(2) * x^2 + P2(3) * x + P2(4) ; % Modified p(t)
error = 1e-10 ;

while abs(f(c)) > error
    if f(c) < 0 && f(a) < 0
        a = c ;       % Updating a
    else
        b = c ;       % Updating b
    end
    c = (a + b) / 2 ; % Updating c
end

format long 
tn = c ;
tn = round(tn, 6) % Value of tn rounded to 6 decimal places

% The value of t at p(t) = H - D is approximately 3.333649
% Therefore, for the model parameters provided, the camera should 
% be triggered as close to 3.333649 seconds as possible. 
%%%%%%%%%%%%%%%%%%%%%%%%% END 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%% TASK 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if watertouchsearch
jumperHeight = 1.75;

[wt_len, wt_spring, wt_res_time, wt_acc] = WaterTouchParams(jumperHeight, H);

wt_dvdt = @(y, v) g - C .* abs(v) .* v - max(0, wt_spring/m .*(y-wt_len));
[wt_position, wt_velocity] = RK4Coupled(wt_dvdt, timeSpan, h, 0, 0, false);


% no touch
f=figure('Position',[100 100 900 250]);
plot(timeSpan, position)
yline(H, 'b','River')
yline(0, 'k', 'Jump Point')
yline(DECK, 'k', 'Deck')
xline(stopVal, 'k', {'Stop (s): ',num2str(stopVal)}, ...
    'LabelOrientation', 'aligned', 'LabelHorizontalAlignment', 'left')
title('Distance from Water - Given Parameters')
subtitle('*Y Axis Reversed*')
m=timeSpan(position==max(position));
line([m m], [H position(position==max(position))], 'Color', 'r');
text(m+0.2,(H+position(position==max(position)))/2,{'Distance from water (m):', num2str(H-position(position==max(position)))},'Color','k')
xlabel('Time(s)')
ylabel({'Relative Postion','from Platform (m)'})
ylim([0 80])
xlim([0 20])
set(gca, 'YDir','reverse')
saveas(f, ['fig6','.png'])

f=figure('Position',[100 100 900 250]);
plot(timeSpan, wt_position+jumperHeight)
yline(H, 'b','River')
yline(0, 'k', 'Jump Point')
yline(DECK, 'k', 'Deck')
xline(stopVal, 'k', ['Stop (s): ',num2str(stopVal)], ...
    'LabelOrientation', 'aligned', 'LabelHorizontalAlignment', 'left')
title('Postion Adjusted By Jumper Height vs Time')
subtitle(['*Y Axis Reversed* Max Gs: ', num2str(wt_acc)])
xlabel('Time(s)')
ylabel('Relative Height from Platform (m)')
ylim([0 80])
xlim([0 stopVal+5])
set(gca, 'YDir','reverse')

disp(['Water touch Spring Constant: ', num2str(wt_spring),'N/m'])
disp(['Water touch Rope Length: ', num2str(wt_len), 'm'])
saveas(f, ['fig7','.png'])
end
%%%%%%%%%%%%%%%%%%%%%%%%% END 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
