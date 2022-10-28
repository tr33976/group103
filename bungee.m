clc
clear
format long g
close all;

% Base parameters
h = 0.0001;
timeSpan = 0:h:80;
maxBounces = 10;
defaultparams = true; 
watertouchsearch = false;

%%%%%% DEFAULT PARAMETERS - DONT CHANGE %%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Velocity ODE - Captures parameters set above when run
% Make sure to run the whole script together to avoid
% unrecognised changes
dvdt = @(y, v) g - C .* abs(v) .* v - max(0, K .*(y-L));

%%%%%%%%%%%%%%%%%%%%%%%%% TASK 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of position and velocity vs time
% ODE evaluated using RK4 method
[position, velocity] = RK4Coupled(dvdt, timeSpan, h, 0, 0);

%%%%% PLOTS %%%%%
% Plot of position vs time
f = figure('Position',[100 100 700 200]) ;
plot(timeSpan, position)
title('Position (m from platform) vs Time')
xlabel('Time(s)')
ylabel({'Relative Postion','from Platform (m)'})
saveas(f, ['fig1','.png'])

% Plot of position vs time with bounce annotations
f = figure('Position',[100 100 900 250]) ;
minimaIDX = islocalmin(position) ;
plot(timeSpan(minimaIDX),position(minimaIDX),'g*')
hold on
iter = find(minimaIDX == 1, maxBounces) ;
stopVal = 0 ;
for ii = 1:length(find(minimaIDX == 1, 10))
    item = iter(ii) ;
    text(timeSpan(item), position(item) -5, num2str(ii),'Color', 'k')
    if ii == length(find(minimaIDX == 1, 10))
        stopVal = timeSpan(item) ;
    end
end
hold on
plot(timeSpan, position)
yline(H, 'b','River')
yline(0, 'k', 'Jump Point')
yline(cam, 'k', 'Deck')
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
f = figure('Position',[100 100 900 275]) ;
plot(timeSpan, velocity)
title('Velocity of the Jumper vs Time')
xlabel('Time(s)')
ylabel('Velocity (m/s)')
[max_velocity, max_velocity_idx] =  max(abs(velocity)) ;
hold on                      
plot(timeSpan(max_velocity_idx),max_velocity, '*') 
text(timeSpan(max_velocity_idx)+0.5, max_velocity+0.5, ['Max Velocity of ', num2str(max_velocity), 'm/s reached ' ...
    'at ', num2str(round(timeSpan(max_velocity_idx),3)), ' seconds'])
xlim([0, 60])
ylim([-20, 25])
saveas(f, ['fig3','.png'])
fprintf('Max Velocity: %.3fm/s\n', max_velocity);
%%%%%%%%%%%%%%%%%%%%%%%%% END 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% TASK 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of acceleration values with numerical differentiation
acceleration = CentralDifferentiation(velocity, h) ;

%%%%% PLOTS %%%%%
f = figure('Position',[100 100 900 275]) ;
plot(timeSpan, acceleration)
title('Acceleration of the Jumper vs Time')
xlabel('Time(s)')
ylabel('Acceleration (m/s^2)')
xlim([0, 80])
[max_acceleration, max_acceleration_idx] = max(acceleration) ;
max_acceleration = round(max_acceleration,3) ;
[max_acceleration_abs, max_acceleration_idx_abs] = max(abs(acceleration)) ;
max_acceleration_abs = round(max_acceleration_abs) ;
hold on
plot(timeSpan(max_acceleration_idx), max_acceleration, '*')
plot(timeSpan(max_acceleration_idx_abs), -max_acceleration_abs, '*')
text(timeSpan(max_acceleration_idx)+0.5, max_acceleration+0.5, ['Max Acc ', num2str(max_acceleration/9.8), 'G: ' ...
    , num2str(round(timeSpan(max_acceleration_idx),2)), 's'])
text(timeSpan(max_acceleration_idx_abs)+0.5, -max_acceleration_abs+0.5, ['Max Neg Acc ', num2str(max_acceleration_abs/9.8), 'G: ' ...
    , num2str(round(timeSpan(max_acceleration_idx_abs),2)), 's'])
xlim([0, 60])
saveas(f, ['fig4','.png'])
fprintf('Max Absolute Gs: %.3fG\n', abs(max_acceleration_abs/9.8)) ;
%%%%%%%%%%%%%%%%%%%%%%%%% END 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% TASK 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
endTime = 60 ;

% Evaluating integral to estimate distance
% Numerical integration performed with Simpson's rule
distTravelled = SimpsonMethod(endTime, timeSpan, velocity, h );

fprintf('Total distance travelled: %.2fm\n', distTravelled) ;
%%%%%%%%%%%%%%%%%%%%%%%%% END 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% TASK 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finding four points to be used for polynomial interpolation 
point_coords = ones(2, 4) ;

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

point_coords = [x0 x1 x2 x3; y0 y1 y2 y3] ;

% Lagrange polynomial interpolation 
P = Lagrange(point_coords) ;

% Modified version of p(t) 
P_mod = @(x) P(x) - 43 ;

% Symbolic expression of p(t) (No calculations performed, just for visualisation)
syms x
P_sym = vpa(expand(P_mod(x)), 4) ; % Symbolic expression of p(t)

% Estimated value of t at p(t) = H - D
tn = Bisection(x1, x2, P_mod, 1e-10) ;

tn = round(tn, 6);
fprintf('Trigger time: %.3fm\n', tn) ;

%%%%% PLOTS %%%%%
f = figure('Position',[100 100 1000 350]) ;
t = tiledlayout(1,2, "TileSpacing",'loose') ;
title(t, "Jumper Camera Crossing", 'FontSize',15, 'FontWeight','bold') ;
nexttile ;
plot(timeSpan, position) ;
hold on ;
fplot(P1) ;
plot(tn, cam, 'k.', 'MarkerSize',18) ;
yline(cam, 'k', 'Deck') ;
xlim([0, tn+10]) ;
ylim([-10, max(position)+10]) ;
xlabel('Time(s)') ;
ylabel({'Relative Postion','from Platform (m)'}) ;
legend('Jumper', 'Interplotating Polynomial') ;
set(gca, 'YDir','reverse') ;
nexttile ;
plot(timeSpan, position) ;
hold on ;
fplot(P1) ;
plot(tn, cam, 'k.', 'MarkerSize',18) ;
yline(cam, 'k', 'Deck') ;
xlabel('Time(s)') ;
text(tn+0.009,cam,{'Deck Crossing', [num2str(round(tn,3)), '(secs)']},'Color','k') ;
xlim([tn-0.1, tn+0.1]) ;
ylim([cam-0.1, cam+0.1]) ;
set(gca, 'YDir','reverse') ;
saveas(f, ['fig5','.png']) ;
%%%%%%%%%%%%%%%%%%%%%%%%% END 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%% TASK 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if watertouchsearch
    jumperHeight = 1.75 ;
    
    % Optimal parameters for water touch
    [wt_len, wt_spring, wt_res_time, wt_acc] = WaterTouchParams(jumperHeight, H);
    
    % Re-evaluation of numerical solution with new parameters
    wt_dvdt = @(y, v) g - C .* abs(v) .* v - max(0, wt_spring/m .*(y-wt_len));
    [wt_position, wt_velocity] = RK4Coupled(wt_dvdt, timeSpan, h, 0, 0);
    
    %%%%% PLOTS %%%%%
    % Without water touch
    f = figure('Position',[100 100 900 250]);
    plot(timeSpan, position)
    yline(H, 'b','River')
    yline(0, 'k', 'Jump Point')
    yline(cam, 'k', 'Deck')
    xline(stopVal, 'k', {'Stop (s): ',num2str(stopVal)}, ...
        'LabelOrientation', 'aligned', 'LabelHorizontalAlignment', 'left')
    title('Distance from Water - Given Parameters')
    subtitle('*Y Axis Reversed*')
    m = timeSpan(position == max(position)) ;
    line([m m], [H position(position == max(position))], 'Color', 'r') ;
    text(m+0.2,(H+position(position == max(position)))/2, {'Distance from water (m):', num2str(H-position(position == max(position)))}, 'Color', 'k')
    xlabel('Time(s)')
    ylabel({'Relative Postion','from Platform (m)'})
    ylim([0 80])
    xlim([0 20])
    set(gca, 'YDir','reverse')
    saveas(f, ['fig6','.png'])
    
    % Water touch
    f = figure('Position',[100 100 900 250]) ;
    plot(timeSpan, wt_position+jumperHeight)
    yline(H, 'b','River')
    yline(0, 'k', 'Jump Point')
    yline(cam, 'k', 'Deck')
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
