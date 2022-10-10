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
%to do: Implement function for divided diff interpolation
% div difference is easier to get the polynomial back as a 
% anon function for plotting
%subset timspan and position by the first bounce

%Interpolating function for 4 points either side of deck crossing

%iterate root finding method for exact crossing value to maybe 1E-5 error
%iterate root finding method for exact crossing value to low error

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