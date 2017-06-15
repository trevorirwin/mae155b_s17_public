%%%%%%%%%%%%%
% Changes made @(06/14/17, 8:35pm) by Trevor Irwin:
% 
% Created biodiesel variant of script.
%
% Changes made @(06/14/17, 2:45pm) by Justin Laughlin:
% 
% Pa = 89137 [Pa] -> 89888 [Pa]
% (https://www.mide.com/pages/air-pressure-at-altitude-calculator)
% h = 3504 ft, T = 95 F
% 
% Reworded "Engine Pressure" to "Chamber Pressure"
% 
% gamma = 1.40 -> 1.22
% (http://www.braeunig.us/space/comb.htm)
% Chamber pressure ~ 25*Patm
% 
% Increased fontsize for axes, title, legend, and tick marks
% 
% Increased linesize
% 
%%%%%%%%%%%%%

%
% Cleanup
%
close all;
clear all;
clc;
%
% Misc Parameters
%
fsbig   = 16;   % larger fontsize for title/axes
fsmed   = 12;   % fontsize for legend/ticks
lwbig   = 2.5;  % bolder linewidth
lwsm    = 1;    % thinner linwdith
%
% Load data
%
load('bio-test2-counter-v1.mat');
load('bio-test2-analog-v1.mat')
adata = adata_bio;
cdata = cdata_bio;
%
% Offset T = 0 to start of data.
%
test_offset = 286.59; % [s] Time at which test started
adata(:,1) = adata(:,1) - test_offset;
cdata(:,1) = cdata(:,1) - test_offset;
%
% Chamber Pressure Plot
%
fig = figure('name', 'Chamber Pressure');
plot(adata(:,1), adata(:,3),'linewidth',lwbig);
hold on;
xlabel('Time [s]','fontsize',fsbig);
ylabel('Chamber Pressure [psi]','fontsize',fsbig);
title('Chamber Pressure vs Time (Biodiesel Test 2)','fontsize',fsbig);
set(gca, 'FontSize',fsmed)
%
% Flow Rate Plot
%
fig = figure('name', 'Mass Flow Rates');
plot(cdata(:,1), cdata(:,9),'linewidth',lwbig);
hold on;
plot(cdata(:,1), cdata(:,10),'linewidth',lwbig);
xlabel('Time [s]','fontsize',fsbig);
ylabel('Flow Rate [GPM]','fontsize',fsbig);
title('Flow Rate vs Time (Biodiesel Test 2)','fontsize',fsbig);
lgd = legend('Fuel', 'LOX');
lgd.FontSize = fsmed;
set(gca, 'FontSize',fsmed)
%
% Create dummy longer counter data
%
crow = 1;
cdata_long = zeros(size(adata));
for arow = 1:length(adata)
    if adata(arow,1) < cdata(crow,1)
        cdata_long(arow,:) = cdata(crow,:);
    else
        crow = crow + 1;
        cdata_long(arow,:) = cdata(crow,:);
    end
end
%
% Thrust and ISP Calculation
%
Pc = adata(:,3);
Pc = Pc*6894.76; % [PSI -> Pa] Chamber pressure

A       = 0.0013380618;     % [m^2] Throat area
Ae      = 0.00684321212;    % [m^2] Exit Area
Pa      = 89888;            % [Pa] Take ambient pressure at 3504' elevation
Pe      = Pa;               % [Pa] Assume exit pressure = ambient
gam     = 1.22;             % Ratio of specific heats
rho_bio = 880;              % [kg/m^3] Assume constant density Biodiesel
rho_lox = 1141;             % [kg/m^3] Assume constant density LOX
g0      = 9.80665;          % [m/s^2] Standard gravity

T = A.*Pc.*(sqrt(((2.*(gam.^2))/(gam-1)).*((2./(gam+1)).^((gam+1)./(gam-1))).*...
    (1-(Pe./Pc).^((gam-1)./gam)))+((Pe./Pc)-(Pa./Pc))*(Ae./A)); % [N] thrust

T = real(T); % Sqrt gives some bogus values, take only real.

T_imp = real(T.*0.224809); % [N -> lbf]

q_dot_bio = cdata_long(:,9)*6.30902e-5; % [gal/min -> m^3/s]
q_dot_lox = cdata_long(:,10)*6.30902e-5; % [gal/min -> m^3/s]

m_dot_bio = q_dot_bio.*rho_bio; % [kg/s]
m_dot_lox = q_dot_lox.*rho_lox; % [kg/s]

m_dot_t = m_dot_bio + m_dot_lox; % [kg/s]

Isp = T./(m_dot_t.*g0); %[s]

%
% Plot Thrust and ISP
%
fig = figure('name', 'Thrust');
plot(adata(:,1), T_imp,'linewidth',lwbig);
hold on;
xlabel('Time [s]','fontsize',fsbig);
ylabel('Thrust [lbf]','fontsize',fsbig);
title('Thrust vs Time (Biodiesel Test 2)','fontsize',fsbig);
set(gca, 'FontSize',fsmed)

fig = figure('name', 'Specific Impulse');
i = find(adata(:,1) < 7); % Limit to 7s because it blows up when flow rate drops
plot(adata(i,1), Isp(i),'linewidth',lwbig);
hold on;
xlabel('Time [s]','fontsize',fsbig);
ylabel('Specific Impulse [s]','fontsize',fsbig);
title('Specific Impulse vs Time (Biodiesel Test 2)','fontsize',fsbig);
set(gca, 'FontSize',fsmed)

%
% Plot Thrust and Flow Rate
%
fig = figure('name', 'Thrust & Flow Rate');
yyaxis left;
plot(adata(:,1), T_imp,'linewidth',lwbig);
hold on;
xlabel('Time [s]','fontsize',fsbig);
ylabel('Thrust [lbf]','fontsize',fsbig);
title('Thrust & Flow Rate vs Time (Biodiesel Test 2)','fontsize',fsbig);
yyaxis right;
plot(adata(:,1), m_dot_t,'linewidth',lwbig);
ylabel('Mass Flow Rate (Total) [kg/s]','fontsize',fsbig);
lgd = legend('Thrust', 'Mass Flow Rate');
lgd.FontSize = fsmed;
set(gca, 'FontSize',fsmed)