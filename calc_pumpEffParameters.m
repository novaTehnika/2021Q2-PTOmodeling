% calc_pumpEffParameters.m script m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 12/13/2021
%
% PURPOSE/DESCRIPTION:
% This script calculates the coefficienct for the Wilson pump and motor 
% efficiency model, given one data point for flow rate and one data point 
% for torque. An assumption is made for the volume fraction coefficient
% (V_r) and for mechanical efficienct at a low pressure (500kPa).
%
% FILE DEPENDENCY:
%
% UPDATES:
% 12/13/2021 - Created.
%
% Copyright (C) 2022  Jeremy W. Simmons II
% 
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program. If not, see <https://www.gnu.org/licenses/>.
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

mu = 1.089e-3;
beta = 2.2e9;

D = 444e-6/2/pi; % [m^3/rad]
P = 6e6; % [Pa]
Q = 44.6/3600; % [m^3/h -> m^3/s]
wq = 2*pi*1700/60; % [rpm -> rad/s]

V_r = 1.103;
C_s = mu*wq/P*(1-P/beta*(V_r+1)-Q/wq/D);


PP1 = 105e3; % [W]
P1 = 8e6; % [Pa]
P2 = 500e3; % [Pa]
wm = 2*pi*1700/60; % [rpm -> rad/s]

eta_m1 = wm*P1*D/PP1;

A = [mu*wm/P1, 1; mu*wm/P2, 1];
b = [1/eta_m1 - 1; 1/0.8/eta_m1 - 1];
x = A\b;
C_v = x(1);
C_f = x(2);


%%
p = linspace(1e6,8e6,20);
w = linspace((200)*2*pi/60,wm,20);
for j = 1:length(p)
    for k = 1:length(w)
        eta_v(j,k) = 1 - C_s*(p(j)/(mu*abs(w(k))))...
            - (p(j)/beta)*(V_r + 1);         % Pump volumetric efficiency
        eta_m(j,k) = 1/(1 + C_v*mu*abs(w(k))/...
            p(j) + C_f);                      % Pump mechanical efficiency
    end
end


figure
surf(w,p*1e-6,eta_v)
xlabel('Angular Velocity (rad/s)')
ylabel('Pressure (MPa)')
zlabel('Efficieny')
title('Volumetric Efficiency in Pumping Mode')
xlim([0 w(end)])
ylim(1e-6*[0 p(end)])

figure
surf(w,p*1e-6,eta_m)
xlabel('Angular Velocity (rad/s)')
ylabel('Pressure (MPa)')
zlabel('Efficieny')
title('Mechanical Efficiency in Pumping Mode')
xlim([0 w(end)])
ylim(1e-6*[0 p(end)])


%%
p = linspace(1e6,8e6,20);
w = linspace((200)*2*pi/60,wm,20);
for j = 1:length(p)
    for k = 1:length(w)
        eta_v(j,k) = 1/(1 + C_s*(p(j)/(mu*abs(w(k))))...
            + (p(j)/beta)*(V_r + 1));         % Pump volumetric efficiency
        eta_m(j,k) = (1 - C_v*mu*abs(w(k))/...
            p(j) - C_f);                      % Pump mechanical efficiency
    end
end


figure
surf(w,p*1e-6,eta_v)
xlabel('Angular Velocity (rad/s)')
ylabel('Pressure (MPa)')
zlabel('Efficieny')
title('Volumetric Efficiency in Motoring Mode')
xlim([0 w(end)])
ylim(1e-6*[0 p(end)])

figure
surf(w,p*1e-6,eta_m)
xlabel('Angular Velocity (rad/s)')
ylabel('Pressure (MPa)')
zlabel('Efficieny')
title('Mechanical Efficiency in Motoring Mode')
xlim([0 w(end)])
ylim(1e-6*[0 p(end)])
