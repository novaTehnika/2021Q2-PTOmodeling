% run_linearPTO.m script m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 08/23/2022
%
% PURPOSE/DESCRIPTION:
% This script serves as a shell for running a single simulation
% using the model contained in sys_coulombPTO.m and run by solved by 
% sim_linearPTO.m.
% The parameter initiallization fuction are called within this
% script before the sim_linearPTO.m script is called.
%
% FILE DEPENDENCY:
% sys_linearPTO.m
% sim_linearPTO.m
% parameters_linearPTO.m
%
% UPDATES:
% 08/23/2022 - Created from run_coulombPTO.m.
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
% clc
addpath('WEC model')
addpath(['WEC model' filesep 'WECdata'])
addpath('Linear damping PTO') 
addpath('Sea States')
addpath('Solvers')
%% %%%%%%%%%%%%   SIMULATION PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation timeframe
par.Tramp = 250; % [s] excitation force ramp period
par.tstart = 0; %[s] start time of simulation
par.tend = 2000; %[s] end time of simulation

% Solver parameters
par.odeSolverRelTol = 1e-4; % Rel. error tolerance parameter for ODE solver
par.odeSolverAbsTol = 1e-4; % Abs. error tolerance parameter for ODE solver
par.MaxStep = 1e-2;

% Sea State and Wave construction parameters
par.wave.Hs = 1.75;
par.wave.Tp = 7;
par.WEC.nw = 1000; % num. of frequency components for harmonic superposition 
par.wave.rngSeedPhase = 3; % seed for the random number generator

% load parameters
par = parameters_linearPTO(par,...
    'nemohResults_vantHoff2009_20180802.mat','vantHoffTFCoeff.mat');

% Define initial conditions
y0 = [  0, ...
        0, ...
        zeros(1,par.WEC.ny_rad)];

%% Special modifications to base parameters
par.Cpto = 50e6; % [Nms/rad] PTO reaction torque
% par.imprt.WEC.I_inf = 24366095;
% par.WEC.I_inf = par.imprt.WEC.I_inf;
%% %%%%%%%%%%%%   COLLECT DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% run simulation
tic
out = sim_linearPTO(y0,par);
dur = toc

nw = par.WEC.nw
% par.tend
% par.MaxStep
% par.wave.rngSeedPhase

itVec = find(out.t>=par.tstart);
PP = mean(-out.T_pto(itVec).*out.theta_dot(itVec))
RMSelev = sqrt(mean(out.waveElev(itVec).^2));
RMS_Te = sqrt(mean(out.T_wave(itVec).^2));
RMS_theta = sqrt(mean(out.theta(itVec).^2));
RMS_PTOtorque = sqrt(mean(out.T_pto(itVec).^2));
A = [PP, dur];
return
%% %%%%%%%%%%%%   PLOTTING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(out.t,out.waveElev)
xlabel('time (s)')
ylabel('elevation (m)')
title('Wave Elevation')

%%
figure
plot(out.t,out.theta)
xlabel('time (s)')
ylabel('position (rad)')
% ylim([-pi/2 pi/2])
% xlim([0 2])

%%
figure
xlabel('time (s)')

yyaxis left
hold on
plot(out.t,out.theta_dot)
ylabel('angular velocity (rad/s)')
% ylim(10*[-pi/2 pi/2])

yyaxis right
hold on
plot(out.t,1e-6*out.T_pto)
ylabel('torque, PTO (MNm)')
% ylim(2*1e-6*[-par.Tcoulomb par.Tcoulomb])
% xlim([0 2])

%%
figure
hold on
plot(out.t,1e-6*out.T_wave)
plot(out.t,1e-6*out.T_pto)
plot(out.t,1e-6*out.T_rad)
ylabel('torque (MNm)')

mean(-out.T_pto(:).*out.theta_dot(:))
mean(out.power.P_WEC)
