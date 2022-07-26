% run_coulombPTOcompressiblePump.m script m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 06/23/2021
%
% PURPOSE/DESCRIPTION:
% This script serves as a shell for running a single simulation
% using the model contained in sys_coulombPTOcompressiblePump.m and run by 
% solved by sim_coulombPTOcompressiblePump.m.
% The parameter initiallization fuction are called within this
% script before the sim_coulombPTOcompressiblePump.m script is called. 
%
% FILE DEPENDENCY:
% sys_coulombPTOcompressiblePump.m
% sim_coulombPTOcompressiblePump.m
% parameters_coulombPTOcompressiblePump.m
% stateIndex_coulombPTOcompressiblePump.m
%
% UPDATES:
% 06/23/2021 - Created.
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
addpath('WEC model\WECdata') 
addpath('Coulomb damping PTO with compressible pump') 
%% %%%%%%%%%%%%   SIMULATION PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation timeframe
par.tstart = 0; %[s] start time of simulation
par.tend = 3000; %[s] end time of simulation

% Solver parameters
par.odeSolverRelTol = 1e-4; % Rel. error tolerance parameter for ODE solver
par.odeSolverAbsTol = 1e-4; % Abs. error tolerance parameter for ODE solver
par.MaxStep = 5e-5;

% Sea State and Wave construction parameters
par.wave.Hs = 2.75;
par.wave.Tp = 12;
par.WEC.nw = 100; % num. of frequency components for harmonic superposition 
par.wave.rngSeedPhase = 3; % seed for the random number generator

% load parameters
par = parameters_coulombPTOcompressiblePump(par);

% Define indicies of state vector
stateIndex_coulombPTOcompressiblePump

% Initial Conditions
y0(iyp_wpA) = par.p_lout;
y0(iyp_wpB) = par.p_lout;
y0(iyWEC) = [0, ...
             0, ...
             zeros(1,par.WEC.ny_rad)];

%% Special modifications to base parameters
% par.Dwp = 1; % [m^3/rad] pump displacement

%% %%%%%%%%%%%%   COLLECT DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% run simulation
tic
out = sim_coulombPTOcompressiblePump(y0,par);
toc

%% %%%%%%%%%%%%   PLOTTING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(out.t,out.waveElev)
xlabel('time (s)')
ylabel('elevation (m)')
title('Wave Elevation')

figure
plot(out.t,out.theta)
xlabel('time (s)')
ylabel('position (rad)')
% ylim([-pi/2 pi/2])
% xlim([0 2])

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
% xlim([0 2])

figure
ax(1) = subplot(3,1,1);
hold on
plot(out.t,1e-6*out.T_pto)
plot(out.t,1e-6*out.T_rad)
ylabel('torque (MNm)')
legend('T_{pto}','T_{rad}')

ax(2) = subplot(3,1,2);
hold on
plot(out.t,1e-6*par.p_lout*ones(out.nt,1))
plot(out.t,1e-6*par.p_hin*ones(out.nt,1))
plot(out.t,1e-6*out.p_wpA)
plot(out.t,1e-6*out.p_wpB)
ylabel('Pressure (MPa)')
legend('p_{lout}','p_{hin}','p_{wpA}','p_{wpB}')

ax(3) = subplot(3,1,3);
hold on
plot(out.t,1e3*out.q_wpP)
plot(out.t,1e3*out.q_wpT)
ylabel('Flow rate (Lpm)')
legend('q_{wpP}','q_{wpT}')

linkaxes(ax,'x');

mean(-out.T_pto(:).*out.theta_dot(:))
