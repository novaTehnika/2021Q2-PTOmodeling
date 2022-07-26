% run_coulombPTOshuntValve.m script m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 06/23/2021
%
% PURPOSE/DESCRIPTION:
% This script serves as a shell for performing parameter variation studies
% using simulations executed by the simPTO_VXXxXX.m function m-file.
% Parameter and function initiallization fuctions are called within this
% script before the sim_coulombPTOcompressiblePump_VXXxXX script is called. 
%
% FILE DEPENDENCY:
% sys_coulombPTOcompressiblePump.m
% sim_coulombPTOcompressiblePump.m
% parameters_coulombPTOcompressiblePump.m
% stateIndex_coulombPTOcompressiblePump.m
% WEC model/flapModel.m
%
% UPDATES:
% 12/30/2021 - Created from
% study_coulombPTOcompressiblePump.m.
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
addpath(['WEC model',filesep,'WECdata']) 
addpath('Coulomb damping PTO with shunt valve') 
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
par = parameters_coulombPTOshuntValve(par);

% Define indicies of state vector
stateIndex_coulombPTOshuntValve

% Initial Conditions
y0(iyp_wpA) = par.p_lout;
y0(iyp_wpB) = par.p_lout;
y0(iyWEC) = [0, ...
             0, ...
             zeros(1,par.WEC.ny_rad)];

%% Special modifications to base parameters
% Tcoulomb = 1e6; % [Nm] target PTO reaction torque
% par.Dwp = 1; % [m^3/rad] pump displacement
par.control.duty = 0.8;

%% %%%%%%%%%%%%   COLLECT DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% run simulation
tic
out = sim_coulombPTOshuntValve(y0,par);
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
ax(1) = subplot(4,1,1);
hold on
plot(out.t,1e-6*out.T_pto)
plot(out.t,1e-6*out.T_rad)
ylabel('torque (MNm)')
legend('T_{pto}','T_{rad}')

ax(2) = subplot(4,1,2);
hold on
plot(out.t,1e-6*par.p_lout*ones(out.nt,1))
plot(out.t,1e-6*par.p_hin*ones(out.nt,1))
plot(out.t,1e-6*out.p_wpA)
plot(out.t,1e-6*out.p_wpB)
ylabel('Pressure (MPa)')
legend('p_{lout}','p_{hin}','p_{wpA}','p_{wpB}')

ax(3) = subplot(4,1,3);
hold on
plot(out.t,1e3*out.q_wpP)
plot(out.t,1e3*out.q_wpT)
% plot(out.t,1e3*out.q_sv)
ylabel('Flow rate (Lpm)')
legend('q_{wpP}','q_{wpT}','q_{sv}')

ax(4) = subplot(4,1,4);
hold on
plot(out.t,1e3*out.q_sv)
ylabel('Flow rate (Lpm)')
legend('q_{sv}')

linkaxes(ax,'x');

mean(-out.T_pto(:).*out.theta_dot(:))
