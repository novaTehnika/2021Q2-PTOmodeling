% run_seriesPTO.m script m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 07/12/2021
%
% PURPOSE/DESCRIPTION:
% This script serves as a shell for running a single simulation
% using the model contained in sys_seriesPTO.m and run by solved by 
% sim_seriesPTO.m.
% The parameter initiallization fuction are called within this
% script before the sim_seriesPTO.m script is called.
%
% FILE DEPENDENCY:
% sys_seriesPTO.m
% sim_seriesPTO.m
% parameters_seriesPTO.m
% stateIndex_seriesPTO.m
% initialConditionDefault_seriesPTO.m
%
% UPDATES:
% 07/12/2021 - Created.
% 07/08/2022 - Added intial conditions and start time (in parameter struc) 
% as inputs to sim_seriesPTO.m
% 08/02/2022 - added ramp period to par struct and included tstart so that
% ramp ends at t=0;
% 05/12/2023 - added specification for down sampling with check for being
% multiple of step size for fixed time step solvers
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
addpath('Series-type PTO')
addpath('Components')
addpath(['Components' filesep 'Pipeline'])
addpath('Sea States')
addpath('Solvers')
%% %%%%%%%%%%%%   SIMULATION PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation timeframe
par.Tramp = 250; % [s] excitation force ramp period
par.tstart = 0; %[s] start time of simulation
par.tend = 2000; %[s] end time of simulation

% Solver parameters
% par.odeSolverRelTol = 1e-4; % Rel. error tolerance parameter for ODE solver
% par.odeSolverAbsTol = 1e-4; % Abs. error tolerance parameter for ODE solver
par.MaxStep = 1e-4;             % [s] for fixed time solver, this is the step size for solver
par.downSampledStepSize = 1e-2; % [s] specifies time step for data output
if mod(par.downSampledStepSize,par.MaxStep)
    warning('down-sampled time step is not an integer multiple of the maximum step size')
end

% Sea State and Wave construction parameters
par.wave.Hs = 2.75;
par.wave.Tp = 12;
par.WEC.nw = 1000; % num. of frequency components for harmonic superposition 
par.wave.rngSeedPhase = 3; % seed for the random number generator

% load parameters
par = parameters_seriesPTO(par,...
    'nemohResults_vantHoff2009_20180802.mat','vantHoffTFCoeff.mat');

% Define initial conditions
stateIndex_seriesPTO % load state indices, provides 'iy_...'
initialConditionDefault_seriesPTO % default ICs, provides 'y0'

%% Special modifications to base parameters
par.Sro = 5000; % [m^3]
par.D_WEC = 0.3;         % [m^3/rad] flap pump displacement
par.control.p_hout_nom = 7e6; % [Pa]

%% %%%%%%%%%%%%   COLLECT DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
out = sim_seriesPTO(y0,par);
toc

%% %%%%%%%%%%%%   PLOTTING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% WEC pos., vel., and Torque
bottomEdge = 1;
leftEdge = 2;
width = 7.5; % one column: 3+9/16, two column: 7.5
height = 6;
fontSize = 8;
lineWidth = 1;

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];


ax(1) = subplot(3,1,1);
plot(out.t,out.waveElev)
xlabel('time (s)')
ylabel('elevation (m)')
title('Wave Elevation')

ax(2) = subplot(3,1,2);
xlabel('time (s)')
hold on

yyaxis left
plot(out.t,out.theta)
xlabel('time (s)')
ylabel('position (rad)')
% ylim([-pi/2 pi/2])

yyaxis right
plot(out.t,out.theta_dot)
ylabel('angular velocity (rad/s)')
% ylim(10*[-pi/2 pi/2])

ax(3) = subplot(3,1,3);
hold on
plot(out.t,1e-6*out.T_wave)
plot(out.t,1e-6*out.T_pto)
plot(out.t,1e-6*out.T_rad)
plot(out.t,1e-6*out.T_hydroStatic)
ylabel('Torque (MNm)')

legend('T_{w}','T_{PTO}','T_{r}','T_{h}')
      
linkaxes(ax,'x');

sgtitle('WEC Behaviour')

%% Pressure
bottomEdge = 1;
leftEdge = 3;
width = 7.5; % one column: 3+9/16, two column: 7.5
height = 6;
fontSize = 8;
lineWidth = 1;

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

ax(1) = subplot(2,1,1);
plot(out.t,out.p_hin)
hold on
plot(out.t,out.p_hout)
plot(out.t,out.p_RO)
xlabel('Time (s)')
ylabel('Pressure (Pa)')
legend('p_{hin}','p_{hout}','p_{RO}')

ax(2) = subplot(2,1,2);
plot(out.t,out.p_lin)
hold on
plot(out.t,out.p_lout)
xlabel('Time (s)')
ylabel('Pressure (Pa)')
legend('p_{lin}','p_{lout}')

linkaxes(ax,'x');

sgtitle('Pressures')

%% Flow rates
bottomEdge = 1;
leftEdge = 3;
width = 7.5; % one column: 3+9/16, two column: 7.5
height = 6;
fontSize = 8;
lineWidth = 1;

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

ax(1) = subplot(3,1,1);
plot(out.t,out.q_pm)
hold on
plot(out.t,out.q_perm)
xlabel('Time (s)')
ylabel('Flow rate (m^3/s)')
legend('q_{pm}','q_{perm}')

ax(2) = subplot(3,1,2);
plot(out.t,out.q_wp(:))
hold on
plot(out.t,out.qHP(:,1))
plot(out.t,out.qHP(:,end))
xlabel('Time (s)')
ylabel('Flow rate (m^3/s)')
legend('q_{w}','q_{h,in}','q_{h,out}')

ax(3) = subplot(3,1,3);
plot(out.t,out.q_wp(:))
hold on
plot(out.t,out.qLP(:,1))
plot(out.t,out.qLP(:,end))
xlabel('Time (s)')
ylabel('Flow rate (m^3/s)')
legend('q_{w}','q_{l,out}','q_{l,in}')

linkaxes(ax,'x')

sgtitle('Flow rates')


%% Controller behavior
bottomEdge = 1;
leftEdge = 3;
width = 7.5; % one column: 3+9/16, two column: 7.5
height = 8;
fontSize = 8;
lineWidth = 1;

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

ax(1) = subplot(4,1,1);
hold on
plot(out.t,1e-6*out.p_hout)
plot(out.t,1e-6*out.control.p_filt)
plot(out.t,1e-6*out.p_RO)
xlabel('Time (s)')
ylabel('Pressure (MPa)')
legend('p_{hout}','p_{filt}','p_{RO}')

ax(2) = subplot(4,1,2);
hold on
yyaxis left
plot(out.t,60/(2*pi)*out.control.w_pm_nom)
plot(out.t,60/(2*pi)*out.w_pm)
ylabel('Shaft speed (rpm)')
yyaxis right
plot(out.t,1e-3*out.Tgen)
ylabel('Torque (kNm)')
legend('nominal','actual','Generator')

ax(3) = subplot(4,1,3);
hold on
plot(out.t,1e3*out.q_pm)
plot(out.t,1e3*out.q_wp)
ylabel('Flow rate (Lpm)')
legend('q_{pm}','q_{w}')

ax(4) = subplot(4,1,4);
hold on
yyaxis left
plot(out.t,out.control.errInt_p_filt)
ylabel('Error integral')
yyaxis right
plot(out.t,out.control.errInt_w_pm)
legend('pressure control','speed control')
ylabel('Error integral')
xlabel('Time (s)')

linkaxes(ax,'x')

sgtitle('Controller behaviour')
