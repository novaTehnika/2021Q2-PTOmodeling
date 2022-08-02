% study_coulombPTO_dampingStudy_singleSS.m script m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 08/02/2021
%
% PURPOSE/DESCRIPTION:
% This script serves as a shell for performing parameter variation studies
% using the model contained in sys_coulombPTO.m and run by solved by 
% sim_coulombPTO.m. Specifically this script varies the Coulomb damping
% torque for a given sea state to produce a mean power curve.
% Parameter and function initiallization fuctions are called within this
% script before the sim_coulombPTO.m script is called. 
%
% FILE DEPENDENCY:
% sys_coulombPTO.m
% sim_coulombPTO.m
% parameters_coulombPTO.m
%
% UPDATES:
% 08/02/2021 - Created.
% 07/08/2022 - added initial conditions as arguement to sim_coulombPTO()
% and simulation start time as parameter.
% 08/02/2022 - added ramp period to par struct and included tstart so that
% ramp ends at t=0 and modified average power and WEC velocity calculations
% to exclude ramp period.
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
addpath('Coulomb damping PTO') 
%% %%%%%%%%%%%%   SIMULATION PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation timeframe
par.tramp = 100; % [s] excitation force ramp period
par.tstart = 0-par.tramp; %[s] start time of simulation
par.tend = 3000; %[s] end time of simulation

% Solver parameters
par.odeSolverRelTol = 1e-9; % Rel. error tolerance parameter for ODE solver
par.odeSolverAbsTol = 1e-9; % Abs. error tolerance parameter for ODE solver
par.MaxStep = 1e-2;

% Sea State and Wave construction parameters
par.wave.Hs = 2.75;
par.wave.Tp = 12;
par.WEC.nw = 100; % num. of frequency components for harmonic superposition 
par.wave.rngSeedPhase = 3; % seed for the random number generator

% load parameters
par = parameters_coulombPTO(par);

% Define initial conditions
y0 = [  0, ...
        0, ...
        zeros(1,par.WEC.ny_rad)];

%% %%%%%%%%%%%%   Study Variables  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nVar = 31;
% Tcoulomb = 1e6*linspace(1,10,nVar1);% [Nm] PTO reaction torque
Tcoulomb = 1e6*logspace(log10(0.1),log10(20),nVar);% [Nm] PTO reaction torque

saveSimData = 1; % save simulation data (1) or just output variables (0)

%% %%%%%%%%%%%%   COLLECT DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parfor iVar = 1:nVar
    param = par; % store individual parameter strucs for parallel compute

    % change design parameter
    param.Tcoulomb = Tcoulomb(iVar); 
    
    % run simulation
    tic
    out = sim_coulombPTO(y0,param);
    toc
    
    % Calculate metrics
    t_vec = find(out.t>=par.tstart);
    PP(iVar) = mean(-out.T_pto(t_vec).*out.theta_dot(t_vec));
    theta_dot_ave(iVar) = mean(abs(out.theta_dot(t_vec)));

    if saveSimData
        simOut(iVar) = out;
    end

end
%% %%%%%%%%%%%%   PLOTTING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
xlabel('Torque (MPa)')
title('WEC Power Absorption, Coulomb damping')
yyaxis left
hold on
plot(1e-6*Tcoulomb,1e-3*PP)
ylabel('Power (kW)')
% ylim(10*[-pi/2 pi/2])

yyaxis right
hold on
plot(1e-6*Tcoulomb,theta_dot_ave)
ylabel('Speed, mean (rad/s)')
% ylim()
% xlim([0 2])
