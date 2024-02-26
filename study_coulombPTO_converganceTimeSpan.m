% study_coulombPTO_converganceTimeSpan.m script m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 02/26/2024
%
% PURPOSE/DESCRIPTION:
% This script serves as a shell for performing parameter variation studies
% using the model contained in sys_coulombPTO.m and run by solved by
% sim_coulombPTO.m. Specifically this script varies the number of frequency
% componenets used to construct the wave elev and excitation force for the
% WEC model.
% Parameter and function initiallization fuctions are called within this
% script before the sim_coulombPTO.m script is called. 
%
% FILE DEPENDENCY:
% sys_coulombPTO.m
% sim_coulombPTO.m
% parameters_coulombPTO.m
%
% UPDATES:
% 02/26/2024 - Created from study_coulombPTO_dampingStudy.m
%
% Copyright (C) 2024  Jeremy W. Simmons II
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
addpath('Coulomb damping PTO')
addpath('Sea States')
addpath('Solvers')
%% %%%%%%%%%%%%   SIMULATION PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation timeframe
par.Tramp = 250; % [s] excitation force ramp period
par.tstart = 0; %[s] start time of simulation
par.tend = 2000; %[s] end time of simulation

% Solver parameters
par.odeSolverRelTol = 1e-9; % Rel. error tolerance parameter for ODE solver
par.odeSolverAbsTol = 1e-9; % Abs. error tolerance parameter for ODE solver
par.MaxStep = 1e-2;

% Sea State and Wave construction parameters
par.wave.Hs = 2.75;
par.wave.Tp = 12;
par.WEC.nw = 1000; % num. of frequency components for harmonic superposition 
par.wave.rngSeedPhase = 3; % seed for the random number generator

% load parameters
par = parameters_coulombPTO(par,...
    'nemohResults_vantHoff2009_20180802.mat','vantHoffTFCoeff.mat');

% Define initial conditions
y0 = [  0, ...
        0, ...
        zeros(1,par.WEC.ny_rad)];

%% %%%%%%%%%%%%   Study Variables  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.Tcoulomb = 4e6;

nVar = 11;
% Tcoulomb = 1e6*linspace(1,10,nVar1);% [Nm] PTO reaction torque
Tspan = logspace(log10(100),log10(4000),nVar);% [Nm] PTO reaction torque

saveSimData = 0; % save simulation data (1) or just output variables (0)

%% %%%%%%%%%%%%   COLLECT DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parfor iVar = 1:nVar
    param = par; % store individual parameter strucs for parallel compute

    % change design parameter
    param.tend = Tspan(iVar); 

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
xlabel('Torque (MNm)')
title('WEC Power Absorption, Coulomb Damping')
yyaxis left
hold on
plot(Tspan,1e-3*PP)
ylabel('Power (kW)')
% ylim(10*[-pi/2 pi/2])

yyaxis right
hold on
plot(Tspan,theta_dot_ave)
ylabel('Speed, mean (rad/s)')
% ylim()
% xlim([0 2])
