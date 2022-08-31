function par = parameters_coulombPTOcompressiblePump(par,filenameCoeff,filenameRadSS)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters_coulombPTOcompressiblePump.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 12/29/2021
%
% PURPOSE/DESCRIPTION:
% This function loads default parameters into the "par" structure. This
% includes calling a simular function to load parameters for the WEC model.
% File names for the hydrodynamic data for the WEC are passed through that
% function. Not all parameters can be modified after this function runs
% without issue. Changing parameters in any design study script should be
% check for affecting the other parameters.
%
% This function is for a model of a WEC with simple Coulomb damping for a
% PTO that includes the effects of compressibilty of the fluid.
%
% FILE DEPENDENCY:
% parameters_WECmodel.m
%
% UPDATES:
% 12/29/2021 - created.
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
    
     % fluid and entrianed gas properties
    par.rho = 1023; % [kg/m3] density of air
    par.mu = 9.4e-4; % [Pa-s]  Dynamic (absolute) viscosity
    par.beta = 2.2e9; % [Pa]  Bulk Modulus of air free fluid
    par.p_vap = 0.037e5; % [Pa] vapour pressure of seawater
    par.R = 0.0001; % [-] fraction  Baseline fraction of air by volume entrained in the hydraulic fluid at atm
    par.p_o = 101.3e3; % [Pa]  Atmospheric pressure (reference)
    par.gamma = 1.4; % [-]ratio of specific heats for air
    
    % WEC parameters
    par = parameters_WECmodel(par,filenameCoeff,filenameRadSS);
    
    % PTO damping
    par.p_lout = 10e5; % [Pa] Pressure at WEC-driven pump inlet
    par.p_hin = 6e6; % [Pa] Pressure at WEC-driven pump outlet

    % WEC-driven pumpp
    par.D_wp = 1000e-3; % [m^3/rad] kinematic displacement
    par.eta_wp = 0.9; % [-] mechanical efficiency, pumping chamber pressure as input
    par.thetaMax = pi/2; % [rad] Maximum angular travel
    par.V_wpdead = 0;
    
    % Check valves
     % Valve flow coefficient
    par.kv_cvMin = 0; % [m^3/Pa^(1/2)] Minimum valve coefficient
    par.kv_cvTx = 2e-3; % [m^3/Pa^(1/2)] Maximum valve coefficient, pump inlet to pump chamber
    par.kv_cvxP = 2e-3; % [m^3/Pa^(1/2)] Maximum valve coefficient, pump chamber to pump outlet

     % Stroke schedule
      % pump inlet to pump chamber
    par.p_crackTx = 1e5; % [Pa] Cracking pressure
    par.dp_strokeTx = 2e5; % [Pa] pressure difference between cracking and max stroke

      % pump chamber to pump outlet
    par.p_crackxP = 1e5; % [Pa] Cracking pressure
    par.dp_strokexP = 2e5; % [Pa] pressure difference between cracking and max stroke

end