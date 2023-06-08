function out = sim_coulombPTO(y0,par)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sim_coulombPTO.m function m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 06/29/2021
%
% PURPOSE/DESCRIPTION:
% This script simulates a simple wave energy PTO with coulomb
% damping.
%
% FILE DEPENDENCY:
% sys_coulombPTO.m
% WEC model/flapModel.m
% ode1.m
%
% UPDATES:
% 06/29/2021 - Created.
% 12/5/2021 - Assignment to the struc "out" moved to replace intermediate 
% assignments.
% 07/08/2022 - added initial conditions as arguement and simulation start 
% time as parameter.
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


%% %%%%%%%%%%%%   SOLUTION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Solve for states of the dynamic system
 % Set-up solver
    tspan = [par.tstart-par.Tramp par.tend];   % time interval

 % Solver options
    options = odeset('RelTol',par.odeSolverRelTol,...
                     'AbsTol',par.odeSolverAbsTol,...
                     'MaxOrder',3);
    options.MaxStep = par.MaxStep; % min(0.1*par.WEC.w(:));

 % Run solver
%     [t,y] = ode23(@(t,y) sys(t,y,par), tspan, y0, options);
%     [t,y] = ode15s(@(t,y) sys(t,y,par), tspan, y0, options);
    
    dt = par.MaxStep;
    [t, y] = ode1(@(t,y) sys(t,y,par)',tspan(1),dt,tspan(2),y0);
    
%% %%%%%%%%%%%%   POST-PROCESS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameters
    out.par = par;

    % Select desired time indices
    itVec = find(t >= par.tstart);
    
    % Extract system states from simulation results
    out.t = t(itVec);
    out.y = y(itVec,:);
    
    out.theta = y(itVec,1); % [rad] position of the WEC
    out.theta_dot = y(itVec,2); % [rad/s] angular velocity of the WEC
    
    % Post-process non-state variables and state derivatives
    syspost = @(t,y,par) sysPost(t,y,par);

    nt_ramp = itVec(1)-1;
    parfor it = 1:length(itVec)

        [dydt(it,:), nonState(it)] = ...
                            syspost(t(it+nt_ramp),y(it+nt_ramp,:),par); 
        
        % Move WEC torque results up a level in nonState stucture so that
        % they can be used like arrays in assiging to the output structure
        temp(it).T_hydroStatic = nonState(it).torqueWEC.hydroStatic;
        temp(it).T_wave = nonState(it).torqueWEC.wave;
        temp(it).T_rad = nonState(it).torqueWEC.radiation;  
                    
    end

    % State derivatives
    out.dydt = dydt;
    
     % System input: wave elevation
    out.waveElev = [nonState(:).waveElev]';
    
     % forces on WEC
    out.T_pto = [nonState(:).T_pto]';
    out.T_hydroStatic = [temp(:).T_hydroStatic]';
    out.T_wave = [temp(:).T_wave]';
    out.T_rad = [temp(:).T_rad]';
        
%% Post-process analysis
    % Energy analysis
    out.power.P_WEC = -out.T_pto.*out.theta_dot;

    % Energy Balance
    
    % Mass Balance


%% %%%%%%%%%%%%   FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function dydt = sys(t,y,par)
    	[dydt, ~ ] = sys_coulombPTO(t,y,par);
    end

    function [dydt, nonState] = sysPost(t,y,par)
    	[dydt, nonState] = sys_coulombPTO(t,y,par);
    end

end
