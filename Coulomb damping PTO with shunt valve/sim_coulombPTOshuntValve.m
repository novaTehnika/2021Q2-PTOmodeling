function out = sim_coulombPTOshuntValve(y0,par)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sim_coulombPTOshuntValve.m function m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 12/29/2021
%
% PURPOSE/DESCRIPTION:
% This script simulates a simple wave energy PTO with coulomb
% damping and PWM shunting. The WEC-driven pump is modeled with 
% compressible fluid filled pumping chambers and a check valve rectifier 
% that is in comunication with constant pressure sources. The shunt valve 
% creates a flow path between the two pumping chambers. The flow 
% coefficient of the shunt valve follows a trapazoidal profile. The model 
% includes the hydrodynamic WEC model.
%
% FILE DEPENDENCY:
% sys_coulombPTOshuntValve.m
% stateIndex_coulombPTOshuntValve.m
% WEC model/flapModel.m
% ode1.m
%
% UPDATES:
% 12/29/2021 - created from sim_coulombPTOcompressiblePump.m.
% 07/08/2022 - replaced state indice specification in file with seperate
% script, stateIndex_coulombPTOshuntValve.m, and added initial condition
% and start time as inputs.
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
stateIndex_coulombPTOshuntValve % load state indices

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
    y = ode1(@(t,y) sys(t,y,par)',tspan(1),dt,tspan(2),y0);
    t = tspan(1) : dt : tspan(2);
%% %%%%%%%%%%%%   POST-PROCESS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameters
    out.par = par;

    % determine number of time steps
    out.nt = length(t); % number of time steps

    % Extract system states from simulation results
    out.t = t;
    out.y = y;
    
    out.p_wpA = y(:,iyp_wpA);
    out.p_wpB = y(:,iyp_wpB);
    
    out.theta = y(:,iytheta); % [rad] position of the WEC
    out.theta_dot = y(:,iytheta_dot); % [rad/s] angular velocity of the WEC
    
    % Post-process non-state variables and state derivatives
    syspost = @(t,y,par) sysPost(t,y,par);
    
    parfor it = 1:out.nt 
        [dydt(it,:), nonState(it)] = syspost(t(it),y(it,:),par); 
        
        % Move WEC torque results up a level in nonState stucture so that
        % they can be used like arrays in assiging to the output structure
        temp(it).T_hydroStatic = nonState(it).torqueWEC.hydroStatic
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

     % Flow rates
    out.q_wpAP = [nonState(:).q_wpAP]';
    out.q_wpBP = [nonState(:).q_wpBP]';
    out.q_wpTA = [nonState(:).q_wpTA]';
    out.q_wpTB = [nonState(:).q_wpTB]';
    out.q_wpP = [nonState(:).q_wpP]'; % net flow out of the WEC-driven pump
    out.q_wpT = [nonState(:).q_wpT]'; % net flow into the WEC-driven pump
    out.q_sv = [nonState(:).q_sv]'; % flow rate through shunt valve, A to B

        
%% Post-process analysis
    %% Energy analysis
    out.power.P_WEC = -out.T_pto.*out.theta_dot;
    out.power.P_wp_mech = (1-out.par.eta_wp).*out.power.P_WEC;
    out.power.P_cvxP = out.q_wpAP.*(out.p_wpA - par.p_hin) ...
                     + out.q_wpBP.*(out.p_wpB - par.p_hin);
    out.power.P_cvTx = out.q_wpTA.*(par.p_lout - out.p_wpA) ...
                     + out.q_wpTB.*(par.p_lout - out.p_wpB);
    out.power.P_cv = out.power.P_cvxP + out.power.P_cvTx;
    out.power.P_sv = out.q_sv.*(out.p_wpA - out.p_wpB);

    %% Energy Balance
    
    %% Mass Balance
    

%% %%%%%%%%%%%%   FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function dydt = sys(t,y,par)
    	[dydt, ~ ] = sys_coulombPTOshuntValve(t,y,par);
    end

    function [dydt, nonState] = sysPost(t,y,par)
    	[dydt, nonState] = sys_coulombPTOshuntValve(t,y,par);
    end

end
