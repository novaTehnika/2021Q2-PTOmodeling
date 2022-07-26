function out = sim_seriesPTO(par)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sim_seriesPTO.m function m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 06/29/2021
%
% PURPOSE/DESCRIPTION:
% This script simulates a wave energy PTO with a long pipeline in the
% high-pressure branch and low-pressure branch. The flow through the
% WEC-driven pump is the input to the system and is determined from a power
% spectrum representing waves. The dynamics of the WEC are ignored and the
% the power spectrum is used directly to calculate the flow rate of the
% WEC-driven pump. A variety of pipeline models are implemented for
% comparison. Accumulators are modeled as linear capacitive elements. The
% low-pressure source node, at the inlet of the low-pressure pipeline, is
% modeled as ... The load on the system is modeled as a linear resistance;
% the resistance of the load is calculated from the mean flow of the
% WEC-driven pump (a predetermined value) and a specified nominal mean
% load pressure. (OR maybe just a specified load resistance)
%
% FILE DEPENDENCY:
% sys_seriesPTO.m
% stateIndex_seriesPTO.m
%
% UPDATES:
% 06/29/2021 - First version creation.
% 07/08/2022 - replaced state indice specification in file with seperate
% script, stateIndex_seriesPTO.m.
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
stateIndex_seriesPTO % load state indices

%% %%%%%%%%%%%%   SOLUTION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Solve for states of the dynamic system
 % Set-up solver
    tspan = [par.tstart par.tend];   % time interval
        
 % Solver options
    options = odeset('RelTol',par.odeSolverRelTol,...
                     'AbsTol',par.odeSolverAbsTol,...
                     'MaxOrder',3);
    options.MaxStep = par.MaxStep;

 % Run solver
%     [t,y] = ode15s(@(t,y) sys(t,y,par), tspan, y0, options);
    
    dt = par.MaxStep;
    y = ode1(@(t,y) sys(t,y,par)',tspan(1),dt,tspan(2),y0);
    t = tspan(1) : dt : tspan(2);
%% %%%%%%%%%%%%   POST-PROCESS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameters
    out.par = par;
    
    % determine number of time steps
    nt = length(t); % number of time steps
    
    % Extract system states from simulation results
    out.t = t;
    out.y = y;
    
    out.p_lin = y(:,iyp_lin);
    out.p_lout = y(:,iyp_lout);
    out.p_hin = y(:,iyp_hin); % [Pa] pressure at inlet of HP pipeline
    out.p_hout = y(:,iyp_hout);
    out.p_RO = y(:,iyp_RO); % [Pa] RO feed pressure
    
    out.w_pm = y(:,iyw_pm); % [rad/s] angular velocity of the pump/motor
    
    out.control.p_filt = y(:,iyp_filt);
    out.control.errInt_p_filt = y(:,iy_errInt_p_filt);
    out.control.errInt_w_pm = y(:,iy_errInt_w_pm);
    
    out.theta = y(:,iytheta); % [rad] position of the WEC
    out.theta_dot = y(:,iytheta_dot); % [rad/s] angular velocity of the WEC
    
    out.yLP = y(:,iyLPPL);
    out.qLP = y(:,iyqLP);
    out.pLP = y(:,iypLP);
    
    out.yLP = y(:,iyHPPL);
    out.qHP = y(:,iyqHP);
    out.pHP = y(:,iypHP);
            
    % Post-process non-state variables and state derivatives
    syspost = @(t,y,par) sysPost(t,y,par);

    parfor it = 1:nt 
        [dydt(it,:), nonState(it), control(it)] = syspost(t(it),y(it,:),par);
        
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
          
     % Control signals
    out.control.w_pm_nom = [control(:).w_pm_nom]';
   
     % torque on pump/motor and generator shafts
    out.Tgen = [control(:).Tgen]';
    out.Tpm = [nonState(:).Tpm]';
    
     % WEC-driven pump flow
    out.q_wp = [nonState(:).q_w]';
    
     % pump/motor flow
    out.q_pm = [nonState(:).q_pm]';
    
     % RO performance
    out.q_perm = [nonState(:).q_perm]';

%% Post-process analysis
    %% Energy analysis
    out.power.P_WEC = -out.T_pto.*out.theta_dot;


    %% Energy Balance
    
    %% Mass Balance    

    
%% %%%%%%%%%%%%   FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function dydt = sys(t,y,par)
    	[dydt, ~ , ~ ] = sys_prototype(t,y,par);
    end

    function [dydt, nonState, control] = sysPost(t,y,par)
    	[dydt, nonState, control] = sys_prototype(t,y,par);
    end

    function [dydt, nonState, control] = sys_prototype(t,y,par)
        [dydt, nonState, control] = sys_seriesPTO(t,y,par);
    end
end
