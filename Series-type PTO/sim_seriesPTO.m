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
% WEC model/flapModel.m
% ode1.m
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
    tspan = [par.tstart-par.Tramp par.tend];   % time interval
        
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
    itVec = find(out.t >= par.tstart);
    
    % Extract system states from simulation results
    out.t = t(itVec);
    out.y = y(itVec,);
    
    out.p_lin = y(itVec,iyp_lin);
    out.p_lout = y(itVec,iyp_lout);
    out.p_hin = y(itVec,iyp_hin); % [Pa] pressure at inlet of HP pipeline
    out.p_hout = y(itVec,iyp_hout);
    out.p_RO = y(itVec,iyp_RO); % [Pa] RO feed pressure
    
    out.w_pm = y(itVec,iyw_pm); % [rad/s] angular velocity of the pump/motor
    
    out.control.p_filt = y(itVec,iyp_filt);
    out.control.errInt_p_filt = y(itVec,iy_errInt_p_filt);
    out.control.errInt_w_pm = y(itVec,iy_errInt_w_pm);
    
    out.theta = y(itVec,iytheta); % [rad] position of the WEC
    out.theta_dot = y(itVec,iytheta_dot); % [rad/s] angular velocity of the WEC
    
    out.yLP = y(itVec,iyLPPL);
    out.qLP = y(itVec,iyqLP);
    out.pLP = y(itVec,iypLP);
    
    out.yLP = y(itVec,iyHPPL);
    out.qHP = y(itVec,iyqHP);
    out.pHP = y(itVec,iypHP);
            
    % Post-process non-state variables and state derivatives
    syspost = @(t,y,par) sysPost(t,y,par);

    parfor it = itVec 
        [dydt(it,:), nonState(it), control(it)] = syspost(t(it),y(it,:),par);
        
        % Move WEC torque results up a level in nonState stucture so that
        % they can be used like arrays in assiging to the output structure
        temp(it).T_hydroStatic = nonState(it).torqueWEC.hydroStatic
        temp(it).T_wave = nonState(it).torqueWEC.wave;
        temp(it).T_rad = nonState(it).torqueWEC.radiation;    
    end
    
    % State derivatives
    out.dydt = dydt(itVec);
    
     % System input: wave elevation
    out.waveElev = [nonState(itVec).waveElev]';
    
     % forces on WEC
    out.T_pto = [nonState(itVec).T_pto]';
    out.T_hydroStatic = [temp(itVec).T_hydroStatic]';
    out.T_wave = [temp(itVec).T_wave]';
    out.T_rad = [temp(itVec).T_rad]';
          
     % Control signals
    out.control.w_pm_nom = [control(itVec).w_pm_nom]';
   
     % torque on pump/motor and generator shafts
    out.Tgen = [control(itVec).Tgen]';
    out.Tpm = [nonState(itVec).Tpm]';
    
     % WEC-driven pump flow
    out.q_wp = [nonState(itVec).q_w]';
    
     % pump/motor flow
    out.q_pm = [nonState(itVec).q_pm]';
    
     % RO performance
    out.q_perm = [nonState(itVec).q_perm]';

     % ERU flow rate
    out.q_brine = [nonState(itVec).q_brine]';
    out.q_ERUfeed = [nonState(itVec).q_ERUfeed]';

     % Charge Pump
    out.q_c = [nonState(itVec).q_c]';

     % Pressure relief valves
    out.q_linPRV = [nonState(itVec).q_linPRV]';
    out.q_loutPRV = [nonState(itVec).q_loutPRV]';
    out.q_hinPRV = [nonState(itVec).q_hinPRV]';
    out.q_houtPRV = [nonState(itVec).q_houtPRV]';
    out.q_ROPRV = [nonState(itVec).q_ROPRV]';

%% Post-process analysis

    %% Energy analysis
    % WEC-drive pump
    out.power.P_WEC = -out.T_pto.*out.theta_dot;
    out.power.P_wpLoss = out.T_pto.*out.theta_dot ...
                       - out.q_wp.*(out.p_hin-out.p_lout);

    % Pump/motor and generator
    out.power.P_pmLoss = (out.p_hout - out.p_RO).*out.q_pm ...
                       - out.Tpm.*out.w_pm;
    out.power.P_gen = ((out.Tgen.*out.w_pm > 0).*par.eta_g ...
                    + (out.Tgen.*out.w_pm < 0)./par.eta_g) ...
                    .*out.Tgen.*out.w_pm;
    out.power.P_genLoss = out.Tgen.*out.w_pm - out.power.P_gen;

    % Charge pump
    out.power.P_cElec = 1/(par.eta_c*par.eta_m)*out.q_c.*(out.p_lin-out.p_o);
    out.power.P_cLoss = out.power.P_c - out.q_c.*(out.p_lin-out.p_o);

    % ERU
    PbalERU = 1/(par.eta_ERUv*par.eta_ERUm)*out.q_ERUfeed.*(out.p_RO - out.p_lin) ...
            - par.eta_ERUv*par.eta_ERUm*out.q_brine.*(out.p_RO - out.p_o);
    out.power.P_eruElec = ((PbalERU > 0)./par.eta_m ...
                        + (PbalERU < 0).*par.eta_m) ...
                        .*PbalERU;
    out.power.P_ERUloss = out.power.P_eruElec - PbalERU;

    % Pipelines
    out.power.P_pLineLossLP = [nonState.pLinePPfricLP]';
    out.power.P_pLineLossHP = [nonState.pLinePPfricHP]';

    % Pressure relief valves
    out.power.P_linPRV = out.q_linPRV.*(out.p_lin-out.p_o);
    out.power.P_loutPRV = out.q_loutPRV.*(out.p_lout-out.p_o);
    out.power.P_hinPRV = out.q_hinPRV.*(out.p_hin-out.p_lout);
    out.power.P_houtPRV = out.q_houtPRV.*(out.p_hout-out.p_lin);
    out.power.P_ROPRV = out.q_ROPRV.*(out.p_RO-out.p_lin);

    % Electrical Energy Storage
    out.power.P_battery = out.power.P_gen ...
                        - out.power.P_cElec - out.power.P_eruElec;
    out.power.deltaE_battery = trapz(out.t,out.power.P_battery);

    %% Energy Balance
    deltaE_lin = par.pc_lin*par.Vc_lin...
                *log(out.p_lin(end)/out.p_lin(1));
    deltaE_lout = par.pc_lout*par.Vc_lout...
                *log(out.p_lout(end)/out.p_lout(1));
    deltaE_hin = par.pc_hin*par.Vc_hin...
                *log(out.p_hin(end)/out.p_hin(1));
    deltaE_hout = par.pc_hout*par.Vc_hout...
                *log(out.p_hout(end)/out.p_hout(1));
    deltaE_RO = par.pc_RO*par.Vc_RO...
                *log(out.p_RO(end)/out.p_RO(1));

    P_bnds = out.power.P_WEC - out.q_perm.*(out.p_RO - par.p_o) ...
            - out.power.P_wpLoss - out.power.P_pmLoss - out.power.P_genLoss ...
            - out.power.P_cLoss - out.power.P_ERUloss ...
            - out.power.P_pLineLossLP - out.power.P_pLineLossHP ...
            - out.power.P_linPRV - out.power.P_loutPRV ...
            - out.power.P_hinPRV - out.power.P_houtPRV - out.power.P_ROPRV;


    out.Ebal = deltaE_lin + deltaE_lout + deltaE_hin + deltaE_hout + deltaE_RO ...
           + out.power.deltaE_battery + trapz(out.t,P_bnds);

    %% Mass Balance
    deltaV_lin = (1./out.p_lin(1) - 1./out.p_lin(end)) ...
                 *par.pc_lin*par.Vc_lin;
    deltaV_lout = (1./out.p_lout(1) - 1./out.p_lout(end)) ...
                 *par.pc_lout*par.Vc_lout;
    deltaV_hin = (1./out.p_hin(1) - 1./out.p_hin(end)) ...
                 *par.pc_hin*par.Vc_hin;
    deltaV_hout = (1./out.p_hout(1) - 1./out.p_hout(end)) ...
                 *par.pc_hout*par.Vc_hout;
    deltaV_RO = (1./out.p_RO(1) - 1./out.p_RO(end)) ...
                 *par.pc_RO*par.Vc_RO;

    out.deltaV_total = deltaV_lin + deltaV_lout ...
                     + deltaV_hin + deltaV_hout + deltaV_RO;

    out.Vbal = out.deltaV_total + trapz(out.t,out.q_perm) ...
               - trapz(out.t,out.q_c);

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
