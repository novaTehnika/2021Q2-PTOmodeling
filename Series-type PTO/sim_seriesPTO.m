function out = sim_seriesPTO(y0,par)
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
% 05/12/2023 - added specification for down sampling in the solver
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
% Define indicies of state vector
iyp_lin = [];
iyp_lout = [];
iyp_hin = [];
iyp_hout = [];
iyp_RO = [];

iyw_pm = [];

iyp_filt = [];
iy_errInt_p_filt = [];
iy_errInt_w_pm = [];
iycontrol = [];

iytheta = [];
iytheta_dot = [];
iyrad = [];

iyLPPL = [];
iyqLP = [];
iypLP = [];
iyHPPL = [];
iyqHP = [];
iypHP = [];

stateIndex_seriesPTO % load state indices

%% %%%%%%%%%%%%   SOLUTION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Solve for states of the dynamic system
 % Set-up solver
    tspan = [par.tstart-par.Tramp par.tend];   % time interval
    
 % Solver options
%     options = odeset('RelTol',par.odeSolverRelTol,...
%                      'AbsTol',par.odeSolverAbsTol,...
%                      'MaxOrder',3);
%     options.MaxStep = par.MaxStep;
    dt = par.MaxStep;
    downSampleRate = floor(par.downSampledStepSize/dt);

 % Run solver
    [t, y] = ode1(@(t,y) sys(t,y,par)', ...
                                tspan(1),dt,tspan(2),y0,downSampleRate);

%% %%%%%%%%%%%%   POST-PROCESS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameters
    out.par = par;
    
    % Select desired time indices
    itVec = find(t >= par.tstart);
    
    % Extract system states from simulation results post ramp
    out.t = t(itVec);
    out.y = y(itVec,:);
    
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

    nt_ramp = itVec(1)-1;
    parfor it = 1:length(itVec)
        
        [dydt(it,:), nonState(it), control(it)] = ...
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
    out.q_brine = [nonState(:).q_brine]';
    out.q_feed = out.q_perm + out.q_brine;

     % ERU flow rate
    out.q_ERUfeed = [nonState(:).q_ERUfeed]';

     % Charge Pump
    out.q_c = [nonState(:).q_c]';

     % Pressure relief valves
    out.q_linPRV = [nonState(:).q_linPRV]';
    out.q_loutPRV = [nonState(:).q_loutPRV]';
    out.q_hinPRV = [nonState(:).q_hinPRV]';
    out.q_houtPRV = [nonState(:).q_houtPRV]';
    out.q_ROPRV = [nonState(:).q_ROPRV]';

%% Post-process analysis

    %% Energy analysis
    % WEC-drive pump
    out.power.P_WEC = -out.T_pto.*out.theta_dot;
    out.power.P_wpLoss = out.power.P_WEC ...
                       - out.q_wp.*(out.p_hin-out.p_lout);

    % Pump/motor and generator
    out.power.P_pmLoss = (out.p_RO - out.p_hout).*out.q_pm ...
                       - out.Tpm.*out.w_pm;
    out.power.P_gen = -((-out.Tgen.*out.w_pm > 0).*par.eta_g ...
                    + (-out.Tgen.*out.w_pm < 0)./par.eta_g) ...
                    .*out.Tgen.*out.w_pm;
    out.power.P_genLoss = -out.Tgen.*out.w_pm - out.power.P_gen;
    out.power.deltaE_pmgenFlywheel = 0.5*out.par.Jpm*(out.w_pm(end).^2-out.w_pm(1).^2);

    % Charge pump
    out.power.P_cElec = 1/(par.eta_c*par.eta_m)*out.q_c.*(out.p_lin-out.par.p_o);
    out.power.P_cLoss = out.power.P_cElec - out.q_c.*(out.p_lin-out.par.p_o);

    % ERU
    P_ERUfeed = out.q_ERUfeed.*(out.p_RO - out.p_lin);
    P_ERUbrine = out.q_brine.*(out.p_RO - out.par.p_o);
    PbalERU = 1/(par.eta_ERUv*par.eta_ERUm)*P_ERUfeed ...
            - par.eta_ERUv*par.eta_ERUm*P_ERUbrine;
    out.power.P_ERULoss = PbalERU + P_ERUbrine - P_ERUfeed;
    out.power.P_ERUelec = ((PbalERU > 0)./par.eta_m ...
                        + (PbalERU < 0).*par.eta_m) ...
                        .*PbalERU;
    out.power.P_ERUelecLoss = out.power.P_ERUelec - PbalERU;

    % Pipelines
    out.power.P_pLineLossLP = [nonState.pLinePPfricLP]';
    out.power.P_pLineLossHP = [nonState.pLinePPfricHP]';

    % Pressure relief valves
    out.power.P_linPRV = out.q_linPRV.*(out.p_lin-out.par.p_o);
    out.power.P_loutPRV = out.q_loutPRV.*(out.p_lout-out.par.p_o);
    out.power.P_hinPRV = out.q_hinPRV.*(out.p_hin-out.p_lout);
    out.power.P_houtPRV = out.q_houtPRV.*(out.p_hout-out.p_lin);
    out.power.P_ROPRV = out.q_ROPRV.*(out.p_RO-out.p_lin);

    % Electrical Energy Storage
    out.power.P_battery = out.power.P_gen ...
                        - out.power.P_cElec - out.power.P_ERUelec;
    out.power.deltaE_battery = trapz(out.t,out.power.P_battery);

    %% Energy Balance
    % Change in available potential energy in accumulators
    dp = 1e2;
     % Low-pressure onshore
    cap = @(p) capAccum(p,par.pc_lin,par.Vc_lin,par.f,par) + lineCap(p,1,par);
    deltaE_lin = deltaE_NI(out.p_lin(1),out.p_lin(end),cap,dp);

     % Low-pressure offshore
    cap = @(p) capAccum(p,par.pc_lout,par.Vc_lout,par.f,par) + lineCap(p,1,par);
    deltaE_lout = deltaE_NI(out.p_lout(1),out.p_lout(end),cap,dp);

     % High-pressure offshore
    cap = @(p) capAccum(p,par.pc_hin,par.Vc_hin,par.f,par) + lineCap(p,2,par);
    deltaE_hin = deltaE_NI(out.p_hin(1),out.p_hin(end),cap,dp);

     % High-pressure onshore
    cap = @(p) capAccum(p,par.pc_hout,par.Vc_hout,par.f,par) + lineCap(p,2,par);
    deltaE_hout = deltaE_NI(out.p_hout(1),out.p_hout(end),cap,dp);

     % RO feed inlet
    cap = @(p) capAccum(p,par.pc_RO,par.Vc_RO,par.f,par);
    deltaE_RO = deltaE_NI(out.p_RO(1),out.p_RO(end),cap,dp);

    % Change in available potential and kinetic energy in pipelines
     % Low-pressure
    cap = @(p) lineCap(p,1,par);
    deltaE_LP = 0;
    for i = 1:par.n_seg(1)-1
        deltaE_LP = deltaE_LP + deltaE_NI(out.pLP(1,i),out.pLP(end,i),cap,dp);
    end
    deltaE_LP = deltaE_LP + 0.5*par.I(1)*sum(out.qLP(end,:).^2 - out.qLP(1,:).^2);

     % High-pressure
    cap = @(p) lineCap(p,2,par);
    deltaE_HP = 0;
    for i = 1:par.n_seg(1)-1
        deltaE_HP = deltaE_HP + deltaE_NI(out.pHP(1,i),out.pHP(end,i),cap,dp);
    end
    deltaE_HP = deltaE_HP + 0.5*par.I(1)*sum(out.qLP(end,:).^2 - out.qLP(1,:).^2);

    deltaE_HP_trap = trapz(out.t,(out.qHP(:,1:end-1)-out.qHP(:,2:end)).*(out.pHP(:,:)-out.par.p_o)) ...
                   + 0.5*par.I(2)*sum(out.qHP(end,:).^2 - out.qHP(1,:).^2);

    % Total change in stored energy in the system
    deltaE_sys = deltaE_LP + deltaE_HP ...
        + deltaE_lin + deltaE_lout + deltaE_hin + deltaE_hout + deltaE_RO ...
        + out.power.deltaE_battery + out.power.deltaE_pmgenFlywheel;

    % Power flow at boundaries
    P_in = out.power.P_WEC;
    P_out = out.q_perm.*(out.p_RO - out.par.p_o);
    P_loss = out.power.P_wpLoss + out.power.P_pmLoss + out.power.P_genLoss ...
            + out.power.P_cLoss + out.power.P_ERULoss + out.power.P_ERUelecLoss...
            + out.power.P_pLineLossLP + out.power.P_pLineLossHP ...
            + out.power.P_linPRV + out.power.P_loutPRV ...
            + out.power.P_hinPRV + out.power.P_houtPRV + out.power.P_ROPRV;
    P_bnds = P_in - P_out - P_loss;

    % Total balance of energy: energy added minus change in energy stored
    out.Ebal = trapz(out.t,P_bnds) - deltaE_sys;
    out.Ebal_error = out.Ebal/trapz(out.t,P_in);

    %% Mass Balance
    % Change in fluid volume in accumulators
    dp = 1e2;
     % Low-pressure onshore
    cap = @(p) capAccum(p,par.pc_lin,par.Vc_lin,par.f,par) + lineCap(p,1,par);
    deltaV_lin = deltaV_NI(out.p_lin(1),out.p_lin(end),cap,dp);

     % Low-pressure offshore
    cap = @(p) capAccum(p,par.pc_lout,par.Vc_lout,par.f,par) + lineCap(p,1,par);
    deltaV_lout = deltaV_NI(out.p_lout(1),out.p_lout(end),cap,dp);

     % High-pressure offshore
    cap = @(p) capAccum(p,par.pc_hin,par.Vc_hin,par.f,par) + lineCap(p,2,par);
    deltaV_hin = deltaV_NI(out.p_hin(1),out.p_hin(end),cap,dp);

     % High-pressure onshore
    cap = @(p) capAccum(p,par.pc_hout,par.Vc_hout,par.f,par) + lineCap(p,2,par);
    deltaV_hout = deltaV_NI(out.p_hout(1),out.p_hout(end),cap,dp);

     % RO feed inlet
    cap = @(p) capAccum(p,par.pc_RO,par.Vc_RO,par.f,par);
    deltaV_RO = deltaV_NI(out.p_RO(1),out.p_RO(end),cap,dp);

    % Change in stored fluid in pipelines
     % Low-pressure
    cap = @(p) lineCap(p,1,par);
    deltaV_LP = 0;
    for i = 1:par.n_seg(1)-1
        deltaV_LP = deltaV_LP + deltaV_NI(out.pLP(1,i),out.pLP(end,i),cap,dp);
    end
     % High-pressure
    cap = @(p) lineCap(p,2,par);
    deltaV_HP = 0;
    for i = 1:par.n_seg(1)-1
        deltaV_HP = deltaV_HP + deltaV_NI(out.pHP(1,i),out.pHP(end,i),cap,dp);
    end

    % Total change in stored volume in the system
    out.deltaV_total = deltaV_LP + deltaV_HP + deltaV_lin + deltaV_lout ...
                     + deltaV_hin + deltaV_hout + deltaV_RO;

    % Flow at boundaries
    qbnds = out.q_c ...
            - (out.q_perm + out.q_brine + out.q_linPRV + out.q_loutPRV);

    % Total balance of flow: flow in minus change in volume stored
    out.Vbal = trapz(out.t,qbnds) - out.deltaV_total;
    out.Vbal_error = out.Vbal/trapz(out.t,out.q_wp);

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
