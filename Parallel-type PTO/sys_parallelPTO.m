function [dydt, nonState, control] = sys_parallelPTO(t,y,par)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sys_parallelPTO.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 6/8/2023
%
% PURPOSE/DESCRIPTION:
% Calculate the state derivatives for a wave energy PTO with a parallel-
% type architecture. The model includes the hydrodynamic WEC
% model found in the 2021Q2 Hydrodynamic modeling project.
%
% FILE DEPENDENCY: 
% pipelineNPi.m
% nonStateVarsPTOwPL.m
%
% UPDATES:
% 6/8/2023 - created from sys_seriesPTO.m.
%
% Copyright (C) 2023  Jeremy W. Simmons II
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

stateIndex_parallelPTO

iWEC = [iytheta iytheta_dot iyrad];


% Calculate control input and the change in controller states (if any)
[dydt_control, control] = controller(t,y,par);

% Calculate non-state variables like coeffs., forces, and flow rates
nonState = nonStateVars(t,y,par);

% Calculate the hydrodynamics WEC state derivatives and output nonstate
% variables like forces and wave elevation
[dydt_WEC, nonState.torqueWEC, nonState.waveElev] = ...
                    flapModel(t,y(iWEC),nonState.T_pto,par);

% Calculate the pipeline state derivatives
 % low-pressure pipeline
[q_inLP, q_outLP, dydt_pLineLP, pLsoln] = ...
        pipelineNPi(t,y(iyLPPL),y(iyp_lin),y(iyp_lout),par,1,1);
nonState.pLinePPfricLP = pLsoln.PPfric;

 % high-pressure pipeline
[q_inHP, q_outHP, dydt_pLineHP, pLsoln] = ...
            pipelineNPi(t,y(iyHPPL),y(iyp_hin),y(iyp_hout),par,2,1);
nonState.pLinePPfricHP = pLsoln.PPfric;

% State derivatives
dydt = zeros(11 + par.WEC.ny_rad,1);

dydt(iyp_lin) = 1/nonState.C_lin*(nonState.q_c - q_inLP ...
                + nonState.q_pm - nonState.q_ERUfeed ...
                - nonState.q_linPRV + nonState.q_houtPRV);
dydt(iyp_lout) = 1/nonState.C_lout*(q_outLP - nonState.q_w ...
                - nonState.q_loutPRV + nonState.q_hinPRV);
dydt(iyp_hin) = 1/nonState.C_hin*(nonState.q_w - q_inHP ...
                - nonState.q_hinPRV);
dydt(iyp_hout) = 1/nonState.C_hout*(q_outHP - nonState.q_pm ...
                - nonState.q_houtPRV ...
                - nonState.q_feed + nonState.q_ERUfeed);

dydt(iyw_pm) = 1/par.Jpm*(nonState.Tpm - control.Tgen);

dydt(iycontrol) = dydt_control;

dydt(iytheta) = dydt_WEC(1); % angular velocity
dydt(iytheta_dot) = dydt_WEC(2); % angular acceleration
dydt(iyrad) = dydt_WEC(3:end); % radiation damping states for WEC model

dydt(iyLPPL) = dydt_pLineLP;
dydt(iyHPPL) = dydt_pLineHP;
%% %%%%%%%%%%%%   FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [dydt_control, control] = controller(t,y,par)
        %% PI control of P_Hout using w_pm
         % Error
        err_p_filt = y(iyp_filt) - par.control.p_hout_nom;
         % Feedforward
        w_ff = 0*par.control.w_pm_ctrl.min;
         % Control signal
        w_pm_nom = w_ff ...
            + (par.control.w_pm_ctrl.kp*err_p_filt ...
            + par.control.w_pm_ctrl.ki*y(iy_errInt_p_filt));
        control.w_pm_nom = (w_pm_nom <= par.control.w_pm_ctrl.max ...
                          & w_pm_nom >= par.control.w_pm_ctrl.min) ...
                          * w_pm_nom ...
                          + (w_pm_nom > par.control.w_pm_ctrl.max) ...
                          * par.control.w_pm_ctrl.max ...
                          + (w_pm_nom < par.control.w_pm_ctrl.min) ...
                          * par.control.w_pm_ctrl.min;
                      
        %% PI control of w_pm using T_gen
         % Error
        err_w_pm = y(iyw_pm) - control.w_pm_nom;
         % Feedforward
        T_ff = par.D_pm*(par.control.p_hout_nom - par.control.p_lin_nom);
         % Control signal
        control.Tgen = T_ff ...
            + (par.control.Tgen_ctrl.kp*err_w_pm ...
            + par.control.Tgen_ctrl.ki*y(iy_errInt_w_pm));
        
         % deriviatives for filtered signal and error integrals (w/
         % anti-wind-up)
        dydt_control = [    % filtered data
                        (y(iyp_hout) - y(iyp_filt))...
                        /par.control.tau_pfilt;

                            % error integral for pressure control
                        (y(iy_errInt_p_filt) < 1e12 ...
                        & y(iy_errInt_p_filt) > -1e12) ...
                        *err_p_filt;
                        
                            % error integral for speed control
                        (y(iy_errInt_w_pm) < 1e0 ...
                        & y(iy_errInt_w_pm) > -1e0)...
                        *err_w_pm];
    end

    function nonState = nonStateVars(t,y,par)
        % Accumulator capacitance
        nonState.C_lin = capAccum(y(iyp_lin),par.pc_lin,par.Vc_lin,par.f,par) ...
                        + lineCap(y(iyp_lin),1,par);
        nonState.C_lout = capAccum(y(iyp_lout),par.pc_lout,par.Vc_lout,par.f,par) ...
                        + lineCap(y(iyp_lout),1,par);
        nonState.C_hin = capAccum(y(iyp_hin),par.pc_hin,par.Vc_hin,par.f,par) ...
                        + lineCap(y(iyp_hin),2,par);
        nonState.C_hout = capAccum(y(iyp_hout),par.pc_hout,par.Vc_hout,par.f,par) ...
                        + lineCap(y(iyp_hout),2,par);

        % WEC-driven pump
        delta_p_wp = y(iyp_hin)-y(iyp_lout);
        WECpumpPumping = y(iytheta_dot)*delta_p_wp < 0;
        WECpumpMotoring = y(iytheta_dot)*delta_p_wp >= 0;
        nonState.q_w = (WECpumpPumping*par.eta_v_WEC ...
                    + WECpumpMotoring/par.eta_v_WEC)...
                    * par.D_WEC*abs(y(iytheta_dot));

        nonState.T_pto = -sign(y(iytheta_dot))...
                    * (WECpumpPumping/par.eta_m_WEC ...
                    + WECpumpMotoring*par.eta_m_WEC)...
                    * (abs(y(iytheta_dot)) > 5e-3)...
                    * par.D_WEC*delta_p_wp;

        % house power pump/motor
        delta_p_pm = y(iyp_lin) - y(iyp_hout);
        nonState.pmPumping = y(iyw_pm)*(delta_p_pm) >= 0;
        nonState.pmMotoring = y(iyw_pm)*(delta_p_pm) < 0;

        volLoss_pm = par.APP.C_s*(delta_p_pm/(par.mu*abs(y(iyw_pm)))) ...
                           + (delta_p_pm/par.beta)*(par.APP.V_r + 1);
        nonState.eta_v_pm = nonState.pmPumping*(1 - volLoss_pm) ...
                          + nonState.pmMotoring/(1 + volLoss_pm);
        nonState.q_pm = (nonState.pmPumping*nonState.eta_v_pm ...
                        + nonState.pmMotoring/nonState.eta_v_pm) ...
                        *par.D_pm*y(iyw_pm);

        mechLoss_pm = par.APP.C_v*par.mu*abs(y(iyw_pm))/delta_p_pm + par.APP.C_f;
        nonState.eta_m_pm = nonState.pmPumping/(1 + mechLoss_pm) ...
                          + nonState.pmMotoring*(1 - mechLoss_pm);
        nonState.Tpm = (nonState.pmPumping/nonState.eta_m_pm ...
                        + nonState.pmMotoring*nonState.eta_m_pm)...
                        *par.D_pm*delta_p_pm;

        % Reverse osmosis module
        nonState.q_perm = par.Sro*par.Aperm*(y(iyp_hout) - par.p_perm - par.p_osm);
        nonState.q_feed = 1/par.Y*nonState.q_perm;

        % ERU
        nonState.q_brine = nonState.q_feed - nonState.q_perm;
        nonState.q_ERUfeed = (par.eta_ERUv)^2*nonState.q_brine;

        % Charge Pump
        dP_SO = (y(iyp_lin) - par.p_o) - par.cn*par.w_c^2; % difference between shut-off pressure and current pressure differential
        nonState.q_c = (dP_SO < 0)* sqrt(dP_SO/par.cq);

        % Pressure relief valves
        nonState.q_linPRV = prv(y(iyp_lin),par.linPRV.p_crack,par.linPRV.C);
        nonState.q_loutPRV = prv(y(iyp_lout),par.loutPRV.p_crack,par.loutPRV.C);
        nonState.q_hinPRV = prv(y(iyp_hin),par.hinPRV.p_crack,par.hinPRV.C);
        nonState.q_houtPRV = prv(y(iyp_hout),par.houtPRV.p_crack,par.houtPRV.C);

    end

end
