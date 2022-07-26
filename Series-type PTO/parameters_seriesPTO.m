function par = parameters_seriesPTO(par)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters_seriesPTO.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 7/6/2021
%
% PURPOSE/DESCRIPTION:
% 
%
% FILE DEPENDENCY: NA
%
% UPDATES:
% 7/6/2021 - created.
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
    par.mu = 9.4e-4; % [Pa-s]  Kinematic (absolute) viscosity
    par.beta = 2.2e9; % [Pa]  Bulk Modulus of air free fluid
    par.p_vap = 0.037e5; % [Pa] vapour pressure of seawater
    par.R = 0.0001; % [-] fraction  Baseline fraction of air by volume entrained in the hydraulic fluid at atm
    par.p_o = 101.3e3; % [Pa]  Atmospheric pressure (reference)
    par.gamma = 1.4; % [-]ratio of specific heats for air
    
    % WEC parameters
    par = parameters_WECmodel(par);

    % WEC-pump
     % pumping chamber
    par.D_WEC = 0.442;         % [m^3/rad] flap pump displacement
    par.eta_v_WEC = 1;
    par.eta_m_WEC = 0.9;                % Flap pump mechanical efficiency
     % check valves
    par.A_check_h = 2.3e-3;           % [m^2]  Check valve area - high pressure (c3,c4)
    par.A_check_l = 2.8e-3;        % [m^2]  Check valve area - low pressure (c1,c2)
    par.P_crack = 101.3e3;          % [Pa]  Cracking pressure of check valve 20.7kPa = 3psi; 68.9kPa = 10psi

    % RO module parameters
    par.p_perm = par.p_o;
    par.p_osm = 2.7e6;
    par.Sro = 10000; % [m^3]
    par.Aperm = 2.57e-12; % [m^3/(N-s)] permeabiity coefficient (Yu and Jenne,2018)
%     par.Aperm = 0.0116e-6*.75;
    par.Y = 0.4;
    
    % power control unit
      % motor
    par.D_pm = 444e-6; % [m^3/rev]  Motor displacement
    par.w_pm_max = 1700/60*2*pi; % [rad/s] maximum speed of motor
    par.w_pm_min = 1/60*2*pi; % [rad/s] minimum speed of motor

       % efficiency model coeffs. (McCandlish and Dory model)
        % Axial piston motor (data from Danfoss APP 1.5)
    par.APM.C_s1 = 2.09e-9;
    par.APM.C_s2 = 2.57e-12;
    par.APM.V_r = 1.103;
    par.APM.C_v1 = 535760;
    par.APM.C_v2 = 8205.8;
    par.APM.C_f1 = 0.136;
    par.APM.C_f2 = -5.9e-4;
    
        % Axial piston pump (data from Danfoss APP 1.5)
    par.APP.C_s1 = 1.14e-9;
    par.APP.C_s2 = 1.42e-12;
    par.APP.V_r = 1.103;
    par.APP.C_v1 = -1798400;
    par.APP.C_v2 = 9672.7;
    par.APP.C_f1 = 0.1634;
    par.APP.C_f2 = -3.5e-4;
    
      % generator
       % rotor inertia estimated by NEMA MG-1 (14.46)
    n_poles = 3; % [qty.] number of poles of induction motor
    genPowerRating = 100; % [hp] power rating of induction motor
    par.Jpm = ( 0.02*2^n_poles*genPowerRating^(1.35-.05*n_poles/2) ) ...
                * 0.0421401101; % [lb ft^2 --> kg m^2] rotor inertia
    
    % Accumulators
     % LPA at inlet to LP pipelinerge pressure
    par.Vc_lout = 1000e-3; % [m^3] gas volume at charge pressure
    par.pc_lout = 0.2e6; % [Pa] charge pressure
     % LPA at outlet from LP pipeline
    par.Vc_lin = 1000e-3; % [m^3] gas volume at charge pressure
    par.pc_lin = 0.2e6; % [Pa] charge pressure
     % HPA at inlet to HP pipeline
    par.Vc_hin = 3000e-3; % [m^3] gas volume at charge pressure
    par.pc_hin = 4e6; % [Pa] charge pressure
     % HPA at outlet from HP pipeline and inlet to PCU
    par.Vc_hout = 2000e-3; % [m^3] gas volume at charge pressure
    par.pc_hout = 4e6; % [Pa] charge pressure
     % HPA at inlet to RO module and outlet from PCU
    par.Vc_RO = 3000e-3; % [m^3] gas volume at charge pressure
    par.pc_RO = 4e6; % [Pa] charge pressure

     
    % Pipelines    
     % LP pipeline
    par.L_line(1) = 500; % [m] length of LP pipeline
    par.d_line(1) = 0.2; % [m] diameter of LP pipeline
    par.A_line(1) = pi/4*par.d_line(1); % crosssectional flow area
    par.n_seg(1) = 2; % minimum of 2
    par.I(1) = par.rho*(par.L_line(1)/par.n_seg(1))/par.A_line(1);
     % HP pipeline
    par.L_line(2) = 500; % [m] length of HP pipeline
    par.d_line(2) = 0.1; % [m] diameter of HP pipeline
    par.A_line(2) = pi/4*par.d_line(2); % crosssectional flow area
    par.n_seg(2) = 2; % minimum of 2
    par.I(2) = par.rho*(par.L_line(1)/par.n_seg(2))/par.A_line(1);
    

    % Contoller Parameters
    par.control.p_hout_nom = 7.64e6; % [Pa]
    par.control.p_RO_nom = 4e6; % [Pa] (not actually used for as control input)
    par.control.p_RO_max = 8.2e6; % [Pa]
    par.control.p_RO_min = max(3e6,par.pc_RO); % [Pa]
    pumpSpeed = @(Pfeed) par.Sro*par.Aperm/par.Y...
                            *(Pfeed - par.p_osm - par.p_o)...
                            /(2*pi*par.D_pm); % [rad/s]
    w_pm_max_RO = pumpSpeed(par.control.p_RO_max);
    w_pm_min_RO = pumpSpeed(par.control.p_RO_min);
    par.control.w_pm_ctrl.nom = pumpSpeed(par.control.p_RO_nom);
    par.control.w_pm_ctrl.max = min(w_pm_max_RO,par.w_pm_max);
    par.control.w_pm_ctrl.min = max(w_pm_min_RO,par.w_pm_min);
    par.control.w_pm_ctrl.kp = 1e-4;
    par.control.w_pm_ctrl.ki = 6e-6;
    par.control.tau_pfilt = 0.1; % [s] time constant for LPF for pressure signal
    
     % PI control of w_pm using T_gen
    par.control.Tgen_ctrl.kp = 1e3;
    par.control.Tgen_ctrl.ki = 2e5;
    
 % Charging system (Intake & Boost pump)
    par.cn = 7;
    par.cq = -1e6;
    par.w_c = (3600)*2*pi/60; % [rpm -> rad/s]

 % Charging system (Intake & Boost pump)
    par.eta_boostPump = 0.7;  % pumping efficiency of pressure boost pump
    par.A_check_boost = 30e-3;
    par.omega_boost = 2600*2*pi/60; % [rpm -> rad/s] Boost pump shaft speed
    par.P_tank = .65e6; 
%% 
par.L_intake1 = 2; %320; %length for inlet to boost pump
par.L_intake2 = 660; %340; %length from boost pump to plant
par.d_intake = 0.11;  % Intake line diameter

 % discharge pipeline
par.L_discharge = 20; % [m] length for discharge pipeline
par.d_discharge = par.d_intake;          % [m] Discharge line ID
     
    %% %%%%%%%%%%%%   FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     function par = parameters_WECmodel(par)
%         par = parameters_WECmodel_V02x00(par);
%     end

end