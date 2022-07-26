function [dydt, nonState] = sys_coulombPTOcompressiblePump(t,y,par)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sys_coulombPTOcompressiblePump.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 12/14/2021
%
% PURPOSE/DESCRIPTION:
% Calculate the state derivatives for a simple wave energy PTO with coulomb
% damping. The WEC-driven pump is modeled with compressible fluid filled 
% pumping chambers and a check valve rectifier that is in comunication with
% constant pressure sources. The model includes the hydrodynamic WEC
% model.
%
% FILE DEPENDENCY: 
%
% UPDATES:
% 12/14/2021 - created from sys_coulombPTO.m and 
%              sys_seriesPTO.m.
% 07/08/2022 - replaced state indice specification in file with seperate
% script, stateIndex_coulombPTOcompressiblePump.m.
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
stateIndex_coulombPTOcompressiblePump

% recover nonstate variables
nonState = nonStateVars(t,y,par); 

% Calculate the hydrodynamics WEC state derivatives and output nonstate 
% variables like forces and wave elevation
[dydt_WEC, nonState.torqueWEC, nonState.waveElev] = ...
                    flapModel(t,y(iyWEC),nonState.T_pto,par);

% State derivatives
dydt = zeros(iyWEC(end),1);
dydt(iyp_wpA) = 1/nonState.CwpA*(nonState.q_wpA ...
                                + par.D_wp*y(iytheta_dot));
dydt(iyp_wpB) = 1/nonState.CwpB*(nonState.q_wpB ...
                                - par.D_wp*y(iytheta_dot));
dydt(iytheta) = dydt_WEC(1); % angular velocity
dydt(iytheta_dot) = dydt_WEC(2); % angular acceleration
dydt(iyrad) = dydt_WEC(3:end); % radiation damping states for WEC model

%% %%%%%%%%%%%%   FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     function [dydt, torqueFlap, waveElev] = flapModel(t,y, ptoTorque, par)
%         [dydt, torqueFlap, waveElev] = ...
%             flapModel_V01x00(t,y, ptoTorque, par);
%     end

    function nonState = nonStateVars(t,y,par)
        % compressibility of pumping chambers
        nonState.CwpA = capPump(y(iyp_wpA), ...
                        par.D_wp*(par.thetaMax - y(iytheta)) + par.V_wpdead);
        nonState.CwpB = capPump(y(iyp_wpB), ...
                        par.D_wp*(par.thetaMax + y(iytheta)) + par.V_wpdead);
        
        % Flow through the check valves
         % pressure differential (used multiple time)
        dp_TA = par.p_lout - y(iyp_wpA);
        dp_AP = y(iyp_wpA) - par.p_hin;
        dp_TB = par.p_lout - y(iyp_wpB);
        dp_BP = y(iyp_wpB) - par.p_hin;
        
         % valve coefficient, a function of the pressure difference
        kv_TA = kvcv(dp_TA,par.kv_cvMin,par.kv_cvTx,par.p_crackTx,par.dp_strokeTx);
        kv_AP = kvcv(dp_AP,par.kv_cvMin,par.kv_cvxP,par.p_crackxP,par.dp_strokexP);
        kv_TB = kvcv(dp_TB,par.kv_cvMin,par.kv_cvTx,par.p_crackTx,par.dp_strokeTx);
        kv_BP = kvcv(dp_BP,par.kv_cvMin,par.kv_cvxP,par.p_crackxP,par.dp_strokexP);
        
         % Flow rate calculation
        nonState.q_wpTA = orifice(dp_TA,kv_TA);
        nonState.q_wpAP = orifice(dp_AP,kv_AP);
        nonState.q_wpTB = orifice(dp_TB,kv_TB);
        nonState.q_wpBP = orifice(dp_BP,kv_BP);
        
        % net flow into of the pumping chambers
        nonState.q_wpA = nonState.q_wpTA - nonState.q_wpAP;
        nonState.q_wpB = nonState.q_wpTB - nonState.q_wpBP;
        
        % net flow out of the WEC-driven pump
        nonState.q_wpP = nonState.q_wpAP + nonState.q_wpBP;
        
        % net flow into the WEC-driven pump
        nonState.q_wpT = nonState.q_wpTA + nonState.q_wpTB;
        
        % Reaction torque of WEC-driven pump
        nonState.T_pto = 1/par.eta_wp*par.D_wp*(y(iyp_wpB) - y(iyp_wpA));
        
        function C = capPump(p,V)
            % Description: Calculate teh capacitance based on teh
            % effeective bulk modulus (the fluid having entrained gas) and
            % the volume of the chamber.
            
            % calculate effective bulk modulus 
             % via Cho method (as arranged in Yudell, 2017)
            beta_eff = par.beta* ...
            (((p/par.p_o)^(1/par.gamma)*exp((par.p_o-p)/par.beta)+par.R)...
            /(par.R/par.gamma*par.beta/p+(p/par.p_o)^(1/par.gamma)...
            *exp((par.p_o-p)/par.beta))); 

             % isothermal bulk modulus
            %  beta_eff = par.beta/(1 + par.beta*(par.R*par.p_o/p^2));

            % calculate capacitance
            C = V/beta_eff;
        end
        
        function kv = kvcv(delta_p,kv_min,kv_max,p_crack,delta_p_stroke)
            % Description: calculate the flow coefficient of a valve
            % assuming that it varies linearly with the pressure 
            % differential, is at a minimum below the cracking pressure and
            % some margin from the cracking pressure where the valve is
            % fully stoked
            kv = min(kv_max,max(kv_min,...
                kv_max*(delta_p - p_crack)/(delta_p_stroke)));
        end
        
        function q = orifice(delta_p,kv)
            % Description: Calculate the flow rate through an orifice
            % assuming a knife edge orifice.
            q = kv*sqrt(abs(delta_p))*sign(delta_p);
        end
        
    end
end
