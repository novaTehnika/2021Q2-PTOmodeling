function [dydt, nonState] = sys(t,y,par)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sys_coulombPTO.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 6/23/2021
%
% PURPOSE/DESCRIPTION:
% Calculate the state derivatives for a simple wave energy PTO with coulomb
% damping. For use as an intermediate step in modeling PTO archetectures in
% the 2021Q2 PTO modeling project. The model includes the hydrodynamic WEC
% model found in the 2021Q2 Hydrodynamic modeling project (but really
% developed in 2018) and assumes a constant PTO torque magnitide.
%
% FILE DEPENDENCY: 
% pipeline_VXXxXX.m
% nonStateVarsPTOwPL_VXXxXX
%
% UPDATES:
% 6/23/2021 - V01x00 - created from sysPTOwPL_V01x00.m to be
%           used in simPTO_V02x00.m.
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
iytheta = 1;
iytheta_dot = 2;
iyrad = (1:par.WEC.ny_rad) + iytheta_dot; % state vector indices for radiation damping states for WEC model
iWEC = [iytheta iytheta_dot iyrad];

nonState = nonStateVars(t,y,par); % recover nonstate variables

[dydt_WEC, nonState.torqueWEC, nonState.waveElev] = ...
                    flapModel(t,y(iWEC),nonState.T_pto,par);

dydt = zeros(2+par.WEC.ny_rad,1);
dydt(iytheta) = dydt_WEC(1); % angular velocity
dydt(iytheta_dot) = dydt_WEC(2); % angular acceleration
dydt(iyrad) = dydt_WEC(3:end); % radiation damping states for WEC model

%% %%%%%%%%%%%%   FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function nonState = nonStateVars(t,y,par)
        nonState.T_pto = (abs(y(iytheta_dot)) > 1e-2)* ...
                        -sign(y(iytheta_dot))*par.Tcoulomb;
    end
end
