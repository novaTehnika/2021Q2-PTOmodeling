% initialConditionDefault_seriesPTO.m script m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 07/08/2022
%
% PURPOSE/DESCRIPTION:
% This script loads default intial conditions for the seriesPTO model.
%
% FILE DEPENDENCY:
%
% UPDATES:
% 07/08/2022 - Created.
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

y0(1,iyp_lin) = 1e6; % [Pa]
y0(1,iyp_lout) = 1e6; % [Pa]
y0(1,iyp_hin) = par.control.p_hout_nom; % [Pa]
y0(iyp_hout) = par.control.p_hout_nom; % [Pa]
y0(iyp_RO) = par.control.p_RO_nom; % [Pa]

y0(iyw_pm) = par.control.w_pm_ctrl.nom;

y0(iyp_filt) = par.control.p_hout_nom;
y0(iy_errInt_p_filt) = 0;
y0(iy_errInt_w_pm) = 0;

y0(iytheta) = 0.1*pi/2;
y0(iytheta_dot) = .1;
y0(iyrad) = zeros(1,length(iyrad));

y0(iyLPPL) = y0(iyp_lout)*mod(2:2*par.n_seg(1),2)';
y0(iyHPPL) = y0(iyp_hin)*mod(2:2*par.n_seg(2),2)';