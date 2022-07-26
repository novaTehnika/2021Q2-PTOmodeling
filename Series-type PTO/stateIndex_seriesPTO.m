% stateIndex_seriesPTO.m script m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 07/08/2022
%
% PURPOSE/DESCRIPTION:
% This script loads the state indices for the seriesPTO model
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

iyp_lin = 1;
iyp_lout = 2;
iyp_hin = 3;
iyp_hout = 4;
iyp_RO = 5;

iyw_pm = 6;

iyp_filt = 7; 
iy_errInt_p_filt = 8; 
iy_errInt_w_pm = 9;
iycontrol = [iyp_filt; iy_errInt_p_filt; iy_errInt_w_pm];

iytheta = 10;
iytheta_dot = 11;
iyrad = (1:par.WEC.ny_rad) + iytheta_dot; % state vector indices for radiation damping states for WEC model

iyLPPL = (1:2*par.n_seg(1)-1) + iyrad(end);
iyqLP = (1:2:2*par.n_seg(1)-1) + (iyLPPL(1)-1);
iypLP = (2:2:2*par.n_seg(1)-1) + (iyLPPL(1)-1);

iyHPPL = (1:2*par.n_seg(2)-1) + iyLPPL(end);
iyqHP = (1:2:2*par.n_seg(2)-1) + (iyHPPL(1)-1);
iypHP = (2:2:2*par.n_seg(2)-1) + (iyHPPL(1)-1);