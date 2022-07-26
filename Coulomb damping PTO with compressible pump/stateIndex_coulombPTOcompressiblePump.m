% stateIndex_coulombPTOcompressiblePump.m script m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 07/08/2022
%
% PURPOSE/DESCRIPTION:
% This script loads the state indices for the coulombPTOcompressiblePump
% model
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

iyp_wpA = 1;
iyp_wpB = 2;
iytheta = 3;
iytheta_dot = 4;
iyrad = (1:par.WEC.ny_rad) + iytheta_dot; % state vector indices for radiation damping states in WEC model
iyWEC = [iytheta iytheta_dot iyrad];