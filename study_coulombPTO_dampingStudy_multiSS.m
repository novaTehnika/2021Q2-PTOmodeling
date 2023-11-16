% study_coulombPTO_dampingStudy_multiSS.m script m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 08/11/2021
%
% PURPOSE/DESCRIPTION:
% This script serves as a shell for performing parameter variation studies
% using the model contained in sys_coulombPTO.m and run by solved by 
% sim_coulombPTO.m. Specifically this script varies the Coulomb damping
% torque for a given set of sea states to produce a mean power curves for 
% each. Parameter and function initiallization fuctions are called within 
% this script before the sim_coulombPTO.m script is called. 
%
% FILE DEPENDENCY:
% sys_coulombPTO.m
% sim_coulombPTO.m
% parameters_coulombPTO.m
%
% UPDATES:
% 08/11/2021 - (V01x00) - Created from 
% studyPTO_coulombPTO_dampingStudy_singleSS.m.
% 07/08/2022 - added initial conditions as arguement to sim_coulombPTO()
% and simulation start time as parameter.
% 08/02/2022 - added ramp period to par struct and included tstart so that
% ramp ends at t=0 and modified average power and WEC velocity calculations
% to exclude ramp period.
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
clear
% clc
addpath('WEC model') 
addpath(['WEC model' filesep 'WECdata']) 
addpath('Coulomb damping PTO') 
addpath('Sea States')
addpath('Solvers')
%% %%%%%%%%%%%%   SIMULATION PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation timeframe
par.Tramp = 250; % [s] excitation force ramp period
par.tstart = 0; %[s] start time of simulation
par.tend = 2000; %[s] end time of simulation

% Solver parameters
par.odeSolverRelTol = 1e-9; % Rel. error tolerance parameter for ODE solver
par.odeSolverAbsTol = 1e-9; % Abs. error tolerance parameter for ODE solver
par.MaxStep = 1e-2;

% Wave construction parameters
par.WEC.nw = 1000; % num. of frequency components for harmonic superposition 
par.wave.rngSeedPhase = 3; % seed for the random number generator

%% %%%%%%%%%%%%   Study Variables  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch 1
    case 1
        Hs = [2.34 2.64 5.36 2.05 5.84 3.25];
        Tp = [7.31 9.86 11.52 12.71 15.23 16.50];
    case 2
        load('SSdata_HumboltBay_1D.mat')
end
nSS = length(Tp);

nVar1 = 121;
% Tcoulomb = 1e6*linspace(1,10,nVar1);% [Nm] PTO reaction torque
Tcoulomb = 1e6*logspace(log10(0.1),log10(15),nVar1);% [Nm] PTO reaction torque

saveSimData = 0; % save simulation data (1) or just output variables (0)

%% %%%%%%%%%%%%   COLLECT DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iSS = 1:nSS
    par.wave.Hs = Hs(iSS);
    par.wave.Tp = Tp(iSS);
    
    % load base parameters
    par = parameters_coulombPTO(par,...
    'nemohResults_vantHoff2009_20180802.mat','vantHoffTFCoeff.mat');
    
    % Define initial conditions
    y0 = [  0, ...
            0, ...
            zeros(1,par.WEC.ny_rad)];

    parfor iVar1 = 1:nVar1
        param = par; % store individual parameter strucs for parallel compute

        % change design parameter
        param.Tcoulomb = Tcoulomb(iVar1); 

        % run simulation
        tic
        out = sim_coulombPTO(y0,param);
        toc

        % Calculate metrics
        t_vec = find(out.t>=par.tstart);
        if max(out.theta) > pi ||  min(out.theta) < -pi
            PP(iSS,iVar1) = 0;
            theta_dot_ave(iSS,iVar1) = 0;
        else 
            PP(iSS,iVar1) = mean(-out.T_pto(t_vec).*out.theta_dot(t_vec));
            theta_dot_ave(iSS,iVar1) = mean(abs(out.theta_dot(t_vec)));
        end

        if saveSimData
            simOut(iSS,iVar1) = out;
        end

    end
end


%% %%%%%%%%%%%%   SAVE DATA   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = ['data_coulombPTO_dampingStudy_',date];
files = ls;
nfiles = size(files,1);
k = 1;
for j = 1:nfiles
    if strfind(files(j,:),filename)
        k = k+1;
    end
end
save([filename,'_',num2str(k)])

return

%% %%%%%%%%%%%%   PLOTTING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tarHs = 5.25;
% tarTe = 9.5;
% temp = find(Hs == tarHs);
% ISS = temp(Te(temp)== tarTe)
% Hs(iSS)
% Te(iSS)
%%
iSS = 7
figure
xlabel('Torque (MNm)')
title(['WEC Power Absorption, Coulomb Damping',newline,...
            'Sea State ',num2str(iSS)])
yyaxis left
hold on
plot(1e-6*Tcoulomb,1e-3*PP(iSS,:))
ylabel('Power (kW)')
% ylim(10*[-pi/2 pi/2])

yyaxis right
hold on
plot(1e-6*Tcoulomb,theta_dot_ave(iSS,:))
ylabel('Speed, mean (rad/s)')
% ylim()
% xlim([0 2])

%%
figure
hold on
title(['WEC Power Absorption, Coulomb Damping'])
% for iSS = 1:114
%     if weight(iSS) >= 0.5 && weight(iSS) < 1
%         p = plot(1e-6*Tcoulomb,1e-3*PP(iSS,:),'color',0.9*[1 1 1],'linewidth',1)
%         p.Annotation.LegendInformation.IconDisplayStyle = 'off';
%     end
% end
for iSS = 1:114
    if weight(iSS) >= 1 && weight(iSS) < 2
        p = plot(1e-6*Tcoulomb,1e-3*PP(iSS,:),'--','color',0.5*[1 1 1],'linewidth',1);
        p.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end
for iSS = 1:114
    if weight(iSS) >= 2
        p = plot(1e-6*Tcoulomb,1e-3*PP(iSS,:),'color',0*[1 1 1],'linewidth',1);
        p.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end
yLim = ylim;

plot(1e-6*Tcoulomb,4e-11*Tcoulomb.^2,'r','LineWidth',2)

plot([-99 -100],[-99 -100],'--','color',0.5*[1 1 1],'linewidth',1);
plot([-99 -100],[-99 -100],'color',0*[1 1 1],'linewidth',1);

legend('quadratic function','sea state probability > 1%','sea state probability > 2%')

ylim(yLim)
xlim([0 20])
xlabel('Torque (MNm)')
ylabel('Power (kW)')
