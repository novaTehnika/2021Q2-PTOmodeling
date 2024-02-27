addpath('WEC model') 
addpath(['WEC model' filesep 'WECdata'])

nw = floor(logspace(log10(1),log10(1e4),100));
for inw = 1:numel(nw)
    [~, A_constructed(inw), ~, ~,A_ref(inw), ~] = test(nw(inw));
end

%% Plotting
black = [0 0 0];
maroon = [122 0 25]/256;
gold = [255 204 51]/256;
blue = [0 75 135]/256;
orange = [226 100 55]/256;
green = [63 150 87]/256;
pink = [255 192 203]/256;
blue1 = [0 150 255]/256;
purple = [128 0 128]/256;
green1 = [0 255 150]/256;
color = [maroon; gold; blue; orange; green; pink; blue1; purple; green1];

linestyles = {'-', '--', '-.', ':','-', '--', '-.', ':',};

supTitleFontSize = 9;
subTitleFontSize = 9;
axFontSize = 8;
bottomEdge = 1;
leftEdge = 3;
width = 4;
height = 3;
lineWidth = 0.5;

%%
width = 4;
height = 3;

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

semilogx(nw,A_constructed,'--','color',maroon,'LineWidth',lineWidth)
hold on
semilogx(nw,A_ref,'k-','LineWidth',lineWidth)

yLim = ylim;
ylim([0 yLim(2)]);

grid on

title(['Wave Power Spectral Density Function Convergence:',newline, ...
        'Number of Wave Frequency Components'], ...
        'Interpreter','latex','FontSize',supTitleFontSize,'fontname','Times')
xlabel('no. freq. components', ...
        'Interpreter','latex','FontSize',supTitleFontSize,'fontname','Times')
ylabel('integral of PSD', ...
        'Interpreter','latex','FontSize',supTitleFontSize,'fontname','Times')

ax = gca;
ax.FontName = 'times';
ax.FontSize = axFontSize;

leg = legend('ideal','equal energy method');
leg.FontSize = axFontSize;
leg.FontName = 'Times';
set(leg, 'Location', 'best')

%%

width = 4;
height = 3;

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

loglog(nw,(A_ref-A_constructed)./A_ref,'k-','LineWidth',lineWidth)

grid on

title(['Wave Power Spectral Density Function Convergence:',newline, ...
        'Number of Wave Frequency Components'], ...
        'Interpreter','latex','FontSize',supTitleFontSize,'fontname','Times')
xlabel('no. freq. components', ...
        'Interpreter','latex','FontSize',supTitleFontSize,'fontname','Times')
ylabel('error of integral', ...
        'Interpreter','latex','FontSize',supTitleFontSize,'fontname','Times')

ax = gca;
ax.FontName = 'times';
ax.FontSize = axFontSize;



%%
[w_constructed, ~, S_w_constructed, w_ref, ~, S_w_ref] = test(20);

width = 4;
height = 3;

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

semilogx(w_ref,S_w_ref,'k-')
hold on
scatter(w_constructed,S_w_constructed,'kx')

grid on

title(['Wave Power Spectral Density Function Construction:',newline, ...
        'Equal Energies Method'], ...
        'Interpreter','latex','FontSize',supTitleFontSize,'fontname','Times')
xlabel('frequency (rad/s)', ...
        'Interpreter','latex','FontSize',supTitleFontSize,'fontname','Times')
ylabel('power spectral density', ...
        'Interpreter','latex','FontSize',supTitleFontSize,'fontname','Times')

ax = gca;
ax.FontName = 'times';
ax.FontSize = axFontSize;

leg = legend('function','discretization');
leg.FontSize = axFontSize;
leg.FontName = 'Times';
set(leg, 'Location', 'best')

%% %%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%
function [w_constructed, A_constructed, S_w_constructed, w, A_ref, S_w] = test(nw)
    load('nemohResults_vantHoff2009_20180802.mat','OSWC_vantHoff2009')
    par.imprt.WEC.w = OSWC_vantHoff2009.w;          % [rad/s] frequency 
    par.wave.Hs = 2.75;
    par.wave.Tp = 12;
    par.WEC.nw = nw;

    % determine new freqencies based on equal energy method
     % determine the area under the curve of the wave spectrum for
     % discretization
      % high density grid of frequencies evenly spaced
    w = linspace(par.imprt.WEC.w(1),par.imprt.WEC.w(end),1e6); 
      % wave spectrum calculated for high-density even grid
    S_w = PiersonSpec(w,par);
      % area for each bin
    a = trapz(w,S_w)/(par.WEC.nw+1);
    
      % calculate the cumulative area under the curve
    A = cumtrapz(w,S_w);
    
     % determine the frequencies for an equal area grid by finding upper and
     % lower bounds for the bins
    par.WEC.w = zeros(par.WEC.nw,1);
    par.WEC.dw = zeros(par.WEC.nw,1);
    w_lb = zeros(par.WEC.nw,1); % lower bound for each bin
    w_ub = zeros(par.WEC.nw,1); % upper bound for each bin
      % min freq. set to a significant value
    wmin = w(find(A >= 0.01*a,1,'first')); 
      % max freq. set to a significant value
    wmax = w(find(A <= (A(end) - 0.01*a),1,'last')); % max freq. set to a significant value
    for iw = 1:par.WEC.nw
        % set lower bound
        if iw == 1
            w_lb(iw) = wmin;
        else
            w_lb(iw) = w_ub(iw-1);
        end
        % find upper bound to satisfy area increment
        if iw == par.WEC.nw
            w_ub(iw) = wmax;
        else
           w_ub(iw) = w(find(A >= a*iw,1,'first'));
        end
	    par.WEC.dw(iw) = w_ub(iw) - w_lb(iw);
        par.WEC.w(iw) =  (w_ub(iw) + w_lb(iw))/2;
    end
    
    % calculate the wave spectrum for the equal area frequencies
    par.wave.S_w = PiersonSpec(par.WEC.w,par);

    % outputs
    w_constructed = par.WEC.w;
    S_w_constructed = par.wave.S_w;
    A_constructed = sum(par.wave.S_w.*par.WEC.dw);
    A_ref = A(end);
end

function S_w = PiersonSpec(w,par)
    % Based on Falnes (2002) "Ocean Waves and Oscillating Systems:..."
    S_w = 10*pi^5*par.wave.Hs^2/par.wave.Tp^4./w.^5 ... 
        .*exp(-20*pi^4/par.wave.Tp^4./w.^4)/(2*pi);
end