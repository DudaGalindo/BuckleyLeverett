%% Problem description
% Analytical solution of Buckley–Leverett equation using Welge method
% Brooks–Corey type relative permeability curves are used
function [index, dfw, Sw, Xsw, Swf, nts] = Welge_Solution(L, A, phi, k, muo, muw, Sor, Swc, no, nw)
%% 1 Parameters initialization
% 1.1 rock properties
theta = pi*0.0; % angle: x-direction vs horizontal[rad]
% 1.2 fluid properties
rho_o = 1.e3; % oil density [kg/m^3]
rho_w = 1.0e3; % water density [kg/m^3]
dlt_rho = rho_w-rho_o; % density difference [kg/m^3]

% 1.3 relative permeability parameters
kro_max = 1.; % maximum permeability of oil phase
krw_max = 1.; % maximum permeability of water phase

% 1.4 other initial parameters
dlt_Sw = 1e-3; % constant saturation step
Sw = ((Swc+eps):dlt_Sw:(1-Sor))'; % water saturation vector
qt = 1.0; % constant water injection rate [m^3/s]
g = 9.8067; % gravity acceleration constant [m/s^2]
time0 = 0; % initial calculation time [s]
timef = 0.2; % final calculation time [s]
nts = 4; % time steps for calculating
nsw = size(Sw,1); % number of water saturation vector
format = '%4.2e'; % precision format for plotting legend
%% 2 Relative permeabilities, fractional flow function and its derivative
% 2.1 initialization of function handles
kr_w = @(S) krw_max.*(( S-Swc)/(1-Swc-Sor)).^nw; % relative permeability of water
kr_o = @(S) kro_max.*((1-S-Sor)/(1-Swc-Sor)).^no; % relative permeability of oil
mob = @(S) (kr_o(S)./muo).*(muw./kr_w(S)); % mobility function
f_w = @(S) 1./(1+mob(S))-A*k.*kr_o(S)./... % fractional flow function fw
    (qt*muo)*dlt_rho*g*sin(theta)./(1+mob(S));
df_w = @(S) ((nw.*mob(S)./(S-Swc))+(no.*mob(S))... % derivate of fw
./(1-S-Sor))./(1+mob(S)).^2+A*k.*kr_o(S)*no*dlt_rho*g*sin(theta)./...
((1-S-Sor)*qt*muo.*(1+mob(S)))+A*k.*kr_o(S)*dlt_rho*g*sin(theta)./...
(qt*muo.*(1+mob(S)).^2).*(((nw.*mob(S)./(S-Swc))+(no.*mob(S))./...
(1-S-Sor))./(1+mob(S)).^2);
% 2.2 relative permeability and fractional flow vectors
krw = kr_w(Sw); % relative permeability vector of water
kro = kr_o(Sw); % relative permeability vector of oil
fw = f_w (Sw); % fractional flow vector
dfw = df_w(Sw); % fractional flow derivate vector
%% 3 Calculate advance frontal water saturation
% 3.1 Evaluate advance frontal water saturation Swf
[index, Swf, dfwf, bf] = calSatFront(Sw, fw, dfw);
% 3.2 Calculate the time when Swf reaches at production well
krw_swf = kr_w(Swf);
kro_swf = kr_o(Swf);
fw_swf = f_w (Swf);
dfw_swf = df_w(Swf);
t_pw = (A*phi*L)/(qt*dfw_swf);

%% 4 Calculate water saturation profile
t = linspace(time0, timef, nts);
dfwt = dfw(end:-1:index); % inverted sequence of dfw
Swt = Sw(end:-1:index); % inverted sequence of Sw
for ti = 1:nts
    for i = 1:(nsw-index+1)
        Xsw(i,ti) = qt*t(ti)/(A*phi)*dfwt(i);
    end
end
%---------------------------------------------------
%% 5 Plot results
% 5.1 Plot the relative permeability curves
%{
h_fig1 = figure(4);
set(h_fig1, 'color', 'w', 'NumberTitle', 'off', 'Name', 'Relative permeability curves: kr');

plot(Sw, krw, '-b', Sw, kro, '-r', 'LineWidth', 2.0);
axis([0.0 1.0 0.0 1.0]);
axis square;
set(gca, 'Fontname', 'Times New Roman','FontSize',10);
title('Relative permeability curves')
xlabel('{\it S}_w');
ylabel('{\it k}_{r {\it \beta}}');
set(gca, 'YTick', 0:0.2:1);
h_legend1 = legend('water phase', 'oil phase');
set(h_legend1,'Box', 'on', 'Location', 'best');

% 5.2 Plot the fractional flow function curve
h_fig2 = figure(2);
set(h_fig2, 'color', 'w', 'NumberTitle', 'off', 'Name', 'Fractional Flow Curve fw');
plot(Sw, fw, 'b', 'LineWidth', 2.0);
hold on;
plot(Sw, (dfwf.*Sw+bf), '-r', 'LineWidth', 2.0);
plot(Swc, 0, 'ro', 'Markersize', 8);
plot(Swf, fw(index), 'ro', 'Markersize', 8);
plot([Swf Swf], [0 fw(index)], '--r', 'LineWidth', 1.5);
plot([0 Swf], [fw(index) fw(index)], '--r', 'LineWidth', 1.5);
axis([0.0 1.0 0.0 1.0]);
axis square;
set(gca, 'Fontname', 'Times New Roman','FontSize',10);
title('fractional flow function curve')
xlabel('{\it S}_w');
ylabel('{\it f}_w');
set(gca, 'YTick', 0:0.2:1);

% 5.3 Plot the derivate curve of fractional flow function
h_fig3 = figure(3);
set(h_fig3, 'color', 'w', 'NumberTitle', 'off', 'Name', 'Derivatives of Fractional Flow Curve dfw/dSw');
plot(Sw, dfw, '-b', 'LineWidth', 2.0);
set(gca, 'XLim', [0 1]);
set(gca, 'Fontname', 'Times New Roman','FontSize',10);
title('(c) derivatives of fractional flow function')
xlabel('{\it S}_w');
ylabel('d {\it f}_w / d {\it S}_w');
axis square;

% 5.4 Plot the water saturation profiles
h_fig4 = figure(4);
set(h_fig4, 'color', 'w', 'NumberTitle', 'off', 'Name', 'Water Saturation Profiles');

SwTime= [];
for ti = 1:nts
    if ti == nts
        plot(Xsw(:,ti), Swt, '-r', 'LineWidth', 2.0);
        hold on;
        plot([0 Xsw(1,ti)], [1-Sor 1-Sor], '-r', 'LineWidth', 2.0);
        SwTime = [ SwTime; ['Time = ' num2str(t(ti)/86400, format) ' day'] ];
        plot([Xsw(end,ti) Xsw(end,ti)], [Swf Swc], '-r', 'LineWidth', 2.0);
        plot([Xsw(end,ti) L], [Swc Swc], '-r', 'LineWidth', 2.0);
    end
end
set(gca,'YLim', [0 1], 'YTick', 0:0.2:1, 'XLim', [0 L]);
set(gca, 'Fontname', 'Times New Roman','FontSize',10);
title('(d) water saturation profiles')
xlabel('{\it x} (m)');
ylabel('{\it S}_w');
h_legend2 = legend(SwTime);
set(h_legend2,'Box', 'on', 'Location', 'best');
axis square;
%}
end
