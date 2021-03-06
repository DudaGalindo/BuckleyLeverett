%% Input data file 
L = 1.0; % domain length [m]
A = 1.0; % area of cross-section [m^2]
phi = 1;
theta = pi*0.0; % angle: x-direction vs horizontal[rad]
% 1.2 fluid properties
mi_o = 1.0e-3; % oil phase viscosity [Pa*s]
mi_w = 1.0e-3; % water phase viscosity [Pa*s]

% 1.3 relative permeability parameters
Sor = 0.0; % residual oil saturation
Swc = 0.0; % connate water saturation
n_o = 2.00; % exponent of oil phase
n_w = 2.00; % exponent of water phase

% 1.4 other initial parameters

t0 = 0.; % initial calculation time [s]
t_final = .5; % final calculation time [s]

%% Run Test - mesh with 8 elements

solver = 'FOUM';
M = Mesh;
M.n_el = 8;
M.L = L;
Sw = ones(1,M.n_el).*Swc;
Sw(1) = 1-Sor;
Sw(end) = Swc;
x8 = linspace(0,L,8);
t = t0;

[Sw8] = main_solver(M, phi, mi_o, mi_w, Sw, Swc, Sor, n_o, n_w, t, t_final, solver);

%% Run Test - mesh with 16 elements
 
solver = 'FOUM';
M = Mesh;
M.n_el = 16;
M.L = L;
Sw = ones(1,M.n_el).*Swc;
Sw(1) = 1-Sor;
Sw(end) = Swc;
x16 = linspace(0,L,16);
t = t0;
[Sw16] = main_solver(M, phi, mi_o, mi_w, Sw, Swc, Sor, n_o, n_w, t, t_final, solver);

 
%% Run Test - mesh with 32 elements
 
solver = 'FOUM';
M = Mesh;
M.n_el = 32;
M.L = L;
Sw = ones(1,M.n_el).*Swc;
Sw(1) = 1-Sor;
Sw(end) = Swc;
x32 = linspace(0,L,32);
t = t0;
[Sw32] = main_solver(M, phi, mi_o, mi_w, Sw, Swc, Sor, n_o, n_w, t, t_final, solver);
 
%% Run Test - mesh with 64 elements
 
solver = 'FOUM';
M = Mesh;
M.n_el =64;
M.L = L;
Sw = ones(1,M.n_el).*Swc;
Sw(1) = 1-Sor;
Sw(end) = Swc;
x64 = linspace(0,L,64);
t = t0;
[Sw64] = main_solver(M, phi, mi_o, mi_w, Sw, Swc, Sor, n_o, n_w, t, t_final, solver);

%% Run Test - mesh with 128 elements
 
solver = 'FOUM';
M = Mesh;
M.n_el =128;
M.L = L;
Sw = ones(1,M.n_el).*Swc;
Sw(1) = 1-Sor;
Sw(end) = Swc;
x128 = linspace(0,L,128);
t=t0;
[Sw128] = main_solver(M, phi, mi_o, mi_w, Sw, Swc, Sor, n_o, n_w, t, t_final, solver);

%% Run Test - mesh with 256 elements
 
solver = 'FOUM';
M = Mesh;
M.n_el = 256;
M.L = L;
Sw = ones(1,M.n_el).*Swc;
Sw(1) = 1-Sor;
Sw(end) = Swc;
x256 = linspace(0,L,256);
t = t0;
[Sw256] = main_solver(M, phi, mi_o, mi_w, Sw, Swc, Sor, n_o, n_w, t, t_final, solver);

%% Run Test - mesh with 512 elements
 
solver = 'FOUM';
M = Mesh;
M.n_el =512;
M.L = L;
Sw = ones(1,M.n_el).*Swc;
Sw(1) = 1-Sor;
Sw(end) = Swc;
x512 = linspace(0,L,512);
t = t0;
[Sw512] = main_solver(M, phi, mi_o, mi_w, Sw, Swc, Sor, n_o, n_w, t, t_final, solver);


%% Run Test - mesh with 1024 elements
 
solver = 'FOUM';
M = Mesh;
M.n_el =1024;
M.L = L;
Sw = ones(1,M.n_el).*Swc;
Sw(1) = 1-Sor;
Sw(end) = Swc;
x1024 = linspace(0,L,1024);
t = t0;
[Sw1024] = main_solver(M, phi, mi_o, mi_w, Sw, Swc, Sor, n_o, n_w, t, t_final, solver);

%% Analytical Solution using Welge's Approach
[Sw, Xsw, Swf, nts] = Welge_Solution();
 %% Plotting Area
 fig1 = figure(1);
 set(fig1, 'color', 'w')
 hold on
 title('Water Saturation Profile - FOUM')
 xlabel('Distance')
 ylabel('Water Saturation')
 plot(x8,Sw8,'r')
 hold on
 plot(x16, Sw16, 'k')
 plot(x32, Sw32, 'b')
 plot(x64, Sw64, 'm')
 plot(x128, Sw128, 'Color', [0.01 0.012 0.782])
 plot(x256, Sw256, 'Color', [0.4 0.15 0.2])
 plot(x512, Sw512, 'Color', [0.02 0.8 0.2])
 plot(x1024, Sw1024, 'Color', [0.67 0.6 0.07])
 hold on
 plot(Xsw(:,nts), Swt, '-r', 'LineWidth', 2.0);
 hold on
 plot([0 Xsw(1,nts)], [1-Sor 1-Sor], '-r', 'LineWidth', 2.0);
 plot([Xsw(end,nts) Xsw(end,nts)], [Swf Swc], '-r', 'LineWidth', 2.0);
        plot([Xsw(end,nts) L], [Swc Swc], '-r', 'LineWidth', 2.0);
 legend('8 elements', '16 elements', '32 elements', '64 elements', ...
     '128 elements', '256 elements', '512 elements', '1024 elements', ...
     'Welge Analytical Solution')
 grid minor

% 5.4 Plot the water saturation profiles
%{
h_fig4 = figure(4);
set(h_fig4, 'color', 'w', 'NumberTitle', 'off', 'Name', 'Water Saturation Profiles');

SwTime= [];
ti = nts;
plot(Xsw(:,ti), Swt, '-r', 'LineWidth', 2.0);
hold on;
plot([0 Xsw(1,ti)], [1-Sor 1-Sor], '-r', 'LineWidth', 2.0);
SwTime = [ SwTime; ['Time = ' num2str(t(ti)/86400, format) ' day'] ];
plot([Xsw(end,ti) Xsw(end,ti)], [Swf Swc], '-r', 'LineWidth', 2.0);
plot([Xsw(end,ti) L], [Swc Swc], '-r', 'LineWidth', 2.0);

set(gca,'YLim', [0 1], 'YTick', 0:0.2:1, 'XLim', [0 L]);
set(gca, 'Fontname', 'Times New Roman','FontSize',10);
title('(d) water saturation profiles')
xlabel('{\it x} (m)');
ylabel('{\it S}_w');
h_legend2 = legend(SwTime);
set(h_legend2,'Box', 'on', 'Location', 'best');
axis square;
%}