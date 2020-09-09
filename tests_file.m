% Arquivo criado para fazer a simula��o de todos os casos e plotar os
% resultados.

%% Input data file 
L = 1.0; % domain length [m]
A = 1.0; % area of cross-section [m^2]
phi = 1; % porosity
mi_o = 1.0e-3; % oil phase viscosity [Pa*s]
mi_w = 1.0e-3; % water phase viscosity [Pa*s]

% Relative permeability parameters
Sor = 0.1; % residual oil saturation
Swc = 0.1; % connate water saturation
n_o = 2.00; % exponent of oil phase
n_w = 2.00; % exponent of water phase

% Other initial parameters
t0 = 0.; % initial calculation time [s]
t_final = 0.2; % final calculation time [s]

%% Run Test - mesh with 8 elements
M = Mesh;
M.n_el = 8;
M.L = L;
Sw = ones(1,M.n_el).*Swc;
Sw(1) = 1-Sor;
Sw(end) = Swc;
x8 = linspace(0,L,8);
t = t0;
solver = 'FOUM';
[Sw8] = main_solver(M, phi, mi_o, mi_w, Sw, Swc, Sor, n_o, n_w, t, t_final, solver);
solver = 'LLF';
[Sw8_MUSCL_LLF] = main_solver(M, phi, mi_o, mi_w, Sw, Swc, Sor, n_o, n_w, t, t_final, solver);
solver = 'Roe';
[Sw8_MUSCL_Roe] = main_solver(M, phi, mi_o, mi_w, Sw, Swc, Sor, n_o, n_w, t, t_final, solver);

%% Run Test - mesh with 16 elements
 
M = Mesh;
M.n_el = 16;
M.L = L;
Sw = ones(1,M.n_el).*Swc;
Sw(1) = 1-Sor;
Sw(end) = Swc;
x16 = linspace(0,L,16);
t = t0;
solver = 'FOUM';
[Sw16] = main_solver(M, phi, mi_o, mi_w, Sw, Swc, Sor, n_o, n_w, t, t_final, solver);
solver = 'LLF';
[Sw16_MUSCL_LLF] = main_solver(M, phi, mi_o, mi_w, Sw, Swc, Sor, n_o, n_w, t, t_final, solver);
solver = 'Roe';
[Sw16_MUSCL_Roe] = main_solver(M, phi, mi_o, mi_w, Sw, Swc, Sor, n_o, n_w, t, t_final, solver);

 
%% Run Test - mesh with 32 elements
 
M = Mesh;
M.n_el = 32;
M.L = L;
Sw = ones(1,M.n_el).*Swc;
Sw(1) = 1-Sor;
Sw(end) = Swc;
x32 = linspace(0,L,32);
t = t0;
solver = 'FOUM';
[Sw32] = main_solver(M, phi, mi_o, mi_w, Sw, Swc, Sor, n_o, n_w, t, t_final, solver);
solver = 'LLF';
[Sw32_MUSCL_LLF] = main_solver(M, phi, mi_o, mi_w, Sw, Swc, Sor, n_o, n_w, t, t_final, solver);
solver = 'Roe';
[Sw32_MUSCL_Roe] = main_solver(M, phi, mi_o, mi_w, Sw, Swc, Sor, n_o, n_w, t, t_final, solver);
 
%% Run Test - mesh with 64 elements
 
M = Mesh;
M.n_el =64;
M.L = L;
Sw = ones(1,M.n_el).*Swc;
Sw(1) = 1-Sor;
Sw(end) = Swc;
x64 = linspace(0,L,64);
t = t0;
solver = 'FOUM';
[Sw64] = main_solver(M, phi, mi_o, mi_w, Sw, Swc, Sor, n_o, n_w, t, t_final, solver);
solver = 'LLF';
[Sw64_MUSCL_LLF] = main_solver(M, phi, mi_o, mi_w, Sw, Swc, Sor, n_o, n_w, t, t_final, solver);
solver = 'Roe';
[Sw64_MUSCL_Roe] = main_solver(M, phi, mi_o, mi_w, Sw, Swc, Sor, n_o, n_w, t, t_final, solver);

%% Run Test - mesh with 128 elements

M = Mesh;
M.n_el =128;
M.L = L;
Sw = ones(1,M.n_el).*Swc;
Sw(1) = 1-Sor;
Sw(end) = Swc;
x128 = linspace(0,L,128);
t=t0;
solver = 'FOUM';
[Sw128] = main_solver(M, phi, mi_o, mi_w, Sw, Swc, Sor, n_o, n_w, t, t_final, solver);
solver = 'LLF';
[Sw128_MUSCL_LLF] = main_solver(M, phi, mi_o, mi_w, Sw, Swc, Sor, n_o, n_w, t, t_final, solver);
solver = 'Roe';
[Sw128_MUSCL_Roe] = main_solver(M, phi, mi_o, mi_w, Sw, Swc, Sor, n_o, n_w, t, t_final, solver);

%% Run Test - mesh with 256 elements

M = Mesh;
M.n_el = 256;
M.L = L;
Sw = ones(1,M.n_el).*Swc;
Sw(1) = 1-Sor;
Sw(end) = Swc;
x256 = linspace(0,L,256);
t = t0;
solver = 'FOUM';
[Sw256] = main_solver(M, phi, mi_o, mi_w, Sw, Swc, Sor, n_o, n_w, t, t_final, solver);
solver = 'LLF';
[Sw256_MUSCL_LLF] = main_solver(M, phi, mi_o, mi_w, Sw, Swc, Sor, n_o, n_w, t, t_final, solver);
solver = 'Roe';
[Sw256_MUSCL_Roe] = main_solver(M, phi, mi_o, mi_w, Sw, Swc, Sor, n_o, n_w, t, t_final, solver);

%% Run Test - mesh with 512 elements
 

M = Mesh;
M.n_el =512;
M.L = L;
Sw = ones(1,M.n_el).*Swc;
Sw(1) = 1-Sor;
Sw(end) = Swc;
x512 = linspace(0,L,512);
t = t0;
solver = 'FOUM';
[Sw512] = main_solver(M, phi, mi_o, mi_w, Sw, Swc, Sor, n_o, n_w, t, t_final, solver);
solver = 'LLF';
[Sw512_MUSCL_LLF] = main_solver(M, phi, mi_o, mi_w, Sw, Swc, Sor, n_o, n_w, t, t_final, solver);
solver = 'Roe';
[Sw512_MUSCL_Roe] = main_solver(M, phi, mi_o, mi_w, Sw, Swc, Sor, n_o, n_w, t, t_final, solver);


%% Run Test - mesh with 1024 elements
 

M = Mesh;
M.n_el =1024;
M.L = L;
Sw = ones(1,M.n_el).*Swc;
Sw(1) = 1-Sor;
Sw(end) = Swc;
x1024 = linspace(0,L,1024);
t = t0;
solver = 'FOUM';
[Sw1024] = main_solver(M, phi, mi_o, mi_w, Sw, Swc, Sor, n_o, n_w, t, t_final, solver);
solver = 'LLF';
[Sw1024_MUSCL_LLF] = main_solver(M, phi, mi_o, mi_w, Sw, Swc, Sor, n_o, n_w, t, t_final, solver);
solver = 'Roe';
[Sw1024_MUSCL_Roe] = main_solver(M, phi, mi_o, mi_w, Sw, Swc, Sor, n_o, n_w, t, t_final, solver);

%% Analytical Solution using Welge's Approach
k = 1; % permeabilidade absoluta
[index, dfw, Sw, Xsw, Swf, nts] = Welge_Solution(L, A, phi, k, mi_o, mi_w, Sor, Swc, n_o, n_w); %chamando funcao
t = linspace(t0, t_final, nts);
dfwt = dfw(end:-1:index); % inverted sequence of dfw
Swt = Sw(end:-1:index); % inverted sequence of Sw
ti = nts;
qt = 1.0;
nsw = size(Sw,1);
for i = 1:(nsw-index+1)
    Xsw(i,ti) = qt*t(ti)/(A*phi)*dfwt(i);
end

 %% Plotting Area
 
 % Plot results for FOUM only
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

 % Plot results for MUSCL-LLF only
fig2 = figure(2);
 set(fig2, 'color', 'w')
 hold on
 title('Water Saturation Profile - MUSCL + LLF')
 xlabel('Distance')
 ylabel('Water Saturation')
 plot(x8,Sw8_MUSCL_LLF,'r')
 hold on
 plot(x16, Sw16_MUSCL_LLF, 'k')
 plot(x32, Sw32_MUSCL_LLF, 'b')
 plot(x64, Sw64_MUSCL_LLF, 'm')
 plot(x128, Sw128_MUSCL_LLF, 'Color', [0.01 0.012 0.782])
 plot(x256, Sw256_MUSCL_LLF, 'Color', [0.4 0.15 0.2])
 plot(x512, Sw512_MUSCL_LLF, 'Color', [0.02 0.8 0.2])
 plot(x1024, Sw1024_MUSCL_LLF, 'Color', [0.67 0.6 0.07])
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
 
 % Plot results for MUSCL-Roe only
 fig3 = figure(3);
 set(fig3, 'color', 'w')
 hold on
 title('Water Saturation Profile - MUSCL + Roe')
 xlabel('Distance')
 ylabel('Water Saturation')
 plot(x8,Sw8_MUSCL_Roe,'r')
 hold on
 plot(x16, Sw16_MUSCL_Roe, 'k')
 plot(x32, Sw32_MUSCL_Roe, 'b')
 plot(x64, Sw64_MUSCL_Roe, 'm')
 plot(x128, Sw128_MUSCL_Roe, 'Color', [0.01 0.012 0.782])
 plot(x256, Sw256_MUSCL_Roe, 'Color', [0.4 0.15 0.2])
 plot(x512, Sw512_MUSCL_Roe, 'Color', [0.02 0.8 0.2])
 plot(x1024, Sw1024_MUSCL_Roe, 'Color', [0.67 0.6 0.07])
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
 
 % Plot results for comparison between methods
 fig4 = figure(4);
 set(fig4, 'color', 'w')
 hold on
 title('Water Saturation Profile')
 xlabel('Distance')
 ylabel('Water Saturation')
 hold on
 plot(x256, Sw256, 'Color', [.5 0.7 0.2])
 plot(x256, Sw256_MUSCL_LLF, 'Color', [0.4 0.15 0.2])
 hold on
 plot(Xsw(:,nts), Swt, '-r', 'LineWidth', 2.0);
 hold on
 plot([0 Xsw(1,nts)], [1-Sor 1-Sor], '-r', 'LineWidth', 2.0);
 plot([Xsw(end,nts) Xsw(end,nts)], [Swf Swc], '-r', 'LineWidth', 2.0);
        plot([Xsw(end,nts) L], [Swc Swc], '-r', 'LineWidth', 2.0);
 legend('256 elements - UPW','256 elements - MUSCL',  ...
     'Welge Analytical Solution')

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