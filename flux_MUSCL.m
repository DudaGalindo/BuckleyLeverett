function [dfw_vols] = flux_MUSCL(M, Sw, mi_w, mi_o, Sor, Swc, n_o, n_w, solver)
%% Obten��o dos gradientes 
conec = M.faces_conec;
dSw_face = Sw(conec(:,2)) -Sw(conec(:,1));
dSw_vols(M.contour_vols_ID) = 0;
for i=2:M.n_el-1
    dSw_vols(i) = Sw(i+1) - Sw(i-1); 
end
dSw_face_neig = zeros(length(conec(:,1)),2);
dSw_face_neig(:,1) = dSw_vols(conec(:,1)) - dSw_face;
dSw_face_neig(:,2) = dSw_vols(conec(:,2)) - dSw_face;

%% C�lculo do r:
r = zeros(length(conec(:,2)),2);
r(:,1) = dSw_face./transpose(dSw_face_neig(:,1));
r(:,2) = dSw_face./transpose(dSw_face_neig(:,2));
r(dSw_face_neig==0)=0;


%% C�lculo da fun��o limitadora (usando aqui Van Leer):
phi = (r + abs(r))./(r + 1);
phi(:,2) = -phi(:,2);
%phi = (r.^2 + r)./(r.^2 + 1);
%phi(:,2) = -phi(:,2);
phi(r==-1)=0;

%% Extrapola��o da fun��o � esquerda e � direita da face:
Sw_face = Sw(conec) + phi/2 .* dSw_face_neig;

%% C�lculo da permeabilidade nos dois lados da face:
kro_face = zeros(length(conec(:,2)),2);
krw_face = zeros(length(conec(:,2)),2);
for i=1:2
    [kro_face(:,i), krw_face(:,i)] = Relative_Perm(Sw_face(:,i), Sor, Swc, n_o, n_w);
end

%% C�lculo das mobilidades nos dois lados da face
lamb_w_face = krw_face/mi_w;
lamb_o_face = kro_face/mi_o;
fw_face = lamb_w_face./(lamb_o_face + lamb_w_face);

dfw_dSw_face = zeros(length(conec(:,2)),3);
delta = 0.001;
Swm_face = sum(Sw_face,2)*0.5;
for i=1:2
    Sw_face_plus = Sw_face(:,i);
    Sw_face_minus = Sw_face(:,i);
    Sw_face_plus = Sw_face_plus + delta/2;
    Sw_face_plus(Sw_face_plus>1) = 1.;
    Sw_face_minus = Sw_face_minus - delta/2;
    Sw_face_minus(Sw_face_minus<0) = 0;
    [kro_face_plus, krw_face_plus] = Relative_Perm(Sw_face_plus, Sor, Swc, n_o, n_w);
    [kro_face_minus, krw_face_minus] = Relative_Perm(Sw_face_minus, Sor, Swc, n_o, n_w);
    lamb_w_face_plus = krw_face_plus/mi_w;
    lamb_o_face_plus = kro_face_plus/mi_o;
    lamb_w_face_minus = krw_face_minus/mi_w;
    lamb_o_face_minus = kro_face_minus/mi_o;
    fw_face_plus = lamb_w_face_plus./(lamb_o_face_plus + lamb_w_face_plus);
    fw_face_minus = lamb_w_face_minus./(lamb_o_face_minus + lamb_w_face_minus);
    dfw_dSw_face(:,i) = (fw_face_plus - fw_face_minus)./(Sw_face_plus - Sw_face_minus);    
end

Swm_face_minus = Swm_face;
Swm_face_plus = Swm_face;
Swm_face_plus = Swm_face_plus + delta/2;
Swm_face_plus(Swm_face_plus>1) = 1.;
Swm_face_minus = Swm_face_minus - delta/2;
Swm_face_minus(Swm_face_minus<0) = 0;

[krom_face_plus, krwm_face_plus] = Relative_Perm(Swm_face_plus, Sor, Swc, n_o, n_w);
[krom_face_minus, krwm_face_minus] = Relative_Perm(Swm_face_minus, Sor, Swc, n_o, n_w);
lambm_w_face_plus = krwm_face_plus/mi_w;
lambm_o_face_plus = krom_face_plus/mi_o;
lambm_w_face_minus = krwm_face_minus/mi_w;
lambm_o_face_minus = krom_face_minus/mi_o;
fwm_face_plus = lambm_w_face_plus./(lambm_o_face_plus + lambm_w_face_plus);
fwm_face_minus = lambm_w_face_minus./(lambm_o_face_minus + lambm_w_face_minus);
dfw_dSw_face(:,3) = (fwm_face_plus - fwm_face_minus)./(Swm_face_plus - Swm_face_minus);

alpha = max(abs(dfw_dSw_face),[],2);
%alpha_HG(abs(Sw_face(:,2) - Sw_face(:,1))<1e-2)=max(abs(dfw_dSw_face(abs(Sw_face(:,2) - Sw_face(:,1))<1e-2)),[],2);
%% C�lculo do Fluxo por LLF:
if strcmp(solver, 'LLF')
    fw_face = 1/2*(sum(fw_face,2) - alpha.*(Sw_face(:,2) - Sw_face(:,1)));
end

%% C�lculo do fluxo por Roe
if strcmp(solver, 'Roe')
    % c�lculo do alpha por Rankine-Hugoniot = A (matriz de Roe).
    alpha_RH = (fw_face(:,2) - fw_face(:,1))./(Sw_face(:,2) - Sw_face(:,1));
    alpha_RH((Sw_face(:,2) - Sw_face(:,1)==0)) = 0;
    A = abs(alpha_RH);
    fw_face = 1/2*(sum(fw_face,2) - A.*(Sw_face(:,2) - Sw_face(:,1)));
end

% Abaixo � mostrado um artificio usando matriz esparsa pra calcular o
% fluxo fracional no volume (dfw_vols) como o balan�o do que entra 
% menos o que sai pelas faces. 
% Esse artificio foi utilizado aqui considerando que podem 
% haver problemas 2-D e 3-D tamb�m.
% Esse artificio tamb�m � usado pra calcular o delta de
% satura��o (dSw_vols) no volume. 
cph = 1;
lines = cat(1,repmat(cph,length(conec(:,1)),1), repmat(cph, length(conec(:,2)),1));
cols = cat(1,conec(:,1),conec(:,2));
data = cat(1,-fw_face, fw_face);
dfw_vols = sparse(lines,cols,data);


end