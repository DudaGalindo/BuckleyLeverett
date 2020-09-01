function [dfwdSw, dfw_vols] = flux_MUSCL(M, Sw, mi_w, mi_o, Sor, Swc, n_o, n_w)
%% Obten��o dos gradientes 
conec = M.faces_conec;
dSw_face = Sw(conec(:,2)) -Sw(conec(:,1));
dSw_face_neig = zeros(length(conec(:,2)),2);
dSw_face_neig(2:end,1) = dSw_face(1:end-1);
dSw_face_neig(1:end-1,2) = dSw_face(2:end);
dSw_face_neig(end,2) = 0;
dSw_face_neig(1,1) = 0;

%% C�lculo do r:
r = zeros(length(conec(:,2)),2);
r(:,1) = dSw_face./transpose(dSw_face_neig(:,1));
r(:,2) = dSw_face./transpose(dSw_face_neig(:,2));
r(dSw_face_neig==0)=0;

%% C�lculo da fun��o limitadora (usando aqui Van Leer):
phi = (r + abs(r))./(r+1);
phi(:,2) = -phi(:,2);

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

dfw_dSw_face = zeros(length(conec(:,2)),2);
delta = 0.001;
for i=1:2
    Sw_face_plus = Sw_face(:,i);
    Sw_face_minus = Sw_face(:,i);
    Sw_face_plus = Sw_face_plus + delta/2;
    Sw_face_minus = Sw_face_minus - delta/2;
    [kro_face_plus, krw_face_plus] = Relative_Perm(Sw_face_plus, Sor, Swc, n_o, n_w);
    [kro_face_minus, krw_face_minus] = Relative_Perm(Sw_face_minus, Sor, Swc, n_o, n_w);
    lamb_w_face_plus = krw_face_plus/mi_w;
    lamb_o_face_plus = kro_face_plus/mi_o;
    lamb_w_face_minus = krw_face_minus/mi_w;
    lamb_o_face_minus = kro_face_minus/mi_o;
    fw_face_plus = lamb_w_face_plus./(lamb_o_face_plus + lamb_w_face_plus);
    fw_face_minus = lamb_w_face_minus./(lamb_o_face_minus + lamb_w_face_minus);
    dfw_dSw_face(:,i) = (fw_face_plus - fw_face_minus)/delta;
end
alpha = max(abs(dfw_dSw_face),[],2);

%% C�lculo do Fluxo por LLF:
fw_face = 1/2*(sum(fw_face,2) - alpha.*(Sw_face(:,2) - Sw_face(:,1)));

cph = 1;
lines = cat(1,repmat(cph,length(conec(:,1)),1), repmat(cph, length(conec(:,2)),1));
cols = cat(1,conec(:,1),conec(:,2));
data = cat(1,-fw_face, fw_face);
dfw_vols = sparse(lines,cols,data);

lines = cat(1,repmat(cph,length(conec(:,1)),1), repmat(cph, length(conec(:,2)),1));
cols = cat(1,conec(:,1),conec(:,2));
data = cat(1,-Sw_face(:,2), Sw_face(:,1));
dSw_vols = sparse(lines,cols,data);
dSw_vols(1) = Sw_face(1,1) - Sw(1);
dSw_vols(end) = Sw_face(end,2) - Sw(end);

dfwdSw = dfw_vols./dSw_vols;
dfwdSw(dSw_vols==0)=0;

end