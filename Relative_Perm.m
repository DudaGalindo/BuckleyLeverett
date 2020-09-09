%% Função criada para o cálculo da permeabilidade relativa
% É considerado que kro_0 = krw_0 = 1.

function [kro, krw] = Relative_Perm(Sw, Sor, Swc, n_o, n_w)
    kro = ((1 - Sw - Swc)/(1 - Swc - Sor)).^n_o;
    krw = ((Sw - Swc)/(1 - Swc - Sor)).^n_w;
end