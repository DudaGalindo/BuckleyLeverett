%% Fun��o principal
function [Sw] = main_solver(M, phi, mi_o, mi_w, Sw, Swc, Sor, n_o, n_w, t, t_final, solver)
    v = 1; %velocidade total unit�ria - simplificacao
    while t<t_final
        [kro, krw] = Relative_Perm(Sw, Sor, Swc, n_o, n_w); %calculo da permeabilidade relativa
        lamb_o = kro./mi_o; %mobilidade da fase �leo
        lamb_w = krw./mi_w; %mobilidade da fase �gua
        fw = lamb_w./(lamb_o + lamb_w); %fluxo fracional da fase �gua
        
        % Abaixo obtemos dfw_vols - o fluxo resultante no volume, de acordo
        %com o fluxo calculado nas faces que comp�em o seu contorno.
        if strcmp(solver,'FOUM')
            [dfw_vols] = flux_FOUM(M, fw); 
        else
             [dfw_vols] = flux_MUSCL(M, Sw, mi_w, mi_o, Sor, Swc, n_o, n_w, solver);
        end
        dt = Update_Time(M, v, t, t_final); %calculo do passo de tempo
        dx = M.L/M.n_el; % espa�amento dos blocos da malha
        Sw(2:end) = Sw(2:end) + dt * v/phi * dfw_vols(2:end)/dx; %c�lculo da satura��o
        t = t + dt; %atualizando o tempo decorrido de simula��o 
    end
end