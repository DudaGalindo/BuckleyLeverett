%% Função principal
function [Sw] = main_solver(M, phi, mi_o, mi_w, Sw, Swc, Sor, n_o, n_w, t, t_final, solver)
    v = 1; %velocidade total unitária - simplificacao
    while t<t_final
        [kro, krw] = Relative_Perm(Sw, Sor, Swc, n_o, n_w); %calculo da permeabilidade relativa
        lamb_o = kro./mi_o; %mobilidade da fase óleo
        lamb_w = krw./mi_w; %mobilidade da fase água
        fw = lamb_w./(lamb_o + lamb_w); %fluxo fracional da fase água
        
        % Abaixo obtemos dfw_vols - o fluxo resultante no volume, de acordo
        %com o fluxo calculado nas faces que compõem o seu contorno.
        if strcmp(solver,'FOUM')
            [dfw_vols] = flux_FOUM(M, fw); 
        else
             [dfw_vols] = flux_MUSCL(M, Sw, mi_w, mi_o, Sor, Swc, n_o, n_w, solver);
        end
        dt = Update_Time(M, v, t, t_final); %calculo do passo de tempo
        dx = M.L/M.n_el; % espaçamento dos blocos da malha
        Sw(2:end) = Sw(2:end) + dt * v/phi * dfw_vols(2:end)/dx; %cálculo da saturação
        t = t + dt; %atualizando o tempo decorrido de simulação 
    end
end