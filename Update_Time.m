%% Fun��o para calcular o passo de tempo
% OBS: Na pr�tica est� sendo utilizado aqui um passo de tempo fixo, uma vez
% que v � constante, vide (CARVALHO,2005).
function [dt] = Update_Time(M, v, t, t_final)
    CFL = 0.05;
    dx = M.L/M.n_el;
    dt = CFL *(abs(dx./(v)));
    if t + dt > t_final
        dt = t_final - t;
    end
end