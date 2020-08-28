function [Sw] = main_solver(M, phi, mi_o, mi_w, Sw, Swc, Sor, n_o, n_w, t, t_final, solver)
    v = 1; %velocidade total unitária - facilitar as contas
    while t<t_final
        [kro, krw] = Relative_Perm(Sw, Sor, Swc, n_o, n_w);
        lamb_o = kro./mi_o;
        lamb_w = krw./mi_w;
        fw = lamb_w./(lamb_o + lamb_w);
        if solver == 'FOUM'
            [dfw_dSw, dSw] = flux_FOUM(M, fw, Sw);
        else
             [dfw_dSw, dSw] = flux_FOUM(M, fw, Sw);
        end
        dt = Update_Time(M, dfw_dSw, v, t, t_final);
        vols_internal = setdiff(M.all_vols_ID, M.contour_vols_ID);
        dx = M.L/M.n_el;
        Sw(vols_internal) = Sw(vols_internal) + dt * v/phi * dfw_dSw(vols_internal) .* dSw(vols_internal)/dx;
        t = t + dt;
    end
end