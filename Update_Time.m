function [dt] = Update_Time(M, dfw_dSw, v, t, t_final, type)
    CFL = 0.05;
    if strcmp(type,'MUSCL')
        CFL = 0.05;
    end
    dx = M.L/M.n_el;
    dt = CFL *(abs(dx./(v)));
    if t + dt > t_final
        dt = t_final - t;
    end
end