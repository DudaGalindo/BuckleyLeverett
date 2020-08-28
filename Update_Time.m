function [dt] = Update_Time(M, dfw_dSw, v, t, t_final)
    CFL = 0.5;
    dx = M.L/M.n_el;
    dt = CFL * min(abs(dx./(v*dfw_dSw)));
    if t + dt > t_final
        dt = t_final - t;
    end
end