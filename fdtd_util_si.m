function txt = fdtd_util_si(v)

    if (v >= 1e12)
        txt = sprintf("%g T", v/1e12);
    elseif (v >= 1e9)
        txt = sprintf("%g G", v/1e9);
    elseif (v >= 1e6)
        txt = sprintf("%g M", v/1e6);
    elseif (v >= 1e3)
        txt = sprintf("%g k", v/1e3);
    elseif (v >= 1)
        txt = sprintf("%g", v);
    elseif (v >= 1e-3)
        txt = sprintf("%g m", v*1e3);
    elseif (v >= 1e-6)
        txt = sprintf("%g u", v*1e6);
    elseif (v >= 1e-9)
        txt = sprintf("%g n", v*1e9);
    elseif (v >= 1e-12)
        txt = sprintf("%g p", v*1e12);
    else
        txt = sprintf("%g f", v*1e15);
    end

end

