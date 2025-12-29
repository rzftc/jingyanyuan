function rho = calculate_Rho_local(n_AC, deltaP_AC, n_EV, deltaP_EV)
    AC_series = n_AC .* deltaP_AC;
    EV_series = n_EV .* deltaP_EV;
    if length(AC_series) > 2
        rho = corr(AC_series, EV_series, 'Type', 'Spearman');
    else
        rho = 0;
    end
end