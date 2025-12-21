function SDCI = calculate_SDCI_local(n_AC, n_EV, deltaP_AC, deltaP_EV)
    AC = n_AC .* deltaP_AC;
    EV = n_EV .* deltaP_EV;
    min_vals = min(AC, EV);
    max_vals = max(AC, EV);
    total_max = sum(max_vals);
    if total_max == 0
        SDCI = 0;
    else
        SDCI = sum(min_vals) / total_max;
    end
end