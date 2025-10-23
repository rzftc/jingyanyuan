% Filename: calculate_dispatch_cost_and_mc.m
function [cost, marginal_cost] = calculate_dispatch_cost_and_mc(delta_P, c1, c2)
    % 计算给定调节功率下的成本和边际成本（基于二次成本函数）
    % Inputs:
    %   delta_P: 调节功率 (kW)
    %   c1: 线性成本系数 (元/kW)
    %   c2: 二次成本系数 (元/kW^2)
    % Outputs:
    %   cost: 总调节成本 (元)
    %   marginal_cost: 边际成本 (元/kW)

    cost = c1 * delta_P + c2 * (delta_P^2);
    marginal_cost = c1 + 2 * c2 * delta_P;
end