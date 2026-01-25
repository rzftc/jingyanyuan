虚拟电厂调节能力分析与优化系统 (VPP Regulatory Capacity Analysis System)

这是一个基于 MATLAB 的综合仿真与分析平台，旨在评估、预测和优化由 空调 (AC) 和 电动汽车 (EV) 构成的虚拟电厂 (VPP) 的电网调节能力。该项目涵盖了从底层物理建模、用户行为不确定性分析、大规模蒙特卡洛模拟到上层风险约束优化调度的全过程。

📁 目录结构与功能模块

代码库按照功能逻辑划分为多个子文件夹，核心模块说明如下：

基础数据与建模

0input_data/: 数据生成器。包含生成 AC/EV 基础参数（额定功率、电池容量、热力学参数）和行为模式（出行时间、设定温度）的脚本，支持 Excel 导出。

核心脚本: generateACParameters.m, generateEVData.m

1initialize/: 初始化接口。负责读取 Excel 配置文件并实例化 MATLAB 结构体对象。

2AC/: 空调物理模型。基于一阶等效热参数模型 (ETP)，包含状态方程更新、聚合参数 (A, B, C) 计算及 PI 控制逻辑。

3EV/: 电动汽车物理模型。包含 SOC 状态演化、基线充电功率计算、充放电死区及灵活性潜力评估逻辑。

5userUncertainties/: 用户行为模型。包含基于激励电价的用户响应概率模型 (Price-Incentive Response)。

预测与分析核心

9predict/: 规模预测模块。

集成多种算法：Bass 扩散模型、Logistic 模型、ARIMA、灰色预测 GM(1,1)、随机森林 (Random Forest)、支持向量回归 (SVR)。

用于预测未来年份（至 2040 年）省/节点级的设备保有量。

10spaceRelation/: 空间特性分析。基于地理信息分析设备在不同功能区（居民区、工作区等）的分布特征。

6SensAnalyses/: 灵敏度分析。评估激励价格变化对聚合调节潜力的影响。

仿真与优化调度

7main/: 仿真主程序。

AC_main_Stateful_Simulation.m: 执行 AC 集群的状态化时序仿真，生成调节潜力基线。

AC_Result_Plotter*.m: 包含丰富的绘图工具，用于可视化仿真结果。

4analysis/ & 4.1analysis/: 确定性优化求解。

包含 贪心算法 (Greedy)、遗传算法 (GA)、混合整数线性规划 (MILP)。

引入 PTDF (功率传输分布因子) 矩阵处理网络潮流约束。

考虑成本最小化与互补性指标 (SDCI, Rho) 的多目标优化。

11cvar/: 风险约束随机优化 (核心)。

CVaR 模型: 基于条件风险价值的随机规划，量化调节缺额风险。

鲁棒优化: 包含与传统鲁棒优化的对比分析 (run_scenario_H_robust_comparison.m)。

场景生成: 基于蒙特卡洛模拟生成大量不确定性场景。

🚀 核心算法与技术亮点

物理驱动的聚合建模:

AC: 实现了基于状态一致性的聚合模型，将成千上万台空调的异构参数映射为统一的 S(t+1) = A*S(t) + B*P(t) + C 状态空间方程。

EV: 考虑出行链约束和电池物理特性的能量边界模型。

多层级优化架构:

上层 (GA): 优化参与调节的设备组合，平衡系统经济性与调节指标（SDCI/Rho）。

下层 (MILP/Greedy): 在满足网络潮流 (PTDF) 和物理约束的前提下，进行精确的功率分配。

风险管理:

引入 CVaR (Conditional Value at Risk) 理论，在调度策略中显式考虑极端场景下的违约风险，平衡经济性与可靠性。

🛠️ 环境要求与安装

基础环境

MATLAB: 推荐 R2020b 或更高版本。

必需工具箱 (Toolboxes)

为确保所有脚本正常运行，请安装以下工具箱：

Optimization Toolbox: 用于 linprog, intlinprog, quadprog 等求解器。

Global Optimization Toolbox: 用于 ga (遗传算法)。

Statistics and Machine Learning Toolbox: 用于 TreeBagger (随机森林), fitrsvm (SVR) 及相关统计函数。

Econometrics Toolbox: 用于 arima (时间序列预测)。

可选求解器

IBM CPLEX: 代码中部分高级调度脚本（如 solve_deterministic_dispatch.m 中的 MILP 问题）预留了 CPLEX 接口。如果未安装，代码通常会回退到 MATLAB 内置求解器或需要手动调整选项。

📈 使用指南与工作流程 (Workflow)

建议按照以下逻辑顺序运行项目代码：

步骤 1: 数据准备与预测

运行 0input_data/generateExampleExcel_real_24.m 生成初始的设备参数 Excel 文件。

(可选) 运行 9predict/predict_node_level_new.m 预测未来的设备规模增长趋势。

步骤 2: 物理仿真与潜力评估

运行 7main/经研院代码/AC_main_Stateful_Simulation.m。

该脚本将读取 Excel 数据，执行时序仿真，计算聚合调节潜力（上/下调边界）。

结果将保存为 .mat 文件，供后续优化使用。

步骤 3: 场景生成与风险建模

运行 11cvar/main_scenario_generation_soc.m。

基于蒙特卡洛方法生成大量随机场景，提取可靠调节域。

步骤 4: 优化调度

VPP 调度: 运行 4.1analysis/main_vpp_optimizer.m。

加载分块仿真结果，执行考虑网络约束的优化调度。

风险分析: 运行 11cvar/run_scenario_H_robust_comparison.m。

对比 CVaR 随机优化与传统鲁棒优化的经济性和安全性。

📊 关键指标说明

指标

全称

说明

SDCI

Signed Direct Complementarity Index

同向互补性指数。衡量 AC 和 EV 在同一调节方向（如同为削峰）上的能力互补程度。

Rho

Spearman's Rank Correlation

斯皮尔曼秩相关系数。评估 AC 和 EV 调节潜力在时间序列上的相关性趋势。

CVaR

Conditional Value at Risk

条件风险价值。用于量化尾部风险（极端场景下的调节缺额）。

📝 注意事项

并行计算: 仿真脚本广泛使用了 parfor。建议在运行前使用 parpool 开启 MATLAB 并行池，以显著提高计算速度。

路径设置: 请确保所有子文件夹（如 1initialize, 2AC 等）都已添加到 MATLAB 的搜索路径中，或者在根目录下运行脚本。

数据文件: 示例 Excel 文件通常生成在 0input_data 目录下，主程序读取时请确认文件名一致。