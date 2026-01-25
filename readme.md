# 用户侧柔性资源调节能力评估与虚拟电厂(VPP)仿真系统

## 项目简介

本项目是一个基于 MATLAB 开发的综合仿真平台，旨在评估和优化用户侧分布式资源（主要是空调系统 AC 和电动汽车 EV）的调节潜力。系统集成了资源建模、不确定性分析、空间特性分析、风险约束优化（CVaR）以及规模预测等模块，用于支持虚拟电厂（VPP）的聚合调度策略研究。

主要功能包括：

- **资源建模**：精细化的空调（热力学模型）和电动汽车（出行与充电行为）负荷建模。
- **调节能力评估**：计算不同激励和价格信号下的聚合调节功率（上调/下调）。
- **优化调度**：基于遗传算法（GA）和贪心算法的日前/实时调度优化。
- **风险管理**：引入条件风险价值（CVaR）评估调度策略在不确定性下的风险。
- **空间分析**：分析资源在不同功能区域的空间分布特性。
- **趋势预测**：利用机器学习（随机森林、SVR）和统计模型预测负荷及设施规模增长。

## 目录结构说明

代码库主要包含以下模块目录：

### 1. 数据与初始化

- **`0input_data/`**: 负责生成仿真所需的原始数据。
   - 包含生成空调参数 (`generateACParameters.m`) 和电动汽车出行数据 (`generateEVData.m`) 的脚本。

- **`1initialize/`**: 初始化脚本。
   - 从Excel读取数据并初始化对象结构体 (`initializeACsFromExcel.m`, `initializeEVsFromExcel.m`)。
   - 生成小时级调节信号。

### 2. 核心资源模型

- **`2AC/` (空调模块)**:
   - 包含单体空调的热力学仿真、状态更新、基线功率计算及调节潜力计算 (`calculateACAdjustmentPotentia.m`)。
   - 支持蒙特卡洛模拟 (`run_AC_simulation_MC.m`)。

- **`3EV/` (电动汽车模块)**:
   - 模拟EV的充电行为、SOC变化及响应激励的潜力 (`calculateEVAdjustmentPotentia.m`)。
   - 包含不同充电策略的基线模拟。

### 3. 优化与分析算法

- **`4analysis/`**: 核心优化算法库。
   - **算法**: 包含遗传算法 (GA) (`run_hourly_GA_PTDF.m`) 和贪心算法 (`solve_hourly_dispatch_greedy.m`) 的实现。
   - **功能**: 求解小时级调度、计算SDCI指标、处理非线性约束和目标函数。
   - **电网约束**: 包含基于PTDF（功率传输分布因子）的潮流约束计算。

- **`4.1analysis/`**: VPP层面的高级分析。
   - 多目标优化 (`main_vpp_optimizer_multi_objective.m`) 和结果度量计算。

- **`11cvar/` (风险分析)**:
   - 基于CVaR（条件风险价值）的风险约束优化。
   - 包含IEEE 30节点系统的案例分析 (`case_ieee30.m`)。
   - 处理随机场景生成与鲁棒性比较。

### 4. 辅助分析模块

- **`10spaceRelation/` (空间分析)**:
   - 分析资源在不同功能区域（如居住区、商业区）的空间分布概率 (`analyzeACSpatialCharacteristics.m`)。

- **`5userUncertainties/` (用户不确定性)**:
   - 模拟用户参与度和对激励信号的响应不确定性 (`calculateParticipation.m`)。

- **`6SensAnalyses/` (灵敏度分析)**:
   - 分析调节能力对各种参数变化的灵敏度并绘图。

- **`9predict/` (预测模块)**:
   - 多种预测算法实现：ARIMA, Gompertz, Grey Prediction (GM1.1), Random Forest, SVR。
   - 用于预测节点级别的负荷或设施规模增长。

### 5. 主程序与可视化

- **`7main/`**: 项目的入口脚本和结果绘图。
   - `AC_Error_Analysis_Main.m`: 误差分析主程序。
   - `AC_Result_Plotter.m`: 仿真结果可视化工具。
   - 包含多个场景的运行脚本 (如 `AC_main_1_inc_pi.m`)。

## 环境依赖

- **MATLAB**: 推荐 R2020b 或更高版本。
- **工具箱**:
   - Optimization Toolbox (用于优化求解)
   - Statistics and Machine Learning Toolbox (用于预测和蒙特卡洛模拟)
   - Global Optimization Toolbox (如果使用了特定的GA函数)
   - Parallel Computing Toolbox (可选，用于加速蒙特卡洛模拟)

## 快速开始

1. **数据准备**:
运行 `0input_data` 下的脚本生成基础参数文件，或者确保 Excel 输入文件已存在。
2. **运行单次仿真**:
进入 `7main` 目录，运行 `AC_main_1_inc_pi.m` (或其他 `_main` 结尾的脚本) 来执行一次标准的空调调节潜力仿真。
3. **运行风险分析**:
进入 `11cvar` 目录，运行 `main_scenario_generation_soc.m` 生成场景，随后运行 `run_AC_simulation_MC_soc.m` 进行蒙特卡洛模拟。
4. **查看结果**:
使用 `7main/AC_Result_Plotter.m` 或各个模块自带的 `plot_*.m` 脚本查看生成的图表。

## 核心算法简介

- **蒙特卡洛模拟 (MCS)**: 用于处理用户行为（如EV到达时间、AC设定温度）的随机性，生成大量场景以评估聚合潜力的概率分布。
- **遗传算法 (GA)**: 用于解决非凸、非线性的调度优化问题，特别是在考虑复杂的用户响应模型时。
- **CVaR 风险度量**: 在优化目标中加入 Conditional Value at Risk，以平衡调度收益与极端情况下的违约风险。

## 注意事项

- 代码中包含部分硬编码的路径（如Excel文件读取），请根据实际运行环境修改 `1initialize` 中的文件路径。
- `copyMFiles.m` 是一个辅助工具脚本，用于整理代码文件，非核心逻辑。

*本文档由自动化工具根据代码库内容生成。*