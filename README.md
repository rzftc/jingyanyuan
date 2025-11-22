# 虚拟电厂(VPP)调节能力分析与优化系统

这是一个基于 MATLAB 开发的综合仿真与分析框架，旨在评估、预测和优化包含规模化电动汽车（EV）和空调（AC）负荷的虚拟电厂（VPP）的调节能力。

该系统涵盖了从底层设备物理建模、用户激励响应、时序仿真、空间特性分析到上层电网调度优化的全流程。

## 📋 项目简介

本项目主要解决以下问题：
1.  **调节潜力评估**：基于物理模型（AC热力学模型、EV电池模型）计算大规模分布式资源的聚合调节潜力。
2.  **激励响应模拟**：模拟用户在不同电价激励下的参与意愿和负荷弹性。
3.  **优化调度**：在满足电网需求和网络潮流约束（PTDF）的前提下，最小化调度成本并优化资源互补性（SDCI）。
4.  **规模预测**：利用多种统计和机器学习模型预测未来区域内的 EV 和 AC 规模。

## 🛠️ 环境要求

* **MATLAB 版本**: R2021b 或更高版本（推荐）
* **必需工具箱**:
    * Optimization Toolbox (用于 `intlinprog`, `fmincon` 等)
    * Global Optimization Toolbox (用于 `ga` 遗传算法)
    * Statistics and Machine Learning Toolbox (用于 `TreeBagger`, `fitrsvm` 等预测模型)
    * Parallel Computing Toolbox (推荐，用于加速大规模 `parfor` 仿真)
    * Econometrics Toolbox (用于 `arima` 预测)

## 📂 目录结构说明

核心代码位于 `10regulatorycapacity` 文件夹下，各子模块功能如下：

| 目录/文件 | 功能描述 |
| :--- | :--- |
| **`0input_data/`** | **数据生成与存储**。包含生成 EV/AC 初始参数（Excel）的脚本，如 `generateEVData.m`。 |
| **`1initialize/`** | **初始化模块**。负责从 Excel 读取数据并构建 MATLAB 结构体，如 `initializeACsFromExcel.m`。 |
| **`2AC/`** | **空调核心算法**。包含单体 AC 的热力学模型、基线功率计算及 SOC 状态更新逻辑。 |
| **`3EV/`** | **电动汽车核心算法**。包含单体 EV 的充电行为模拟、基线计算及灵活性边界计算。 |
| **`4analysis/`** & **`4.1analysis/`** | **优化与指标计算**。包含分层优化算法（GA+MILP）、网络约束（PTDF）处理及互补性指标（SDCI, Rho）计算。 |
| **`5userUncertainties/`** | **用户行为建模**。计算基于电价的用户参与度概率及响应不确定性。 |
| **`6SensAnalyses/`** | **灵敏度分析**。分析激励价格变化对调节潜力的影响梯度。 |
| **`7main/`** | **仿真主程序入口**。包含大规模设备的时序仿真、分块处理逻辑及绘图脚本。 |
| **`9predict/`** | **规模预测模块**。包含 Bass、Logistic、ARIMA、随机森林、SVR 等多种预测模型。 |
| **`10spaceRelation/`** | **空间分析**。分析负荷的空间分布特性，处理经纬度与功能区映射。 |
| **`11stackburg/`** | **博弈模型**。涉及 Stackelberg 博弈相关的调度计算。 |

## 🚀 快速开始

### 1. 生成数据 (可选)
如果 `0input_data` 中没有 `.xlsx` 数据文件，运行以下脚本生成：
```matlab
cd 0input_data
generateExampleExcel_real_48(200, 100, 0.5); % 生成示例 AC 和 EV 数据
2. 运行调节潜力仿真 (Main Simulation)进入 7main 目录，运行主仿真脚本。该脚本会计算基线负荷和调节潜力上下限。推荐脚本: ac_ev_simulation_block_abselute_hour.m支持分块处理（Chunk Processing），防止内存溢出。使用绝对时间轴（例如 06:00 - 次日 06:00）。自动保存结果到 chunk_results/ 文件夹。3. 运行优化调度 (Optimization)仿真完成后，使用 4.1analysis 中的脚本加载结果并进行电网调度优化。推荐脚本: main_vpp_optimizer_multi_objective.m 或 main_vpp_optimizer.m加载步骤 2 生成的分块数据。执行基于遗传算法（GA）和贪心算法的分层优化。考虑网络潮流约束（PTDF）和互补性指标。4. 运行规模预测 (Prediction)如果需要进行未来的保有量预测：进入 9predict 目录。运行 test2.m 或 predict_node_level.m 进行多模型对比预测（Bass, Logistic, Random Forest 等）。🧩 关键模块详解A. 设备物理模型 (2AC, 3EV)空调 (AC): 采用二阶热阻-热容模型。状态变量为室内温度，通过占空比或变频控制调节功率。关键函数: calculateACABC_single (状态转移矩阵), calculateACAdjustmentPotentia.电动汽车 (EV): 考虑到达时间、离开时间、目标电量和最大充电功率。引入“虚拟 SOC”概念描述调节能力的紧迫度。关键函数: EVbaseP_ChargeUntilFull, calculateEVAdjustmentPotentia_new.B. 分块仿真架构 (7main)为了处理数万甚至数十万台设备，系统采用分块处理机制：读取: 每次从 Excel 读取 chunkSize (如 10000) 行数据。计算: 并行计算 (parfor) 该块内设备的基线和调节潜力。保存: 将该块的计算结果（上调潜力、下调潜力、SOC轨迹）保存为独立的 .mat 文件。聚合: 后续分析脚本 (4.1analysis/load_simulation_data_from_chunks.m) 动态加载并聚合这些文件。C. 优化调度策略 (4analysis)采用了分层优化架构：上层 (GA): 优化每小时参与调度的设备数量 ($n_{AC}, n_{EV}$)。目标是满足电网峰值需求并优化互补性指标。下层 (MILP/Greedy): 在已知参与数量的前提下，选择具体的设备并分配功率。目标函数: 最小化调节成本。约束: 功率平衡、设备物理约束、电网潮流约束 (基于 PTDF 矩阵)。关键指标: SDCI (源荷互补指数) 和 Spearman Rho (相关系数)。D. 预测模型 (9predict)集成了多种预测方法以应对不同数据质量：Bass / Logistic: 适用于长期趋势预测。Grey Prediction (GM1,1): 适用于少数据样本。ARIMA: 时间序列分析。Random Forest / SVR: 引入外部因子（如GDP、政策）的回归预测。📊 结果可视化系统包含丰富的绘图工具，位于 7main 和 4.1analysis：plot_vpp_analysis.m: 绘制优化前后的调度对比图，展示 SDCI 和 Rho 指标的改善。plot_sense.m: 绘制不同激励电价下的调节潜力灵敏度曲线。AC_Result_Plotter.m: 专门用于绘制空调集群的状态演变（SOC、温度、功率）。注意：本代码库中包含大量的 .m 脚本，部分为历史版本或测试脚本（如 old_ver 文件夹）。请优先使用本文档推荐的主程序入口。