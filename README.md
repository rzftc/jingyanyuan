# 代码库结构与功能概览（虚拟电厂 AC/EV 调节潜力仿真框架）

> 本概览文档面向**快速上手与维护**：梳理目录结构、核心脚本职责、典型运行路径与产出物。代码主要使用 **MATLAB (.m)**。

---

## 一、顶层结构

```
jingyanyuan-main/
├─ README.md                         # 仓库内原始说明（详细背景/公式/示例）
└─ 10regulatorycapacity/            # 主体代码（数据→建模→仿真→分析→预测→调度）
```

---

## 二、核心目录与职责

> 下列各目录互相衔接，建议按“**数据输入 → 初始化 → 个体模型 → 主仿真 → 分析/可视化 → 规模预测/空间分析/应急调度**”的链条使用。

### 1) `10regulatorycapacity/0input_data` —— 输入数据生成/示例
- **用途**：生成示例 Excel/CSV 输入，或按 24/48 步长构造 EV 负荷样本。
- **代表脚本**
  - `generateExampleExcel_real_24.m` / `generateExampleExcel_real_48.m`：生成 **真实结构**示例 Excel。
  - `generateEVData.m` / `generateEVData_48.m` / `generateEV_48.m`：按指定粒度生成/整理 EV 数据。

### 2) `10regulatorycapacity/1initialize` —— 初始化装置清单（AC/EV）
- **用途**：从 Excel/CSV 读入原始清单与参数，构造仿真所需结构体。
- **代表脚本**
  - `initializeACsFromExcel.m`：初始化空调（AC）对象列表与参数。
  - `initializeEVsFromExcel.m`：初始化电动汽车（EV）对象列表与充电参数。

### 3) `10regulatorycapacity/2AC` —— AC 个体模型与可调节潜力
- **用途**：计算单台 AC 的基线、S 指标（灵活度）、三参数 ABC、以及在激励下的功率调整潜力。
- **代表脚本**
  - `ACbaseP_single.m`：单体 AC **基线功率**计算。
  - `calculateACS_single.m` / `calculateACABC_single.m`：单体 **S 指标**与 **ABC** 参数。
  - `calculateACAdjustmentPotentia.m`：AC 的**向上/向下调节潜力**评估。
  - `incentiveTempAC.m`：AC 在**温控/价格激励**下的响应。

### 4) `10regulatorycapacity/3EV` —— EV 个体模型与可调节潜力
- **用途**：描述 EV 充电行为、SOC（含虚拟/修正）与灵活性窗口，评估在不同激励机制下的调节潜力。
- **代表脚本**
  - `EVbaseP_single.m` / `EVbaseP_single_longstep.m` / `EVbaseP_ChargeUntilFull.m`：EV **基线充电曲线**（不同策略/步长）。
  - `calculateEVS_single.m` / `calculateEVABC_single.m`：单体 **S/ABC** 指标。
  - `calculateEVAdjustmentPotentia*.m`：EV **向上/向下可调潜力**（含新版/对比版）。
  - 其它：`copyEVStruct.m`、`plotIncentiveResponse.m` 等辅助/可视化。

### 5) `10regulatorycapacity/4analysis` —— 互补性/跟踪/优化分析
- **用途**：统计学与经济性分析、AC-EV 互补性、聚合功率跟踪与**遗传算法(GA)** 优化调度。
- **代表脚本**
  - `calculateSpearmanRho.m` / `calculateSDCI.m`：互补性指标（相关性 rho、SDCI）。
  - `run_hourly_GA_plus_Greedy.m`、`solve_total_time_dispatch_ga.m`：**GA + 贪心**的调度求解器。
  - `objective_function_ga*.m`、`nonlinear_constraints_ga*.m`、`hourly_ga_fitness_constrained*.m`：GA 目标/约束/适应度。
  - `eco_test*.m`：经济性测试与 PTDF 变体。

### 6) `10regulatorycapacity/5userUncertainties` —— 用户侧不确定性/参与度
- **用途**：外生价格/激励曲线、设备参与概率、向上/向下响应差异等的不确定性建模。
- **代表脚本**
  - `calculateParticipation.m`：**参与度**与可用性估计。
  - `incentiveTempEV.m` / `incentiveTempEV_updown.m` / `priceTemp.m`：**激励-响应**模板与价格场景。

### 7) `10regulatorycapacity/6SensAnalyses` —— 敏感性分析
- **用途**：对关键参数（价格、舒适度参数、SOC 修正系数、窗口门限等）进行扫描与梯度分析。
- **代表脚本**
  - `analyzeSensitivity.m` / `plotSensitivityCurves.m`：批量扫描+绘图。
  - `computeGradients.m`：对目标/约束的参数敏感度估计。

### 8) `10regulatorycapacity/7main` —— 主仿真入口（聚合 AC/EV）
- **用途**：按**不同离散步长/求解策略**进行 AC+EV 聚合仿真与绘图。
- **代表脚本**
  - `AC_main.m`：仅 AC 的主流程示例。
  - `ac_ev_simulation_*.m`：AC+EV **联合仿真主程序**（block / improve / new / slow 等版本）。
  - `main_diff_delt_48*.m`：对比 **不同时间步长**（48步/日）策略与精度。
  - `plot_*.m`：仿真结果可视化（充电功率、跟踪误差、r/dt 对比等）。

### 9) `10regulatorycapacity/9predict` —— 规模/基础设施预测
- **用途**：给出 AC/EV **装机规模、充电桩基础设施、节点级扩展**等预测。
- **代表脚本**
  - `ev_scale_prediction.m` / `ac_scale_prediction.m`：装机规模预测。
  - `ev_charging_infra_prediction.m`：充电基础设施规模化需求。
  - `predict_node_level.m`：**节点级**预测。

### 10) `10regulatorycapacity/10spaceRelation` —— 空间特性与功能区分析
- **用途**：统计 **不同功能区** 的设定温度/充电分布，分析 AC/EV 的空间异质性，生成空间分析输入。
- **代表脚本**
  - `classifyLocationByFunctionalArea.m`：根据规则**划分功能区**（如住宅/办公/商业等）。
  - `aggregateSetpointCountsByArea.m` / `calculateSetpointProbabilitiesByArea.m`：**设定点统计/概率**。
  - `analyzeACSpatialCharacteristics.m` / `analyzeEVSpatialCharacteristics.m`：AC/EV 的**空间特性**分析。
  - `generateSpatialAnalysisInputFiles.m` / `loadSpatialAnalysisInputs.m`：空间分析输入/读取。

### 11) `10regulatorycapacity/11stackburg` —— 应急调度与边际成本
- **用途**：在应急/紧急窗口下的**虚拟电厂(VPP)出清**与**边际成本**分析。
- **代表脚本**
  - `dispatch_vpp_emergency_horizon.m`：VPP **应急出清**调度求解。
  - `calculate_dispatch_cost_and_mc.m`：**调度成本与边际成本**计算。

---

## 三、典型工作流（建议顺序）

1. **准备输入**（或复现实验）  
   - 使用 `0input_data` 生成示例 Excel（24/48 步长），或将真实数据整理为相同字段。

2. **初始化对象**  
   - 运行 `1initialize/initializeACsFromExcel.m` 与 `initializeEVsFromExcel.m`，得到 AC/EV 结构体数组。

3. **核对个体模型**  
   - 通过 `2AC`、`3EV` 中的基线/指标脚本，抽样验证 **S/ABC/潜力** 计算是否符合业务认知。

4. **聚合主仿真**  
   - 选择 `7main/ac_ev_simulation_*.m` 之一执行 AC+EV 联合仿真；如仅 AC，使用 `AC_main.m`。

5. **分析与可视化**  
   - 调用 `4analysis` 的互补性/跟踪/优化脚本，或 `6SensAnalyses` 对关键参数做敏感性分析。  
   - 如涉**空间维度**，使用 `10spaceRelation` 进行功能区统计与空间特性评估。

6. **规模与应急**  
   - 使用 `9predict` 做规模/基础设施/节点级预测；在**应急场景**下，使用 `11stackburg` 做应急调度与 MC 评估。

---

## 四、输入/输出约定（摘要）

- **输入**：
  - Excel/CSV（设备清单、地理/功能区、价格/激励、出行与到离网信息、容量与 SOC 参数等）。
  - 统一的**时间步长**（常见：24 或 48 步/日），脚本名已标注步长版本。

- **输出**：
  - 聚合功率轨迹、跟踪误差、S/ABC/SDCI/rho 等指标表与图。
  - 预测结果（规模/桩数/节点级）、空间分布统计。
  - 优化/调度结果（GA 出清、应急出清、成本与边际成本）。

---

## 五、运行环境与依赖（建议）

- **MATLAB R2020a+**（实测更高版本兼容性更好）。
- 可能用到的工具箱：
  - **Optimization Toolbox / Global Optimization Toolbox**（GA 相关函数）。
  - **Statistics and Machine Learning Toolbox**（互补性/统计分析）。
- 操作步骤：克隆仓库 → 打开 `10regulatorycapacity` 为工作目录 → 依序运行「典型工作流」。

---

## 六、命名/版本提示

- 带有 `*_48`、`*_24` 的文件表示**时间步长版本**；`*_new`/`*_improve`/`*_slow` 表示不同实现/性能权衡。
- `calculate*` / `EVbaseP*` / `ACbaseP*`：通常为**算法/模型计算**；`plot*` 为可视化；`*test*` 为示例/自测脚本。
- GA/出清相关脚本集中在 `4analysis` 与 `11stackburg`。

---

## 七、快速上手示例（伪代码顺序）

```matlab
% 1) 生成或准备输入
generateExampleExcel_real_48;

% 2) 初始化
acList = initializeACsFromExcel('input_ac.xlsx');
evList = initializeEVsFromExcel('input_ev.xlsx');

% 3) 联合仿真（选择一种主程序）
ac_ev_simulation_new(acList, evList, /* 参数集 */);

% 4) 分析/可视化
rho = calculateSpearmanRho( /* 轨迹 */ );
sdci = calculateSDCI( /* 轨迹 */ );
run_hourly_GA_plus_Greedy( /* 目标/约束/场景 */ );
```

---

### 维护建议
- 新增算法时请在对应子目录新建 `*_new.m` 并保留旧版（便于 A/B 与回归测试）。
- 与数据/步长强相关的脚本在文件名中**显式标注**（如 `_48`），避免误用。
- 在 README/本概览中同步**入口脚本**与**示例数据字段**说明。

—— 完 ——
