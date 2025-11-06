# 代码库结构与功能概览（虚拟电厂 AC/EV 调节潜力仿真框架）

> 本概览文档面向**快速上手与维护**：梳理目录结构、核心脚本职责、典型运行路径与产出物。代码主要使用 **MATLAB (.m)**，并包含少量 **Mathematica (.wl)** 脚本用于特定分析。

---

## 一、顶层结构

```
jingyanyuan/
├─ README.md                # 本概览文件
└─ 10regulatorycapacity/     # 主体代码（数据→建模→仿真→分析→预测→调度）
```

---

## 二、核心目录与职责 (`10regulatorycapacity/`)

> 下列各目录互相衔接，建议按“**数据输入 → 初始化 → 个体模型 → 主仿真 → 分析/可视化 → 规模预测/空间分析/应急调度**”的链条使用。

### 1) `0input_data` —— 输入数据生成/示例
- **用途**：生成示例 Excel/CSV 输入，或按 24/48 步长构造 EV 负荷样本。
- **代表脚本**: `generateExampleExcel_real_24.m`, `generateEVData.m`

### 2) `1initialize` —— 初始化装置清单（AC/EV）
- **用途**：从 Excel/CSV 读入原始清单与参数，构造仿真所需结构体。
- **代表脚本**: `initializeACsFromExcel.m`, `initializeEVsFromExcel.m`

### 3) `2AC` —— AC 个体模型与可调节潜力
- **用途**：计算单台 AC 的基线、S 指标（灵活度）、三参数 ABC、以及在激励下的功率调整潜力。
- **代表脚本**: `ACbaseP_single.m`, `calculateACS_single.m`, `calculateACAdjustmentPotentia.m`

### 4) `3EV` —— EV 个体模型与可调节潜力
- **用途**：描述 EV 充电行为、SOC（含虚拟/修正）与灵活性窗口，评估在不同激励机制下的调节潜力。
- **代表脚本**: `EVbaseP_single.m`, `calculateEVS_single.m`, `calculateEVAdjustmentPotentia*.m`

### 5) `4analysis` & `4.1analysis` —— 互补性/跟踪/优化分析
- **用途**：进行统计学与经济性分析、AC-EV 互补性、聚合功率跟踪与**遗传算法(GA)** 优化调度。`4.1analysis` 包含额外的分析脚本或变体。
- **代表脚本**: `calculateSpearmanRho.m`, `run_hourly_GA_plus_Greedy.m`, `objective_function_ga*.m`

### 6) `5userUncertainties` —— 用户侧不确定性/参与度
- **用途**：外生价格/激励曲线、设备参与概率、向上/向下响应差异等的不确定性建模。
- **代表脚本**: `calculateParticipation.m`, `incentiveTempEV.m`

### 7) `6SensAnalyses` —— 敏感性分析
- **用途**：对关键参数（价格、舒适度、SOC 修正系数等）进行扫描与梯度分析。
- **代表脚本**: `analyzeSensitivity.m`, `plotSensitivityCurves.m`

### 8) `7main` —— 主仿真入口（聚合 AC/EV）
- **用途**：按**不同离散步长/求解策略**进行 AC+EV 聚合仿真与绘图。
- **代表脚本**: `ac_ev_simulation_*.m` (联合仿真主程序), `AC_main.m` (仅AC), `main_diff_delt_48*.m` (不同步长对比)

### 9) `9predict` —— 规模/基础设施预测
- **用途**：给出 AC/EV **装机规模、充电桩基础设施、节点级扩展**等预测。
- **代表脚本**: `ev_scale_prediction.m`, `predict_node_level.m`

### 10) `10spaceRelation` —— 空间特性与功能区分析
- **用途**：统计 **不同功能区** 的设定温度/充电分布，分析 AC/EV 的空间异质性。
- **代表脚本**: `classifyLocationByFunctionalArea.m`, `analyzeACSpatialCharacteristics.m`

### 11) `11stackburg` —— 应急调度与边际成本
- **用途**：在应急/紧急窗口下的**虚拟电厂(VPP)出清**与**边际成本**分析。
- **代表脚本**: `dispatch_vpp_emergency_horizon.m`, `calculate_dispatch_cost_and_mc.m`

### 辅助脚本
- `copyMFiles.m`: 一个工具脚本，用于按需复制或整理目录中的 `.m` 文件。

---

## 三、典型工作流（建议顺序）

1.  **准备输入** (`0input_data`)：生成或准备 Excel 输入数据。
2.  **初始化对象** (`1initialize`)：运行 `initialize*` 脚本，加载设备列表。
3.  **核对个体模型** (`2AC`, `3EV`)：抽样验证单个设备的行为模型。
4.  **聚合主仿真** (`7main`)：选择一个 `ac_ev_simulation_*.m` 脚本执行联合仿真。
5.  **分析与可视化** (`4analysis`, `6SensAnalyses`, `10spaceRelation`)：进行互补性、敏感性或空间特性分析。
6.  **规模与应急** (`9predict`, `11stackburg`)：执行规模预测或应急调度仿真。

---

## 四、输入/输出约定（摘要）

-   **输入**：Excel/CSV 文件，包含设备清单、地理/功能区、价格/激励、出行与到离网信息等。注意脚本名中标注的**时间步长**（如 `_24` 或 `_48`）。
-   **输出**：聚合功率轨迹、跟踪误差、各类指标（S/ABC/SDCI/rho）、预测结果、优化调度方案等图表和数据。

---

## 五、运行环境与依赖

-   **MATLAB R2020a+**
-   **Optimization Toolbox / Global Optimization Toolbox** (用于遗传算法)
-   **Statistics and Machine Learning Toolbox** (用于统计分析)

---

## 六、快速上手示例（伪代码）

```matlab
% 1) 生成或准备输入 (目录: 0input_data)
generateExampleExcel_real_48;

% 2) 初始化 (目录: 1initialize)
acList = initializeACsFromExcel('input_ac.xlsx');
evList = initializeEVsFromExcel('input_ev.xlsx');

% 3) 联合仿真 (目录: 7main)
ac_ev_simulation_new(acList, evList, /* ...其他参数... */);

% 4) 分析 (目录: 4analysis)
rho = calculateSpearmanRho(aggregatedPowerAC, aggregatedPowerEV);
run_hourly_GA_plus_Greedy(/* ...目标/约束/场景... */);
```

---

### 维护建议
- 新增算法时，请在对应子目录新建 `*_new.m` 并保留旧版，便于回归测试。
- 与数据/步长强相关的脚本，在文件名中**显式标注**（如 `_48`），避免误用。
- 更新主要功能或目录结构后，请同步更新本文档。

—— 完 ——