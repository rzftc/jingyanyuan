# 虚拟电厂（AC/EV）调节潜力仿真框架 README

> 本仓库围绕“空调（AC）与电动汽车（EV）”聚合体的调节潜力建模、仿真、可视化与调度优化，包含数据生成、个体建模、互补性指标、敏感性分析、空间特性分析、规模预测，以及应急调度等模块。

---

## 1. 项目结构

```
10regulatorycapacity/
  0input_data/                 # 示例/合成输入数据生成脚本（AC/EV/Excel模板）
  1initialize/                 # 从Excel初始化 AC/EV 结构体
  2AC/                         # AC 个体建模与调节能力计算
  3EV/                         # EV 个体建模与调节能力计算
  4analysis/                   # 互补性指标、优化与调度求解
    old_ver/                   # 旧版目标函数/实验脚本
  5userUncertainties/          # 参与率、电价激励、温度/价格模板
  6SensAnalyses/               # 敏感性分析
  7main/                       # 主要仿真与绘图入口脚本
  9predict/                    # AC/EV 规模与基础设施预测
  10spaceRelation/             # 空间特性（功能区）分析
  combineMFilesToTXT.m         # 合并.m至单个txt
  keepMFilesOnly.m             # 清理非 .m 文件的小工具
```

**关键入口脚本（7main）**：
- `AC_main.m`：仅 AC 的快速仿真与绘图
- `ac_ev_simulation_new.m`：全功能 AC+EV 框架（并行处理、分块/个体结果）
- `ac_ev_sinulation_slow.m`：逐步同步的“慢速”版本（便于状态追踪）
- `main_diff_delt_48.m` / `main_diff_delt_48_updown.m`：不同时间步/价格数组的批量仿真（生成结果供绘图）
- `plot_pch.m`、`plot_diffdt.m`、`plot_diffr.m`、`plot_sense.m`：可视化脚本

---

## 2. 环境依赖

- MATLAB R2021a 或更高版本（推荐）
- 工具箱（按模块需要）：
  - Optimization Toolbox（`intlinprog` 等，调度/混合整数优化）
  - Global Optimization Toolbox（`ga` 遗传算法，跨时域二进制调度）
  - Parallel Computing Toolbox（`parfor` 并行加速）
  - Statistics and Machine Learning Toolbox（如 `tiedrank`，计算秩相关）

> 若仅运行基础仿真与绘图，可不安装全部工具箱；但运行 GA/整数规划/并行部分时需相应工具箱。

---

## 3. 快速开始（Quick Start）

### 3.1 准备输入数据（Excel/合成数据）
- 进入 `10regulatorycapacity/0input_data/`。
- 可用脚本：
  - `generateExampleExcel_real_48.m` / `generateExampleExcel_real_24.m`：生成示例 Excel（48×/24×时间分辨率）。
  - `generateEV_48(num_devices_ev, area_type, selected_models)`：按区域（`residential`/`commercial`）与车型集合（或 `all`）生成 EV 数据与对应 Excel。

**示例调用**：
```matlab
% 生成48步长的示例数据/Excel
%（运行前将当前目录切到 0input_data 或把其加入路径）
% generateExampleExcel_real_48;

% 仅生成 EV 示例（500辆，居民区，使用默认车型集合）
% generateEV_48(500, 'residential');
```

> 生成的 Excel/CSV 模板供 `1initialize/initializeACsFromExcel.m` 与 `initializeEVsFromExcel.m` 读取。若模板文件名与入口脚本默认名不一致，请在入口脚本中手动调整（见下）。

### 3.2 运行入口脚本

**A) 仅 AC（快速）**
```matlab
% 在 MATLAB 当前目录切到 7main
AC_main   % 或在命令行运行 AC_main.m
```
关键参数：`T_total`（仿真总时长，小时）、`dt`（步长，小时）、`base_price`（基准电价）；脚本会从 Excel 初始化 AC 阵列，并在多种激励档位下评估上/下调潜力与曲线。

**B) 全功能 AC + EV（推荐）**
```matlab
% 7main 目录
ac_ev_simulation_new
```
此脚本：
- 逐时（或分块）并行遍历 AC/EV 个体，更新状态，计算上/下调潜力；
- 输出聚合结果（`aggregate_results.mat`）与可选的个体级结果（便于后续优化/绘图）。

**C) 同步/慢速版本**
```matlab
% 7main 目录
ac_ev_sinulation_slow
```
此版本额外保存 `aggregate_results_slow_synced.mat`，便于与逐步状态对齐的后续分析。

**D) 生成绘图所需批量结果**
```matlab
% 7main 目录
main_diff_delt_48        % 或 main_diff_delt_48_updown
```
运行后会保存包含 `results_3D` 结构体的 `.mat` 文件（不同时间步/价格档位）。

### 3.3 可视化（示例）
- `plot_pch.m`：双轴显示“EV总充电功率 & 在线车辆数”的时变曲线。

用法要点：在脚本顶部设置数据文件名（如 `mat_filename = 'dt_5m.mat'`），以及要绘制的价格索引 `p_idx_for_plot`，脚本会自动调整坐标轴并以 400 DPI 导出 PNG。

---

## 4. 核心模型与指标

### 4.1 个体层建模
- **AC**：基线功率 `ACbaseP_single`、状态 `calculateACS_single`、调节潜力 `calculateACAdjustmentPotentia`。
- **EV**：状态更新 `calculateEVS_single`、调节潜力 `calculateEVAdjustmentPotentia` / `*_new`、辅助 `calculateEVABC_single`。

### 4.2 互补性与相关性
- **SDCI（同向互补性指数）**：`calculateSDCI.m`
- **ρ（Spearman 秩相关）**：`calculateSpearmanRho.m`（内部使用 `tiedrank` 处理并列次序）

### 4.3 调度优化（分析/求解）
- **小时内最小成本调度**：`solve_hourly_dispatch.m`（混合整数规划，`intlinprog`）
- **跨时域二进制调度（GA）**：`solve_total_time_dispatch_ga.m`（遗传算法，位串编码）
- **应急指令分解（逐时）**：`11stackburg/dispatch_vpp_emergency_horizon.m`（给定上/下调需求序列与 AC/EV 成本-功率边界，返回逐时分配与系统边际成本）

---

## 5. 空间特性分析（可选）
位于 `10spaceRelation/`，包含：
- `classifyLocationByFunctionalArea.m`、`aggregateSetpointCountsByArea.m` 等工具
- `analyzeACSpatialCharacteristics.m` / `analyzeEVSpatialCharacteristics.m`：按功能区输出设定点/概率分布图
- `generateSpatialAnalysisInputFiles.m`、`loadSpatialAnalysisInputs.m`：从 Excel 组织/加载空间分析输入

---

## 6. 规模与基础设施预测（可选）
位于 `9predict/`，包含：
- `ev_scale_prediction.m` / `ac_scale_prediction.m`：高/低情景下的年度保有量（或规模）预测，支持 Logistic / 回归等方法
- `ev_charging_infra_prediction.m`：充电基础设施的推演
- `predicttest.m`：示例脚本与调试输出

---

## 7. 实用小工具

- **只保留 .m 文件**：`keepMFilesOnly.m`
  - 运行后弹窗选择目录，默认“移入回收站”（Windows），并提供预览/二次确认；支持可选的“直接删除/清空空目录”。
  - **使用方式**：将脚本置于 MATLAB 路径，命令行输入 `keepMFilesOnly`，按指引执行。

- **合并 .m 到单个 txt**：`combineMFilesToTXT.m`
  - 递归收集子目录 `*.m`，写入一个 UTF-8 文本（含相对路径标头）。
  - **示例**：
    ```matlab
    combineMFilesToTXT(pwd, fullfile(pwd, 'all_m_code.txt'))
    ```

---

## 8. 常见问题（FAQ）

- **Excel 文件名不匹配**：入口脚本默认 `AC_template.xlsx` 等，如使用生成器脚本导出的其他命名（如 `AC_template1.xlsx`/`EV_template2.xlsx`），请在入口脚本中修改文件名或将文件重命名到默认名。
- **缺少工具箱**：
  - GA 相关报错 → 安装 Global Optimization Toolbox；
  - `intlinprog` 报错 → 安装 Optimization Toolbox；
  - `parfor` 报错 → 安装 Parallel Computing Toolbox；
  - 秩相关计算异常 → 安装 Statistics and Machine Learning Toolbox。
- **图像未生成**：确认 `plot_*` 脚本中的数据文件路径/名称与前序仿真输出一致；脚本会导出 `*.png` 到当前目录。

---

## 9. 复现实验建议
1. 使用 `0input_data` 生成可复现实例 Excel。
2. 运行 `ac_ev_simulation_new.m`（或 `AC_main.m`）得到聚合/个体结果。
3. 运行 `main_diff_delt_48*.m` 产出 `results_3D`，再用 `plot_pch.m` 等完成可视化。
4. 若需要调度优化：先保存个体上/下调矩阵，再运行 `solve_hourly_dispatch.m`（分钟级）或 `solve_total_time_dispatch_ga.m`（跨时域）。
5. 进行空间/规模预测可选步骤，视项目需要启用。

---

## 10. 备注
- 仓库未包含统一 License 文件；如需开源发布或对外共享，请补充许可与版权说明。
- 若需要 Docker/CI 或学术引用模板，可在 issue 中提出，本 README 可继续扩展。
