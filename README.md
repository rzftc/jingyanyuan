# 虚拟电厂（AC/EV）调节潜力仿真框架 README

> 本仓库围绕“空调（AC）与电动汽车（EV）”聚合体的调节潜力建模、仿真、可视化与调度优化，包含数据生成、个体建模、互补性指标、敏感性分析、空间特性分析、规模预测，以及应急调度等模块。

## 1. 项目概述
本项目提供虚拟电厂（含交流负荷AC与电动汽车EV）调节潜力的仿真框架，支持多时间尺度仿真、聚合控制策略验证及调节潜力分析，可模拟不同激励机制下的EV响应特性与SOC（状态-of-charge）动态变化，新增在线车辆统计与灵活性窗口计算功能。

## 2. 核心功能
- 电动汽车（EV）充电行为与调节潜力仿真
- 多时间尺度（长/短步长）控制策略实现
- 虚拟SOC（原始/修正）计算与动态修正
- 聚合功率跟踪与调节潜力（向上/向下）分析
- 激励响应机制模拟（含参与度计算与灵活性窗口）
- 在线车辆统计与状态分布分析
- 多场景参数敏感性分析（支持激励电价范围扫描）

## 3. 代码结构
### 3.1 主要模块
- `EV/`：电动汽车相关仿真代码
  - `6main/`：主程序入口
    - `main_v5_S_modified.m`：支持修正SOC计算的EV仿真主程序（最新版）
    - `main_v4_100.m`：基于Excel数据接口的完整版主程序
  - `4EVupdate/`：EV状态更新函数
    - `updateLockState.m`：EV闭锁状态（LockON/LockOFF/ON/OFF）动态切换
    - `calculateVirtualSOC_upgrade.m`：优化版虚拟SOC计算（含原始/修正SOC）
    - `calculateVirtualSOC_upgrade_ori.m`：原始版虚拟SOC计算（保留对比用）
  - `999wastecode/`：历史版本代码
    - `main_v2.m`：基础版本主程序（遵循论文模型结构）
    - `main.m`：早期版本主程序
    - `main_acc.m`：初始化向量化版本主程序
    - `dataconvert.m`：数据转换辅助脚本
- `10regulatorycapacity/`：调节潜力分析模块
  - `7main/main_diff_delt_48_updown.m`：全功能虚拟电厂调节潜力分析系统（支持激励电价扫描）
  - 实用工具脚本：
    - `keepMFilesOnly.m`：仅保留指定目录下的.m文件
    - `combineMFilesToTXT.m`：合并所有.m文件到单个TXT
  - `999wastecode/`：
    - `dtw.m`：动态时间规整算法实现
- `0inputdata/`：输入数据模板
  - `residential_all_models.xlsx`：EV参数模板（最新版主程序默认读取）
  - `evdata.xlsx`：早期版本EV参数文件

### 3.2 关键主程序说明
| 程序路径 | 功能描述 | 核心特性 |
|---------|---------|---------|
| `EV/6main/main_v5_S_modified.m` | 修正SOC仿真主程序 | 支持短时间步状态更新、修正SOC计算、个体EV状态跟踪 |
| `EV/6main/main_v4_100.m` | Excel数据接口完整版 | 基于Excel参数输入，支持100辆EV仿真，含基础可视化 |
| `10regulatorycapacity/7main/main_diff_delt_48_updown.m` | 调节潜力分析系统 | 激励电价范围扫描、在线车辆统计、灵活性窗口计算 |
| `EV/999wastecode/main_v2.m` | 基础版本主程序 | 保留论文模型结构，支持基础功率跟踪与SOC分析 |

## 4. 核心算法说明
### 4.1 虚拟SOC计算
#### 4.1.1 原始SOC（`S_original`）
反映实际充电量与期望充电量的偏差：
```
S_original = -(E_actual - E_exp) / (C * r)
```
其中：
- `E_actual`：实际累计充电量
- `E_exp`：期望累计充电量（基于匀速充电模型）
- `C`：电池容量，`r`：调节系数

#### 4.1.2 修正SOC（`S_modified`）
融合原始SOC与充电紧急程度指标，公式为：
```
S_modified = alpha1 * S_original + alpha2 * I_value
```
其中紧急程度指标`I_value`计算：
```
I_value = tanh(kappa*(rho - 0.5) - gamma*(P_current/P_N)*t_elapsed)
```
- `rho`：剩余充电时间比率（`tau_rem / 剩余可用时间`）
- `tau_rem`：剩余充电时间（动态更新）
- `kappa`：tanh函数陡峭度（默认4）
- `gamma`：时间衰减项强度（默认0.05）
- `t_elapsed`：入网后经过时间

### 4.2 状态更新机制（`updateLockState.m`）
EV状态（LockON/LockOFF/ON/OFF）动态切换逻辑：
1. **离网判断**：当前时间超出离网时间（`t > t_dep`）→ 强制LockOFF
2. **电量满足**：实际电量≥目标电量（`E_actual ≥ E_tar`）→ LockOFF
3. **时间不足**：剩余时间无法完成充电需求 → LockON
4. **正常状态**：根据当前功率方向切换子状态（ON充电/OFF空闲）
5. **边界保护**：修正SOC触及[-1,1]边界时强制切换子状态（防止振荡）

### 4.3 激励响应机制
- **参与度计算**：基于激励电价与基准电价的对比动态生成参与概率
- **灵活性窗口**：
  ```
  E_reg_min = max(E_tar_original - deltaE_down, E_in)  % 下调调节边界
  E_reg_max = min(E_tar_original + deltaE_up, C_EV)   % 上调调节边界
  ```
  其中`deltaE_up`/`deltaE_down`由激励电价与功率约束计算得出

### 4.4 多时间尺度控制
- **长步长（`dt_long`）**：聚合控制决策、基准功率计算、λ*优化
- **短步长（`dt_short`）**：EV状态实时更新、功率响应模拟、SOC动态修正

## 5. 使用步骤
1. **数据准备**：
   - 手动准备：使用`0inputdata/residential_all_models.xlsx`填写EV参数
   - 自动生成：运行`generateEVParameters_real(excelFile, 100, 0.6)`生成100辆EV数据
2. **选择主程序**：
   - 基础仿真：`EV/6main/main_v5_S_modified.m`
   - 调节潜力分析：`10regulatorycapacity/7main/main_diff_delt_48_updown.m`
3. **参数配置**：根据需求修改时间参数（`dt_short`/`dt_long`）、激励电价范围等
4. **运行程序**：直接在MATLAB中运行主程序，自动生成仿真结果
5. **结果查看**：程序自动生成可视化图表，结果数据存储于`results`结构体

## 6. 参数说明
| 关键参数 | 说明 | 示例值 |
|---------|------|-------|
| `dt_short` | 短时间步长（分钟） | 0.6分钟 |
| `dt_long` | 长时间步长（分钟） | 60分钟 |
| `t_sim` | 仿真总时长（分钟） | 1440分钟（24小时） |
| `simulation_start_hour` | 仿真开始时间（小时） | 6（早上6点） |
| `simulation_end_hour` | 仿真结束时间（小时） | 30（次日早上6点） |
| `t_adj` | 调节时长（小时） | 1小时 |
| `p_incentive_range` | 激励电价扫描范围（元） | 0~50（步长2） |
| `alpha1/alpha2` | SOC修正权重系数 | 0.8/0.2 |
| `kappa/gamma` | 紧急程度指标参数 | 4/0.05 |

## 7. 实用小工具
- **只保留 .m 文件**：`keepMFilesOnly.m`
  - 功能：保留所选文件夹（含子文件夹）中的.m文件，删除其余所有文件
  - 特性：支持Windows回收站、预览待删除文件、二次确认；可选清理空目录
  - **使用方式**：将脚本置于MATLAB路径，命令行输入`keepMFilesOnly`，按指引执行
  - 可调参数：`forceDelete`（强制删除）、`removeEmptyDirs`（清理空目录）等

- **合并 .m 到单个 txt**：`combineMFilesToTXT.m`
  - 功能：递归收集子目录所有`.m`文件，按相对路径排序后写入UTF-8编码文本
  - 特性：自动过滤隐藏文件夹（Windows）、保留文件结构标头
  - **示例**：
    ```matlab
    combineMFilesToTXT(pwd, fullfile(pwd, 'all_m_code.txt'))
    ```

## 8. 常见问题（FAQ）
- **Excel 文件名不匹配**：最新主程序默认读取`residential_all_models.xlsx`，若使用其他命名需在入口脚本中修改。早期版本可能使用`evdata.xlsx`。
- **缺少工具箱**：
  - GA 相关报错 → 安装 Global Optimization Toolbox；
  - `intlinprog` 报错 → 安装 Optimization Toolbox；
  - `parfor` 报错 → 安装 Parallel Computing Toolbox；
  - 秩相关计算异常 → 安装 Statistics and Machine Learning Toolbox。
- **图像未生成**：确认`plot_*`脚本中数据路径与仿真输出一致，脚本默认导出`*.png`到当前目录。
- **SOC计算异常**：检查`C`（电池容量）和`r`（调节系数）是否为非零值，避免除零错误。

## 9. 可视化输出
仿真程序自动生成关键分析图表，包括：
- 功率跟踪与λ*动态对比（长/短步长数据）
- 聚合SOC与个体SOC（原始/修正）变化分析
- 特定EV（如第10辆）的SOC与当前功率对比
- 调节潜力（向上/向下）随激励电价的变化曲线
- 在线车辆数量统计与状态分布

## 10. 备注
- 仓库未包含统一 License 文件；如需开源发布或对外共享，请补充许可与版权说明。
- 若需要 Docker/CI 配置或学术引用模板，可在 issue 中提出，本 README 可继续扩展。
- 历史版本代码（`999wastecode/`）仅作参考，建议使用`6main/`和`10regulatorycapacity/7main/`下的最新程序。


### 下载说明
目前代码以源代码形式提供，可通过以下方式获取：
1. 克隆代码仓库（若提供版本控制）：
   ```bash
   git clone [仓库URL]
   ```
2. 直接下载代码压缩包（请联系项目维护者获取最新版本）

建议使用 MATLAB 2020b 及以上版本运行，确保已安装所需工具箱（见FAQ部分）。