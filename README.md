# Regulatory Capacity Analysis System

基于 MATLAB 开发的监管能力分析与预测框架。

## 目录结构说明

项目核心代码位于 `10regulatorycapacity` 文件夹中，模块划分如下：

| 模块路径 | 功能描述 |
| :--- | :--- |
| **0input_data** | 输入数据存储与预处理 |
| **1initialize** | 系统参数与环境初始化 |
| **2AC / 3EV** | 核心计算组件 |
| **4analysis** | 基础分析模块 |
| **5userUncertainties** | 用户不确定性因素建模 |
| **6SensAnalyses** | 灵敏度分析 |
| **7main** | **主程序入口** |
| **9predict** | 预测模型 |
| **10spaceRelation** | 空间关系计算 |
| **11stackburg** | Stackelberg 博弈模型 |

## 快速开始

1. 确保安装 MATLAB 环境。
2. 进入 `10regulatorycapacity/7main` 目录。
3. 运行相应的主脚本开始分析。