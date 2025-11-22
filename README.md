# 监管能力分析系统

一个基于 MATLAB 的框架，用于监管能力的分析与预测。

## 目录结构

核心代码位于 `10regulatorycapacity` 文件夹，模块简介：

| 模块路径             | 功能                |
|-------------------|-------------------|
| **0input_data**  | 数据存储与预处理       |
| **1initialize**  | 参数和环境初始化       |
| **2AC / 3EV**    | 核心计算模块         |
| **4analysis**    | 基础分析功能         |
| **5userUncertainties** | 用户不确定性建模      |
| **6SensAnalyses** | 灵敏度分析          |
| **7main**         | 主程序入口           |
| **9predict**      | 预测模块            |
| **10spaceRelation** | 空间关系计算         |
| **11stackburg**   | Stackelberg 博弈模型 |

## 快速开始

1. 安装 MATLAB。
2. 进入 `10regulatorycapacity/7main` 目录。
3. 运行主脚本进行分析。