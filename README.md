# 🔬 DML-Monte-Carlo: 蒙特卡洛模拟驱动下的双重机器学习模型优选方法 毕业论文   

![R Version](https://img.shields.io/badge/R-4.2+-blue.svg)
![License](https://img.shields.io/badge/License-MIT-green.svg)
![Build](https://img.shields.io/badge/Build-Passing-brightgreen.svg)
![Status](https://img.shields.io/badge/Status-Academic_&_Engineering_Ready-orange.svg)

本开源项目提供了一套从零硬编码（Hard-coded）的**双重机器学习 (Double Machine Learning, DML)** 仿真与评估框架。项目致力于解决高维观测数据因果推断中的“维度诅咒”与“正则化偏误”问题。通过整合 Neyman 正交化与严格的交叉拟合（Cross-fitting）技术，并提出并验证了高度优化的 **Stacking-LM (线性元学习器融合)** 架构，实现了在极端非线性场景下的无偏因果效应估计。

## 📌 项目背景与解决痛点

在传统的计量经济学中，高维混淆变量会导致严重的遗漏变量偏误；而直接将预测型机器学习模型（如深度学习、集成树模型）引入因果推断，又会因为模型极强的拟合能力而引发**正则化偏误 (Regularization Bias)** 与**过拟合偏误 (Overfitting-induced Bias)**。

本项目严格遵循 Chernozhukov et al. (2018) 的经典 DML 理论，通过 R 语言自主构建了包含 5 种数据生成过程 (DGP) 和 5 类底层机器学习算法的超大规模蒙特卡洛仿真（Monte Carlo Simulation）系统，科学界定了不同辅助模型（Nuisance Parameters Estimators）在因果推断中的能力边界。

---

## 🌟 核心架构与算法创新

### 1. Stacking-LM “复杂+简单”融合架构 (最优解)
本项目摒弃了单一基学习器的局限性，自研了定制化的 Stacking 模型融合架构：
* **基学习器层 (复杂特征提取)：** 采用随机森林 (RF)、加入 $L_2$ 正则化的神经网络 (NN)、梯度提升树 (GBM) 和袋装法 (Bagging)，极致捕捉高维非线性混淆机制。
* **元学习器层 (方差约束)：** 采用正则化线性模型 (LM) 作为次级学习器。
* **架构优势：** 通过对照实验 (`lm.R` vs `rf.R`) 科学证明：在 DML 的正交残差阶段，使用简单线性模型做融合能够有效实施二次正则化，防止使用复杂元学习器（如 Stacking-RF）引发的二次过拟合，从而严格保障了 DML 要求的 $N^{-1/4}$ 渐近收敛速率。

### 2. DML1 与 DML2 双估计器底层实现
代码在底层计算因果效应 $\theta$ 时，同时支持两种核心估计逻辑：
* **DML1：** 在交叉拟合的每一折独立计算 $\theta$ 后求均值。
* **DML2：** 全局池化 (Pooling) 所有折的正交残差后进行统一回归。仿真数据证明了 DML2 在有限小样本下具有更小的方差和更优的稳定性。

### 3. 面向极端非线性场景的极限工程防御
为应对深层非线性 DGP（情形5）导致的梯度爆炸或数值发散，系统内嵌了极其严密的工程防线：
* **稳健推断底座：** 在最终估计 $\theta$ 时，优先调用高崩溃点 (Breakdown point) 的 **MM-稳健回归估计 (`MASS::rlm`)** 抑制离群残差的杠杆效应；并严格使用 **HC1 异方差稳健标准误 (`sandwich::vcovHC`)** 构建置信区间，确保统计推断的绝对严谨。
* **数据洗脱与截断：** 对倾向得分机制实施严格边界保护（$m_0(X) \in [0.1, 0.9]$），并基于 IQR 构建自适应分位数盖帽法 (Adaptive Outlier Capping)，从根本上消除了除零错误与残差极值。
* **深层特征标准化：** 针对神经网络模块，强制在每一折（Fold）内部进行 `Center` 与 `Scale` 独立缩放，防止数据泄露并保证激活函数的非饱和状态。

---

## 📊 核心量化发现 (Simulation Results)

基于 $N \in [100, 10000]$ 的大跨度样本与 5 种复杂度递增的 DGP 场景，蒙特卡洛模拟揭示了以下核心结论：

1. **绝对无偏的收敛能力：** 在大样本 ($N=10000$) 场景下，Stacking-LM 架构的绝对偏差 (Absolute Bias) 完美逼近 **0.025**，实现了无偏的因果处理效应估计。
2. **碾压级的方差控制：** 在包含指数与对数混合生成的极端非线性场景（情形 5）中，Stacking-LM 架构有效抵御了噪声干扰，其 **RMSE 较 Stacking-RF 大幅降低了 32%**。
3. **自适应动态权重分配：** Stacking-LM 能够根据不同 DGP 复杂度自动调节基学习器权重。在线性场景中 Bagging 占主导，而在深层非线性场景中，RF 与 GBM 的权重自动提升。

---

### 📈 1. 模型贡献度分析 (Model Contribution)
> **分析洞察：** 如下对比所示，Stacking-LM（线性次级学习器）展现了更稳定的权重分配，有效避免了 Stacking-RF 在复杂场景下由于“黑盒叠加”导致的权重坍塌。

| 权重占比百分比图 | 动态贡献曲线 (LM) | 动态贡献曲线 (RF) |
| :--- | :--- | :--- |
| ![Contribution Percent](https://github.com/user-attachments/assets/51d4f29c-3c74-4ffa-9342-224a61440f35) | ![LM Contribution](https://github.com/user-attachments/assets/19261d6f-5c83-4c0c-8791-34d5507a3e0c) | ![RF Contribution](https://github.com/user-attachments/assets/12239860-b161-4e72-86cc-efb9eca76c20) |

---

### 🎯 2. 估计值的无偏性与稳定性 (Bias & Mean Comparison)
> **分析洞察：** 在全样本量跨度下，Stacking-LM 的估计均值（Mean）更迅速地收敛至真实值 1.0，且偏差（Bias）波动远小于 Stacking-RF。

#### 🔹 偏差 (Bias) 对比
* **左图 (Stacking-LM):** 展现了极佳的收敛斜率。
* **右图 (Stacking-RF):** 在小样本下存在明显的过拟合偏误。

![LM Bias](https://github.com/user-attachments/assets/30d25ba6-2bc5-4a29-ba40-29fec9f22c60) 
![RF Bias](https://github.com/user-attachments/assets/c68637ca-8085-456c-a5e0-af411e33ce64)

#### 🔹 均值 (Mean) 对比
![LM Mean](https://github.com/user-attachments/assets/8ea92b41-31f6-4c02-a601-f9b29f1b5bc7)
![RF Mean](https://github.com/user-attachments/assets/da796cec-3f93-423c-b9bf-8e723e13ec6e)

---

### 📉 3. 预测性能与 RMSE 比较 (Performance/RMSE)
> **分析洞察：** 均方根误差 (RMSE) 是衡量模型优劣的核心指标。实验证明，线性元学习器 (LM) 能够有效起到正则化作用，将非线性误差压制在最低水平。

| Stacking-LM 性能 (最优) | Stacking-RF 性能 (对照) |
| :--- | :--- |
| ![LM RMSE](https://github.com/user-attachments/assets/d4224876-7f49-479b-98ec-74980945200d) | ![RF RMSE](https://github.com/user-attachments/assets/aa360c2f-c7b8-4f82-b445-5286960aba54) |

---

## 📩 获取完整实验数据
本项目生成了超过 20 项学术级资产（包括完整的 Excel 结果表与 HC1 稳健标准误汇总）。如需获取详细原始数据或进行学术交流，请通过以下方式联系作者：

* **📧 Email:** 3035402815@qq.com
* **💬 WeChat:** zc17709337386

---
👤 **Developer:** Chi Zhang | 🏆 **本项目为作者本科优秀毕业设计核心代码库（专业 Top 4）。**




---

## 🛠️ 代码库目录结构

```text
.
├── BAGGAGE+ 5.R        # 基于袋装树 (Bagged Trees) 的 DML 实现基准
├── BOOSTING 5+.R       # 基于梯度提升 (GBM) 的 DML 实现 (含内部交叉验证调参)
├── nn gai.R            # 基于神经网络 (NN) 的 DML 实现 (含严格特征缩放与L2衰减)
├── rf.R                # 基于随机森林 (RF) 的基准对比 (严格控制K-Fold样本隔离)
├── lm.R                # 🌟核心：Stacking-LM 融合架构主程序与可视化系统
└── README.md           # 项目工程文档
```

## 🚀 快速开始与复现

### 环境配置
请确保您的运行环境为 R (>= 4.2.0)。脚本内置了自动依赖检查机制，将自动安装并加载核心算法包（包括 `caret`, `nnet`, `gbm`, `randomForest`, `sandwich`, `MASS`, `ggplot2` 等）。

### 运行完整仿真流水线
要复现最终的 Stacking-LM 仿真并生成学术对比图表，请直接执行核心脚本：

```R
# 运行主融合架构程序
source("lm.R")
```

*注：程序运行后，将自动生成超过 20 项学术级资产，包括长/宽格式的 CSV/Excel 结果表、HC1 稳健标准误汇总表、以及由 `ggplot2` 绘制的 RMSE、Bias 随样本量衰减的高分辨率对比折线图和学习器权重热力图。*

---
👤 **Developer:** Chi Zhang | 🏆 *本项目为作者本科优秀毕业设计核心代码库（最终答辩成绩：专业 Top 4）
# =====================================================================
# 版权所有 (c) 2024-2026 张弛（Chi Zhang）。保留所有权利。
# 标题：蒙特卡洛模拟驱动下的双重机器学习模型优选方法 毕业论文   
# 作者：张弛（北京化工大学）
# 备注：本代码库属于一项专有的毕业论文内容。
# 严禁未经授权的学术出版、商业用途，或对本代码架构及 DGP（数据生成过程）设计的剽窃行为。
# =====================================================================
