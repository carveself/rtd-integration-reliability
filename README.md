# rtd-integration-reliability

开放量子输运的能量积分可靠性研究：怎样区分目标达标、误差指标可信和计算量减少。

这是可运行、含原始数据与失败记录的纯仿真研究项目，适合个人电脑复现。当前成果支持数值研究报告，**尚未建立新方法论文的创新性与可靠性主张**。

TCAD（Technology Computer-Aided Design，工艺与器件计算机辅助设计）用于用数值模型研究半导体器件。这里采用 NEGF（Non-Equilibrium Green's Function，非平衡格林函数）计算开放系统输运，以 RTD（Resonant Tunneling Diode，共振隧穿二极管）中的窄透射峰检查能量积分。模型恒质量、预设线性偏压降、相干弹性输运，没有自洽 Poisson equation（泊松方程，用电荷密度计算电势）的反馈，也没有实验标定。

## 最新结果：选择性谱引导

24 个新案例中，选择性方法均达到 0.1% 电流目标，计入粗探测后比普通自适应基线少约 **12.31%** 实轴输运点评估。但旧的 308 项回归仍出现 **3 次误差低估**，最坏低估比为 **16.74**，所以不能作为可靠替代算法发布。

完整谱引导在旧回归中没有观察到低估，但能量点评估为原基线约 **1.91 倍**；极点引导为约 **2.03 倍**。它们是经验检查结果，并非数学认证。

![效率与可靠性对照](phase4/results/efficiency_reliability.png)

[第四阶段方法与结果](phase4/README.md) · [未解决的失败记录](phase4/docs/FAILURE_AUDIT.md) · [第三阶段导师视角自审](docs/ADVISOR_REVIEW.md)

## 安装与复现

在 Python 3.12.14、NumPy 2.3.5、SciPy 1.18.1、Matplotlib 3.11.2、mpmath 1.4.1 实际运行。无需 GPU 或商业 TCAD 软件。

```bash
git clone https://github.com/carveself/rtd-integration-reliability.git
cd rtd-integration-reliability
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
python -W error tests/test_phase3.py
python -W error phase4/tests/test_selective.py
python phase4/run_phase4.py
python phase4/analyze_phase4.py
```

Windows 使用 `.venv\Scripts\activate` 激活环境，并在相应终端设置线程变量。17 项测试已通过。新研究运行会覆盖 `phase4/results` 与 `phase4/results_regression` 的同名文件；建议在副本中复现。只重画现有结果可直接运行 `python phase4/analyze_phase4.py`。

## 研究记录

| 阶段 | 问题与结论 | 入口 |
|---|---|---|
| Pilot | 粗能量网格会改变固定样本变异和有限集设计选择；已有自适应算法是强反证 | `archives/rtd-research-pilot.zip` |
| Phase 2 | 目标驱动细分的净节省有限，且经验误差指标出现 5 次低估 | `archives/rtd-research-phase2.zip` |
| Phase 3 | 谱/极点引导消除这 5 次已知低估，但计算量接近翻倍 | [第三阶段完整入口](PHASE3_README.md) |
| Phase 4 | 选择性引导节省新案例计算，却留下 3 次旧低估 | [第四阶段](phase4/README.md) |

前几阶段原始 ZIP 均保留，历史文档中“尚未上传”等陈述对应当时状态。当前仓库展开提供第三、第四阶段代码和数据；旧阶段 ZIP 可独立解压复现。教学入门项目为 [mini-tcad-negf-rtd](https://github.com/carveself/mini-tcad-negf-rtd)，其原文件没有被覆盖。

## 如何判断结果可信到什么程度

- 候选算法停止后才查看参考值；谱定位、粗采样和复数求根的成本均计入。
- 独立带状求解、QUADPACK、细网格求积及局部任意精度计算交叉检查参考值，原警告与失败均保留。
- 新案例与旧回归分开，后续压力测试也单列。固定高度扰动是合成情景，不是实测制造分布。
- CV（Coefficient of Variation，变异系数）研究使用固定样本；局部误差覆盖与统计/设计目标达标不能混为一谈。
- 计时为单机器短重复的描述性结果，不声称普遍加速、完整漏峰保证或真实器件预测精度。

算法已有相关先例，见 [文献审计](docs/LITERATURE.md) 与 [第四阶段说明](phase4/README.md)。本项目由工具辅助制作，复现者应明确区分参考实现和自己独立完成的推导、修改及实验。导师视角文档是内部模拟审阅，不代表真实导师意见。

GitHub 网页上传提交按文件组保留说明；本地开发提交记录见 `COMMIT_HISTORY.txt`。文件校验清单为 `MANIFEST.sha256`。
