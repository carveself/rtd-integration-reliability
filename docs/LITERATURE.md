# 文献与创新边界

本轮进一步阅读了 Gao 等作者稿的 II.2 节，而非仅依赖摘要：

- Gao et al., *Efficient self-consistent quantum transport simulator for quantum devices*, Journal of Applied Physics 115, 133707 (2014), DOI 10.1063/1.4870288。[作者稿全文](https://arxiv.org/html/1403.7564)，[原始投稿记录](https://arxiv.org/abs/1403.7564)。该节明确将封闭系统本征能量用于自适应能量网格，并采用最小/最大间距及几何增长控制。故本征能量引导不是本项目的新贡献。
- Baumgartner et al., *Adaptive Energy Integration of Non-Equilibrium Green's Functions*, NSTI-Nanotech 2007, vol.3, pp.145–148。[作者机构原文](https://www.iue.tuwien.ac.at/pdf/ib_2007/hashed_links/5d4PPhrpGwMLGJjaQJrCY_us.pdf)。上一阶段已阅读，表明 NEGF 自适应积分已有长期研究。
- Agarwal et al., *Adaptive quadrature for sharply spiked integrands*, Journal of Computational Electronics 9, 252–255 (2010), DOI 10.1007/s10825-010-0338-3。[作者机构记录](https://docs.lib.purdue.edu/nanopub/722/)。仍未取得并完整核验算法全文；不能声称本项目胜过它或已忠实复刻。

本轮实现的是实部接触边界矩阵铅笔、Schur 补复数极点提示及共用 Gauss 求积，不是 Gao 等的完整 CBR–Poisson 求解器。`spectral` 只能称为文献启发的强基线。新增极点定位也有广泛的散射理论背景；本轮未完成极点自适应方法的系统新颖性检索，因此不作首创声明。

当前结果支持写一份透明的数值实验报告；不足以支持“首次提出共振引导”“认证误差界”“优于现有最佳算法”“高保真器件预测”等投稿主张。若继续争取方法论文，至少需要完整复现适当的公开先例、证明新增步骤的增量价值，再使用未参与开发的新问题验证。
