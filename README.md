# CT_linear —— 不稳定性模求解（从 CT_pinch_fixed 中**只**抽出的那一部分）

V0.87 本项目是一个CT系数计算中间程序，程序开发重度使用Deepseek Harness。使用DIII-D cyclone base ITG/TEM 进行benchmark。


**只保留"解不稳定性模"所需的代码**（局域线性色散关系的求解 → (ω_r, γ)(k)），
**不含** C_T 后处理（`ct_ct` / `ct_ct_exact` / `ct_mode_deriv` / `ct_setgrad`）、
不含完整版主程序 `CT_toroidal.m` 与出图函数 `ct_plot.m`、也不含依赖 C_T 的
`ct_selftest.m` / `ct_selfcheck.m` / `matcheck.js`。


---

## 1　解的是什么

```text
D(ω) = [Γ₀(b_i) − 𝒥_i]/τ + [Γ₀(b_e) − 𝒥_e] = 0

𝒥_s = ⟨ W_s (ω − ω_*s^T)/(ω − k∥v∥ − ω_Ds) ⟩ ,   W_s = J₀²(k⊥ρ_s x⊥) × Maxwell 权重
ω_Ds = −ω_D0s g_s(θ)(2x∥² + x⊥²) ,   ω_*s^T = ω_*s[1 + η_s(x² − 3/2)]
```

* x⊥ 用 Gauss-Legendre 求积；**x∥ 解析积分**（二次分母因子化 + 部分分式 + 分支正确的 f 函数，`ct_xpar.m` / `ct_fz.m`）。
* 电子通道 = 通行电子（Z 函数闭式，可切 `'adiabatic'`）+ 捕获电子（`'kinetic'` 或文档 §9.4 的 `'bounce'` 弹跳平均）。
* 根用复平面绕数 + 自适应加密 + 多起点 Newton 求（`ct_roots.m`），并支持多尺度盒与解析 ∂D/∂ω。
* 支的判定：**Re ω < 0 → ITG**，**Re ω > 0 → TEM**，各取 γ 最大者。

## 2　文件清单

**求解核心（复制自完整版，未改）**

| 文件 | 作用 |
|---|---|
| `ct_read_profile.m` / `ct_write_profile.m` | 剖面 txt 的唯一读写入口（严格报错） |
| `ct_geometry.m` / `ct_geom_fp.m` / `ct_geom_diff.m` | 与 k 无关的几何/单位/派生量 + 几何指纹 |
| `ct_configure.m` | 配置单一入口：预设归一化 + 枚举校验 + 指纹比对 |
| `CT_toroidal_defaults.m` | 全部 cfg 字段（C_T 相关字段在本文件夹里**不被读取**） |
| `ct_model.m` | 逐 k 模型：k⊥、k∥、ω_*、ω_D0、b_s、Γ₀、J₀²、FBW、ν_eff、弹跳平均表 |
| `ct_dispersion.m` / `ct_chanresp.m` / `ct_channel.m` | 组装 D(ω) 与单物种响应 𝒥_s（含三个偏导数） |
| `ct_xpar.m` / `ct_fz.m` / `ct_zfun.m` | x∥ 解析积分、f 函数、等离子体色散函数 Z |
| `ct_pole_integral.m` / `ct_ellip_ke.m` / `ct_gaussleg.m` | 单极点约化积分、椭圆积分、GL 节点 |
| `ct_trapped_avg.m` | 捕获电子弹跳平均表（含 ŝ、α 矩） |
| `ct_roots.m` / `ct_dcount.m` | 复平面高精度求根 + 求值计数器 |
| `ct_balloon.m` / `ct_theta_avg.m` | θ 结构的两条路线（气球模本征值 / 角度平均，可选） |

**新增脚本（本文件夹独有）**

| 文件 | 作用 |
|---|---|
| `CT_linear.m` | **主求解器**：剖面 → 几何 → 逐 k 求根 → 分 ITG/TEM 支 → 汇总（无 C_T） |
| `plot_linear.m` | 标准四联图：γ(k)、ω_r(k)、γ/\|ω_r\|(k)、k∥(k)；可选叠加 GX 参考 |
| `plot_linear_roots.m` | 单个 k 的**复平面根谱** + 该 k 的根的强弱 |
| `linear_demo.m` | 一键演示：`'cbc'`（默认）/ `'kpar'` / `'alpha'` / `'tem'`，自动出图到 `figs\` |
| `linear_selftest.m` | 只看线性部分的 9 项精简自检（见 §4） |
| `sync_core.cjs` | 从完整版重新同步"求解核心"文件（见 §5） |
| **`derivation_verified.tex` / `.pdf`** | **【最新·权威版】逐步推导＋物理原理＋代码逐条核对**（22 页，XeLaTeX+ctex）：每一步都给 ①物理原理 ②推导 ③`文件.m:行`＋表达式原文的代码核对 ④校验项；含 §10 逐条核对总表（28 行）与 §11 失效边界、附录（自检↔公式、历史修正 9 条） |
| **`推导与代码核对.md` / `.pdf`** | 上者的 Markdown 姊妹版（同构同内容，Unicode 公式，便于检索/diff；Edge 打印版） |
| `公式推导与模型处理.md` / `.html` / `.pdf` | 较早的一版推导文档（12–13 页）：公式与模型处理，含 kinetic vs pole 对照节 |
| `formulas_derivation.tex` / `.pdf` | 上者对应的纯 LaTeX 排版版（15 页） |
| `build_doc.cjs` | Markdown→HTML→PDF 的转换器（`node build_doc.cjs 输入.md 输出.html`，再用 Edge headless 打印） |
| `profiles\profile_cbc.txt` | Cyclone base case 剖面（r/a=0.5, q=1.4, ŝ=0.8, R/L_T=6.9, R/L_n=2.2） |
| `profiles\profile_example.txt` | 包自带示例剖面 |

## 3　怎么用

```matlab
addpath('G:\DSH_plasma\CT_pinch\CT_linear');
cd('G:\DSH_plasma\CT_pinch\CT_linear');

% —— 一键演示（CBC 绝热电子扫描 + 出图到 figs\）——
out = linear_demo('cbc');
linear_demo('kpar');     % 三种 k∥ 模型对比(高 k 端的 γ 截止由平行动力学决定)
linear_demo('alpha');    % s-α 的 α 扫描
linear_demo('tem');      % TEM 算例(离子近乎被动 + 动力学电子)
```



```matlab
% —— 手动用 ——
cfg = CT_toroidal_defaults();
cfg.profile_file = 'profiles\profile_cbc.txt';
cfg.k_unit   = 'rho_i';        % 与 GX 的 ky 定义一致(k_θρ_i)
cfg.kmin = 0.05; cfg.kmax = 0.8; cfg.nk = 16;
cfg.ft_model = 'fixed'; cfg.ft_fixed = 0;   % 绝热电子
cfg.passing_e = 'adiabatic';
cfg.kpar_preset = 'cbc_benchmark';          % k∥ = balloon_rms + θ_rms=1
cfg = ct_configure(cfg);
out = CT_linear(cfg);

plot_linear(out);                       % 四联图
plot_linear(out, struct('gx_ref',true, 'file','figs\my_scan.png'));  % 叠加 GX 参考并保存
plot_linear_roots(out, 8);              % 第 8 个 k 的复平面根谱
plot_linear_roots(out, 0.283, struct('file','figs\spec.png'));       % 按 k_yρ_s 找最近点

disp(out.summary);                      % 汇总表(cell; 可直接写文件)
```

`CT_linear` 的输出结构（供二次开发）：

| 字段 | 内容 |
|---|---|
| `out.k_in` / `out.k` | 用户单位波数 / k_yρ_s |
| `out.ITG(i)` / `out.TEM(i)` | 每个 k 的两支记录：`k, branch, omega_r, gamma, omega_r_n, gamma_n, ratio, kpar, kperp, wstar_i, wstar_e, wD0_e, Dval, src, exists` |
| `out.kpar(k)` | 该 k 的 k∥（来自 `ct_model`，与支无关） |
| `out.roots{k}` | 该 k 搜索盒内的**全部根**（复数，SI）—— 用于根谱图 / 检查是否漏根 |
| `out.nD(k)` / `out.nexpand(k)` / `out.nroots(k)` | 求值次数 / 自动扩盒次数 / 盒内根数 |
| `out.balloon` | 气球模本征值（仅当 `cfg.balloon_k` 非空） |
| `out.summary` | 汇总表（cell） |

坐标轴单位：`plot_linear` 默认按 `out.cfg.k_unit` 取横轴（`'rho_i'` 时画 k_θρ_i），
可用 `opt.kunit` 覆盖；纵轴单位由 `cfg.omega_unit`（默认 `c_s/a`）决定。

## 4　自检（`linear_selftest.m`，9 项）

| 项 | 判据 | 实测（CBC，ky=0.4） |
|---|---|---|
| L1 | Σ_p w_p = 1 且 Σ w_p x⊥² = 1 | 1.0e−14 |
| L2 | den² 族恒等式 ∂M₀/∂ω = −M₀b（数值差分独立复核） | 2.0e−10 |
| L3 | k∥ → 0 的连续性 |D(0) − D(ε)|/|D| | 6.9e−12 |
| L4 | CBC 绝热电子存在 γ>0 且 Re ω<0 的 ITG 根 | γ = 0.18074, ω = −0.26280 |
| L5 | 根的残差 |D(z)| ≈ 0 | 1.6e−15 |
| L6 | `root_grid` 的 `fast` 与 `high` 给同一根 | 相对差 0 |
| L7 | `no_curv` ⇒ g_i = g_e = 0 | 0 |
| L8 | 角度平均的小域极限 = 局域根 | 2.2e−5 |
| L9 | 场矩阵（气球模）在局域根处取极小、Ω 翻倍升 >5 倍 | 3.89e−4 → 5.22e−1（升 1344 倍） |


## 5　与完整版的关系 / 同步

```text
node sync_core.cjs          # 从 ..\CT_pinch_fixed\matlab_toroidal_CT 重新复制"求解核心"文件
node sync_core.cjs --check  # 只比较, 不复制(列出哪些文件与完整版不同)
```

* 新增/修改功能请改**完整版**，然后同步过来（保持"核心代码逐字相同"这条约束）。
* 本文件夹的 `CT_linear.m` / `plot_linear*.m` / `linear_demo.m` / `linear_selftest.m`
  是**本文件夹独有**的脚本，不会被同步覆盖。
* `CT_toroidal_defaults.m` 里仍有 C_T 相关字段（`ct_exact` / `cD_model` / `ms_*` …），
  在本文件夹里它们**不被任何代码读取**，保留只是为了与完整版保持同一份配置定义。

## 6　已知道的限制（与完整版同一套）

* 局域单 θ 模型：k∥ 由"几何因子 + 等效角宽"给出，不是解出来的；θ 结构的影响见
  `ct_balloon.m` / `ct_theta_avg.m` 与完整版的 `benchmark_cbc\ballooning.md`（气球模本征值
  目前**物理根未能可靠收敛**，模块会在筛不出物理根时告警而不是给数）。
* 静电、无 E×B 剪切、无径向模结构、无电磁效应。
* 捕获电子通道默认 `'kinetic'`；与全 GK 对照（尤其 TEM/动力学电子）建议 `trap_model='bounce'`
  + `fbw_trap_model='banana_avg'`（见 `linear_demo('tem')`）。
