function cfg = CT_toroidal_defaults()
%CT_TOROIDAL_DEFAULTS   CT_toroidal 的默认输入配置
%=========================================================================================
% 【请注意】本文件必须与 CT_toroidal.m 放在同一目录(或把该目录加入 MATLAB 路径)。
% MATLAB 的局部函数在文件外不可见, 因此默认配置必须是独立的 .m 文件。
%
% 用法
%   cfg = CT_toroidal_defaults();               % 取这套默认值(只有配置, 无剖面)
%   cfg.profile_file = 'profile_example.txt';   % ★ 必填: 指定剖面 txt
%   out = CT_toroidal(cfg);                     % 跑主程序
%   (直接 >> CT_toroidal 会报错: 未提供剖面文件 —— 这是有意的)
%
% =============================== 参数总表 ===============================
%
% ---------- 0. 剖面来源（**唯一**: 标准 txt 文件） ----------
%   cfg.profile_file   ''    ★ 剖面文件名(相对当前目录或程序目录, 也可绝对路径)
%   本程序**只接受标准 txt 剖面文件**, 格式为每行 "名字 = 数值"('# ' 之后为注释):
%       必需: rho,R0,B0,a,q,shear,A_ion,theta,n_e,T_e  + 四项 R/L_X(或物理梯度)
%       可选: r  kappa  n_i  T_i
%   结构体 / cfg 字段两种旧方式**已删除**; 若在 cfg 中写剖面字段会直接报错。
%   完整字段表、格式规则与全部报错情形见 ct_read_profile.m 头部注释;
%   生成剖面文件可用 ct_write_profile.m; 示例见 profile_example.txt。
%
% ---------- A. 波数扫描 ----------
%   kmin   0.0    k_y rho_s   扫描下限(k=0 为退化点, 自动跳过)
%   kmax   5.0    k_y rho_s   扫描上限   ★ 用户要求的 0-5 范围
%   nk     26     -           扫描点数(步长 0.2)
%   k_unit 'rho_s' -          波数单位: 'rho_s'(k_y rho_s, 默认) | 'rho_i'(k_theta rho_i)
%                             换算 k_theta rho_i = k_y rho_s * sqrt(2 T_i/T_e)
%
% ---------- C. 物理模型 ----------%   disp_model 'kinetic' -  色散关系模型
%                          'kinetic' = **环形推广**, 完整二维速度空间数值积分(默认)
%                          'pole'    = **单极点约化**, 略去 k∥v∥, 积分塌缩为单极点
%   ct_model   'both'    -  【未使用/死配置】早期版本的 C_T 提取开关, **代码不再读取**
%                          (自检 [T3d] 与核查报告 P7 已确认; 保留仅为兼容旧 cfg)
%   passing_e  'kinetic' -  通行电子: 'kinetic'(含平行朗道阻尼, Z 函数闭式) | 'adiabatic'
%   ft_model   'eps'     -  捕获份额: 'eps' = sqrt(2eps/(1+eps)) | 'fixed'
%   ft_fixed   0.5       -  ft_model='fixed' 时使用
%   kpar_model 'shear'   -  'shear': k∥ = k_theta*s_hat*theta ; 'qR': k∥ = kpar_factor/(q R0)
%   kpar_factor 1.0      -  kpar_model='qR' 时的因子(建议 0.3-1.0)
%   trap_geom  'local'   -  捕获电子进动几何因子 g_e: 'local'(= g_i) 或直接给数值
%   fbw_model  'banana'  -  有限香蕉宽度稳定化: 'banana'(默认) | 'banana_avg' | 'off'
%                           'banana'     J_0²(k⊥rho_b) 取热香蕉宽度; 在 k⊥rho_b = 2.405
%                                        处恰为零(关掉捕获电子响应) => TEM 支在
%                                        k_y rho_s ≈ 0.55 后消失。与 weiland_CT 参考表一致。
%                           'banana_avg' 对热速度分布取平均 ⟨J_0²⟩ = e^{-b}I_0(b),
%                                        b = (k⊥rho_b)²/2: 恒正、单调衰减, 无假零点
%                                        => TEM 支延伸到 k_y rho_s ≈ 0.85(见程序说明 §6.7)
%                           'off'        关闭(高 k 的 gamma 会无界增长, 仅作对比用)
%   fbw_factor 1.0       -  rho_b = fbw_factor*q*rho_i/sqrt(eps)
%   nu_eff_model 'zero'  -  捕获电子的**碰撞(去捕获)展宽** [v4]('zero' | 'collisional' | 'tglf' | 'fixed')
%                           ω -> ω + i·nu_eff 只作用于捕获电子通道; 'zero' = 与旧版逐位一致。
%                           'collisional' 用 nu_ei/eps (与 full_gk_ct 同公式):
%                             nu_ei = 2.91e-6*n_e[cm^-3]*Zeff*lnL/Te[eV]^1.5 [s^-1],  lnL 见下
%                           这是与 TGLF 的 detrapping 项对应的耗散: 本剖面 nu_ei/eps ≈ 2.4e5 s^-1
%                           而 |ω_*e|(k_y rho_s=0.2) ≈ 8.7e4 s^-1, 故碰撞足以把 TEM 的 γ 压低约 2 倍
%                           (见程序说明 §6.9; 自检 [T19])。
%                           'tglf' 用 TGLF 自己的 detrapping 率 cb1 = 0.163*sqrt(k∥v_te/√2·nu_ei
%                           ·(1+0.82 Zeff)) —— **随 k 变化, 且 k∥=0(θ=0)时恒为 0**(程序会告警)。
%                           与 TGLF 对比时可优先试 'tglf'; 若 θ=0 则它失效, 此时 'collisional'
%                           往往过强(本剖面 0.388 c_s/a 会把根全部杀光), 需用 'fixed' 给中等值。
%   nu_eff_factor 1.0    -  'collisional' 时为 nu_ei/eps 的倍数; 'fixed' 时直接是 nu_eff [s^-1]
%   Zeff 1.0, lnL []     -  库仑对数留空时用 24 - ln(sqrt(n_e[cm^-3])/Te[eV])
%   polar_term 'Gam0'    -  极化项: 'Gam0' = Γ_0(b) (本文档 PDF 约定, 默认)
%                                     | 'J0' = <J_0(b)> (与 weiland_CT 兼容)
%   no_curv    false     -  true = 关闭曲率漂移(仅自检/平板基准用)
%   omega_unit 'cs/a'    -  **ω_r 与 γ 的输出单位**(默认 'cs/a'):
%                             'cs/a'   : ω a / c_s ,  c_s = sqrt(T_e/m_i)  (= ρ_s Ω_ci / a)
%                             'krad/s' : ω / 1e3
%                             'rad/s'  : SI
%                           注: C_T 只依赖频率之比 u_* = |ω_r|/(|ω_D0e||c_D|), 与单位无关;
%                               out 结构体里 omega_r/gamma 始终是 SI(rad/s),%                               omega_r_n/gamma_n 是按 omega_unit 归一化的值。
%
% ---------- D. 数值 ----------
%   N_perp  64     x⊥ 方向 Gauss-Legendre 节点数(**唯一剩下的速度空间网格**)
%                  x∥ 方向已由 ct_xpar.m 解析积分, 因此不再需要 N_par / U_par。
%   N_par   —      【已废弃】旧版 x∥ 方向的节点数; 现在只作为占位保留, 不再被读取
%   U_par   —      【已废弃】旧版 x∥ 积分上限, 同上
%   sigma   —      【已废弃】旧版"朗道解析延拓的路径下移量"; 该做法对二次分母是错的
%                  (等价于文档② §2.5 的"整函数 Z"陷阱), 已删除。D(ω) 现在与它无关,
%                  自检 [T14] 直接检验 σ 无关性。
%   nre     48     求根粗网格(实部)单元数
%   nim     24     求根粗网格(虚部)单元数
%   root_im0 1e-3  搜索盒虚部下界(以 max(|w*_i|,|w*_e|) 为单位)
%   root_im1 4.0   搜索盒虚部上界
%   root_re  5.0   搜索盒实部半宽
%
% ---------- E. C_T 与输出 ----------
%   cD_model 'thermal' -  |c_D| 的取法: 'thermal' = (4/3)|g_e| (由固定能量角平均导出, 默认)
%                                        'fixed'   = cfg.cD_abs (兼容 weiland_CT 的 2.0)
%   cD_abs   2.0       -  cD_model='fixed' 时使用
%
%   ---- [v6] C_T 的两种口径(见 ct_ct_exact.m 与 PDF 式 (3)) ----
%   CT_pole     = 3/2 - u_*                        单极点约化 + 冻结模(旧口径)
%   CT_explicit = -Im[∂𝒥/∂Y]/Im[∂𝒥/∂X]            精确核 + 冻结模
%   CT_exact    = -Im[d𝒥/dY]/Im[d𝒥/dX]            精确核 + 模移动项(PDF 式 (3))
%                 d/dX ≡ ∂/∂X + (∂/∂ω)(dω/dX) 是沿本征值支的总导数
%   三者关系: CT_pole ──换用精确核──> CT_explicit ──打开模移动项──> CT_exact
%
%   ct_exact   true     -  是否计算含模移动项的 CT_exact(需要 disp_model='kinetic')
%                          关闭时逐 k 只算 CT_pole, 程序明显更快
%   ct_species 'auto'   -  CT_exact 用哪个物种的粒子通量:
%                          'auto' = ITG 支(ω_r<0)取离子、TEM 支(ω_r>0)取电子
%                          'i' / 'e' = 强制取离子/电子通道
%   ms_frozen  false    -  true = **冻结模基准**(∂ω/∂g ≡ 0): 此时 CT_exact ≡ CT_explicit,
%                          用于把"单极点约化的误差"和"模移动项"分开(自检 [T22d])
%   ms_step    1e-5     -  求 dω/dX_s、dω/dY_s 的**相对**步长
%                          (h = ms_step·max(|X_s|,|Y_s|,|ω|), 中心差分)
%                          太大 → 可能跳到别的根; 太小 → 差分被舍入吃掉
%   ms_track_tol 0.25   -  不跳支判据: 扰动后 |Δω| ≤ ms_track_tol·max(|ω_*s|,|ω|)
%   xi       2.0       -  谱权重指数 w ∝ (gamma/<k⊥²>)^xi, 用于合成 C_T
%   selftest true      -  是否做自检 [T1]-[T22d] (共 57 项; MATLAB R2023a 实测)
%   plots    true      -  是否出图
%   verbose  true      -  是否逐 k 打印
%=========================================================================================

% ---------------- 0. 剖面来源（**唯一**: 标准 txt 文件） ----------------
% 【重要】本程序只接受一种剖面输入: 标准文本文件(名字 = 数值, # 注释)。
%   结构体 / cfg 字段两种旧方式**已删除**; 若在 cfg 里写剖面字段会直接报错。
%   文件格式、字段表、以及全部报错情形见 ct_read_profile.m 头部注释。
cfg.profile_file = '';          % ★ 必填: 外部剖面文件名, 例 'profile_example.txt'
                                %   为空 -> CT_toroidal 立即报错并给出提示
                                %   不在当前目录时, 程序会自动到本程序目录里找

% ---------------- A. 波数扫描 ----------------
cfg.kmin    = 0.0;
cfg.kmax    = 5.0;              % ★ 0-5 k_y rho_s
cfg.nk      = 26;
cfg.k_unit  = 'rho_s';

% ---------------- C. 物理模型 ----------------
cfg.disp_model = 'kinetic';     % 'kinetic'(环形推广, 默认) | 'pole'(单极点约化)
cfg.ct_model   = 'both';        % 【未使用/死配置】保留仅为兼容旧 cfg
% ---- [v7] 响应函数的分子约定(本次修正的核心; 详见 ../修订说明_CORRECTIONS.md §1) ----
%   'standard'(默认) : 𝒥_s = ⟨W(ω - ω_*s^T)/(ω - k∥v∥ - ω_Ds)⟩
%                      = 标准局域线性 GK 方程的解(Ivanov & Adkins 2023 式 (3)/(12);
%                        空间漂移算子 (v∥b̂+v_d)·∇ 只作用在非绝热部分 h_s 上)。
%   'legacy'         : 分子再多一项 -ω_Ds = v6 及以前的做法(与文档① (2.6) 一致)。
%                      该多出项在 k∥=0、零梯度时把离子通道的贡献精确抹成 0(自检 [T23a] 守门),
%                      并使默认剖面 k=0.2 的根由 (0.0905,0.0206) 变成 (0.1303,0.0740) c_s/a。
%                      仅供与 v6 结果逐位对照, **不可用于物理结论**。
cfg.j_numerator = 'standard';   % 'standard' | 'legacy'
cfg.passing_e  = 'kinetic';     % 'kinetic' | 'adiabatic'
cfg.ft_model   = 'eps';         % 'eps' | 'fixed'
cfg.ft_fixed   = 0.5;
% ---- [v8] 平行波数 k∥ 的物理模型(开关) ----
%   局域模型里 k∥ 决定离子/通行电子的平行朗道共振与通行电子 Z 函数自变量 ζ = ω/(k∥v_te)。
%   θ=0 且 k∥=0 时**完全没有平行阻尼** ⇒ γ 在高 k 端无界增长(与全 GK 参考比值可达 2.4 倍)。
%     'zero'         k∥ = 0                 —— 纯二维; 显式声明用(等价于 v7 在 θ=0 的默认)
%     'balloon'(默认) k∥ = k_θ·ŝ·θ·G        —— 局域气球模关系(θ=0 时仍为 0)
%     'balloon_rms'  k∥ = k_θ·ŝ·θ_rms·G     —— 气球模**包络平均**: θ=0 也有有限 k∥,
%                      推荐用于与全 GK(GENE/GS2/GX/TGLF)对照; 见 cfg.kpar_theta_rms
%     'qR'           k∥ = kpar_factor/(qR)  —— 包络口径 Δk∥ ~ 1/(qR·Δθ) 的等价写法
%                      (kpar_factor ≈ 1/Δθ, 与 k_θ 无关; Δθ = 包络在位相角上的宽度)
%     'shear_legacy' k∥ = k_⊥·ŝ·θ           —— v7 及以前的形式(相当于 G = 1), 量级偏大
%                      1–2 个数量级, 任何 θ≠0 的算例都被过度平行阻尼; 仅作回归对照。
%   几何因子 G(auto = ρ_i/(qR)): 三种口径的推导与 CBC 实测标定见 ct_model.m 文件头
%   (局部 (m,n) 口径给 r/(qR); 归一化单位口径给 ρ_i/(qR); 包络口径给 1/(qR·Δθ))。
cfg.kpar_model  = 'balloon';    % 'zero' | 'balloon' | 'balloon_rms' | 'qR' | 'shear_legacy'
cfg.kpar_factor = 1.0;          % 'qR' 时的倍数(≈ 1/Δθ, Δθ = 气球模包络角宽)
cfg.kpar_geom   = 'auto';       % 几何因子 G: 'auto'(= ρ_i/(qR)) | 'r'(= r/(qR)) | 数值
cfg.kpar_theta_rms = 1.0;       % 'balloon_rms' 的等效包络角宽 [rad](CBC 标定建议 1–3)
cfg.kpar_sign   = +1;           % k∥ 的符号 ±1(影响朗道相位/分支; 默认 +1)
cfg.kpar_preset = '';           % 一键预设: ''(不启用) | 'cbc_benchmark' | 'none'
%   'cbc_benchmark' ≡ kpar_model='balloon_rms' + kpar_theta_rms=1.0 + kpar_geom='auto'
%                   (CBC 标定最优: 平均 γ/GX = 1.015, 见 benchmark_cbc/kpar_models.md)
%   'none'          ≡ kpar_model='zero' (显式声明"无平行动力学")
%   ★ 优先级(v8.2): **用户显式给出的 kpar_model/theta_rms/factor 优先于预设** ——
%     预设只填"仍是出厂默认值"的字段, 所以 `preset + kpar_model='qR'` 得到 qR 而不会被静默改掉。
cfg.trap_geom  = 'local';       % 'local' | 数值

% ---- [v8.4] s-α 几何: 压力梯度参数 α = −q²R dβ/dr (β = 8πp/B²) ----
%   α 是 Shafranov 位移/压力梯度参数(Kim & Han 2017 式 (5) 后; Ma et al. 2023 用它解释
%   捕获电子进动方向反转)。它进入三处:
%     ① 磁漂移极向结构 g(θ) = cosθ + sinθ(ŝθ − α sinθ)        → cfg.lam(G) 与各响应
%     ② 垂直波数 WKB 因子 k⊥² = k_y²[1 + (ŝθ − α sinθ)²]      → b、J0、FBW、k∥
%     ③ 平行波数 k∥ = k_y(ŝθ − α sinθ)/(qR)                   → 'balloon'/'balloon_rms'
%   ★ 三者都在 θ=0 处退化为无 α 的形式 ⇒ 在**默认 θ=0** 的单 θ 模型里, α 主要通过
%     **捕获电子的弹跳平均**(沿整条轨道积分)起作用: B(ξ₀) 里多出 −α⟨sin²θ⟩_b 一项,
%     α>0 使捕获电子进动减小、足够大时**方向反转**(Ma 2023 的共振机理)。
%   典型量级: α ≈ 0.1–3(β ~ 1%、R/L_p ~ 数)。设 α = 0 且 trap_precession='legacy' 时
%   与 v8.3 逐位相同。
cfg.alpha      = 0;             % 数值; 0 = 关闭 s-α 效应(θ≠0 时 g 仍含 ŝθsinθ)
cfg.trap_precession = 'auto';   % 'auto'(α≠0 ⇒ 'full') | 'legacy'(仅 cosθ,cos²θ) | 'full'
cfg.fbw_model  = 'banana';      % 'banana' | 'banana_avg' | 'off'
cfg.fbw_factor = 1.0;
cfg.polar_term = 'Gam0';        % 'Gam0'(PDF 约定) | 'J0'
cfg.no_curv    = false;
cfg.omega_unit = 'cs/a';        % ★ 频率输出单位: 'cs/a'(默认) | 'krad/s' | 'rad/s'

% ---- 碰撞(去捕获)展宽 TEM 的进动共振 [v4] ----
%   ω -> ω + i·nu_eff 只作用于捕获电子通道。这是与 TGLF 的 detrapping 项、
%   以及工作区 full_gk_ct 对应的耗散机制; 原版 nu=0, 使 TEM 的 γ 系统性偏大。
%     'zero'(默认)  nu_eff = 0            —— 与旧版逐位一致(自检 [T19a])
%     'collisional' nu_eff = factor*nu_ei/eps,  nu_ei = 2.91e-6*n_e[cm^-3]*Zeff*lnL/Te[eV]^1.5
%     'fixed'       nu_eff = nu_eff_factor  [s^-1]
cfg.nu_eff_model  = 'zero';     % 'zero' | 'collisional' | 'fixed'
cfg.nu_eff_factor = 1.0;        % 'collisional' 时为 nu_ei/eps 的倍数; 'fixed' 时为绝对值 [s^-1]
cfg.Zeff          = 1.0;        % 有效电荷数(影响 nu_ei)
cfg.lnL           = [];         % 库仑对数; 留空则用 24 - ln(sqrt(n_e[cm^-3])/Te[eV])

% ---- 捕获电子通道的模型 [v5] ----
%   'kinetic'(默认) 旧行为: 全 3D 麦氏速度积分 + **单一尖锐进动共振** + J0²(k⊥rho_b) 因子。
%                   该项在高 k 端依赖那个假零点, 谱形会被"雕刻"(见程序说明 §6.8)。
%   'bounce'        文档① §9.4/(9.8): 捕获区 (能量 x, 投角 ξ0) 上取平均, 进动频率用
%                   **弹跳平均** ⟨ω_De⟩_b = -ω_D0e x² B(ξ0); B 随陷阱参数 κ 变号
%                   (深捕获 +1, 边际捕获 -1) —— 共振被总体分布抹平, 且**没有假零点**。
%                   ⚠ 但捕获通道因此对 k 近乎**标度不变**(ω_*e、ω_D0e、k∥v_te 都 ∝ k),
%                   于是 γ ∝ k、高 k 无截断 ⇒ 相对 TGLF 会明显**高估** ITG;
%                   故弹跳分支必须配轨道宽度因子 cfg.fbw_trap_model='banana_avg'(默认)。
cfg.trap_model = 'kinetic';     % 'kinetic' | 'bounce'
cfg.N_trap_x   = 32;            % 'bounce' 时的能量节点数
cfg.N_trap_xi  = 24;            % 'bounce' 时的投角(ξ0)节点数

% ---- 弹跳平均通道的**轨道宽度(FOW)因子** [v5.1] ----
% 【为什么需要】文档① (9.8) 的平均里没有 J0²/Bessel 因子, 故 'off' 是 (9.8) 的字面实现;
%   可是弹跳平均消去 k∥v∥ 之后, 捕获通道里所有频率都 ∝ k, 模型对 k 标度不变 ⇒ γ ∝ k。
%   唯一能破坏标度不变性、给出高 k 截断的机制是有限轨道宽度 k⊥ρ_b (ρ_b = q ρ_i/√ε)。
%   ★ 另一个必须加它的理由: **程序内部一致性** —— 'kinetic'(默认)分支一直乘着
%     M.FBW = J0²(k⊥ρ_b)(见 ct_dispersion.m), 弹跳分支若完全不要这个因子, 同一个物理
%     效应在两个通道里一个有一个没有, 谱形会随之差出 3 倍以上(实测 k_yρ_s=3: 0.150 vs 0.05)。
%     'off'         1                        —— (9.8) 字面; 与 v5 逐位一致(会高估 ITG, 运行时告警)
%     'banana'      J0²(k⊥ρ_b)               —— v4 约定, 在 k⊥ρ_b = 2.405 处有**假零点**(谱形被雕刻)
%     'banana_avg'(默认) ⟨J0²⟩ = e^{-b_b} I0(b_b) —— 恒正单调衰减, 无假零点
%     'energy'      Σ_a w_x(a) J0²(k⊥ρ_b x_a) —— 再对能量分布平均(抑制最弱, 高 k 仍有根)
%   仅作用于捕获通道; 通行电子仍用 cfg.fbw_model。
%   【与 TGLF 的对照】本剖面(ρ=0.45)取 'bounce'+'banana_avg' 时,
%     γ/(k_yρ_s) = 0.027(k=0.4) / 0.031(k=1.0) / 0.028(k=1.2) / 0.005(k=2.0),
%     与 TGLF 蓝线 0.027/0.032/0.019/0.002 在 k≈0.4-1.0 上一致到 2%, 高 k 也如期衰减。
cfg.fbw_trap_model = 'banana_avg';  % 'off' | 'banana' | 'banana_avg' | 'energy'

% ---------------- D. 数值 ----------------
cfg.N_perp    = 64;             % x⊥ 节点数(唯一的速度空间网格; x∥ 已解析积分)
cfg.N_par     = 64;             % 【已废弃】旧版 x∥ 节点数, 保留仅为兼容旧配置
cfg.U_par     = 5.5;            % 【已废弃】旧版 x∥ 积分上限, 保留仅为兼容旧配置
cfg.sigma     = 0.6;            % 【已废弃】旧版 σ 下移量, 已不再使用(自检 [T14])

% ---- [v8] 色散关系根搜索: 预置精度档 root_grid ----
%   v7 固定 48×24 粗网格 + 30 个 |D| 最小点, 在 Cyclone base case 上**漏掉 15 个 k 点中的 11 个**
%   (benchmark_cbc/README.md §5.1)。v8 增加: 自适应加密(refine 级)、多起点 Newton(nseed)、
%   边界绕数触发自动扩盒(expand)、可选的负虚部(阻尼根)搜索、以及诊断统计(R(1).diag)。
%     'fast'     48×24, refine=1, nseed=16  —— 只用于快速浏览(仍比 v7 稳, 但会漏边际模)
%     'standard' 64×32, refine=2, nseed=24  —— **默认**; 已通过 CBC 15/15 不漏根验证
%     'high'     96×48, refine=3, nseed=40  —— 边际模/近简并根; 成本约 standard 的 2.5 倍
%     'ultra'   128×64, refine=3, nseed=64  —— 交付级验证用; 成本约 standard 的 4.5 倍
%   也可直接给 cfg.nre/nim/root_refine/root_nseed 覆盖(见下)。
cfg.root_grid    = 'standard';  % 'fast' | 'standard' | 'high' | 'ultra'
cfg.nre          = [];          % 粗网格实部单元数(留空 = 由 root_grid 决定)
cfg.nim          = [];          % 粗网格虚部单元数(留空 = 由 root_grid 决定)
cfg.root_refine  = [];          % 自适应加密级数(留空 = 由 root_grid 决定)
cfg.root_nseed   = [];          % |D| 最小点作 Newton 种子的个数(留空 = 由 root_grid 决定)
cfg.root_expand  = 1.6;         % 盒边界出现绕数时自动扩盒的倍数(<=1 关闭)
cfg.root_maxexpand = 2;         % 最多扩盒次数
cfg.root_negim   = false;       % true = 同时搜索下半平面(阻尼根)
cfg.root_imfrac  = 1.0;         % root_negim=true 时下沿 = -imfrac*root_im1*wsc
cfg.root_maxroots = 16;         % 最多返回根数(超出会 warning)
cfg.root_verbose = false;       % 打印每次搜索的诊断(零点单元/边界绕数/扩盒/求值次数)
cfg.root_im0  = 1e-3;
cfg.root_im1  = 4.0;
cfg.root_re   = 5.0;
% ---- [v8.5] 多尺度搜索盒与网格收敛自检 ----
%   动机(实测): α=1.0 的 CBC 变体里 γ=0.041 的**弱低频根**被 ±4ω_*T 的大盒漏掉, 却被
%   ±1.5ω_*T 的小盒抓到 —— 同一个盒的网格不可能同时"覆盖大 |ω|"又"分辨弱根"。
%   因此 root_scales 允许给一串**同心盒的尺度**(以 root_re/root_im 为基准), 依次搜索后把
%   根取并集; root_boxes 允许直接给盒(以 max(|ω_*i|,|ω_*e|) 为单位)。
%     root_scales = []      → 单盒(默认, 与 v8.4 逐位相同; CBC 回归用)
%     root_scales = [0.4 1] → 两个盒: 0.4 倍(细分辨弱根) 与 1 倍(标准覆盖)
%     root_scales = [0.3 1 3] → 三盒(最强鲁棒性, 代价约 3 倍)
cfg.root_scales = [];           % 同心盒尺度列表(空 = 只用 root_re/root_im 的单盒)
cfg.root_boxes  = {};           % 显式盒: {struct('re',[lo hi],'im',[lo hi])} × wsc 单位
cfg.root_converge = false;      % true = 再把网格加密一倍重搜一遍, 统计"新出现的根"(漏根证据)
cfg.root_maxconv  = 1;          % 收敛自检最多做几轮(每轮不得新增根即停)

% ---- [v8.5] 配置的单一入口与"陈旧几何"检查 ----
%   ct_configure.m 做三件事: ① 归一化预设/别名 ② 校验枚举值拼写 ③ 把当前 cfg 与几何指纹
%   G.fp(建 G 时记录的 cfg 字段)比对, 一次性列出**所有**不一致字段并告警(提示重建 G)。
%   quiet_geometry_check = true 时 ct_model **跳过** ③ —— 只给单元测试用(它们故意拿改过的
%   cfg 配同一套几何来测开关行为, 那些"不一致"是被测对象)。
cfg.quiet_geometry_check = false;
% ---- [v8.5] 气球模本征值(ct_balloon.m)在主流程里的挂载点 ----
%   cfg.balloon_k 给出**要做气球模本征值**的 k_yρ_s 列表(空 = 不做, 主流程行为与 v8.4 相同)。
%   本征值本身不进入 out.ITG/out.TEM(那两套记录的字段契约不变), 而是单独打印一张表并列在
%   out.balloon 里。逐 k 的局域根与气球模本征值的差别见 benchmark_cbc/ballooning.md。
cfg.balloon_k    = [];          % 例: [0.1414 0.2828 0.4243]
cfg.balloon_Nth  = 81;          % θ 网格节点数
cfg.balloon_Npar = 12;          % x∥ 节点数
cfg.balloon_thmax= 4*pi;        % θ 域半宽
cfg.balloon_elec = '';          % ''(由 passing_e 决定) | 'adiabatic' | 'kinetic'

% ---------------- E. C_T 与输出 ----------------
cfg.cD_model  = 'thermal';      % 'thermal'(= 4|g_e|/3) | 'fixed'
cfg.cD_abs    = 2.0;

% ---- [v6] 含模移动项的 C_T(PDF 式 (3)) ----
%   逐 k 同时给出 CT_pole(单极点闭式)、CT_explicit(精确核+冻结模)、CT_exact(精确核+模移动项)。
%   ct_exact = false 时跳过 CT_exact 的计算(明显更快), 逐 k 表的该列显示 "—"。
cfg.ct_exact     = true;        % 是否算含模移动项的 C_T(需要 disp_model='kinetic')
cfg.ct_species   = 'auto';      % 'auto'(ITG->离子, TEM->电子) | 'i' | 'e'
cfg.ms_frozen    = false;       % true = 冻结模基准(∂ω/∂g ≡ 0), 此时 CT_exact ≡ CT_explicit
cfg.ms_step      = 1e-5;        % dω/dX_s、dω/dY_s 的相对步长(中心差分)
cfg.ms_track_tol = 0.25;        % 不跳支判据: |Δω| ≤ ms_track_tol·max(|ω_*s|,|ω|)
cfg.ms_verbose   = false;       % 是否逐 k 打印模移动项的诊断(显式/移动/合计)

cfg.xi        = 2.0;
cfg.selftest  = true;
cfg.plots     = true;
cfg.verbose   = true;
end
