function [D, dDdw] = ct_dispersion(M, omega, disp_model, want_deriv)
%CT_DISPERSION  局域线性色散函数 D(ω), 两种模型
%=========================================================================================
% 准中性(n_i = n_e):   D(ω) = [1 - H_e] + (1/τ)[1 - H_i] = 0 ,  τ = T_i/T_e
% 其中  1 - H_s = Γ_0(b_s) - J_s ,  Γ_0(b) = e^{-b} I_0(b) = <J_0²(k⊥ρ_s)> ,  b_s = k⊥²ρ_s²/2
%
% ---------------------------------------------------------------------------------------
% disp_model = 'kinetic'  (默认) —— **环形推广**: 保留完整速度空间积分
%   J_s = < J_0²(k⊥ρ_s) (ω - ω_Ds - ω_*s^T) / (ω - k∥ v_t x∥ - ω_Ds) >
%   ω_*s^T = ω_*s [1 + η_s (x² - 3/2)] ,   ω_Ds = -ω_D0s (2x∥² + x⊥²) g_s
%
%   离子 : 全动理学(平行朗道共振 + 曲率/∇B 漂移共振)
%   电子 : (1-f_t) 通行电子 + f_t 弹跳平均捕获电子(只有进动漂移共振, k∥ = 0)
%          通行电子可切换 'kinetic'(含平行朗道阻尼, Z 函数闭式) / 'adiabatic'(玻尔兹曼)
%          有限香蕉宽度 J_0²(k⊥ρ_b) 作用在捕获电子通道上(高 k 端的主要稳定机制)
%
%   **碰撞(去捕获)[v4]**: 捕获电子通道的进动共振带一个碰撞宽度,
%      ω -> ω + i·ν_eff,  ν_eff 由 cfg.nu_eff_model 决定('zero' 默认/'fixed'/'collisional')。
%      这是与 TGLF(矩矩阵里的 detrapping 项)以及工作区 full_gk_ct 对应的耗散机制;
%      原版完全没有碰撞, 使 TEM 的 γ 系统性偏大(见程序说明 §6.9)。
%
%   **二次分母的处理(关键)**: 分母含 ω_Ds ∝ x∥², 对 x∥ 是二次的。本模型按
%   《DGH问题_...详解.md》§2 的精确路线: 每个 x⊥ 上把 x∥ 积分解析做出
%   (部分分式 + 分支正确的 f 函数 = ct_xpar.m/ct_fz.m), 只对 x⊥ 求积。
%   不使用"下移 σ 的路径"(旧版做法, 对二次分母等价于文档② §2.5 的整函数 Z 陷阱)。
%   与 σ 无关: 自检 [T14] 直接比较不同 cfg.sigma 下的 D(ω)。
%
% ---------------------------------------------------------------------------------------
% disp_model = 'pole' —— **单极点约化** (PDF §7.3–7.5)
%   略去共振分母中的 k∥ v∥, 并把 ω_Ds 写成单变量 ω_Ds ≈ C_s u (u = x²):
%       J_s = 1 + ω_*s * I_s ,   I_s 由 ct_pole_integral(ω, C_s, η_s) 给出
%       C_s = -(4/3) g_s ω_D0s     (由 <2x∥²+x⊥²>_u = (4/3)u 严格得到; 见 ct_ct.m)
%   于是速度空间积分塌缩为**单个极点** u_* = ω/C_s 的留数 + 主值。
%   这是 PDF 第 7 章的核心结果, 也是 C_T = 3/2 - u_* 的来源。
%   【注意】该模型对二维(TEM、外中平面附近)模较可靠; 对 k∥ v_t >> ω_D 的模会高估漂移共振。
%
% ---------------------------------------------------------------------------------------
% [v6] **响应函数的唯一出处**
%   本函数不再自己算 𝒥_s, 而是调用 ct_chanresp.m(它对每个物种给出 𝒥_s 以及
%   ∂𝒥_s/∂X_s、∂𝒥_s/∂Y_s、∂𝒥_s/∂ω)。这样"色散关系用到的响应"与"C_T 里求导的响应"
%   在构造上必然一致 —— 这正是 PDF 式 (3) 的前提(见 ct_ct_exact.m 的文件头)。
%   ⚠ 重构后 D(ω) 的数值必须与旧版**逐位相同**, 自检 [T21a] 直接对照两者的复现性
%   (旧版公式: 离子 (ω-ω_*i)I0 - ω_*i η_i I1 - I2; 电子 (1-f_t)H_ps + f_t(Γ0(b_e) - J_t))。
%
% 用法:  D = ct_dispersion(M, omega)                  % 默认 'kinetic'
%        D = ct_dispersion(M, omega, 'pole')
%        [D, dDdw] = ct_dispersion(M, omega, 'kinetic', true)   % [v8.5] 解析 ∂D/∂ω
%=========================================================================================

if nargin < 3 || isempty(disp_model)
    disp_model = 'kinetic';
    if isfield(M.cfg,'disp_model') && ~isempty(M.cfg.disp_model)
        disp_model = M.cfg.disp_model;
    end
end
if nargin < 4 || isempty(want_deriv), want_deriv = false; end

% 把(可能被调用者覆盖的)模型选择放回 cfg, 供 ct_chanresp 读取
Ms = M;
if ~isfield(Ms, 'cfg') || ~isstruct(Ms.cfg), Ms.cfg = struct(); end
Ms.cfg.disp_model = disp_model;

tau = M.G.tau;

% ---- [v8.5] 解析导数: ct_chanresp 在 want_deriv=true 时给出 𝒥 的 ∂/∂ω(以及 ∂/∂X, ∂/∂Y);
%   单极点(pole)模型没有解析导数 ⇒ 只返回 D, 导数置 NaN(调用者退回数值差分)。
want_deriv = want_deriv && ~strcmpi(disp_model,'pole');
Ji = ct_chanresp(Ms, 'i', omega, want_deriv);
Je = ct_chanresp(Ms, 'e', omega, want_deriv);

D  = (M.adiab(M.b_i) - Ji.J)/tau + (M.adiab(M.b_e) - Je.J);
if nargout > 1 || want_deriv
    if want_deriv && isfield(Ji,'JW') && isfield(Je,'JW') && Ji.ok && Je.ok
        dDdw = -(Ji.JW/tau + Je.JW);        % A(b) 与 ω 无关 ⇒ 导数只来自 𝒥
    else
        dDdw = NaN;
    end
end
if nargout < 2, clear dDdw; end
end
