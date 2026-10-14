function G = ct_geometry(P, cfg)
%CT_GEOMETRY  由外部剖面构建磁面几何与**与 k 无关**的特征量
%=========================================================================================
% 输入
%   P   : ct_read_profile 给出的剖面结构体
%   cfg : 配置(见 CT_toroidal_defaults.m); 本函数只用到 A_ion / ft_model / ft_fixed /
%         trap_geom / no_curv
%
% 输出 G 的字段
%   mi,me          离子/电子质量 [kg]
%   Ti,Te          [J] ;  tau = Ti/Te
%   r, eps         r [m], eps = r/R0
%   Rth, Bth       该极向角处的 R, B  (B ≈ B0 R0/R 的圆截面近似)
%   rho_i, rho_s, rho_e   热回旋半径 / 离子声回旋半径 / 电子热回旋半径 [m]
%   vti, vte       热速度 [m/s]
%   ft             捕获电子份额
%   eta_i, eta_e
%   lam_i, lam_e   磁漂移几何因子 g_s (ω_Ds = -ω_D0s (2x∥²+x⊥²) g_s)
%   rho_b, fbw_on  有限香蕉宽度因子
%   kpar_of_theta  函数句柄: k∥ 由 k_⊥ 与 theta 给出
%
% 频率约定(SI, rad/s; q_s 带符号)
%   ω_*s = -k⊥ T_s/(q_s B L_ns)     ==>  离子 ω_*i < 0 ,  电子 ω_*e > 0
%   ω_D0s =  k⊥ T_s/(q_s B R)       ==>  离子 ω_D0i > 0 , 电子 ω_D0e < 0
%   ω_Ds = -ω_D0s (2x∥² + x⊥²) g_s
%
% 回旋半径的两套定义(务必区分)
%   ρ_s^{热} = sqrt(2 T_s/m_s)/|Ω_s|            (b_s = k⊥²ρ²/2 与 J_0 用它)
%   ρ_s^{声} = sqrt(  T_e/m_i)/Ω_ci = c_s/Ω_ci  (**波数单位 k_y ρ_s 用它**)
%   换算: k_θ ρ_i = k_y ρ_s · sqrt(2 T_i/T_e) = k_y ρ_s · sqrt(2 tau)
%=========================================================================================

QE = 1.602176634e-19;  MP = 1.67262192369e-27;  ME = 9.1093837015e-31;

% ---- 剖面字段的兜底派生 ----
% ct_read_profile 会给出这些派生量; 但本函数也允许被**单独调用**(例如 ct_selftest.m
% 里手工构造的测试剖面), 因此在这里补齐, 避免 "无法识别的字段名称 'Lni'" 之类的报错。
req = {'rho','R0','B0','a','q','shear','theta','n_e','T_e','RL_ne','RL_Te','RL_ni','RL_Ti'};
miss = req(~isfield(P, req));
if ~isempty(miss)
    error('ct_geometry: 剖面缺少必需字段: %s', strjoin(miss, ', '));
end
if ~isfield(P,'A_ion') || isempty(P.A_ion), P.A_ion = 2.0; end
if ~isfield(P,'r')     || isempty(P.r),     P.r = P.rho*P.a; end
if ~isfield(P,'n_e')   || isempty(P.n_e),   error('ct_geometry: 缺少 n_e'); end
if ~isfield(P,'T_e')   || isempty(P.T_e),   error('ct_geometry: 缺少 T_e'); end
if ~isfield(P,'n_i')   || isempty(P.n_i),   P.n_i = P.n_e; end
if ~isfield(P,'T_i')   || isempty(P.T_i),   P.T_i = P.T_e; end
if ~isfield(P,'Lne')   || isempty(P.Lne),   P.Lne = P.R0/P.RL_ne; end
if ~isfield(P,'Lni')   || isempty(P.Lni),   P.Lni = P.R0/P.RL_ni; end
if ~isfield(P,'LTe')   || isempty(P.LTe),   P.LTe = P.R0/P.RL_Te; end
if ~isfield(P,'LTi')   || isempty(P.LTi),   P.LTi = P.R0/P.RL_Ti; end
if ~isfield(P,'eta_i') || isempty(P.eta_i), P.eta_i = P.RL_Ti/P.RL_ni; end
if ~isfield(P,'eta_e') || isempty(P.eta_e), P.eta_e = P.RL_Te/P.RL_ne; end

% 把剖面量直接挂到 G 上, 下游一律用 G.xxx(不再经过 G.P.xxx), 少一层间接、少一类笔误
G.P   = P;
G.RL_ne = P.RL_ne;  G.RL_Te = P.RL_Te;  G.RL_ni = P.RL_ni;  G.RL_Ti = P.RL_Ti;
G.Lne = P.Lne;      G.Lni   = P.Lni;    G.LTe   = P.LTe;    G.LTi   = P.LTi;

G.mi = P.A_ion*MP;
G.me = ME;
G.Ti = P.T_i*QE;                 % [J]
G.Te = P.T_e*QE;                 % [J]
G.tau = G.Ti/G.Te;

G.r   = P.r;
G.eps = P.r/P.R0;
G.th  = P.theta;
G.Rth = P.R0 + P.r*cos(P.theta);
G.Bth = P.B0*P.R0/G.Rth;         % 圆截面近似 B ≈ B0 R0/R

G.rho_i = sqrt(2*G.Ti*G.mi)/(QE*G.Bth);      % 离子热回旋半径
G.rho_e = sqrt(2*G.Te*G.me)/(QE*G.Bth);      % 电子热回旋半径
G.rho_s = sqrt(   G.Te*G.mi)/(QE*G.Bth);     % 离子声回旋半径 (k_y ρ_s 的单位)
G.vti   = sqrt(2*G.Ti/G.mi);
G.vte   = sqrt(2*G.Te/G.me);

% ---- 频率的归一化单位 c_s/a ----
%   c_s = sqrt(T_e/m_i)  —— 与上面 ρ_s = c_s/Ω_ci 用的是**同一个** c_s 约定
%   本程序默认把 ω_r 与 γ 以 c_s/a 为单位输出(见 cfg.omega_unit);
%   而 C_T 只依赖频率之比 u_* = |ω_r|/(|ω_D0e| |c_D|), 与单位选择无关。
G.cs   = sqrt(G.Te/G.mi);                    % 离子声速 [m/s]
G.cs_a = G.cs/P.a;                           % 频率归一化单位 [rad/s]

% ---- 捕获份额 ----
if isfield(cfg,'ft_model') && strcmpi(cfg.ft_model,'fixed')
    G.ft = cfg.ft_fixed;
else
    G.ft = sqrt(2*G.eps/(1+G.eps));          % 标准结果约 sqrt(2 eps)
end

% ---- eta ----
G.eta_i = P.RL_Ti/P.RL_ni;
G.eta_e = P.RL_Te/P.RL_ne;

% ---- [v8.4] s-α 压力参数 α = −q²R dβ/dr (β = 8πp/B²) ----
%   α 是 Shafranov 位移/压力梯度参数, 进入三处(见 Kim & Han 2017 式 (5) 后的定义):
%     ① 磁漂移的极向结构: g(θ) = cosθ + sinθ(ŝθ − α sinθ)
%     ② 垂直波数的 WKB 因子: k⊥² = k_y²[1 + (ŝθ − α sinθ)²]
%     ③ 平行波数:           k∥ = k_y(ŝθ − α sinθ)/(qR)
%   ★ 三者都在 θ=0 处退化为无 α 的形式(g=1, k⊥=k_y, k∥-因子=0) ⇒ 本程序默认 θ=0 时,
%     α 只能通过**捕获电子的弹跳平均**(沿整条轨道积分, 见 ct_trapped_avg.m)起作用。
G.alpha = 0;
if isfield(cfg,'alpha') && ~isempty(cfg.alpha), G.alpha = cfg.alpha; end
G.alpha_on = abs(G.alpha) > 0;
G.shear_eff = P.shear*P.theta - G.alpha*sin(P.theta);     % ŝθ − α sinθ (θ 处的有效剪切)
G.kperp_fac = sqrt(1 + G.shear_eff^2);                    % k⊥/k_y 的 s-α 因子

% ---- 磁漂移几何因子 ----
% 局域气球模因子 g = cosθ + sinθ(ŝθ − α sinθ)  (θ=0 时 = 1; α=0 时退化为 cosθ + ŝθsinθ)
if isfield(cfg,'no_curv') && cfg.no_curv
    G.lam_i = 0;  G.lam_e = 0;
else
    G.lam_i = cos(P.theta) + G.shear_eff*sin(P.theta);
    if isfield(cfg,'trap_geom') && isnumeric(cfg.trap_geom)
        G.lam_e = cfg.trap_geom;             % 允许直接给数值(体现弹跳平均的减弱)
    else
        G.lam_e = G.lam_i;
    end
end

% ---- 有限香蕉宽度(FBW): ρ_b = q ρ_i / sqrt(eps) ----
G.rho_b = P.q*G.rho_i/sqrt(max(G.eps,1e-6));
G.fbw_on = ~(isfield(cfg,'fbw_model') && strcmpi(cfg.fbw_model,'off'));
if isfield(cfg,'fbw_factor'), G.fbw_factor = cfg.fbw_factor; else, G.fbw_factor = 1; end

% ---- 平行波数模型 ----
G.q = P.q;  G.shear = P.shear;  G.R0 = P.R0;
G.kpar_model  = 'shear';   if isfield(cfg,'kpar_model'),  G.kpar_model  = cfg.kpar_model;  end
G.kpar_factor = 1.0;       if isfield(cfg,'kpar_factor'), G.kpar_factor = cfg.kpar_factor; end

% ---- 碰撞(去捕获)频率 [v4] ----
% 【为什么要它】TEM 的主要驱动是捕获电子的进动共振; 碰撞把电子从陷阱里"去捕获"出去,
%   等效于给该共振一个宽度 nu_eff。TGLF 在矩矩阵里显式含 detrapping, 本程序原版完全没有
%   碰撞(nu=0)。这里只算**与 k 无关**的部分; 逐 k 的 nu_eff 在 ct_model.m 里给出:
%       lnL    = 24 - ln( sqrt(n_e[cm^-3]) / T_e[eV] )        (未显式给出时)
%       nu_ei  = 2.91e-6 * n_e[cm^-3] * Zeff * lnL / T_e[eV]^1.5      [s^-1]
%       nu_eff = 0                                     'zero'   (默认=原版行为)
%              = nu_eff_factor * nu_ei / eps           'collisional'  (= full_gk_ct 的取法)
%              = nu_eff_factor * 0.163*sqrt(k∥v_te/√2 * nu_ei * (1+0.82*Zeff))
%                                                      'tglf'   (= TGLF 自己的 cb1, 随 k 变化)
%              = nu_eff_factor [s^-1]                  'fixed'
G.Zeff = 1.0;   if isfield(cfg,'Zeff'), G.Zeff = cfg.Zeff; end
G.ne_cm3 = P.n_e*1e-6;
if isfield(cfg,'lnL') && ~isempty(cfg.lnL)
    G.lnL = cfg.lnL;
else
    G.lnL = 24 - log(sqrt(G.ne_cm3)/P.T_e);
end
G.nu_ei = 2.91e-6*G.ne_cm3*G.Zeff*G.lnL/max(P.T_e,1e-9)^1.5;
G.nu_eff_model  = 'zero';  if isfield(cfg,'nu_eff_model'),  G.nu_eff_model  = cfg.nu_eff_model;  end
G.nu_eff_factor = 1.0;     if isfield(cfg,'nu_eff_factor'), G.nu_eff_factor = cfg.nu_eff_factor; end
% 与 k 无关那部分的口径(逐 k 的 nu_eff 在 ct_model.m 里给出; 'tglf' 依赖 k∥ 故这里记 0)
switch lower(char(G.nu_eff_model))
    case 'fixed'
        G.nu_eff = G.nu_eff_factor;
    case 'collisional'
        G.nu_eff = G.nu_eff_factor*G.nu_ei/max(G.eps,1e-6);
    otherwise
        G.nu_eff = 0;
end
% [v8.5] 几何指纹: 记录生成本 G 的 cfg 字段 + 剖面量, 供 ct_configure 检测"改了 cfg 未重建 G"
G.fp = ct_geom_fp(cfg);
G.fp.shear = P.shear;
G.fp.theta = P.theta;
end
