function M = ct_model(G, cfg, k_rhos)
%CT_MODEL  构建**给定波数** k_y ρ_s 处的局域模型(环形推广)
%=========================================================================================
% 输入
%   G      : ct_geometry 的输出
%   cfg    : 配置
%   k_rhos : k_y ρ_s (离子声回旋半径单位; 本程序默认的波数单位)
%
% 输出 M(供 ct_dispersion.m 使用)关键字段
%   kperp      k⊥ [1/m]  (= k_θ, 强气球模近似 k_r << k_θ)
%   kpar       k∥ [1/m]
%   wstar_i, wstar_e        抗磁频率(带符号)
%   wD0_i,  wD0_e           漂移频率标度(带符号)
%   b_i, b_e                FLR 参数  b_s = k⊥²ρ_s²/2
%   adiab(b)                极化项  Γ_0(b) = e^{-b} I_0(b)
%   xperp, wperp            **仅 x⊥ 方向**的 Gauss-Legendre 节点与权重
%   j0i, j0e                J_0²(k⊥ρ_s x⊥) (与 x∥ 无关, 逐节点)
%   FBW                     有限香蕉宽度因子 J_0²(k⊥ρ_b)
%
% 【速度空间: x∥ 解析积分, 只留 x⊥ 数值求积】
%   环形漂移 ω_Ds = -ω_D0s (2x∥² + x⊥²) g_s 含 x∥², 使共振分母对 x∥ 为**二次**——
%   这就是《DGH问题_磁漂移共振二次分母_详解.md》(文档②) 讨论的困难。本程序按文档②
%   §2 的**精确路线**处理: 对每个 x⊥ 把 x∥ 积分解析地做出(二次式因子化 + 部分分式 +
%   分支正确的 f 函数, 见 ct_xpar.m / ct_fz.m), 只对 x⊥ 做一维 Gauss-Legendre 求积。
%
%   【为什么不再用"下移 σ 的路径"】旧版用一个整体下移 σ 的复路径来做所谓的朗道解析
%   延拓。对**线性**分母(单极点)那是对的; 但对**二次**分母, 两个根恒分居实轴两侧,
%   且 ω 在上半平面移动时没有根穿越实轴, 因此不存在需要延拓的跳跃 —— 下移路径额外
%   拾取条带 (-σ, 0) 内极点的残数, 恰好等于文档② §2.5 判为错误、误差 O(1) 的
%   "把根代进整函数 Z"的值。解析路线既精确又与 σ 无关(自检 [T14] 直接检验这一点)。
%=========================================================================================

g = cfg;
% ---- [v8.5] 配置归一化 + 与几何指纹的一致性检查(单一入口, 见 ct_configure.m) ----
%   旧版在 ct_model 里对 alpha / fbw_model / nu_eff_model 各写一段"cfg 与几何副本不一致"的
%   warning, 既重复又容易漏掉新开关。现在统一交给 ct_configure: 它会一次性列出**所有**与
%   G.fp 不一致的字段(并做枚举值校验/预设归一化), 同一组合只告警一次。
%   注意: 归一化(如 kpar_preset 'cbc_benchmark' → kpar_model 'balloon_rms')在 ct_configure
%   里已经做过, 下面第 141 行的预设分支因此不会重复生效(值已不是出厂默认 'balloon')。
%   cfg.quiet_geometry_check = true 时不查指纹 —— 供单元测试(ct_selftest/ct_selfcheck)使用:
%   它们**故意**拿改过的 cfg 配同一套几何, 那些"不一致"是被测行为本身, 不是脚本 bug。
chk_fp = true;
if isfield(g,'quiet_geometry_check') && ~isempty(g.quiet_geometry_check) && g.quiet_geometry_check
    chk_fp = false;
end
if chk_fp && isfield(G,'fp') && isstruct(G.fp)
    g = ct_configure(g, G);
end
defs = {'kpar_model','shear'; 'kpar_factor',1.0; 'fbw_model','banana'; 'fbw_factor',1.0; ...
        'N_perp',64};
for ii = 1:size(defs,1)
    if ~isfield(g, defs{ii,1}) || isempty(g.(defs{ii,1}))
        g.(defs{ii,1}) = defs{ii,2};
    end
end

% ---- 波数换算: k_y ρ_s -> k_θ ρ_i ----
%   k_θ ρ_i = k_y ρ_s · sqrt(2 T_i/T_e)
M.k_rhos = k_rhos;
M.kperp  = k_rhos * sqrt(2*G.tau) / G.rho_i;      % [1/m] = k_θ (ω_*、ω_D、k∥ 的基准)
% [v8.4] s-α 的 FLR/WKB 因子: k⊥² = k_y²[1+(ŝθ−α sinθ)²] ⇒ k⊥ = k_θ·G.kperp_fac
%   ★ 只用于 FLR/Bessel 类量(b_i, b_e, J0, FBW_trap); **不**用于 ω_*、ω_D0 的幅值
%     (局域 s-α 里 ω_* = k_y T/(eBL_n)、ω_D0 = 2ε_n ω_*, 都只含 k_y)。
%   θ=0 或 α=0 时 G.kperp_fac = 1 ⇒ 与 v8.3 逐位相同。
M.kperp_flr = M.kperp * G.kperp_fac;
% ---- [v8.4] s-α 参数 α: **cfg 优先**, 与几何副本不一致时告警(同 v8.3 的 nu_eff/fbw 政策) ----
%   注: G.alpha / G.kperp_fac / G.shear_eff 是 ct_geometry 从 cfg 抄进去的副本;
%   只改 cfg.alpha 而不重建 G 会让 FLR 因子仍按旧 α 计算 —— 这里至少保证 α(进 k∥ 与
%   捕获电子弹跳平均)用 cfg 的值, 并提醒重建。
alpha_v = G.alpha;
if isfield(g,'alpha') && ~isempty(g.alpha) && g.alpha ~= G.alpha
    alpha_v = g.alpha;
    warning('ct_model:staleGeometry', ...
            ['cfg.alpha = %.4g 与几何副本 G.alpha = %.4g 不一致; k∥ 与捕获电子弹跳平均\n' ...
             '  按 cfg 计算, 但 k⊥(FLR)/漂移因子仍是旧 α。建议重建: G = ct_geometry(P, cfg)。'], ...
            alpha_v, G.alpha);
end
M.alpha = alpha_v;
M.trap_precession = 'legacy';
if isfield(g,'trap_precession') && ~isempty(g.trap_precession)
    M.trap_precession = lower(char(g.trap_precession));
end
if strcmp(M.trap_precession,'auto')          % 自动: α≠0 时启用含 ŝ/α 的完整弹跳平均
    if G.alpha_on, M.trap_precession = 'full'; else, M.trap_precession = 'legacy'; end
end

% =======================================================================================
% 平行波数 k∥ —— [v8] 五种模型 + 显式开关
% ---------------------------------------------------------------------------------------
% 局域模型里 k∥ 决定三件事: ① 离子/通行电子的平行朗道共振 (ω → ω - k∥v∥);
% ② 通行电子 Z 函数闭式的自变量 ζ = ω/(k∥v_te); ③ θ=0 时若 k∥=0 则**完全没有平行阻尼**。
%
% 【气球模的几何因子 G(推导与三种口径)】对 (m,n) 分量, k∥ = (nq - m)/(qR)。带剪切的模在极向角
% θ 处的"局部极向模数"由剪切位移给出: 标准 s-α 度规 k_⊥² = k_θ²(1 + ŝ²θ²) ⟹ 径向波数
% k_r(θ) = k_θ·ŝ·θ(与 k_θ 同量纲), 于是
%       k∥(θ) = k_θ·ŝ·θ·G                [cfg: 'balloon' / 'balloon_rms']
% 【v8.4 s-α】含压力梯度 α 时, 有效剪切因子变成 **ŝθ − α sinθ**(见 ct_geometry.m 的
%   G.shear_eff): 标准 s-α 度规是 k_⊥² = k_θ²[1+(ŝθ−α sinθ)²] ⟹ k_r(θ) = −k_θ(ŝθ−α sinθ),
%   于是 k∥ = k_θ(ŝθ−α sinθ)G。α=0 时退回上式; `shear_legacy` 保持 v7 形式(仅作回归对照)。
% G 的**取值口径有三种**(详见本文件后面的"几何因子 G —— 三种口径与实测标定", 那里有 CBC 数值):
%   (i)  局部 (m,n) 口径 m(θ) = nq + k_r·r ⟹ G = r/(qR)     (CBC: 10.8 /m, 偏大 1-2 个量级)
%   (ii) 归一化单位口径(长度以 ρ_i 归一)    ⟹ G = ρ_i/(qR)   (CBC: 0.14 /m) ★auto 默认
%   (iii) 气球模包络口径 Δk∥ ~ 1/(qRΔθ)     ⟹ 与 k_θ 无关      (CBC: 0.43 /m, 即 'qR' 模型)
% v8 之前程序用的是 G = 1(即 k∥ = k_⊥ŝθ), 是这三种之外的另一端(缺因子)。
% G 可用 cfg.kpar_geom('auto' | 'r' | 数值)覆盖。
%
% 【为什么还要 'balloon_rms'】单 θ 模型在 θ=0 处 k∥ ≡ 0, 于是高 k 端没有平行阻尼、γ 无界
% 增长(与全 GK 参考比值可达 2.4)。真实气球模沿场线有有限宽度, 其**包络平均**给
%       k∥^eff = k_θ·ŝ·θ_rms·G ,        [cfg: 'balloon_rms', cfg.kpar_theta_rms]
% θ_rms 是包络的等效角宽(默认 1.0 rad)。它让 θ=0 的算例也能含平行动力学 —— 与全 GK 对照时
% 推荐用它(见 benchmark_cbc/kpar_models.md 的标定结果)。
%
% 五种模型 (cfg.kpar_model):
%   'zero'          k∥ = 0                      —— 纯二维(无平行阻尼); 显式声明用
%   'balloon'       k∥ = k_θ·ŝ·θ·G             —— 角度处的局域气球模关系(θ=0 时仍为 0)
%   'balloon_rms'   k∥ = k_θ·ŝ·θ_rms·G         —— 包络等效(θ=0 时非零, 推荐与全 GK 对比)
%   'qR'            k∥ = kpar_factor/(qR)      —— 连接长度标度(与 θ 无关)
%   'shear_legacy'  k∥ = k_⊥·ŝ·θ               —— v7 及以前的形式(缺几何因子), 仅作回归对照
% 另有: cfg.kpar_factor(qR 用, 默认 1)、cfg.kpar_geom('auto'|[数值])、
%       cfg.kpar_theta_rms(balloon_rms 用, 默认 1.0)、cfg.kpar_sign(+1/-1, 默认 +1)。
% =======================================================================================
% 【几何因子 G —— 三种口径与实测标定(务必读)】
%   对 (m,n) 分量有 k∥ = (nq - m)/(qR)。把"局部极向模数"换成剪切位移有两种写法:
%     (i)  局部 (m,n) 口径: m(θ) = nq + k_r(θ)·r, k_r = k_θ·ŝθ ⟹ k∥ = k_θŝθ·(r/(qR)); 
%     (ii) 归一化单位口径(文献常见的 k̂∥ = k̂_y ŝθ/(q R̂), 长度以 ρ_i 归一, R̂ = R/ρ_i)
%          ⟹ k∥ = k_θŝθ·(ρ_i/(qR))。
%   而**气球模包络**口径完全不同: 包络在位相角 θ 上宽 Δθ ⟹ 平行波数谱宽
%     Δk∥ ~ 1/(qR·Δθ)  —— **与 k_θ 无关**(这正是 'qR' 模型), 对 Δθ ~ O(1) 即 1/(qR)。
%   ★ 三者在 Cyclone base case 上的量级: (i) 10.8 /m, (ii) 0.137 /m, 'qR' 0.43 /m (ky=0.4, θ=1);
%     而 GX 全回旋动理学参考显示该模的**有效平行阻尼是中等强度**, 需要 k∥v_ti/ω_* ~ O(1),
%     即 k∥ ~ 0.1–0.5 /m ⟹ (ii) 与 'qR' 同量级、物理合理; (i) 与 v7 的 G=1 都偏大 1–2 个量级。
%     (ii) 的 0.137 /m 与 k∥ 标定给出的**实测最优 0.137 /m** 相同(同一设置: balloon_rms + θ_rms=1)。
%   故 **auto = ρ_i/(qR)**; 需要 (i) 或其它取值时用 cfg.kpar_geom 给数值覆盖。
%   标定表: ../benchmark_cbc/cbc_kpar_scan.txt
km = 'balloon';
if isfield(g,'kpar_model') && ~isempty(g.kpar_model)
    km = lower(char(g.kpar_model));
    if strcmpi(km,'shear'), km = 'shear_legacy'; end     % 旧名兼容(= v7 的形式)
end
ksign = 1;  if isfield(g,'kpar_sign') && ~isempty(g.kpar_sign), ksign = sign(g.kpar_sign); end
kfac  = 1.0; if isfield(g,'kpar_factor') && ~isempty(g.kpar_factor), kfac = g.kpar_factor; end
th_rms = 1.0;
if isfield(g,'kpar_theta_rms') && ~isempty(g.kpar_theta_rms), th_rms = g.kpar_theta_rms; end
% ---- [v8] 一键预设: cfg.kpar_preset = 'cbc_benchmark' 套用 CBC 标定过的设置 ----
%   标定结果(benchmark_cbc/cbc_kpar_scan.txt): k∥=0 在高 k 端高估 γ 达 2.4 倍;
%   balloon_rms + θ_rms = 1.0(取 auto 几何因子 ⇒ k∥ = 0.069 /m @ky=0.2, 即 0.137 /m @ky=0.4)
%   给出平均 γ 比值 1.015, 是本模型结构下与全 GK 参考最接近的设置。
%   需要更强/更弱的平行阻尼时改 θ_rms(2 已明显过阻尼, 4 完全无模)。
%
%   【v8.2 优先级修正】预设只填"仍处于出厂默认"的字段: 出厂默认 kpar_model = 'balloon',
%   故只有当用户**没有**把 kpar_model 改成别的值时, 预设才把它升为 'balloon_rms';
%   用户显式给出的 kpar_model / kpar_theta_rms / kpar_factor 一律**优先**。
%   (旧写法无条件 `km = 'balloon_rms'; th_rms = 1.0; kfac = 1.0;`, 于是
%    `cfg.kpar_preset='cbc_benchmark'; cfg.kpar_model='qR';` 会被静默改成 balloon_rms ——
%    由 cbc_tem.m 的 k∥ 敏感性表发现: 四种设置的 γ 竟然逐位相同。)
if isfield(g,'kpar_preset') && ~isempty(g.kpar_preset)
    switch lower(char(g.kpar_preset))
        case 'cbc_benchmark'
            if strcmp(km, 'balloon')            % 出厂默认 ⇒ 未被用户显式指定
                km = 'balloon_rms';
            end
            if th_rms == 1.0 && kfac == 1.0     % 出厂默认 ⇒ 保持标定用的值(1.0/1.0)
                th_rms = 1.0;  kfac = 1.0;
            end
        case {'none','zero'}
            if strcmp(km, 'balloon'), km = 'zero'; end
        otherwise
            error('ct_model:badKparPreset', ...
                  '未识别的 cfg.kpar_preset = ''%s''; 可选: cbc_benchmark | none', char(g.kpar_preset));
    end
end
% 几何因子 G = ρ_i/(qR)(见上); 允许 cfg.kpar_geom 给数值覆盖, 或 'r' 用 r/(qR)
Ggeom = G.rho_i/(G.q*G.R0);
if isfield(g,'kpar_geom') && ~isempty(g.kpar_geom)
    if isnumeric(g.kpar_geom)
        Ggeom = g.kpar_geom;
    elseif strcmpi(char(g.kpar_geom),'r')
        Ggeom = G.r/(G.q*G.R0);
    elseif ~strcmpi(char(g.kpar_geom),'auto')
        error('ct_model:badKparGeom', ...
              'cfg.kpar_geom 只能是 ''auto''(=ρ_i/(qR)) | ''r''(=r/(qR)) | 数值; 当前 ''%s''', ...
              char(g.kpar_geom));
    end
end
switch km
    case 'zero'
        M.kpar = 0;
    case 'balloon'
        % [v8.4] s-α: k∥ ∝ (ŝθ − α sinθ); 取**绝对值** —— 该因子可因 ŝθ = α sinθ 而穿零,
        %   但"有效平行波数"不能穿零/变号(否则会人为去掉平行朗道阻尼, 见 α≈ŝ 的情形)。
        M.kpar = M.kperp * abs(G.shear*G.th - alpha_v*sin(G.th)) * Ggeom;
    case 'balloon_rms'
        M.kpar = M.kperp * abs(G.shear*th_rms - alpha_v*sin(th_rms)) * Ggeom;
    case 'qr'
        M.kpar = kfac/(G.q*G.R0);
    case 'shear_legacy'
        M.kpar = M.kperp * G.shear * G.th;             % v7 旧形式(回归用, 不含 α)
    otherwise
        error('ct_model:badKparModel', ...
              ['未识别的 cfg.kpar_model = ''%s''; 可选: zero | balloon | balloon_rms | qR | ' ...
               'shear_legacy'], km);
end
M.kpar = ksign*M.kpar;
M.kpar_model = km;  M.kpar_geom = Ggeom;  M.kpar_theta_rms = th_rms;

% ---- 抗磁频率(带符号) ----
%   局域 s-α 里 ω_* = k_y T/(e B L_n): **只含 k_y**, 不含 s-α 的 k⊥ 因子
M.wstar_i = -(M.kperp*G.Ti)/(QE_()*G.Bth*G.Lni);
M.wstar_e = +(M.kperp*G.Te)/(QE_()*G.Bth*G.Lne);
% ---- 漂移频率标度(带符号); ω_D0 = 2ε_n ω_* 亦只含 k_y ----
M.wD0_i = +(M.kperp*G.Ti)/(QE_()*G.Bth*G.R0);
M.wD0_e = -(M.kperp*G.Te)/(QE_()*G.Bth*G.R0);

% ---- FLR 参数与极化项 ----
%   b_s = k⊥²ρ_s²/2 ;  Γ_0(b) = e^{-b} I_0(b)  (= <J_0²(k⊥ρ)>)
%   [v8.4] 这里用含 s-α 的 k⊥ = k_y√(1+(ŝθ−αsinθ)²); θ=0 或 α=0 时与 v8.3 相同
M.b_i = M.kperp_flr^2*G.rho_i^2/2;
M.b_e = M.kperp_flr^2*G.rho_e^2/2;
if isfield(g,'polar_term') && strcmpi(g.polar_term,'J0')
    M.adiab = @(b) exp(-b/2).*besseli(0,b/2);        % <J_0(b)> 约定(与 weiland_CT 兼容)
else
    M.adiab = @(b) exp(-b).*besseli(0,b);            % Γ_0(b) 约定(本文档 PDF 的约定, 默认)
end

% ---- 有限香蕉宽度因子 ----
% 【两种取法, 见 README §7.3 / 程序说明 §6.7】
%   'banana'     (默认)  FBW = J_0²(k⊥ρ_b), ρ_b = fbw_factor·qρ_i/√ε 取**热**香蕉宽度
%                        —— 与 weiland_CT 的公开参考表逐点一致([T9]), 但 J_0² 在
%                        k⊥ρ_b = 2.405 处**恰为零**, 把捕获电子的非绝热响应整个关掉,
%                        使 TEM 支(ω_r>0)在 k_yρ_s ≈ 0.55 之后就消失(见 §6.7 的表)。
%   'banana_avg' 对捕获电子的热速度分布取平均(ρ_b ∝ 热速度, 垂直 Maxwellian 权重):
%                        ⟨J_0²(k⊥ρ_b)⟩ = e^{-b_b} I_0(b_b),  b_b = (k⊥ρ_b)²/2
%                        —— 与 adiab(Γ_0) 同一函数形式, **恒正、单调衰减**, 没有假零点;
%                        于是 TEM 支延伸到 k_yρ_s ≈ 0.85, 且 k = 0.5 处 γ_TEM 约 2 倍。
%   'off'        关闭(G_geometry 里 fbw_on = false)
%
% 【v8.3 修正: 开关的 cfg 优先 + "陈旧几何"告警】
%   这些开关在 ct_geometry 里留了副本(G.fbw_on / G.nu_eff_model), ct_model 过去**只读副本**,
%   于是"改了 cfg 但没重建 G"时开关**静默无效** —— 由 benchmark_cbc/cbc_tem.m 的 K1c/K5
%   发现(把 fbw_model 设成 'off' 后 D(ω) 与默认逐位相同)。现在: **cfg 优先**, 且当 cfg 与
%   几何副本不一致时明确告警, 提醒调用者重建 G = ct_geometry(P,cfg)。
fbw_model = '';
if isfield(g,'fbw_model') && ~isempty(g.fbw_model), fbw_model = lower(char(g.fbw_model)); end
G_fbw_on = true;
if isfield(G,'fbw_on') && ~isempty(G.fbw_on), G_fbw_on = logical(G.fbw_on); end
if isempty(fbw_model)
    if G_fbw_on, fbw_model = 'banana'; else, fbw_model = 'off'; end
end
if ~G_fbw_on && ~strcmp(fbw_model,'off')
    warning('ct_model:staleGeometry', ...
            ['几何是按 fbw_model=''off'' 建的(G.fbw_on = false), 而 cfg 现在给 ''%s''; ' ...
             '本次以 cfg 为准。建议改开关后重建几何: G = ct_geometry(P, cfg)。'], fbw_model);
elseif G_fbw_on && strcmp(fbw_model,'off') && isfield(g,'fbw_model') && ~isempty(g.fbw_model)
    warning('ct_model:staleGeometry', ...
            ['cfg.fbw_model=''off'' 是在几何建好之后才设的(G.fbw_on 仍为 true); ' ...
             '本次以 cfg 为准(即关闭 FBW)。建议重建几何: G = ct_geometry(P, cfg)。']);
end
if strcmp(fbw_model,'off')
    M.FBW = 1;
else
    xb = g.fbw_factor*M.kperp_flr*G.rho_b;      % [v8.4] 用含 s-α 的 k⊥
    switch fbw_model
        case 'banana_avg'
            bb = xb^2/2;
            M.FBW = exp(-bb)*besseli(0,bb);
        otherwise
            M.FBW = besselj(0, xb)^2;
    end
end
M.fbw_model_used = fbw_model;

% ---- 捕获电子通道的碰撞(去捕获)展宽 [v4] ----
%   ν_eff 的口径由 cfg.nu_eff_model 决定; ct_dispersion 用它把捕获电子通道写成 ω + i·ν_eff。
%   'tglf' 用 TGLF 自己的 detrapping 率(GACODE src_tglf_eigensolver.f90:464)
%       cb1 = 0.163*sqrt(k∥v_te/√2 · ν_ei · (1+0.82*Zeff))       —— 随 k 变化, 且 k∥=0 时为 0
%   这是与 TGLF 对比时**最该用**的口径: 在本剖面 k∥=1/(qR) 处给 ≈1.3e5 s^-1,
%   比 'collisional'(ν_ei/ε = 2.4e5 s^-1) 小约 1.9 倍 —— 即对 TEM 的抑制更弱。
%   [v8.3] 同 FBW: **cfg.nu_eff_model 优先**, 与几何副本不一致时告警(旧写法只读 G 的副本)。
nu_model = '';
if isfield(g,'nu_eff_model') && ~isempty(g.nu_eff_model), nu_model = lower(char(g.nu_eff_model)); end
if isempty(nu_model) && isfield(G,'nu_eff_model') && ~isempty(G.nu_eff_model)
    nu_model = lower(char(G.nu_eff_model));
end
if isfield(G,'nu_eff_model') && ~isempty(G.nu_eff_model) && ...
   ~isempty(g.nu_eff_model) && ~strcmpi(char(G.nu_eff_model), nu_model)
    warning('ct_model:staleGeometry', ...
            ['cfg.nu_eff_model = ''%s'' 与几何里的副本 ''%s'' 不一致; 本次以 cfg 为准。\n' ...
             '  建议改开关后重建几何: G = ct_geometry(P, cfg)。'], ...
            nu_model, char(G.nu_eff_model));
end
nu_fac = G.nu_eff_factor;
if isfield(g,'nu_eff_factor') && ~isempty(g.nu_eff_factor), nu_fac = g.nu_eff_factor; end
switch nu_model
    case 'fixed'
        M.nu_eff = nu_fac;
    case 'collisional'
        M.nu_eff = nu_fac*G.nu_ei/max(G.eps,1e-6);
    case 'tglf'
        kparvthe = abs(M.kpar)*G.vte/sqrt(2);
        M.nu_eff = nu_fac*0.163*sqrt(max(kparvthe*G.nu_ei*(1 + 0.82*G.Zeff), 0));
    otherwise
        M.nu_eff = 0;
end
M.nu_ei  = G.nu_ei;

% ---- 捕获电子通道的**弹跳平均**模型 [v5] ----
%   cfg.trap_model = 'kinetic'(默认, 旧行为: 全 3D 麦氏积分 + 单一尖锐进动共振 + J0²(k⊥ρ_b))
%                  | 'bounce' (文档① §9.4/(9.8): 捕获区 (x, ξ₀) 均匀投角 + 弹跳平均进动频率)
%   弹跳平均表与 k 无关, 在这里一次算好; 逐 k 的求和见 ct_dispersion.m。
%   [v8.4] 表里加入 **s-α 效应**: 投角因子 B(ξ₀) = (2−λ̂)M₁ + λ̂εM₂c, 其中 M₁/M₂c 含
%     ŝ⟨θsinθ⟩ 与 **α⟨sin²θ⟩** 两族新矩 ⇒ α 通过"沿轨道的弹跳平均"改变捕获电子进动的大小与
%     方向(α>0 使其减小甚至反转 —— 这正是 Ma et al. 2023 的共振反转机制)。
%     cfg.trap_precession: 'auto'(默认; α≠0 ⇒ 'full') | 'legacy'(仅 cosθ, cos²θ, v8.3 行为)
%     | 'full'。
M.trap_on = false;
if isfield(g,'trap_model') && strcmpi(g.trap_model,'bounce')
    Nx = 32;  if isfield(g,'N_trap_x')  && ~isempty(g.N_trap_x),  Nx  = g.N_trap_x;  end
    Ni = 24;  if isfield(g,'N_trap_xi') && ~isempty(g.N_trap_xi), Ni = g.N_trap_xi; end
    if strcmp(M.trap_precession,'legacy')
        M.trap = ct_trapped_avg(G.eps, Nx, Ni);                  % 旧形式(ŝ, α 项不进入)
    else
        M.trap = ct_trapped_avg(G.eps, Nx, Ni, G.shear, alpha_v);
    end
    M.trap_on = true;
end

% ---- 弹跳平均通道的**轨道宽度(FOW)因子** [v5.1] ----
% 【为什么需要】文档① (9.8) 的捕获区平均里没有 J0²/Bessel 因子, 故 'off' 是 (9.8) 的字面实现;
%   但弹跳平均把 k∥v∥ 消掉之后, 捕获通道里所有频率(ω_*e、⟨ω_De⟩_b、k∥v_te)都与 k 成正比,
%   于是模型对 k **标度不变** ⇒ γ ∝ k, 高 k 无截断, 相对 TGLF 会明显**高估** ITG。
%   唯一能破坏标度不变性的机制是有限轨道宽度: k⊥ρ_b, ρ_b = q ρ_i/√ε。
%     'off'        1                                —— 文档① (9.8) 字面
%     'banana'     J0²(k⊥ρ_b)                       —— v4 约定(在 k⊥ρ_b = 2.405 处有**假零点**)
%     'banana_avg' ⟨J0²⟩ = e^{-b_b} I0(b_b)         —— 恒正、单调衰减(推荐, 无假零点)
%     'energy'     Σ_a w_x(a) J0²(k⊥ρ_b x_a)        —— 再对捕获电子的能量分布平均(最弱)
%   注意: 该因子只作用于**捕获**通道; 通行电子仍用 M.FBW(见 ct_dispersion.m)。
M.FBW_trap = 1;
if M.trap_on && isfield(g,'fbw_trap_model')
    switch lower(char(g.fbw_trap_model))
        case 'banana'
            M.FBW_trap = besselj(0, M.kperp_flr*G.rho_b).^2;
        case 'banana_avg'
            bb = (M.kperp_flr*G.rho_b)^2/2;
            M.FBW_trap = exp(-bb)*besseli(0, bb);
        case 'energy'
            s = 0;
            for ia = 1:numel(M.trap.x)
                s = s + M.trap.wx(ia)*besselj(0, M.kperp_flr*G.rho_b*M.trap.x(ia)).^2;
            end
            M.FBW_trap = s;
        otherwise
            M.FBW_trap = 1;
    end
end

% ---- 速度空间节点: 只需 x⊥ 方向 ----
%   x∥ 方向由 ct_xpar.m 解析积分(见文件头), 不再需要 x∥ 网格, 也不需要"下移路径"。
%   仍保留 cfg.N_perp(以 x⊥ 节点数), 因为它决定唯一剩下的数值求积。
[xx, wx] = ct_gaussleg(g.N_perp, 0, 6);      % x⊥ 节点(升序)

M.xperp = xx(:);
M.wperp = 2*xx(:).*exp(-xx(:).^2).*wx(:);    % ∫d²x⊥ e^{-x⊥²}(...) 的径向权重

% J_0²(k⊥ρ_s x⊥)(与 x∥ 无关)
M.j0i = besselj(0, M.kperp_flr*G.rho_i*xx(:)).^2;
M.j0e = besselj(0, M.kperp_flr*G.rho_e*xx(:)).^2;

M.G = G;  M.cfg = g;
end

% ---------------------------------------------------------------------------------------
function q = QE_()
q = 1.602176634e-19;
end
