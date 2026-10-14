function [res, Dg] = ct_balloon(G, cfg, k_rhos, opt)
%CT_BALLOON  [v8.5] **标准气球模本征值问题** —— 把平行导数 ∂/∂θ 真正解出来, 不做角度平均
%=========================================================================================
% 【本函数要解决的问题】
%   ct_model/ct_dispersion 是**单 θ 局域**模型: θ 只是一个参数, 于是平行波数只能用
%   "几何因子 G"(ρ_i/(qR)? r/(qR)? 1/(qRΔθ?) 三选一)加一个等效角宽 θ_rms 来凑。CBC 上实测
%   这个 G 的不确定度直接让 γ 差 1–2 倍(见 ct_model.m 文件头, benchmark_cbc/kpar_models.md)。
%   **标准做法**是回到气球模表示(ballooning representation), 把平行方向当微分算子解出来:
%
%       k∥ = −i (1/(qR)) ∂/∂θ        (k∥ 是**算子**, 不是数)
%       k⊥²(θ) = k_θ²[1 + (ŝθ − α sinθ)²]          (s-α 度规)
%       g_s(θ)  = cosθ + (ŝθ − α sinθ) sinθ          (漂移极向结构)
%
%   本征值问题于是是 θ 上的**二阶微分方程**(Schrödinger 型边值问题), 边界 φ(±θ_max)=0。
%   模的平行结构、有效 k∥、模宽都由方程自己给出 —— 不再需要任何 ad hoc 因子。
%
% 【方程(用 c_s/a 归一的无量纲频率 Ω = ω/(c_s/a))】
%   物种 s 的非绝热响应 u_s(θ,x⊥,x∥) 满足(与 ct_dispersion.m 的记账逐项对应)
%       (Ω − Ω_Ds(θ,x) − i·A_par,s·x∥·∂_θ) u_s = J_0(k⊥(θ)ρ_s x⊥)·φ(θ)        (GK)
%   其中  Ω_Ds = −(ω_D0s/(c_s/a))·g_s(θ)·(2x∥² + x⊥²),
%         A_par,s = (v_ts/c_s)·(a/(qR(θ))),  v_ti = √(2T_i/m_i), v_te = √(2T_e/m_e),
%         x∥ 以各物种热速度为单位;  c_i = 1/τ, c_e = 1。
%   准中性:
%       Σ_s c_s [ Γ_0(b_s(θ))·φ(θ) − Σ_{p,q} w_p w_q J_0(k⊥ρ_s x⊥_p)·u_s ] = 0
%
% 【怎么解 —— 场矩阵行列式 + 复平面绕数求根(而不是"广义矩阵本征值问题")】
%   ① 对**给定**的 Ω 解 GK 方程(θ 上是**带状**线性方程组, 每个速度节点一次):
%        u_{s,p,q} = G_{s,p,q}(Ω)·[J_0(·,p) ∘ φ],  G = (Ω I − diag(Ω_D) + A_par x∥ (iD_θ))^{−1}
%   ② 代进准中性得 **M(Ω) φ = 0**, M(Ω) 是 θ 节点上的 Mth×Mth **场矩阵**
%        M(Ω) = diag(Q_r) − Σ_{s,p,q} c_s w_p w_q (Ω − Ω_*s − Ω_*s η_s χ_pq)
%                           · diag(J_0(·,p)) G_{s,p,q}(Ω) diag(J_0(·,p))
%        Q_r = Σ_s c_s Γ_0(b_s(θ_r)) + [绝热电子: 1]
%   ③ 本征值 = det M(Ω) = 0 的根(实正因子归一), 交给 **ct_roots 的函数句柄接口**
%      (绕数 + 自适应加密 + 多起点 Newton ⇒ P0-②③ 新增的高精度搜索)。
%
%   【为什么不用"把 (u,φ) 拼成矩阵再 eigs"】消去 φ 之后 B = I + 低秩 **接近奇异**, 谱里因此
%   混入大量 |Ω| ≫ ω_* 的伪根(实测 |Ω| ~ 40 c_s/a 而物理模 |Ω| ~ 0.2), 而 eigs 只能按"虚部
%   最大"取 ⇒ 取到的全是伪根。场矩阵行列式 + 有界盒 + 三项物理判据能干净地排除它们。
%
% 【怎么把物理根从伪根里挑出来(三项判据, 缺一不可)】
%   本模块的速度空间是**数值求积**(Gauss-Legendre × Gauss-Hermite), 它把响应变成 Ω 的**有理
%   函数**, 于是复平面上会出现大量"离散化伪根"(每个速度节点贡献一个极点; 这是求积口径的固有
%   性质, 与解析 x∥ 积分的局域模型不同)。三项判据:
%     ① **θ 局域化**: 峰值必须在 θ≈0 附近、两端 |φ| 小(排除**边界**伪根);
%     ② **速度空间平滑**(res.uconc = max|u|²/mean|u|²): 物理响应是平滑 φ 的响应 ⇒ O(10);
%        伪根的 u 集中在个别速度节点上 ⇒ O(10²)以上。这一项排除**内部**伪根;
%     ③ **速度网格一致**: 粗网格搜索得到的根, 必须在**细网格**上做局部 Newton 后仍落在 2% 内。
%   三者都过才标 ★物理根; 一个都没有时本函数**明确告警并不给结论**(见文件末尾),
%   而不是把伪根当结果返回。
%
% 【为什么这才是"标准做法"】
%   ① k∥ 不再可调: 它是算子, 平行结构由方程决定; ct_model 里"G 的三种口径"在这里**不存在**。
%   ② 本函数**测出**本征模的有效平行波数 k∥^rms = √(⟨|∂_θφ|²⟩/⟨|φ|²⟩)/(qR) —— 可用于检验
%      局域模型 k∥ 口径的量级(见 ballooning.md), 而不是拿它当拟合参数。
%   ③ 高 k 端的 γ 截止自然出现(平行结构与 FLR 都随 θ 变化)。
%
% 【本版覆盖范围(务必读)】
%   * 离子: 恒为动理学(与 ct_dispersion 同)。
%   * 电子: opt.elec='adiabatic'(当 cfg.passing_e='adiabatic') | 'kinetic'(通行电子, 含平行
%     朗道共振 + 曲率漂移)。'kinetic' **等价于局域模型里 f_t = 0**。
%     ⚠ **捕获电子(弹跳平均)通道本版未含**(它是 θ 上的另一种非局域性, 见 ct_trapped_avg.m /
%       README §7.3); 碰撞同理未含。做法: 用局域模型比较 f_t=ε 与 f_t=0 的差来**估计**其贡献。
%   * 只有 (k_yρ_s, 剖面) 的**局域**气球模问题: 无径向模结构、无 E×B 剪切、无电磁效应。
%
% 【用法】
%   P = ct_read_profile('profile_cbc.txt');  cfg = ct_configure(CT_toroidal_defaults());
%   G = ct_geometry(P, cfg);
%   [res, Dg] = ct_balloon(G, cfg, 0.2828);                 % k_y ρ_s = 0.2828
%   res(1).gamma_n, res(1).omega_r_n, res(1).kpar_rms      % c_s/a 与 [1/m]
%   plot(res(1).th, abs(res(1).phi))                       % 模的平行结构
%   idx = find([res.physical], 1);                          % ★物理根(判据见上)
%
% opt 字段(均可缺省)
%   Nth      61        θ 网格节点数(含两端; 两端 φ=0 直接消去)
%   thmax    4*pi      θ 域半宽(要 ≫ 模宽, 也**不要过大** —— 见"物理根筛选")
%   Nx       8         搜索用 x⊥ 节点数;  Npar  12  搜索用 x∥ 节点数
%   Nx_fine  16        报告用(细) x⊥ 节点数; Npar_fine 24 报告用 x∥ 节点数
%   par_on   true      false = 去掉 ∂_θ 项 ⇒ 逐 θ 解耦的局域问题并集(交叉检验用)
%   elec     ''        ''(由 cfg.passing_e 决定) | 'adiabatic' | 'kinetic'
%   elec_drift true    'kinetic' 电子是否含曲率漂移(局域模型的通行通道**不含**; 量级 ~1%)
%   kperp_fac true     false = k⊥ 不含 (ŝθ−αsinθ) 因子(交叉检验用)
%   Rloc     false     true = 平行连接长度用局域 R(θ)=R0+r cosθ(默认 R0)
%   re_fac   4         搜索盒实部半宽 = re_fac × Ω_*T;  im_fac 3 虚部上界
%   edge_frac 0.3      判据①: 两端 |φ|/max|φ| 上限;  thrms_max 3.0 判据①: θ_rms 上限
%                      ★ 物理气球模**本来就比较宽**(θ_rms ~ 1.5–3; 平行输运把模在 θ 上摊开),
%                        所以 θ_rms 上限不能设太小 —— 实测 CBC 物理模 θ_rms ≈ 2–3。
%   uconc_max 30       判据②: 速度空间集中度上限(max|u|²/mean|u|²; 物理模 ~7–15, 伪根 ≳ 60)
%   thpeak_max 0.6     判据①: |θ_peak| 上限
%   refitol  0.05      判据③: 粗细网格根位置相对差上限(粗网格本身有 ~4% 求积误差, 见文件头)
%   maxroots 200       ★ 交给 ct_roots 的返回上限: **必须给足**。离散化伪根常常比物理根
%                       更不稳定(实测伪根 γ ~ 0.4–1.5 而物理根 γ ~ 0.07), 若按 γ 截断到几十个,
%                       物理根会**被截掉** —— 这正是本模块一度"找不到物理根"的直接原因
%                       (诊断显示 171 个候选中物理根排在 γ 序 48 名之后)。
%                       诊断(每次 ~10 ms)对 200 个候选只花约 2 s, 不构成代价。
%   nre/nim/refine/nseed/expand/maxexpand/maxroots   交给 ct_roots(见该文件)
%   verbose  true      打印候选根与诊断
%   wprobe   []        ★ 只做**场矩阵行列式探针**: 给定一组 Ω, 返回 |det M(Ω)| 后立即返回
%                      (Dg.fprobe)。用于检验组装(D 的零点位置)而不依赖搜根。
%=========================================================================================

if nargin < 4 || isempty(opt), opt = struct(); end
def = struct('Nth',61, 'thmax',4*pi, 'Nx',8, 'Npar',12, 'Nx_fine',16, 'Npar_fine',24, ...
             'par_on',true, 'elec','', 'elec_drift',true, 'kperp_fac',true, 'Rloc',false, ...
             're_fac',4, 'im_fac',3, 'edge_frac',0.3, 'thrms_max',3.0, 'uconc_max',30, ...
             'thpeak_max',0.6, 'refitol',0.05, 'expand',1.6, 'maxexpand',2, 'maxroots',200, ...
             'nre',64, 'nim',32, 'refine',2, 'nseed',24, 'verbose',true, 'wprobe',[]);
fn = fieldnames(def);
for ii = 1:numel(fn)
    if ~isfield(opt, fn{ii}) || isempty(opt.(fn{ii})), opt.(fn{ii}) = def.(fn{ii}); end
end
if isempty(opt.elec)
    opt.elec = 'adiabatic';
    if isfield(cfg,'passing_e') && ~isempty(cfg.passing_e) && strcmpi(cfg.passing_e,'kinetic')
        opt.elec = 'kinetic';
    end
end
t0 = tic;

% =======================================================================================
% [1] 粗网格(搜索)与细网格(报告)两套场矩阵数据
% =======================================================================================
Dc = bl_build(G, cfg, k_rhos, opt, opt.Nx,      opt.Npar);
if opt.par_on && Dc.Mth < 11
    error('ct_balloon:badGrid', ['opt.par_on=true 时内部 θ 节点数至少 11 —— 4 阶差分在 j=2 与 ' ...
          'j=Nth−1 处要用不引用端点的单侧公式(见 ct_d1_reduced); 诊断用 par_on=false 可取小值。']);
end
use_fine = (opt.Nx_fine ~= opt.Nx) || (opt.Npar_fine ~= opt.Npar);
if use_fine
    Df = bl_build(G, cfg, k_rhos, opt, opt.Nx_fine, opt.Npar_fine);
else
    Df = Dc;
end

% =======================================================================================
% [1b] opt.wprobe: 只求**场矩阵行列式**在给定 Ω 处的值然后返回(供自检/诊断用)
%   ★ 这是检验"组装是否与局域模型一致"最直接的量: 取 Ω = 局域解析根, 则 par_on=false 时
%     |det M(Ω)| 必须 ~0(相差几个量级才算通过)。不依赖任何搜根, 因此不会被"伪根"干扰。
if ~isempty(opt.wprobe)
    fq = zeros(size(opt.wprobe));  rhos = fq;  Mdbg = {};
    for iq = 1:numel(opt.wprobe)
        [fq(iq), ~, ~, Mdbg{iq}, rhos(iq)] = bl_detfun(opt.wprobe(iq), Dc); %#ok<AGROW>
    end
    res = struct([]);
    Dg = struct('Mth',Dc.Mth, 'osc',Dc.osc, 'fprobe',fq, 'wprobe',opt.wprobe, ...
                'rho',rhos, 'Mdbg',{Mdbg}, ...
                'Nx',opt.Nx, 'Npar',opt.Npar, 'nroots',0, 'nphys',0, ...
                'nloc',0, 'nconv',0, 'nD',0, 't',toc(t0));
    return;
end

% =======================================================================================
% [2] 搜索盒: **物理上有界**(|Ω| ~ 几倍 Ω_*T), 排除 |Ω| ≫ ω_* 的离散化伪根
% =======================================================================================
box = struct('re', [-opt.re_fac, opt.re_fac]*Dc.osc, 'im', [1e-3, opt.im_fac]*Dc.osc);

% =======================================================================================
% [3] 求根: det M(Ω) = 0(ct_roots 的函数句柄接口; 用粗网格以控制代价)
% =======================================================================================
oroots = struct('nre',opt.nre, 'nim',opt.nim, 'refine',opt.refine, 'nseed',opt.nseed, ...
                'expand',opt.expand, 'maxexpand',opt.maxexpand, 'maxroots',opt.maxroots, ...
                'verbose',false);
R = ct_roots(@(w) bl_detfun(w, Dc), box, oroots);
nD1 = R(1).diag.nD;

% =======================================================================================
% [4] 细网格局部 Newton 修正 ⇒ 判据③; 并给出模结构 / 有效 k∥ / 速度空间集中度
% =======================================================================================
th = Df.th;  Mth = Df.Mth;  dth = Df.dth;
res = struct('omega',{}, 'gamma',{}, 'omega_r',{}, 'gamma_n',{}, 'omega_r_n',{}, ...
             'omega_c',{}, 'phi',{}, 'th',{}, 'localized',{}, 'converged',{}, ...
             'th_peak',{}, 'th_rms',{}, 'kpar_rms',{}, 'edge_ratio',{}, 'uconc',{}, ...
             'physical',{}, 'rankdef',{}, 'Dloc_peak',{}, 'Dloc_theta0',{});
if ~isempty([R.omega])
    om = [R.omega];
    [~, ord] = sort(imag(om), 'descend');
    wgt = ones(Mth,1)*dth;
    for kk = 1:numel(ord)
        zc = om(ord(kk));
        [zf, okn] = bl_newton(@(w) bl_detfun(w, Df), zc, 12, box);
        if ~okn, zf = zc; end
        rel = abs(zf - zc)/max(abs(zc), realmin);
        [~, phi, rk] = bl_detfun(zf, Df);
        mx = max(abs(phi));
        if mx > 0, phin = phi/mx; else, phin = phi; end
        e2     = sum(wgt.*abs(phin).^2);
        th_rms = sqrt(sum(wgt.*abs(phin).^2.*th.^2)/max(e2,realmin));
        dphi   = Df.D1*phin;
        kp     = sqrt(sum(wgt.*abs(dphi).^2)/max(e2,realmin))/(Df.q*Df.R0);
        [~, jm]= max(abs(phin));
        edge   = max([abs(phin(1:min(3,Mth))); abs(phin(max(1,Mth-2):Mth))]);
        uc     = bl_uconc(zf, Df, phin);
        i = numel(res)+1;
        res(i).omega     = zf*Df.cs_a;
        res(i).gamma     = imag(zf)*Df.cs_a;
        res(i).omega_r   = real(zf)*Df.cs_a;
        res(i).gamma_n   = imag(zf);
        res(i).omega_r_n = real(zf);
        res(i).omega_c   = zc;                       % 粗网格上的位置(供核对)
        res(i).phi       = phin;
        res(i).th        = th;
        res(i).th_peak   = th(jm);
        res(i).th_rms    = th_rms;
        res(i).kpar_rms  = kp;
        res(i).edge_ratio= edge;
        res(i).uconc     = uc;
        res(i).rankdef   = rk;
        res(i).localized = (edge < opt.edge_frac) && (th_rms < opt.thrms_max) && ...
                           (abs(th(jm)) < opt.thpeak_max);
        res(i).converged = okn && (rel < opt.refitol);
        res(i).physical  = res(i).localized && res(i).converged && (uc < opt.uconc_max) && ...
                           (imag(zf) > 0);
        res(i).Dloc_peak   = ct_disp_loc(G, cfg, k_rhos, zf*Df.cs_a, th(jm));
        res(i).Dloc_theta0 = ct_disp_loc(G, cfg, k_rhos, zf*Df.cs_a, 0);
    end
end

% ---- 排序: ★物理根 → 局域化候选 → 其余(各自内部按 γ 降序) ----
%   ★ 必须**先做诊断再排序**: 唯一能区分"物理根"的信息(θ 局域化、速度空间集中度)只有诊断里才有,
%     所以不能按 γ 截断(离散化伪根常常比物理根更不稳定, 实测 γ 可差一个量级)。
if ~isempty(res)
    i1 = find([res.physical]);
    i2 = find(~[res.physical] & [res.localized]);
    i3 = find(~[res.localized]);
    res = res([i1, i2, i3]);
end

Dg = struct('Mth',Mth, 'Nx',opt.Nx, 'Npar',opt.Npar, 'Nx_fine',opt.Nx_fine, ...
            'Npar_fine',opt.Npar_fine, 'Nth',Df.Nth, 'thmax',opt.thmax, 'par_on',opt.par_on, ...
            'elec',opt.elec, 'kperp_fac',opt.kperp_fac, 'elec_drift',opt.elec_drift, ...
            'osc',Dc.osc, 'nD',nD1, 'nroots',numel(res), 'nloc',sum([res.localized]), ...
            'nconv',sum([res.converged]), 'nphys',sum([res.physical]), 't',toc(t0));

if opt.verbose
    fprintf(['    [ct_balloon] k_yρ_s=%.4f: det M(Ω)=0; 搜索网格 θ %d × x⊥ %d × x∥ %d, 报告网格 ' ...
             'x⊥ %d × x∥ %d × 物种 %d; 盒 |ReΩ|≤%.3g 0<ImΩ≤%.3g; 粗网格求值 %d; %.1f s; 物理根 %d/%d\n'], ...
            k_rhos, Dc.Mth, opt.Nx, opt.Npar, opt.Nx_fine, opt.Npar_fine, Dc.Ns, ...
            box.re(2), box.im(2), nD1, Dg.t, Dg.nphys, numel(res));
    fprintf('      %-3s %-12s %-12s %-8s %-8s %-8s %-9s %-8s %-6s %s\n', ...
            '#','γ[c_s/a]','ω_r[c_s/a]','θ_peak','θ_rms','两端比','k∥^rms','速度集中','秩亏','判定');
    for k = 1:min(numel(res), 12)
        tag = '非局域(疑伪根)';
        if res(k).physical, tag = '★物理根';
        elseif res(k).localized
            if res(k).uconc >= opt.uconc_max, tag = '速度空间尖峰(伪根)';
            elseif ~res(k).converged,         tag = '粗细网格不一致'; end
        end
        fprintf('      %-3d %-12.5f %-12.5f %-8.3f %-8.3f %-8.4f %-9.4g %-8.1f %-6d %s\n', ...
                k, res(k).gamma_n, res(k).omega_r_n, res(k).th_peak, res(k).th_rms, ...
                res(k).edge_ratio, res(k).kpar_rms, res(k).uconc, res(k).rankdef, tag);
    end
end
if Dg.nphys == 0
    warning('ct_balloon:noPhysical', ...
        ['k_yρ_s=%.4f: %d 个候选根里没有一个同时满足"θ 局域 + 速度空间平滑 + 粗细网格一致" ' ...
         '⇒ 物理根未被可靠识别。\n' ...
         '  本模块速度空间为**数值求积**, 响应是 Ω 的有理函数, 复平面上因此有大量离散化伪根。\n' ...
         '  可试: 提高 opt.Nx_fine/Npar_fine, 调小 opt.thmax(把不可信的 |θ| 区切掉), 或\n' ...
         '  放宽 opt.edge_frac/uconc_max 并人工检查 φ(θ) 与 u 的分布。**结论请勿直接引用**。'], ...
        k_rhos, numel(res));
end
end

% =======================================================================================
% 局部函数
% =======================================================================================
function D = bl_build(G, cfg, k_rhos, opt, Nx, Npar)
%BL_BUILD  组装场矩阵所需的一切: θ 网格几何 / 速度网格 / 物种参数 / 权重
%   (搜索用粗网格与报告用细网格各调一次; 除速度节点外完全相同。)
M     = ct_model(G, cfg, k_rhos);
alpha = M.alpha;
q = G.q;  f0 = G.cs_a;  tau = G.tau;

Nth = opt.Nth;
if Nth < 3, error('ct_balloon:badGrid', 'opt.Nth 至少 3'); end
th_full = linspace(-opt.thmax, opt.thmax, Nth)';
dth     = th_full(2) - th_full(1);
th      = th_full(2:Nth-1);
Mth     = numel(th);

sheb = G.shear*th - alpha*sin(th);
if opt.kperp_fac, kfac = sqrt(1 + sheb.^2); else, kfac = ones(size(th)); end
kperp = M.kperp*kfac;
G0i = M.adiab((kperp*G.rho_i).^2/2);
G0e = M.adiab((kperp*G.rho_e).^2/2);
lam = cos(th) + sheb.*sin(th);
if opt.Rloc, Rpar = G.R0 + G.r*cos(th); else, Rpar = G.R0*ones(size(th)); end
if opt.par_on
    D1 = ct_d1_reduced(Nth, dth);
else
    D1 = sparse(Mth, Mth);
end

[yp, wyp] = ct_gaussleg(Nx, 0, 6);
wperp = 2*yp(:).*exp(-yp(:).^2).*wyp(:);   wperp = wperp/sum(wperp);
[xq, wq] = ct_gh_norm(Npar);
elec_kin = strcmpi(opt.elec,'kinetic');
aa = G.P.a;

D = struct();
D.Mth = Mth;  D.Nth = Nth;  D.th = th;  D.dth = dth;  D.D1 = D1;  D.T = 1i*D1;
D.q = q;  D.R0 = G.R0;  D.cs_a = f0;
D.Ns = 1 + elec_kin;  D.Nx = Nx;  D.yp = yp(:);  D.wperp = wperp;
D.Rpar = Rpar(:);  D.Rloc = opt.Rloc;  D.elec_kin = elec_kin;
D.J0 = cell(1,D.Ns);  D.G0 = cell(1,D.Ns);  D.gam = cell(1,D.Ns);  D.coef = zeros(1,D.Ns);
D.c = zeros(1,D.Ns);  D.X = zeros(1,D.Ns);  D.Y = zeros(1,D.Ns);
D.c(1) = 1/tau;  D.X(1) = M.wstar_i/f0;  D.Y(1) = D.X(1)*G.eta_i;
D.G0{1} = G0i;   D.J0{1} = besselj(0, kperp*G.rho_i*yp(:)');
D.gam{1} = -(M.wD0_i/f0)*lam;
D.coef(1) = (G.vti/G.cs)*aa/q;
if elec_kin
    D.c(2) = 1.0;  D.X(2) = M.wstar_e/f0;  D.Y(2) = D.X(2)*G.eta_e;
    lme = lam;  if ~opt.elec_drift, lme = zeros(size(lam)); end
    D.G0{2} = G0e;  D.J0{2} = besselj(0, kperp*G.rho_e*yp(:)');
    D.gam{2} = -(M.wD0_e/f0)*lme;
    D.coef(2) = (G.vte/G.cs)*aa/q;
end
D = bl_setvel(D, xq, wq);

% Q_r = Σ_s c_s Γ_0(b_s(θ_r)) + [绝热电子: 1]
%   ★ 绝热电子那一项必须是 **1** 而不是 Γ_0(b_e): 电子的局域响应是 𝒥_e = Γ_0(b_e) − 1
%     (与 ct_dispersion.m 的 𝒥_e = (1−f_t)(Γ_0(b_e)−1+J_ps) + f_t𝒥_t 在 f_t=0、
%      passing_e='adiabatic' 时一致), 代进 (Γ_0(b_e) − 𝒥_e) 后两个 Γ_0(b_e) **精确相消**,
%      只剩 Boltzmann 的 1。漏掉它等于把电子稳定化整个删掉 —— 实测 γ 虚增两个量级。
Qr = zeros(Mth,1);
for s = 1:D.Ns, Qr = Qr + D.c(s)*D.G0{s}; end
if ~elec_kin, Qr = Qr + 1; end
D.Qr = Qr;
D.osc = max(abs(D.X(1))*(1+G.eta_i), abs(M.wstar_e/f0)*(1+G.eta_e));
end

% ---------------------------------------------------------------------------------------
function [f, phi, rankdef, Mraw, rho] = bl_detfun(w, D)
%BL_DETFUN  场矩阵 M(Ω) 的行列式(可选给出零空间向量 φ 与秩亏量)
%   det 按 max|M| 归一: 实正因子既不改变零点也不改变辐角, 但避免大 Mth 时上溢/下溢
%   (Mth = 100 时原始 det 可跨几百个量级, 会让 ct_roots 的绕数判据失效)。
Mth = D.Mth;
M = zeros(Mth, Mth);
for r = 1:Mth, M(r,r) = D.Qr(r); end
Is = speye(Mth);
for s = 1:D.Ns
    Js = D.J0{s};  gs = D.gam{s};
    for p = 1:D.Nx
        jc = Js(:,p);  y2 = D.yp(p)^2;
        ojj = jc*jc.';
        for qi = 1:D.Nq
            xv = D.xq(qi);
            m = (p-1)*D.Nq + qi;
            Mat = spdiags(w - gs*(2*xv^2 + y2), 0, Mth, Mth) + ...
                  spdiags(D.cq{s}(:,qi), 0, Mth, Mth)*D.T;
            G = full(Mat \ Is);          % ★ 必须 full: 见下面 ojj.*G 的说明
            al = D.c(s)*D.wpq(m)*(w - D.X(s) - D.Y(s)*D.chi(m));
            % 【MATLAB 陷阱, 本模块踩过一次】ojj 是稠密, 而 Mat\Is 返回的是**稀疏**矩阵;
            %   `稠密 .* 稀疏` 的结果是**稀疏**的, 于是"只保留稀疏那一侧的非零位置"。
            %   par_on=false 时 Mat 是对角阵 ⇒ G 只有对角非零 ⇒ ojj.*G 把**全部非对角元丢掉了**,
            %   场矩阵退化成对角阵(det 与 Ω 几乎无关)。修法: 先把 G 转成 full。
            M = M - al*(ojj.*G);
        end
    end
end
rho = max(abs(M(:)));
Mraw = M;
if ~(rho > 0), f = 0; phi = zeros(Mth,1); rankdef = NaN; return; end
Ms = M/rho;
f = det(Ms);
phi = zeros(Mth,1);  rankdef = NaN;
if nargout > 1
    [~, Sv, Vv] = svd(Ms);
    sv = diag(Sv);
    phi = Vv(:,end);
    rankdef = sum(sv < 1e-10*max(sv));
end
end

% ---------------------------------------------------------------------------------------
function uc = bl_uconc(w, D, phi)
%BL_UCONC  速度空间集中度 max_m|u_m|² / mean_m|u_m|²(m 遍历 (θ,p,q) 全部节点)
%   物理根的 u(θ,x⊥,x∥) 是平滑 φ 的响应 ⇒ O(1–30); 离散化伪根的 u 集中在个别速度节点上
%   ⇒ O(N_vel) ~ 10²–10³。这是区分"内部伪根"与物理根最干净的判据(θ 局域化判据只管边界)。
Mth = D.Mth;  amax = 0;  asum = 0;  n = 0;
for s = 1:D.Ns
    Js = D.J0{s};  gs = D.gam{s};
    for p = 1:D.Nx
        jc = Js(:,p);  y2 = D.yp(p)^2;  rhs = jc.*phi;
        for qi = 1:D.Nq
            xv = D.xq(qi);
            Mat = spdiags(w - gs*(2*xv^2 + y2), 0, Mth, Mth) + ...
                  spdiags(D.cq{s}(:,qi), 0, Mth, Mth)*D.T;
            u = Mat \ rhs;
            a2 = abs(u).^2;
            amax = max(amax, max(a2));  asum = asum + sum(a2);  n = n + numel(a2);
        end
    end
end
if asum > 0, uc = amax/(asum/n); else, uc = NaN; end
end

% ---------------------------------------------------------------------------------------
function [z, ok] = bl_newton(fun, z0, itmax, box)
%BL_NEWTON  复数 Newton(det 没有解析导数 ⇒ 中心差分); 把粗网格根修正到细网格。
%   ★ 必须带**盒约束**: 否则细网格 Newton 会跑到物理盒之外(实测出现过 |Ω| = 4.7 而盒只到 1.85),
%     于是"报告"的根其实是一个远离物理区的伪根。步长既受盒限制也受单步上限限制。
z = z0;  ok = false;
if nargin < 4 || isempty(box)
    box = struct('re',[-inf inf], 'im',[-inf inf]);
end
steplim = 0.5*max(1, abs(z0));
for it = 1:itmax
    f0 = fun(z);
    if ~isfinite(f0), return; end
    if abs(f0) == 0, ok = true; return; end
    h = 1e-6*max(1e-3, abs(z));
    d = (fun(z + h) - fun(z - h))/(2*h);
    if ~isfinite(d) || abs(d) < 1e-300, return; end
    dz = -f0/d;
    if abs(dz) > steplim, dz = dz*steplim/abs(dz); end
    zn = z + dz;
    if real(zn) < box.re(1) || real(zn) > box.re(2) || ...
       imag(zn) < box.im(1) || imag(zn) > box.im(2)
        return;                                  % 出盒 ⇒ 不算收敛(调用者回退到粗网格根)
    end
    z = zn;
    if abs(dz) < 1e-12*max(1, abs(z)), ok = true; return; end
end
end

% ---------------------------------------------------------------------------------------
function D = bl_setvel(D, xq, wq)
%BL_SETVEL  由 x∥ 节点/权重生成与速度网格有关的量: wD / chi / 权重 wpq / 平行系数 cq
xq = xq(:)';  wq = wq(:)';  Nq = numel(xq);
Nx = D.Nx;  yp = D.yp;  wperp = D.wperp;
chi = (xq.^2)'*ones(1,Nx) + ones(Nq,1)*(yp(:).^2)' - 1.5;
wpq = zeros(1, Nx*Nq);
for p = 1:Nx, wpq((p-1)*Nq + (1:Nq)) = wperp(p)*wq; end
D.Nq = Nq;  D.xq = xq;  D.wq = wq;  D.wpq = wpq;  D.chi = reshape(chi, 1, Nx*Nq);
D.wD = cell(1,D.Ns);  D.cq = cell(1,D.Ns);
for s = 1:D.Ns
    fac2 = reshape(2*(xq(:).^2)*ones(1,Nx) + ones(Nq,1)*(yp(:).^2)', 1, Nx*Nq);
    D.wD{s} = D.gam{s} * fac2;
    if D.Rloc
        D.cq{s} = (D.coef(s)./D.Rpar(:))*xq;
    else
        D.cq{s} = (D.coef(s)/D.Rpar(1))*ones(D.Mth,1)*xq;
    end
end
end

% ---------------------------------------------------------------------------------------
function D = ct_disp_loc(G, cfg, k_rhos, omega, theta)
%CT_DISP_LOC  在指定 θ 处用**局域模型**(ct_dispersion)算 D(ω) —— 纯诊断对照(ω 为 SI)
P2 = G.P;  P2.theta = theta;
G2 = ct_geometry(P2, cfg);
M2 = ct_model(G2, cfg, k_rhos);
D  = ct_dispersion(M2, omega, cfg.disp_model);
end

% ---------------------------------------------------------------------------------------
function D1 = ct_d1_reduced(Nth, h)
%CT_D1_REDUCED  ∂/∂θ 的 4 阶差分矩阵(只含内部节点 2..Nth−1)
%   两端值恒为 0, 因此用**不引用端点**的单侧公式; 三组系数之和都精确为 0 ⇒ D1·(常数)= 0,
%   这正是"无剪切、无曲率时 φ = const 必须是精确解"所要求的。
if Nth < 11, error('ct_d1_reduced: Nth 至少 11'); end
Df = zeros(Nth);
for j = 3:Nth-2
    Df(j, j-2) =  1/(12*h);   Df(j, j-1) = -8/(12*h);
    Df(j, j+1) =  8/(12*h);   Df(j, j+2) = -1/(12*h);
end
Df(2, 2:6)             = [-25, 48, -36, 16, -3]/(12*h);
Df(Nth-1, Nth-5:Nth-1) = [3, -16, 36, -48, 25]/(12*h);
D1 = sparse(Df(2:Nth-1, 2:Nth-1));
end

% ---------------------------------------------------------------------------------------
function [x, w] = ct_gh_norm(N)
%CT_GH_NORM  归一化 Gauss-Hermite 节点与权重: ∫dx e^{−x²}f(x) ≈ Σ_q w_q f(x_q), Σ_q w_q = 1
if N < 2, error('ct_gh_norm: N 至少 2'); end
b = sqrt((1:N-1)/2);
J = diag(b,1) + diag(b,-1);
[Vv, Dd] = eig(J);
x = diag(Dd);
w = (Vv(1,:)').^2;   w = w/sum(w);
[x, ix] = sort(x);   w = w(ix);
end
