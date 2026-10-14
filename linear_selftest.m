function ok = linear_selftest(cfg)
%LINEAR_SELFTEST  只看**不稳定性模求解**这一部分的精简自检(不涉及 C_T)
%=========================================================================================
% 本文件夹不包含 ct_selftest.m(它含大量 C_T 项), 因此这里给出与之互补的最小判据集,
% 每一条都对应"求解路径上的一个环节", 且判据**独立重算**(不复用被测代码的内部结论):
%
%   L1  x⊥ 求积权重归一化:      Σ_p w_p = 1 且 Σ_p w_p x⊥² = 1   (麦氏测度自洽)
%   L2  den² 族恒等式:          ∂M₀/∂ω = −M₀b(数值差分独立复核)
%   L3  k∥ → 0 连续过渡:        |D(k∥=0) − D(k∥→0)| / |D| 应 ≪ 1
%   L4  CBC 绝热电子找到 ITG 根: γ > 0 且 Re ω < 0
%   L5  根的残差:               |D(z)| ≲ 1e-12·|D| 尺度
%   L6  网格档位一致性:          fast 与 high 两档给同一根(相对差 < 1e-3)
%   L7  no_curv ⇒ g_i = g_e = 0 (漂移关闭时几何因子必须为零)
%   L8  角度平均的小域极限:      θ∈±0.01 时 ⟨D⟩ 的根 = 局域根(相对差 < 1e-3)
%   L9  场矩阵探针(气球模):      par_on=false 时 |det M| 在局域根处是尖锐极小
%
% 用法:  ok = linear_selftest;                 % 用 CBC 剖面
%        ok = linear_selftest(cfg);            % 用给定配置(需含 profile_file)
% 输出:  通过项数; 每项打印 [PASS]/[FAIL] 与数值。
%=========================================================================================

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(cfg)
    cfg = CT_toroidal_defaults();
    cfg.profile_file = fullfile(here, 'profiles', 'profile_cbc.txt');
end
cfg.plots = false;  cfg.selftest = false;  cfg.verbose = false;
cfg.k_unit = 'rho_i';
if ~isfield(cfg,'ft_model') || isempty(cfg.ft_model), cfg.ft_model = 'fixed'; end
% ★ 强制"绝热电子"口径: f_t = 0 且通行电子 Boltzmann。
%   【为什么必须强制】L9 的场矩阵探针(ct_balloon)目前**只含绝热或纯通行电子**,
%   若这里用 ft_model='eps'(默认, f_t≈0.55)算出的根去探, 两者不是同一个物理问题,
%   探针的极小会明显变浅(实测 |det M| 由 3.9e-4 变成 3.2e-2, 判据会 FAIL)。
cfg.ft_model = 'fixed';  cfg.ft_fixed = 0;      % 绝热电子(与 CBC 基准同口径)
cfg.passing_e = 'adiabatic';
cfg.quiet_geometry_check = true;
if ~isfield(cfg,'kpar_model') || isempty(cfg.kpar_model), cfg.kpar_model = 'zero'; end
cfg = ct_configure(cfg);

fprintf('\n [线性部分自检] %s\n', cfg.profile_file);
P = ct_read_profile(cfg.profile_file);
G = ct_geometry(P, cfg);
n = 0;  ok = 0;
krs = 0.4/sqrt(2);      % ky = k_θρ_i = 0.4 (与完整版 ct_selfcheck 的 D 组同一算例)

% ---- L1 权重归一化 ----
[yp, wyp] = ct_gaussleg(cfg.N_perp, 0, 6);
wp = 2*yp(:).*exp(-yp(:).^2).*wyp(:);
[n,ok] = chk(n, ok, 'L1  x⊥ 权重: Σw = 1 且 Σ w x⊥² = 1', ...
             max(abs([sum(wp)-1, sum(wp.*yp(:).^2)-1])), 0, 1e-12);

% ---- L2 ∂M₀/∂ω = −M₀b ----
M = ct_model(G, cfg, krs);
w0 = (real(M.wstar_i) + 0.7i*abs(M.wstar_i));
h  = 1e-6*max(1, abs(w0));
E1 = ct_channel(M, 'i', w0 + h, true);
E2 = ct_channel(M, 'i', w0 - h, true);
E0 = ct_channel(M, 'i', w0, true);
dnum = (E1.I0 - E2.I0)/(2*h);
[n,ok] = chk(n, ok, 'L2  den² 族恒等式: ∂M₀/∂ω = −M₀b', ...
             abs(dnum + E0.I0b)/max(abs(dnum), realmin), 0, 1e-5);

% ---- L3 k∥ → 0 连续 ----
c0 = cfg;  c0.kpar_model = 'zero';
c1 = cfg;  c1.kpar_model = 'qR';  c1.kpar_factor = 1e-6;      % 极小 k∥
M0 = ct_model(G, c0, krs);  M1 = ct_model(G, c1, krs);
d0 = ct_dispersion(M0, w0);  d1 = ct_dispersion(M1, w0);
[n,ok] = chk(n, ok, 'L3  k∥ → 0 连续: |D(0) − D(ε)|/|D| ≪ 1', ...
             abs(d0 - d1)/max(abs(d0), realmin), 0, 1e-4);

% ---- L4/L5 找到 ITG 根且残差小 ----
Mb = ct_model(G, cfg, krs);
z = local_root(Mb, G);
[n,ok] = chk(n, ok, 'L4  CBC 绝热电子: 存在 γ>0 且 Re ω<0 的 ITG 根', ...
             double(isfinite(z) && imag(z) > 0 && real(z) < 0), 1, 0);
if isfinite(z)
    fprintf('        （局域根 = %.5f %+.5fi  [c_s/a]）\n', real(z)/G.cs_a, imag(z)/G.cs_a);
end
if isfinite(z)
    Ds = ct_dispersion(Mb, z);
    sc = abs(ct_dispersion(Mb, z*1.2));
    [n,ok] = chk(n, ok, 'L5  根残差: |D(z)| ≲ 1e-12·(远处 |D| 作为尺度)', ...
                 abs(Ds)/max(sc, realmin), 0, 1e-12);
else
    n = n + 1;  fprintf('     [SKIP] L5  未找到根\n');
end

% ---- L6 两档精度给同一个根 ----
z1 = root_with_grid(G, cfg, krs, 'fast');
z2 = root_with_grid(G, cfg, krs, 'high');
[n,ok] = chk(n, ok, 'L6  root_grid = fast 与 high 给同一根(相对差 < 1e-3)', ...
             abs(z1 - z2)/max(abs(z2), realmin), 0, 1e-3);

% ---- L7 no_curv ----
cN = cfg;  cN.no_curv = true;
GN = ct_geometry(P, cN);
[n,ok] = chk(n, ok, 'L7  no_curv ⇒ g_i = g_e = 0', ...
             max(abs([GN.lam_i, GN.lam_e])), 0, 0);

% ---- L8 角度平均小域极限 ----
if isfinite(z)
    ra = ct_theta_avg(G, cfg, krs, struct('Nth',3,'thmax',0.01,'avg','uniform','verbose',false));
    if isempty(ra)
        n = n + 1;  fprintf('     [FAIL] L8  角度平均小域极限: 无根\n');
    else
        [n,ok] = chk(n, ok, 'L8  角度平均(±0.01)的根 = 局域根(相对差 < 1e-3)', ...
                     abs(ra(1).omega - z)/abs(z), 0, 1e-3);
    end
else
    n = n + 1;  fprintf('     [SKIP] L8  未找到局域根\n');
end

% ---- L9 场矩阵探针(气球模) ----
%   ★ 探针频率必须用**归一化** Ω = ω/(c_s/a); θ 节点数需 ≥ 11(Mth=1 时归一化 det ≡ ±1);
%     锐度随速度网格增长(实测 Nx=16/Npar=28 时偏离 1% 升 ~10 倍), 故阈值取 0.3 / 0.1。
if isfinite(z)
    zp = z/G.cs_a;
    [~, Dp] = ct_balloon(G, cfg, krs, struct('par_on',false, 'Nth',11, 'Npar',28, 'Nx',16, ...
                       'verbose',false, 'wprobe', zp*[1, 1.01, 1.05, 2]));
    fp = abs(Dp.fprobe);
    % 判据: ① 局域根处 |det M| 是这些探针点里的**最小值**; ② |det M| 在 Ω 翻倍时至少升 5 倍。
    %   (绝对锐度与 k、速度网格有关 —— 实测 ky=0.4/Nx=16 下 ×2 处升 17 倍, 而 +1% 处只升 1.1 倍,
    %    所以判据不能只看近邻点。)
    [n,ok] = chk(n, ok, 'L9  场矩阵: |det M| 在局域根处取极小且 Ω 翻倍时升 > 5 倍', ...
                 double(fp(1) <= 1.05*min(fp) && fp(4)/fp(1) > 5), 1, 0);
    fprintf('        （L9 数值: |det M| = %.2e / %.2e / %.2e / %.2e (根 / +1%% / +5%% / ×2); 锐度 %.2f / %.1f）\n', ...
            fp(1), fp(2), fp(3), fp(4), fp(4)/fp(1), fp(4)/fp(1));
else
    n = n + 1;  fprintf('     [SKIP] L9  未找到局域根\n');
end

fprintf('\n 线性部分自检: 通过 %d / %d\n\n', ok, n);
if ok < n, warning('linear_selftest:fail', '有 %d 项未通过', n - ok); end
end

% =======================================================================================
function z = local_root(M, G)
%LOCAL_ROOT  最大 γ 的不稳定根(SI); 找不到返回 NaN
wsc = max(abs(M.wstar_i)*(1+G.eta_i), abs(M.wstar_e)*(1+G.eta_e));
box = struct('re',[-8 8]*wsc, 'im',[1e-3 6]*wsc);
R = ct_roots(M, box, struct('nre',96,'nim',48,'refine',2,'nseed',32, 'verbose',false));
z = NaN;
if isempty(R), return; end
om = [R.omega];  om = om(imag(om) > 0);
if isempty(om), return; end
[~, j] = max(imag(om));  z = om(j);
end

% ---------------------------------------------------------------------------------------
function z = root_with_grid(G, cfg, krs, gridname)
c = cfg;  c.root_grid = gridname;
M = ct_model(G, c, krs);
z = local_root(M, G);
if ~isfinite(z), z = NaN; end
end

% ---------------------------------------------------------------------------------------
function [n, ok] = chk(n, ok, name, val, expv, tol)
%CHK  打印一项并按 |val−exp| ≤ tol 判定(与 ct_selfcheck 的同一约定)
n = n + 1;
d = abs(val - expv);
pass = (d <= tol) || (isnan(val) && isnan(expv));
if pass, ok = ok + 1;  s = 'PASS'; else, s = 'FAIL'; end
fprintf('     [%s] %-64s 值 = %-13.6g 期望 = %-10.6g 偏差 = %.2e\n', s, name, val, expv, d);
end
