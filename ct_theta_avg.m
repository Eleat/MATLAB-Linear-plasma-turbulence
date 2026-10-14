function [res, Dg] = ct_theta_avg(G, cfg, k_rhos, opt)
%CT_THETA_AVG  [v8.5] **角度平均**(通量面平均)的色散根: ⟨D(ω;θ)⟩_θ = 0
%=========================================================================================
% 三种"把 θ 结构处理掉"的办法, 由弱到强:
%   ① 单 θ 局域(ct_model/ct_dispersion): θ 取一个值(默认 0), k∥ 用几何因子凑 —— 最省, 但
%      θ 只取一点, 谱形与 ω_r 都有结构性偏差(见 benchmark_cbc/cbc_compare_gx.md §3②)。
%   ② **角度平均(本函数)**: 把 D 在 θ 上按权重平均后求根
%         ⟨D(ω)⟩ = Σ_r w_r D(ω; θ_r) / Σ_r w_r = 0
%      它把 θ 上的驱动/阻尼按权重摊平, 比单点合理; 但**平均把 θ 结构抹掉了** —— 物理上
%      正确的对象是"一个模", 而 D(ω;θ) 的零点是逐 θ 不同的函数, 平均再求根并不等于
%      "模在 θ 上的自洽解"。因此它只是**过渡近似**, 正确定标见下。
%   ③ 气球模本征值(ct_balloon): 保留 k∥ = −i(1/(qR))∂/∂θ, 解二阶本征值问题 —— 标准做法。
%      三者的逐 k 数值对比见 benchmark_cbc/ballooning.md。
%
% 【权重】opt.avg: 'gl'(默认, θ 上的 Gauss-Legendre 节点/权重, 精度高) | 'uniform'
%   θ ∈ [−thmax, thmax]。★ 注意: θ=±θmax 的贡献在权重里非零, 而真实模在那里振幅已很小;
%   这是角度平均固有的(不是实现问题), 也正是它需要与气球模本征值对照的原因。
%
% 【实现】D(ω) 以**函数句柄**交给 ct_roots(v8.5 新增的接口), 并同时给出解析导数的句柄
%   dD/dω = Σ_r w_r ∂D/∂ω(θ_r)/Σ w_r  ⇒ Newton 每次迭代只要 2 次求值(见 ct_roots 说明)。
%
% 【用法】
%   [res, Dg] = ct_theta_avg(G, cfg, 0.2828);                 % k_yρ_s = 0.2828
%   [res, Dg] = ct_theta_avg(G, cfg, 0.2828, struct('avg','uniform','Nth',49));
%   res(1).gamma_n, res(1).omega_r_n      % c_s/a
%
% opt 字段(其余同 ct_roots)
%   Nth     33     θ 节点数
%   thmax   3*pi   θ 域半宽
%   avg     'gl'   'gl' | 'uniform'
%   nre/nim/refine/nseed/expand/maxroots/verbose  见 ct_roots
%=========================================================================================

if nargin < 4 || isempty(opt), opt = struct(); end
def = struct('Nth',33, 'thmax',3*pi, 'avg','gl', 'nre',64, 'nim',32, 'refine',2, ...
             'nseed',24, 'expand',1.6, 'maxexpand',2, 'maxroots',16, 'verbose',false);
fn = fieldnames(def);
for ii = 1:numel(fn)
    if ~isfield(opt, fn{ii}) || isempty(opt.(fn{ii})), opt.(fn{ii}) = def.(fn{ii}); end
end

t0 = tic;
% ---- θ 节点与权重 ----
switch lower(char(opt.avg))
    case 'uniform'
        th = linspace(-opt.thmax, opt.thmax, opt.Nth)';
        wt = ones(size(th));
    otherwise
        [xg, wg] = ct_gaussleg(opt.Nth, -opt.thmax, opt.thmax);
        th = xg(:);  wt = wg(:);
end
wt = wt/sum(wt);

% ---- 逐 θ 的局域模型(几何沿 θ 重建) ----
M = ct_model(G, cfg, k_rhos);                 % 只为拿 ω_* 定搜索盒
Nr = numel(th);
Ms = cell(1, Nr);
for r = 1:Nr
    P2 = G.P;  P2.theta = th(r);
    G2 = ct_geometry(P2, cfg);
    Ms{r} = ct_model(G2, cfg, k_rhos);
end
dm = cfg.disp_model;

% ---- ⟨D(ω)⟩ 与其解析导数(函数句柄接口) ----
Dfun  = @(w) avg_D(Ms, wt, w, dm);
dDfun = @(w) avg_dD(Ms, wt, w, dm);

wsc = max(abs(M.wstar_i)*(1+G.eta_i), abs(M.wstar_e)*(1+G.eta_e));
box = struct('re', [-cfg.root_re, cfg.root_re]*wsc, 'im', [cfg.root_im0, cfg.root_im1]*wsc);
o = struct('nre',opt.nre, 'nim',opt.nim, 'refine',opt.refine, 'nseed',opt.nseed, ...
           'expand',opt.expand, 'maxexpand',opt.maxexpand, 'maxroots',opt.maxroots, ...
           'disp_model',dm, 'verbose',opt.verbose, 'dDfun',dDfun, 'analyt_deriv',true);
R = ct_roots(Dfun, box, o);

res = struct('omega',{}, 'gamma',{}, 'omega_r',{}, 'gamma_n',{}, 'omega_r_n',{}, ...
             'k',{}, 'Dval',{}, 'th',{}, 'wt',{});
if ~isempty(R)
    om = [R.omega];
    % ★ 只保留**增长**根(Im ω > 0): ct_roots 在扩盒时可能把盒下沿带到 Im<0, 于是会带出
    %   阻尼根; 平均值本身没有不稳定根时, 返回空数组比返回一个阻尼根更诚实。
    om = om(imag(om) > 0);
    [~, ord] = sort(imag(om), 'descend');
    for k = 1:numel(ord)
        z = om(ord(k));
        res(k).omega = z;      res(k).gamma = imag(z);  res(k).omega_r = real(z);
        res(k).gamma_n = imag(z)/G.cs_a;   res(k).omega_r_n = real(z)/G.cs_a;
        res(k).k = k_rhos;     res(k).Dval = abs(Dfun(z));
        res(k).th = th;        res(k).wt = wt;
    end
end
Dg = struct('avg',opt.avg, 'Nth',Nr, 'thmax',opt.thmax, 'nD',R(1).diag.nD, ...
            'nroots',numel(R), 'ngrow',numel(res), 't',toc(t0));
if opt.verbose
    fprintf('    [ct_theta_avg] k_yρ_s=%.4f: 平均口径 ''%s''(%d 个 θ 节点, ±%.3fπ), 根 %d 个(增长 %d), 求值 %d, %.1f s\n', ...
            k_rhos, opt.avg, Nr, opt.thmax/pi, numel(R), numel(res), Dg.nD, Dg.t);
    for k = 1:numel(res)
        fprintf('      γ=%.5f  ω=%.5f [c_s/a]  |⟨D⟩|=%.2e\n', res(k).gamma_n, res(k).omega_r_n, res(k).Dval);
    end
    if isempty(res)
        fprintf(['      ⚠ 平均色散关系 ⟨D(ω;θ)⟩_θ = 0 **没有增长根**: 这是角度平均的固有缺陷 ——\n' ...
                 '        它把每个 θ 的**局域**响应(含该 θ 处的全部平行阻尼)按权重相加, 于是\n' ...
                 '        |θ| 大的阻尼区把驱动摊掉。正确做法是解气球模本征值(ct_balloon.m)。\n']);
    end
end
end

% ---------------------------------------------------------------------------------------
function v = avg_D(Ms, wt, w, dm)
%AVG_D  ⟨D(ω;θ)⟩_θ
v = 0;
for r = 1:numel(Ms)
    v = v + wt(r)*ct_dispersion(Ms{r}, w, dm);
end
end

function v = avg_dD(Ms, wt, w, dm)
%AVG_DD  ⟨∂D/∂ω(θ)⟩_θ(解析; pole 模型下 ct_dispersion 返回 NaN ⇒ 由 ct_newton 退回差分)
v = 0;
for r = 1:numel(Ms)
    v = v + wt(r)*ct_disp_deriv1(Ms{r}, w, dm);
end
end

function d = ct_disp_deriv1(M, w, dm)
[~, d] = ct_dispersion(M, w, dm, true);
end
