function out = CT_linear(cfg)
%CT_LINEAR  【不稳定性模专用】解局域线性色散关系 D(ω)=0 ⇒ (ω_r, γ)(k_yρ_s)
%=========================================================================================
% 本函数与完整版 CT_toroidal.m **逐行同源**(剖面→几何→逐 k 模型→根搜索→支判定),
% 唯一区别是:**不做 C_T 后处理**(不调用 ct_ct / ct_ct_exact / ct_mode_deriv, 不产出 C_T
% 相关字段), 因此也不需要 c_D、模移动导数那一套。放在独立文件夹里便于单独交付与阅读。
%
% 解的方程(局域单 θ, 静电, 环形推广):
%       D(ω) = [Γ₀(b_i) − 𝒥_i]/τ + [Γ₀(b_e) − 𝒥_e] = 0
%   𝒥_s = ⟨ W_s (ω − ω_*s^T)/(ω − k∥v∥ − ω_Ds) ⟩,  W_s = J₀²(k⊥ρ_s x⊥)×Maxwell 权重
%   x∥ 方向**解析积分**(ct_xpar), 只对 x⊥ 做 Gauss-Legendre 求积。
%
% 支的判定(与完整版一致):
%   Re ω < 0 → ITG 支(离子抗磁方向);  Re ω > 0 → TEM 支(电子抗磁方向); 各取 γ 最大者。
%
% ------------------------------- 用法 -------------------------------
%   cfg = CT_toroidal_defaults();
%   cfg.profile_file = 'profiles\profile_cbc.txt';
%   out = CT_linear(cfg);
%   plot_linear(out);                      % 出标准图(γ、ω_r、γ/|ω_r|、k∥)
%   plot_linear_roots(out, 5);             % 第 5 个 k 上的复平面根谱
%
% ------------------------------- 输出 out -------------------------------
%   out.P, out.G, out.cfg      剖面 / 几何 / 生效配置
%   out.k_in                   用户单位下的波数(与 cfg.k_unit 一致)
%   out.k                      k_yρ_s(内部统一单位)
%   out.ITG, out.TEM           结构体数组, 每个 k 一项; 字段见 rec_empty.m 的说明
%   out.kpar(k)                该 k 的 k∥(来自 ct_model, 与支无关)
%   out.roots{k}               该 k 搜索盒内的**全部**根(复数, SI) —— 供根谱图用
%   out.nD(k), out.nexpand(k)  根搜索的诊断(求值次数 / 自动扩盒次数)
%   out.balloon                气球模本征值(仅当 cfg.balloon_k 非空)
%   out.omega_unit / omega_norm    显示单位标签与换算因子
%
% ------------------------------- 相关脚本 -------------------------------
%   plot_linear.m         标准四联图
%   plot_linear_roots.m   单 k 复平面根谱
%   linear_demo.m         一键演示(CBC 扫描 / k∥ 对比 / α 对比) + 出图
%   linear_selftest.m     只看线性部分的精简自检
%=========================================================================================

if nargin < 1 || isempty(cfg), cfg = CT_toroidal_defaults(); end
cfg = ct_fill_defaults(cfg);
cfg = ct_configure(cfg);                  % 归一化预设/别名 + 枚举校验
t0 = tic;

fprintf('\n');
fprintf('================================================================================\n');
fprintf(' CT_linear : 不稳定性模求解 D(ω)=0  ->  (ω_r, γ)(k_yρ_s)     [无 C_T 后处理]\n');
fprintf('================================================================================\n');

% ---------------- [1] 剖面与几何 ----------------
pfname = ct_profile_path(cfg);
P = ct_read_profile(pfname);
G = ct_geometry(P, cfg);
fprintf(' [1] 剖面: %s\n', P.src_file);
fprintf('     r/a = %.4f   R0 = %.3f m   a = %.3f m   q = %.3f   ŝ = %.3f   θ = %.3f rad\n', ...
        P.rho, P.R0, P.a, P.q, P.shear, P.theta);
fprintf('     R/L_ne = %.3f  R/L_Te = %.3f   τ = %.4f   α = %.4g   f_t = %.4f\n', ...
        P.RL_ne, P.RL_Te, G.tau, G.alpha, G.ft);
fprintf('     k∥ 模型 = ''%s''   k⊥ 口径 = ρ_i/(qR) = %.6g /m   ρ_i = %.4f mm   c_s/a = %.4g rad/s\n', ...
        cfg.kpar_model, G.rho_i/(G.q*G.R0), G.rho_i*1e3, G.cs_a);
if cfg.verbose
    fprintf('     自检: %s\n', ternary(cfg.selftest, '开(ct_selftest 在完整版里, 本文件夹用 linear_selftest)', '关'));
end

% ---------------- [2] 波数扫描 ----------------
kk_in = linspace(cfg.kmin, cfg.kmax, cfg.nk);
switch lower(cfg.k_unit)
    case 'rho_i', kk = kk_in/sqrt(2*G.tau);      % k_θρ_i -> k_yρ_s
    otherwise,    kk = kk_in;
end
[WNORM, WLAB] = ct_omega_unit(G, cfg);
nK = numel(kk);

fprintf('\n [2] 扫描: %s = %.3f .. %.3f (%d 点)  频率单位 = %s\n', ...
        cfg.k_unit, cfg.kmin, cfg.kmax, cfg.nk, WLAB);
if cfg.verbose
    fprintf('     %-9s %-6s %-13s %-13s %-8s %-10s %-10s %s\n', ...
            'k_yρ_s','支','ω_r','γ','γ/|ω_r|','k∥[1/m]','|D|','根来源');
end

ITG = repmat(rec_empty(0), 1, nK);
TEM = repmat(rec_empty(0), 1, nK);
roots_all = cell(1, nK);
kpar_k = nan(1, nK);  nD = zeros(1, nK);  nexp = zeros(1, nK);  nrt = zeros(1, nK);
seeds = [];
gr = ct_root_grid(cfg);

for ik = 1:nK
    k = kk(ik);
    if k <= 0
        ITG(ik) = rec_empty(k);  TEM(ik) = rec_empty(k);
        roots_all{ik} = [];      continue;
    end
    M = ct_model(G, cfg, k);
    kpar_k(ik) = M.kpar;
    wsc = max(abs(M.wstar_i), abs(M.wstar_e));
    box = struct('re', [-cfg.root_re, cfg.root_re]*wsc, ...
                 'im', [cfg.root_im0, cfg.root_im1]*wsc);
    boxes = {box};
    if ~isempty(cfg.root_boxes)
        boxes = cfg.root_boxes;
    elseif ~isempty(cfg.root_scales)
        sc = cfg.root_scales(:)';
        boxes = cell(1, numel(sc));
        for ib = 1:numel(sc)
            boxes{ib} = struct('re', [-cfg.root_re, cfg.root_re]*wsc*sc(ib), ...
                               'im', [cfg.root_im0, cfg.root_im1]*wsc*sc(ib));
        end
    end
    opt = struct('nre',gr.nre, 'nim',gr.nim, 'refine',gr.refine, 'nseed',gr.nseed, ...
                 'expand',cfg.root_expand, 'maxexpand',cfg.root_maxexpand, ...
                 'negim',cfg.root_negim, 'imfrac',cfg.root_imfrac, ...
                 'maxroots',cfg.root_maxroots, 'seeds',seeds, 'boxes',{boxes}, ...
                 'converge',cfg.root_converge, 'maxconv',cfg.root_maxconv, ...
                 'disp_model',cfg.disp_model, 'verbose',cfg.root_verbose);
    R = ct_roots(M, box, opt);
    nD(ik) = R(1).diag.nD;   nexp(ik) = R(1).diag.nexpand;
    if ~isempty(R)
        om = [R.omega];
        roots_all{ik} = om;   nrt(ik) = numel(om);
        mI = find(real(om) < 0);   mE = find(real(om) > 0);
        if ~isempty(mI)
            [~, j] = max(imag(om(mI)));
            ITG(ik) = rec_fill(M, om(mI(j)), 'ITG', WNORM, R(mI(j)).Dval, R(mI(j)).src);
        end
        if ~isempty(mE)
            [~, j] = max(imag(om(mE)));
            TEM(ik) = rec_fill(M, om(mE(j)), 'TEM', WNORM, R(mE(j)).Dval, R(mE(j)).src);
        end
        seeds = om;                       % 沿 k 连续延拓
    else
        roots_all{ik} = [];  seeds = [];
    end

    if cfg.verbose
        if ITG(ik).exists, rec_print(ITG(ik), WLAB); end
        if TEM(ik).exists, rec_print(TEM(ik), WLAB); end
        if ~ITG(ik).exists && ~TEM(ik).exists
            fprintf('     %-9.3f %-6s %-13s %-13s %-8s %-10.4g %-10s %s\n', ...
                    k, '—', '稳定(盒内无 Re ω 任一侧的不稳定根)', '', '', kpar_k(ik), '', '');
        end
    end
end

% ---------------- [3] 可选: 气球模本征值 ----------------
B = [];
if ~isempty(cfg.balloon_k)
    kb = cfg.balloon_k(:)';
    fprintf('\n [3] 气球模本征值(ct_balloon)  %d 个 k\n', numel(kb));
    fprintf('     %-9s %-13s %-13s %-10s %-9s %s\n', 'k_yρ_s','γ[c_s/a]','ω_r[c_s/a]','k∥^rms','θ_rms','判定');
    B = repmat(struct('k',0,'gamma_n',NaN,'omega_r_n',NaN,'kpar_rms',NaN, ...
                      'th_rms',NaN,'physical',false,'phi',[],'th',[]), 1, numel(kb));
    for ib = 1:numel(kb)
        o = struct('Nth',cfg.balloon_Nth, 'Npar',cfg.balloon_Npar, 'thmax',cfg.balloon_thmax, ...
                   'elec',cfg.balloon_elec, 'verbose',false, 'maxroots',60);
        rb = ct_balloon(G, cfg, kb(ib), o);
        if isempty(rb), continue; end
        ip = find([rb.physical]);
        if isempty(ip), [~, j] = max([rb.gamma]); else, [~, jj] = max([rb(ip).gamma]); j = ip(jj); end
        B(ib).k = kb(ib);  B(ib).gamma_n = rb(j).gamma_n;  B(ib).omega_r_n = rb(j).omega_r_n;
        B(ib).kpar_rms = rb(j).kpar_rms;  B(ib).th_rms = rb(j).th_rms;
        B(ib).physical = rb(j).physical;  B(ib).phi = rb(j).phi;  B(ib).th = rb(j).th;
        fprintf('     %-9.4f %-13.5f %-13.5f %-10.4g %-9.3f %s\n', kb(ib), rb(j).gamma_n, ...
                rb(j).omega_r_n, rb(j).kpar_rms, rb(j).th_rms, ...
                ternary(rb(j).physical, '★物理根', '局域化候选(见 ballooning 文档)'));
    end
end

% ---------------- [4] 打包 ----------------
out = struct();
out.P = P;  out.G = G;  out.cfg = cfg;
out.k_in = kk_in;  out.k = kk;
out.ITG = ITG;  out.TEM = TEM;
out.kpar = kpar_k;  out.roots = roots_all;
out.nD = nD;  out.nexpand = nexp;  out.nroots = nrt;
out.balloon = B;
out.omega_unit = WLAB;  out.omega_norm = WNORM;  out.cs_a = G.cs_a;
out.summary = build_summary(kk, ITG, TEM, kpar_k, WLAB);

fprintf('\n [4] 汇总: ITG 支 %d/%d 个 k 有不稳定根; TEM 支 %d/%d\n', ...
        sum([ITG.exists]), nK, sum([TEM.exists]), nK);
if any([ITG.exists])
    g = [ITG.gamma];  kkI = [ITG.k];
    [gmax, jm] = max(g);
    fprintf('     ITG: γ_max = %.5f %s @ k_yρ_s = %.4f;  峰值 k_yρ_s 区间 γ>0.8γ_max: [%.3f, %.3f]\n', ...
            gmax*WNORM, WLAB, kkI(jm), min(kkI(g > 0.8*gmax)), max(kkI(g > 0.8*gmax)));
end
if any([TEM.exists])
    g = [TEM.gamma];
    fprintf('     TEM: γ_max = %.5f %s @ k_yρ_s = %.4f\n', max(g)*WNORM, WLAB, ...
            [TEM(find(g == max(g), 1)).k]);
end
fprintf('     根搜索求值总数 = %d, 用时 %.1f s\n', sum(nD), toc(t0));
fprintf('\n (出图: plot_linear(out)  |  单 k 根谱: plot_linear_roots(out, k序号) )\n\n');
end

% =======================================================================================
% 局部函数
% =======================================================================================
function cfg = ct_fill_defaults(cfg)
%CT_FILL_DEFAULTS  用 CT_toroidal_defaults 补齐缺失/空的**配置**字段(剖面字段不在这里处理)
d = CT_toroidal_defaults();
fn = fieldnames(d);
for ii = 1:numel(fn)
    if ~isfield(cfg, fn{ii}) || isempty(cfg.(fn{ii})), cfg.(fn{ii}) = d.(fn{ii}); end
end
end

% ---------------------------------------------------------------------------------------
function fname = ct_profile_path(cfg)
%CT_PROFILE_PATH  唯一的剖面来源: cfg.profile_file (标准 txt); 一切不合规立即报错
prof_fields = {'rho','r','R0','B0','a','q','shear','kappa','A_ion','theta', ...
               'n_e','T_e','n_i','T_i','RL_ne','RL_Te','RL_ni','RL_Ti', ...
               'dn_e_dr','dT_e_dr','dn_i_dr','dT_i_dr'};
fname = '';
if isfield(cfg,'profile_file') && ~isempty(cfg.profile_file)
    if ~(ischar(cfg.profile_file) || isstring(cfg.profile_file))
        error('CT_linear:badProfileFileType', 'cfg.profile_file 必须是文件名(字符型)');
    end
    fname = char(cfg.profile_file);
end
if isfield(cfg,'profile') && isstruct(cfg.profile) && ~isempty(fieldnames(cfg.profile))
    error('CT_linear:structProfile', 'cfg.profile 结构体输入已删除; 请写成 txt 用 cfg.profile_file 指定');
end
bad = prof_fields(isfield(cfg, prof_fields));
if ~isempty(bad)
    error('CT_linear:profileFieldInCfg', 'cfg 中出现剖面字段: %s (剖面只能写在 txt 里)', strjoin(bad, ', '));
end
if isempty(fname)
    here = fileparts(mfilename('fullpath'));
    hint = '';
    if exist(fullfile(here,'profiles','profile_example.txt'), 'file') == 2
        hint = sprintf('\n    cfg.profile_file = ''%s'';', fullfile(here,'profiles','profile_example.txt'));
    end
    error('CT_linear:noProfileFile', '未提供剖面文件。请设置 cfg.profile_file, 例如%s', hint);
end
end

% ---------------------------------------------------------------------------------------
function gr = ct_root_grid(cfg)
%CT_ROOT_GRID  把精度档 cfg.root_grid 解析成 (nre, nim, refine, nseed); cfg 显式值优先
presets = { 'fast',     48, 24, 1, 16; ...
            'standard', 64, 32, 2, 24; ...
            'high',     96, 48, 3, 40; ...
            'ultra',   128, 64, 3, 64 };
name = 'standard';
if isfield(cfg,'root_grid') && ~isempty(cfg.root_grid), name = lower(char(cfg.root_grid)); end
k = find(strcmp(presets(:,1), name), 1);
if isempty(k)
    warning('CT_linear:root_grid', '未识别的 root_grid=''%s'', 回退 standard', name);  k = 2;
end
gr = struct('name',presets{k,1}, 'nre',presets{k,2}, 'nim',presets{k,3}, ...
            'refine',presets{k,4}, 'nseed',presets{k,5});
if isfield(cfg,'nre')        && ~isempty(cfg.nre),        gr.nre = cfg.nre;          end
if isfield(cfg,'nim')        && ~isempty(cfg.nim),        gr.nim = cfg.nim;          end
if isfield(cfg,'root_refine')&& ~isempty(cfg.root_refine),gr.refine = cfg.root_refine; end
if isfield(cfg,'root_nseed') && ~isempty(cfg.root_nseed), gr.nseed = cfg.root_nseed;  end
end

% ---------------------------------------------------------------------------------------
function [fac, lab] = ct_omega_unit(G, cfg)
%CT_OMEGA_UNIT  频率显示单位: 'cs/a'(默认) | 'krad/s' | 'rad/s'; 值 = rad/s × fac
u = 'cs/a';
if isfield(cfg,'omega_unit') && ~isempty(cfg.omega_unit), u = lower(char(cfg.omega_unit)); end
switch u
    case 'cs/a',   fac = 1/G.cs_a;  lab = 'c_s/a';
    case 'rad/s',  fac = 1;         lab = 'rad/s';
    case 'krad/s', fac = 1e-3;      lab = 'krad/s';
    otherwise
        warning('CT_linear:omega_unit', '未识别 omega_unit=''%s'', 回退 cs/a', u);
        fac = 1/G.cs_a;  lab = 'c_s/a';
end
end

% ---------------------------------------------------------------------------------------
function r = rec_empty(k)
%REC_EMPTY  逐 k 记录的**字段顺序**(供 rec_fill 对齐; 不要改顺序)
r = struct('k',k, 'branch','', 'omega_r',NaN, 'gamma',NaN, ...
           'omega_r_n',NaN, 'gamma_n',NaN, 'ratio',NaN, ...
           'kpar',NaN, 'kperp',NaN, 'wstar_i',NaN, 'wstar_e',NaN, 'wD0_e',NaN, ...
           'Dval',NaN, 'src','', 'exists',false);
end

function r = rec_fill(M, z, branch, fac, Dval, src)
%REC_FILL  由根填记录: omega_r/gamma 是 SI(rad/s), *_n 是按 cfg.omega_unit 归一的值
r = struct('k',M.k_rhos, 'branch',branch, 'omega_r',real(z), 'gamma',imag(z), ...
           'omega_r_n',real(z)*fac, 'gamma_n',imag(z)*fac, ...
           'ratio',abs(imag(z))/max(abs(real(z)),realmin), ...
           'kpar',M.kpar, 'kperp',M.kperp, 'wstar_i',M.wstar_i, 'wstar_e',M.wstar_e, ...
           'wD0_e',M.wD0_e, 'Dval',abs(Dval), 'src',src, 'exists',true);
end

function rec_print(r, ~)
rat = r.ratio;
kp = ternary(abs(r.kpar) > 0, sprintf('%-.4g', r.kpar), '0(k∥=0)');
fprintf('     %-9.3f %-6s %-13.4f %-13.4f %-8.3f %-10s %-10.1e %s\n', ...
        r.k, r.branch, r.omega_r_n, r.gamma_n, rat, kp, r.Dval, r.src);
end

function S = build_summary(kk, ITG, TEM, kpar_k, lab)
%BUILD_SUMMARY  汇总表(cell, 可直接写到文件或用 fprintf 打印)
S = cell(numel(kk)+2, 9);
S(1,:) = {'k_y rho_s', 'branch', ['omega_r [' lab ']'], ['gamma [' lab ']'], ...
          'gamma/|omega_r|', 'kpar [1/m]', '|D|', 'src', 'exists'};
row = 2;
for ii = 1:numel(kk)
    for r = [ITG(ii), TEM(ii)]
        if r.exists
            S{row,1} = r.k;            S{row,2} = r.branch;
            S{row,3} = r.omega_r_n;    S{row,4} = r.gamma_n;
            S{row,5} = r.ratio;        S{row,6} = r.kpar;
            S{row,7} = r.Dval;         S{row,8} = r.src;   S{row,9} = true;
        else
            S{row,1} = kk(ii);  S{row,2} = '—';  S{row,6} = kpar_k(ii);  S{row,9} = false;
        end
        row = row + 1;
    end
end
S = S(1:row-1, :);
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end
