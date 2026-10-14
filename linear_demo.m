function out = linear_demo(which, savefig)
%LINEAR_DEMO  不稳定性模求解的一键演示(含出图)
%=========================================================================================
% 用法:
%   linear_demo                 % 默认 'cbc': CBC 绝热电子扫描 + 四联图 + 单 k 根谱
%   linear_demo('cbc')          % 同上(约 1–3 min, 取决于 root_grid)
%   linear_demo('kpar')         % 三种 k∥ 模型对比(zero / balloon_rms / qR)
%   linear_demo('alpha')        % s-α 的 α 扫描对比
%   linear_demo('tem')          % TEM 算例(离子梯度关掉, 电子动理学)
%   linear_demo('cbc', false)   % 不出图(只算)
%
% 输出: out(最后一次求解的结果); 图存到 figs\ 下。
% 依赖: 本文件夹的核心 .m 文件与 profiles\ 下的剖面。
%=========================================================================================
if nargin < 1 || isempty(which), which = 'cbc'; end
if nargin < 2 || isempty(savefig), savefig = true; end
which = lower(char(which));
here  = fileparts(mfilename('fullpath'));
addpath(here);
figdir = fullfile(here, 'figs');
if ~exist(figdir,'dir'), mkdir(figdir); end

switch which
    case 'cbc',   out = demo_cbc(here, figdir, savefig);
    case 'kpar',  out = demo_kpar(here, figdir, savefig);
    case 'alpha', out = demo_alpha(here, figdir, savefig);
    case 'tem',   out = demo_tem(here, figdir, savefig);
    otherwise
        error('linear_demo:badTask', '未知任务 ''%s''; 可选: cbc | kpar | alpha | tem', which);
end
end

% =======================================================================================
function base = base_cfg(here, pf, nk)
%BASE_CFG  演示用的公共配置
base = CT_toroidal_defaults();
base.profile_file = fullfile(here, 'profiles', pf);
base.k_unit   = 'rho_i';                 % 与 GX 的 ky 定义一致(k_θρ_i)
base.plots    = false;
base.selftest = false;
base.verbose  = true;
base.kmin = 0.05;  base.kmax = 0.80;  base.nk = nk;
base.root_grid = 'standard';
end

% ---------------------------------------------------------------------------------------
function out = demo_cbc(here, figdir, savefig)
%① CBC 绝热电子: 与 GX(Miller) 参考同口径的扫描
fprintf('\n===== [demo] CBC 绝热电子线性扫描 =====\n');
cfg = base_cfg(here, 'profile_cbc.txt', 16);
cfg.ft_model = 'fixed';  cfg.ft_fixed = 0;      % 绝热电子
cfg.passing_e = 'adiabatic';
cfg.kpar_preset = 'cbc_benchmark';              % balloon_rms + θ_rms = 1
cfg = ct_configure(cfg);
out = CT_linear(cfg);
if savefig
    plot_linear(out, struct('gx_ref',true, 'title','CBC 绝热电子: 不稳定性模 γ(k) 与 ω_r(k)', ...
                            'file', fullfile(figdir,'cbc_adiabatic_scan.png')));
    ik = find([out.ITG.exists], 1, 'last');
    if ~isempty(ik), plot_linear_roots(out, ik, struct('file', fullfile(figdir,'cbc_root_spectrum.png'))); end
end
end

% ---------------------------------------------------------------------------------------
function out = demo_kpar(here, figdir, savefig)
%② 三种 k∥ 模型对比 —— 平行动力学如何决定高 k 端的 γ 截止
fprintf('\n===== [demo] k∥ 模型对比 (CBC 绝热电子) =====\n');
models = {'zero', 'balloon_rms', 'qR'};
outs = cell(1, numel(models));
for i = 1:numel(models)
    cfg = base_cfg(here, 'profile_cbc.txt', 14);
    cfg.ft_model = 'fixed';  cfg.ft_fixed = 0;  cfg.passing_e = 'adiabatic';
    cfg.kpar_model = models{i};
    if strcmpi(models{i},'qR'), cfg.kpar_factor = 1.0; end
    if strcmpi(models{i},'balloon_rms'), cfg.kpar_theta_rms = 1.0; end
    cfg = ct_configure(cfg);
    fprintf(' --- %s ---\n', models{i});
    outs{i} = CT_linear(cfg);
end
out = outs{end};
if savefig
    f = figure('Color','w','Position',[80 80 820 560]);  clf(f);
    hold on;  grid on;  box on;
    col = [0.85 0.20 0.15; 0.10 0.45 0.80; 0.20 0.55 0.30];
    for i = 1:numel(models)
        x = outs{i}.k * sqrt(2*outs{i}.G.tau);
        g = nan(size(x));
        for j = 1:numel(x)
            if outs{i}.ITG(j).exists, g(j) = outs{i}.ITG(j).gamma_n; end
        end
        plot(x, g, '-o', 'Color',col(i,:), 'MarkerFaceColor',col(i,:), 'LineWidth',1.6);
    end
    % GX 参考
    [kx, gy] = gx_ref_table();   kx = kx/sqrt(2);
    plot(kx, gy, '-', 'Color',[0.35 0.35 0.35], 'LineWidth',1.8);
    legend([models, {'GX(Miller, 绝热)'}], 'Location','best');
    xlabel('k_y\rho_s');  ylabel('\gamma  [c_s/a]');
    title('k∥ 模型对线性增长率谱的影响(ITG 支)');
    try, exportgraphics(f, fullfile(figdir,'kpar_models_compare.png'), 'Resolution',200); catch, end
    fprintf(' 图已保存: %s\n', fullfile(figdir,'kpar_models_compare.png'));
end
end

% ---------------------------------------------------------------------------------------
function out = demo_alpha(here, figdir, savefig)
%③ s-α 的 α 扫描(θ 取 0; α 通过漂移/FLR/k∥ 与捕获电子弹跳平均起作用)
fprintf('\n===== [demo] s-α 的 α 扫描 (CBC) =====\n');
alphas = [0 0.5 1.0 1.5 3.0];
outs = cell(1, numel(alphas));
for i = 1:numel(alphas)
    cfg = base_cfg(here, 'profile_cbc.txt', 10);
    cfg.ft_model = 'fixed';  cfg.ft_fixed = 0;  cfg.passing_e = 'adiabatic';
    cfg.kpar_preset = 'cbc_benchmark';
    cfg.alpha = alphas(i);
    cfg = ct_configure(cfg);
    fprintf(' --- alpha = %.2f ---\n', alphas(i));
    outs{i} = CT_linear(cfg);
end
out = outs{end};
if savefig
    f = figure('Color','w','Position',[80 80 900 420]);  clf(f);
    tl = tiledlayout(f, 1, 2, 'Padding','compact','TileSpacing','compact');
    cmap = lines(numel(alphas));
    ax1 = nexttile(tl);  hold(ax1,'on');  grid(ax1,'on');  box(ax1,'on');
    ax2 = nexttile(tl);  hold(ax2,'on');  grid(ax2,'on');  box(ax2,'on');
    for i = 1:numel(alphas)
        x = outs{i}.k * sqrt(2*outs{i}.G.tau);
        g = nan(size(x));  w = nan(size(x));
        for j = 1:numel(x)
            if outs{i}.ITG(j).exists, g(j) = outs{i}.ITG(j).gamma_n; w(j) = outs{i}.ITG(j).omega_r_n; end
        end
        plot(ax1, x, g, '-o', 'Color',cmap(i,:), 'LineWidth',1.5, 'DisplayName',sprintf('\\alpha = %.1f', alphas(i)));
        plot(ax2, x, w, '-o', 'Color',cmap(i,:), 'LineWidth',1.5, 'DisplayName',sprintf('\\alpha = %.1f', alphas(i)));
    end
    xlabel(ax1,'k_y\rho_s');  ylabel(ax1,'\gamma [c_s/a]');  title(ax1,'\gamma(\alpha)');
    xlabel(ax2,'k_y\rho_s');  ylabel(ax2,'\omega_r [c_s/a]'); title(ax2,'\omega_r(\alpha)');
    legend(ax1,'Location','best');
    try, exportgraphics(f, fullfile(figdir,'alpha_scan.png'), 'Resolution',200); catch, end
    fprintf(' 图已保存: %s\n', fullfile(figdir,'alpha_scan.png'));
end
end

% ---------------------------------------------------------------------------------------
function out = demo_tem(here, figdir, savefig)
%④ TEM 算例: 关掉离子梯度(离子被动) + 电子动理学 ⇒ 不稳定根应落在 Re ω > 0 一侧
fprintf('\n===== [demo] TEM 算例(离子被动, 动力学电子) =====\n');
cfg = base_cfg(here, 'profile_cbc.txt', 12);
cfg.ft_model  = 'eps';            % 捕获电子按 ε 算
cfg.trap_model= 'bounce';         % 弹跳平均通道
cfg.fbw_trap_model = 'banana_avg';
cfg.passing_e = 'kinetic';
cfg.kpar_preset = 'cbc_benchmark';
cfg = ct_configure(cfg);
% 离子近乎被动: 见下面的说明(不能直接写 0, 校验要求 R/L > 0)
pf0 = ct_read_profile(cfg.profile_file);
% 把离子梯度压到很小: 离子几乎被动 ⇒ 不稳定根只剩电子通道驱动的 TEM。
% ★ 不能用 0: ct_read_profile 严格校验 R/L_X > 0(这是有意的, 免得"梯度为零"被误当成合法剖面),
%   所以这里写一个**很小的正数** 1e-3(等效于离子驱动可忽略), 而不是 0。
tmp = fullfile(tempdir, 'profile_tem_demo.txt');
ct_write_profile(tmp, pf0, 'TEM demo: 离子近乎被动');
txt = fileread(tmp);
txt = regexprep(txt, 'RL_Ti(\s*)=(\s*)[0-9.eE+-]+',   'RL_Ti$1=$2 1.0e-3');
txt = regexprep(txt, 'RL_ni(\s*)=(\s*)[0-9.eE+-]+',   'RL_ni$1=$2 1.0e-3');
fid = fopen(tmp,'w','n','UTF-8');  fwrite(fid, txt);  fclose(fid);
cfg.profile_file = tmp;
cfg = ct_configure(cfg);
out = CT_linear(cfg);
if savefig
    plot_linear(out, struct('title','TEM 算例(离子被动 + 动力学电子): 不稳定根在 Re \omega > 0 一侧', ...
                            'file', fullfile(figdir,'tem_case.png')));
end
delete(tmp);
end

% =======================================================================================
function [kx, gy] = gx_ref_table()
%GX_REF_TABLE  CBC 绝热电子 GX(Miller) 参考(k_θρ_i, c_s/a) —— 与 plot_linear 内的一致
kx = [0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75];
gy = [0.01749 0.04618 0.08260 0.11998 0.15272 0.17800 0.19466 0.20240 ...
      0.20154 0.19293 0.17811 0.15935 0.13928 0.11943 0.09999];
end
