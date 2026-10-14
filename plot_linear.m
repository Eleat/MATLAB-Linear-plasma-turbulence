function h = plot_linear(out, opt)
%PLOT_LINEAR  不稳定性模结果的标准图(四联图 + 可选 GX 参考叠加)
%=========================================================================================
% 用法:
%   out = CT_linear(cfg);  plot_linear(out);
%   plot_linear(out, struct('file','figs\cbc_scan','gx_ref',true,'title','CBC 绝热电子'));
%
% 四个子图:
%   (1) γ(k_yρ_s)               —— ITG 实心圆实线, TEM 空心方块虚线
%   (2) ω_r(k_yρ_s)             —— 同一支的颜色约定; 水平零线标出 ITG/TEM 的分界
%   (3) γ/|ω_r|                 —— 模的"品质"(γ 与实频的相对大小)
%   (4) k∥(k_yρ_s)              —— 该配置下的平行波数(来自 ct_model, 与支无关)
%
% opt 字段
%   file     ''        非空则保存图片(自动按后缀决定格式; 无后缀按 .png 处理)
%   gx_ref   false     true = 叠加 CBC 绝热电子的 GX(Miller) 参考 γ(本文件夹自带该表)
%   title    ''        总标题
%   kunit    ''        ''(自动: out.cfg.k_unit) | 'rho_s' | 'rho_i' —— 横轴单位
%   ylim_g   []        γ 轴上限(留空自动)
%   figure   true      false = 不新建图(只返回句柄; 用于脚本里批量出图)
%=========================================================================================

if nargin < 2 || isempty(opt), opt = struct(); end
def = struct('file','', 'gx_ref',false, 'title','', 'kunit','', 'ylim_g',[], 'figure',true);
fn = fieldnames(def);
for i = 1:numel(fn)
    if ~isfield(opt, fn{i}) || isempty(opt.(fn{i})), opt.(fn{i}) = def.(fn{i}); end
end
if ~isfield(out,'ITG') || ~isfield(out,'TEM')
    error('plot_linear:badInput', 'out 必须是 CT_linear 的返回值(需含 ITG/TEM 字段)');
end

% ---- 横轴: 默认用产出时的单位(与用户配置一致) ----
ku = out.cfg.k_unit;
if ~isempty(opt.kunit), ku = opt.kunit; end
if strcmpi(ku,'rho_i')
    x = out.k * sqrt(2*out.G.tau);   xlab = 'k_θρ_i';
else
    x = out.k;                       xlab = 'k_yρ_s';
end

% ---- 取出两支曲线(只在 exists 处连线, 缺失点断开) ----
gi = nan(size(x));  oi = gi;   ri = gi;
ge = nan(size(x));  oe = ge;   re = ge;
for i = 1:numel(x)
    if out.ITG(i).exists
        gi(i) = out.ITG(i).gamma_n;  oi(i) = out.ITG(i).omega_r_n;  ri(i) = out.ITG(i).ratio;
    end
    if out.TEM(i).exists
        ge(i) = out.TEM(i).gamma_n;  oe(i) = out.TEM(i).omega_r_n;  re(i) = out.TEM(i).ratio;
    end
end

if opt.figure
    h.fig = figure('Color','w', 'Position',[80 80 1120 760], 'Name','CT_linear');
else
    h.fig = gcf;
end
clf(h.fig);
tl = tiledlayout(h.fig, 2, 2, 'Padding','compact', 'TileSpacing','compact');
if isempty(opt.title)
    ttl = sprintf('不稳定性模:  %s   (k∥ 模型 = %s,  τ = %.3f,  θ = %.3f)', ...
                  out.P.src_file, out.cfg.kpar_model, out.G.tau, out.P.theta);
else
    ttl = opt.title;
end
title(tl, ttl, 'Interpreter','none', 'FontWeight','bold');
h.tl = tl;

cols = struct('itg',[0.85 0.20 0.15], 'tem',[0.10 0.45 0.80]);
lab_itg = 'ITG (Re \omega<0)';   if ~any(isfinite(gi)), lab_itg = [lab_itg ' — 本算例无根']; end
lab_tem = 'TEM (Re \omega>0)';   if ~any(isfinite(ge)), lab_tem = [lab_tem ' — 本算例无根']; end

% ---------- (1) γ ----------
ax1 = nexttile(tl);
hold(ax1,'on');  grid(ax1,'on');  box(ax1,'on');
plot(ax1, x, gi, '-o', 'Color',cols.itg, 'MarkerFaceColor',cols.itg, 'LineWidth',1.6, 'DisplayName',lab_itg);
plot(ax1, x, ge, '--s', 'Color',cols.tem, 'MarkerFaceColor','w', 'LineWidth',1.6, 'DisplayName',lab_tem);
if opt.gx_ref
    [kx, gy] = gx_cbc_adiabatic();          % k_θρ_i, γ[c_s/a]
    if strcmpi(ku,'rho_s'), kx = kx/sqrt(2); end
    plot(ax1, kx, gy, '-', 'Color',[0.35 0.35 0.35], 'LineWidth',1.8, ...
         'DisplayName','GX(Miller, 绝热电子)');
end
xlabel(ax1, xlab);  ylabel(ax1, '\gamma  [c_s/a]');
title(ax1, '增长率 \gamma');  legend(ax1,'Location','best');
if ~isempty(opt.ylim_g), ylim(ax1, opt.ylim_g); end
h.ax_g = ax1;

% ---------- (2) ω_r ----------
ax2 = nexttile(tl);
hold(ax2,'on');  grid(ax2,'on');  box(ax2,'on');
yline(ax2, 0, ':', 'Color',[0.6 0.6 0.6]);
plot(ax2, x, oi, '-o', 'Color',cols.itg, 'MarkerFaceColor',cols.itg, 'LineWidth',1.6);
plot(ax2, x, oe, '--s', 'Color',cols.tem, 'MarkerFaceColor','w', 'LineWidth',1.6);
xlabel(ax2, xlab);  ylabel(ax2, '\omega_r  [c_s/a]');
title(ax2, '实频率 \omega_r (符号即支的方向)');
h.ax_w = ax2;

% ---------- (3) γ/|ω_r| ----------
ax3 = nexttile(tl);
hold(ax3,'on');  grid(ax3,'on');  box(ax3,'on');
plot(ax3, x, ri, '-o', 'Color',cols.itg, 'MarkerFaceColor',cols.itg, 'LineWidth',1.6);
plot(ax3, x, re, '--s', 'Color',cols.tem, 'MarkerFaceColor','w', 'LineWidth',1.6);
xlabel(ax3, xlab);  ylabel(ax3, '\gamma/|\omega_r|');
title(ax3, '驱动强度比 \gamma/|\omega_r|');
h.ax_r = ax3;

% ---------- (4) k∥ ----------
ax4 = nexttile(tl);
hold(ax4,'on');  grid(ax4,'on');  box(ax4,'on');
plot(ax4, x, out.kpar, '-^', 'Color',[0.20 0.55 0.30], 'MarkerFaceColor',[0.20 0.55 0.30], ...
     'LineWidth',1.6);
xlabel(ax4, xlab);  ylabel(ax4, 'k_{||}  [1/m]');
title(ax4, sprintf('平行波数 (模型: %s)', out.cfg.kpar_model));
if all(abs(out.kpar) < 1e-12)
    text(ax4, 0.5, 0.5, 'k∥ ≡ 0 (无平行动力学)', 'Units','normalized', ...
         'HorizontalAlignment','center', 'Color',[0.7 0.2 0.2], 'FontWeight','bold');
end
h.ax_k = ax4;

% ---------- 保存 ----------
if ~isempty(opt.file)
    f = char(opt.file);
    if isempty(regexp(f, '\.(png|jpg|jpeg|pdf|eps|fig)$', 'once')), f = [f '.png']; end
    d = fileparts(f);
    if ~isempty(d) && ~exist(d,'dir'), mkdir(d); end
    try
        exportgraphics(h.fig, f, 'Resolution', 200);
    catch
        saveas(h.fig, f);
    end
    fprintf(' 图已保存: %s\n', f);
    h.file = f;
end
end

% =======================================================================================
function [ky_krhoi, gam] = gx_cbc_adiabatic()
%GX_CBC_ADIABATIC  CBC 绝热电子的 GX(Miller 几何) 参考线性增长率
%   来源: benchmark_cbc/cbc_compare_gx.md 的逐波数表(由 GX 官方基准输出提取)。
%   单位: ky 为 k_θρ_i, γ 为 c_s/a。★ 仅作参照带, 不是本程序的目标函数。
ky_krhoi = [0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75];
gam      = [0.01749 0.04618 0.08260 0.11998 0.15272 0.17800 0.19466 0.20240 ...
            0.20154 0.19293 0.17811 0.15935 0.13928 0.11943 0.09999];
end
