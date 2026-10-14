function h = plot_linear_roots(out, idx, opt)
%PLOT_LINEAR_ROOTS  某个 k 上**复平面根谱**图 + 该 k 的 γ/ω_r 剖面(两联图)
%=========================================================================================
% 用途: 看清"扫描给出的那一个根"在复平面上相对于**盒内全部根**的位置 —— 判断
%   (a) 根搜索是否漏根(盒内还有别的更不稳定的根吗?),
%   (b) 支的分界(Re ω 的符号),
%   (c) 根的强弱(γ 与 |ω| 的比)。
%
% 用法:
%   plot_linear_roots(out, 5);                       % 第 5 个 k
%   plot_linear_roots(out, 0.283, struct('file','figs\spectrum.png'));  % 按 k_yρ_s 找最近点
%
% opt 字段
%   file    ''      非空则保存
%   box     true    是否画出根搜索盒(以该 k 的 ω_* 尺度重建)
%   clim    []      颜色标量(默认用 γ)
%=========================================================================================

if nargin < 2 || isempty(idx), idx = 1; end
if nargin < 3 || isempty(opt), opt = struct(); end
def = struct('file','', 'box',true, 'clim',[]);
fn = fieldnames(def);
for i = 1:numel(fn)
    if ~isfield(opt, fn{i}) || isempty(opt.(fn{i})), opt.(fn{i}) = def.(fn{i}); end
end

% ---- 定位 k ----
if idx >= 1 && idx <= numel(out.k) && (idx == round(idx)) && idx > 0
    ik = idx;
else
    [~, ik] = min(abs(out.k - idx));      % 按 k_yρ_s 找最近
end
k = out.k(ik);
om = out.roots{ik};
fac = out.omega_norm;
if isempty(om)
    warning('plot_linear_roots:noRoots', 'k_yρ_s = %.4f 的盒内没有根(或该点被跳过)', k);
    h = struct();  return;
end

gi = NaN; oi = NaN; ge = NaN; oe = NaN;
if out.ITG(ik).exists, gi = out.ITG(ik).gamma_n;  oi = out.ITG(ik).omega_r_n; end
if out.TEM(ik).exists, ge = out.TEM(ik).gamma_n;  oe = out.TEM(ik).omega_r_n; end

h.fig = figure('Color','w', 'Position',[100 100 1080 460], 'Name',sprintf('roots @ k=%.4f', k));
clf(h.fig);
tl = tiledlayout(h.fig, 1, 2, 'Padding','compact', 'TileSpacing','compact');
title(tl, sprintf('复平面根谱   k_{y}\\rho_s = %.4f   (k∥ 模型 = %s,  盒内 %d 个根)', ...
      k, out.cfg.kpar_model, numel(om)), 'FontWeight','bold');

% ---------- 左: 复平面 ----------
ax1 = nexttile(tl);
hold(ax1,'on');  grid(ax1,'on');  box(ax1,'on');
zr = real(om)*fac;  zi = imag(om)*fac;
scatter(ax1, zr, zi, 46, zi, 'filled', 'MarkerEdgeColor',[0.25 0.25 0.25]);
yline(ax1, 0, '-', 'Color',[0.55 0.55 0.55]);
xline(ax1, 0, ':', 'Color',[0.55 0.55 0.55]);
text(ax1, 0.02, 0.96, 'TEM 侧 (Re ω>0)', 'Units','normalized', 'Color',[0.10 0.45 0.80]);
text(ax1, 0.98, 0.96, 'ITG 侧 (Re ω<0)', 'Units','normalized', 'HorizontalAlignment','right', ...
     'Color',[0.85 0.20 0.15]);
if isfinite(gi)
    plot(ax1, oi, gi, 'o', 'MarkerSize',13, 'LineWidth',2.4, 'Color',[0.85 0.20 0.15]);
    text(ax1, oi, gi, '  ITG(选中)', 'Color',[0.85 0.20 0.15], 'FontWeight','bold');
end
if isfinite(ge)
    plot(ax1, oe, ge, 's', 'MarkerSize',12, 'LineWidth',2.4, 'Color',[0.10 0.45 0.80]);
    text(ax1, oe, ge, '  TEM(选中)', 'Color',[0.10 0.45 0.80], 'FontWeight','bold');
end
if opt.box
    wsc = max(abs(out.ITG(ik).wstar_i), abs(out.ITG(ik).wstar_e));
    if ~isfinite(wsc) || wsc == 0, wsc = max(abs(om)); end
    xr = [-out.cfg.root_re out.cfg.root_re]*wsc*fac;
    yr = [out.cfg.root_im0 out.cfg.root_im1]*wsc*fac;
    rectangle(ax1, 'Position',[xr(1) yr(1) diff(xr) diff(yr)], 'EdgeColor',[0.75 0.75 0.75], ...
              'LineStyle','--', 'LineWidth',1.2);
end
xlabel(ax1, sprintf('Re \\omega  [%s]', out.omega_unit));
ylabel(ax1, sprintf('Im \\omega = \\gamma  [%s]', out.omega_unit));
title(ax1, '盒内全部根(颜色 = γ)');
axis(ax1,'tight');

% ---------- 右: 同一 k 的稳定/不稳定说明 + 根的强弱 ----------
ax2 = nexttile(tl);
hold(ax2,'on');  grid(ax2,'on');  box(ax2,'on');
gm = max(zi);  omax = max(abs(zr));
bar(ax2, 1, gm, 0.45, 'FaceColor',[0.85 0.20 0.15]);
bar(ax2, 2, omax, 0.45, 'FaceColor',[0.45 0.45 0.55]);
if isfinite(gi), plot(ax2, 1, gi, 'ko', 'MarkerSize',9, 'LineWidth',2); end
if isfinite(ge), plot(ax2, 1, ge, 'ks', 'MarkerSize',9, 'LineWidth',2); end
set(ax2, 'XTick',[1 2], 'XTickLabel',{'max γ (盒内)','max |Re ω|'});
ylabel(ax2, sprintf('[%s]', out.omega_unit));  xlim(ax2,[0.4 2.6]);
title(ax2, sprintf('该 k 的根: 数量 %d, 最强 γ = %.4g, γ/|ω| = %.3f', ...
      numel(om), gm, gm/max(omax,realmin)));
txt = sprintf(['k∥(该 k) = %.4g /m\nγ_ITG = %s\nγ_TEM = %s\n求值次数 = %d\n扩盒次数 = %d'], ...
      out.kpar(ik), n2s(gi), n2s(ge), out.nD(ik), out.nexpand(ik));
text(ax2, 0.03, 0.06, txt, 'Units','normalized', 'VerticalAlignment','bottom', ...
     'FontName','Consolas', 'FontSize',9, 'BackgroundColor',[0.96 0.96 0.98], ...
     'EdgeColor',[0.8 0.8 0.85]);
h.ax1 = ax1;  h.ax2 = ax2;

if ~isempty(opt.file)
    f = char(opt.file);
    if isempty(regexp(f, '\.(png|jpg|jpeg|pdf|eps|fig)$', 'once')), f = [f '.png']; end
    d = fileparts(f);  if ~isempty(d) && ~exist(d,'dir'), mkdir(d); end
    try, exportgraphics(h.fig, f, 'Resolution', 200); catch, saveas(h.fig, f); end
    fprintf(' 图已保存: %s\n', f);
end
end

% ---------------------------------------------------------------------------------------
function s = n2s(v)
if isfinite(v), s = sprintf('%.5f', v); else, s = '—'; end
end
