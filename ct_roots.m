function R = ct_roots(M, box, opt)
%CT_ROOTS  [v8] 在复平面矩形盒内求色散函数 D(ω) = 0 的根 —— 高精度版
%=========================================================================================
% v7 的做法是"48×24 粗网格绕数 + 30 个 |D| 最小点作 Newton 种子"。在 Cyclone base case 上
% 实测**漏掉 15 个 k 点中的 11 个**(见 ../benchmark_cbc/README.md), 因此 v8 改为:
%
%   (1) **可配置网格**: opt.nre × opt.nim (默认 64×32, 面积是 v7 默认的 ~1.8 倍);
%   (2) **自适应加密**: 对绕数 ≠ 0 的单元做 opt.refine 级二分加密(默认 2 级),
%       把候选根的位置分辨率提高 4 倍, 并把加密后的单元中心加入种子表;
%   (3) **多起点 Newton**: |D| 最小的 opt.nseed 个点(默认 24) + 加密候选 + 外部种子(沿 k 延拓);
%   (4) **自动扩盒**: 若盒边界上出现非零绕数(根可能被切在盒外), 自动把盒按 opt.expand
%       放大(默认 1.6 倍, 最多 opt.maxexpand 次)后重搜, 并在诊断里记录 nexpand;
%   (5) **可选下半平面**: opt.negim = true 时把搜索盒下沿开到 -imfrac·max(box.im),
%       用于找**阻尼根**(默认 false, 与原行为一致);
%   (6) **诊断输出**: R(1).diag 记录求值次数 nD、零点单元数、边界绕数、种子数、收敛数、
%       是否触发扩盒/截断 —— "没找到根"从此可以区分"盒外 / 网格太粗 / 确实无根"。
%
% 【v8.1 两处代价修正】(2026-10-13, 由"默认算例 26 个 k 跑了 80 min"的实测触发)
%   ① **删掉冗余重搜**: 旧写法 `for iexp = 0:opt.maxexpand` 在"盒内无根且边界无绕数"时
%      仍会把**同一个盒**连搜 maxexpand+1 次(实测 k=0.5: nD 2145→6435, 61 s→184 s)。
%      现在只有真的扩了盒(expanded=true)才进入下一轮。
%   ② **Newton 计入计数 + 停滞判据**: 旧版 Newton 用不带计数的 Dfun, 于是 diag.nD 只反映
%      网格求值, 漏掉真正的大头(每步 3 次求值 × 最多 80 步 × 25 起点 ≈ 6000 次);
%      同时给 Newton 加"连续 stallmax 步 |D| 无实质下降即退出"的停滞判据(stallmax 默认 40,
%      标定见 ct_newton 的说明), 无根时不再空转到迭代上限。
%      实测(无根 k, CBC): 190 s → **约 45 s**; 且两个弱根单元(γ/|ω_*| ≲ 3%)仍能找到。
%
% 输入
%   M    : ct_model 的输出
%   box  : struct('re',[lo hi], 'im',[lo hi])
%   opt  : 结构体(可缺省); 字段与默认值:
%          nre, nim   粗网格单元数 (64, 32)
%          refine     自适应加密级数 (2)
%          nseed      |D| 最小点作 Newton 种子的个数 (24)
%          expand     自动扩盒倍数 (1.6; <=1 表示不扩)
%          maxexpand  最多扩盒次数 (2)
%          negim      true 时搜索盒下沿也覆盖负虚部 (false)
%          imfrac     负虚部下沿 = -imfrac*max(box.im) (1.0)
%          maxroots   最多返回多少个根 (16)
%          maxrefcell 单次最多加密多少个候选单元 (64)
%          tol        去重判据 (1e-6 * 盒尺度)
%          seeds      外部 Newton 初值(如上一 k 的解; 连续延拓)
%          disp_model 传给 ct_dispersion
%          verbose    打印诊断
%          stallmax    Newton 连续多少步无实质改进就退出 (40; 见 ct_newton 的标定说明)
%          stalltol    "实质改进"的门槛: |D| 降不到 (1-stalltol) 倍即算停滞 (1e-4)
%          itmax       Newton 迭代上限 (80)
%
% 输出
%   R : 结构体数组, 字段 omega(复)、Dval、src(来源: winding/refined/seed/external)
%       R(1).diag 记录本次搜索的统计量
%
% 【两个必须记住的 MATLAB 坑】(v8 里都踩过一次)
%   ① cand 必须预分配成"空结构体数组" cand = struct('ctr',{},...); 若写成 cand = [] 再
%      cand(end+1) = struct(...) 会报"无法从 struct 转换为 double"。
%   ② 求值计数用独立的 persistent 计数器 ct_dcount.m, 不用嵌套函数 —— 含嵌套函数的文件里
%      所有函数都必须显式 end, 且嵌套函数不能出现在控制流中, 极易踩语法坑。
%
% 用法:  R = ct_roots(M, struct('re',[-4 4],'im',[1e-3 3]));
%=========================================================================================

if nargin < 3, opt = struct(); end
def = struct('nre',64, 'nim',32, 'refine',2, 'nseed',24, 'expand',1.6, 'maxexpand',2, ...
             'negim',false, 'imfrac',1.0, 'maxroots',16, 'maxrefcell',64, 'tol',[], ...
             'seeds',[], 'disp_model',[], 'verbose',false, ...
             'stallmax',40, 'stalltol',1e-4, 'itmax',80, ...
             'analyt_deriv',true, 'dDfun',[], 'boxes',{{}}, 'converge',false, 'maxconv',1);
fn = fieldnames(def);
for ii = 1:numel(fn)
    if ~isfield(opt, fn{ii}) || isempty(opt.(fn{ii})), opt.(fn{ii}) = def.(fn{ii}); end
end

% ---- [v8.5] D 可以由 **函数句柄** 给出(用于气球模本征值矩阵等自造色散函数), ----
%     也可以由 M 结构体给出(走 ct_dispersion); 二者都支持解析导数(每次迭代 2 次求值 vs 3 次)。
if isa(M, 'function_handle')
    Dfun  = M;
    dDfun = opt.dDfun;                 % 句柄形式下, 若给了解析导数就用它
else
    if isempty(opt.disp_model) && isfield(M.cfg,'disp_model') && ~isempty(M.cfg.disp_model)
        opt.disp_model = M.cfg.disp_model;
    end
    if isempty(opt.disp_model), opt.disp_model = 'kinetic'; end
    mD = opt.disp_model;
    Dfun = @(w) ct_dispersion(M, w, mD);
    if opt.analyt_deriv
        dDfun = @(w) ct_disp_deriv(M, w, mD);
    else
        dDfun = [];
    end
end
ct_dcount('reset');                       % 求值计数器(诊断用; 见 ct_dcount.m)
Dwrap = @(w) ct_dcount(Dfun, w);
if isempty(dDfun)
    dDwrap = [];
else
    dDwrap = @(w) ct_dcount(dDfun, w);
end

% ---- [v8.5] 多尺度盒: opt.boxes 给定时, 依次搜索并把根**取并集** ----
%   为什么需要: 小盒的网格相对更细, 能抓到"弱而低频"的根; 大盒才能容纳大 |ω| 支。
%   实测(α=1.0, CBC): γ=0.041 的弱根被 ±4ω_*T 大盒漏掉、被 ±1.5ω_*T 小盒抓到。
boxes = { box };
if ~isempty(opt.boxes)
    boxes = opt.boxes;
end
if isempty(opt.tol)
    opt.tol = 1e-6*max(diff(box.re), diff(box.im));
end

diag = struct('nD',0, 'box',box, 'nexpand',0, 'nzcell',0, 'nboundary',0, ...
              'nrefcell',0, 'nseed',0, 'ncand',0, 'nacc',0, 'Dscale',NaN, ...
              'truncated',false, 'nbox',0, 'nconv',0, 'nconvnew',0);
Dscale = 1;

roots = [];  srcs = {};
for ib = 1:numel(boxes)
    b0 = boxes{ib};
    if opt.negim                            % 可选: 连下半平面一起搜(阻尼根)
        b0.im = [-opt.imfrac*max(abs(b0.im)), max(b0.im)];
    end
    [rt, sr, dg, Dsc] = search_one_box(Dwrap, dDwrap, b0, opt);
    diag.nbox       = diag.nbox + 1;
    diag.nexpand    = diag.nexpand + dg.nexpand;
    diag.nzcell     = max(diag.nzcell, dg.nzcell);
    diag.nboundary  = max(diag.nboundary, dg.nboundary);
    diag.nrefcell   = max(diag.nrefcell, dg.nrefcell);
    diag.nseed      = max(diag.nseed, dg.nseed);
    diag.ncand      = max(diag.ncand, dg.ncand);
    if isnan(Dscale) || ~(Dsc > 0), Dscale = max(Dscale, Dsc); end
    [roots, srcs] = merge_roots(roots, srcs, rt, sr, opt.tol);

    % ---- [v8.5] 网格收敛自检: 把 nre/nim 加倍再搜一遍, 新出现的根并入并计数 ----
    %   "加网格后多出根"就是"网格太粗"的直接证据(本次 α=1.0 的弱根正是这样被抓回的补充)。
    for icv = 1:opt.maxconv
        if ~opt.converge, break; end
        o2 = opt;  o2.nre = 2*opt.nre;  o2.nim = 2*opt.nim;  o2.converge = false;
        [rt2, sr2] = search_one_box(Dwrap, dDwrap, b0, o2);
        n0 = numel(roots);
        [roots, srcs] = merge_roots(roots, srcs, rt2, sr2, opt.tol);
        diag.nconv = diag.nconv + 1;
        diag.nconvnew = diag.nconvnew + (numel(roots) - n0);
        if numel(roots) == n0, break; end        % 已收敛
    end
end


% ---------------- 打包输出 ----------------
if isempty(roots)
    R = struct('omega', {}, 'Dval', {}, 'src', {});
    d = diag;  d.nD = ct_dcount();  d.Dscale = Dscale;
    R(1).diag = d;                                     %#ok<AGROW>
    if opt.verbose
        fprintf(['    [ct_roots] 盒 [%g,%g]x[%g,%g]i 无根: 零点单元 %d(边界 %d), 扩盒 %d, ' ...
                 '求值 %d\n'], box.re, box.im, diag.nzcell, diag.nboundary, diag.nexpand, ct_dcount());
    end
    return;
end

[~, idx] = sort(imag(roots), 'descend');
roots = roots(idx);  srcs = srcs(idx);
if numel(roots) > opt.maxroots
    warning('ct_roots:maxroots', ...
            ['盒内共 %d 个根, 超过 maxroots = %d, 已按 γ 降序截断; ' ...
             '请提高 opt.maxroots 或缩小搜索盒后再下结论。'], numel(roots), opt.maxroots);
    roots = roots(1:opt.maxroots);  srcs = srcs(1:opt.maxroots);
    diag.truncated = true;
end

R = struct('omega', num2cell(roots));
for ii = 1:numel(roots)
    R(ii).Dval = abs(Dfun(roots(ii)));
    R(ii).src  = srcs{ii};
end
diag.nD = ct_dcount();  diag.Dscale = Dscale;  diag.nacc = numel(roots);  diag.box = box;
R(1).diag = diag;
if opt.verbose
    fprintf(['    [ct_roots] 盒 [%g,%g]x[%g,%g]i: 零点单元 %d(边界 %d), 加密 %d, 种子 %d, ' ...
             '收敛 %d, 扩盒 %d, 求值 %d\n'], box.re, box.im, diag.nzcell, diag.nboundary, ...
             diag.nrefcell, diag.nseed, numel(roots), diag.nexpand, ct_dcount());
end
end

% ---------------------------------------------------------------------------------------
function s = seed_src(ii, cand, nfixed)
%SEED_SRC  根的来源标签(cand 项为 winding/refined, 其后依次为 seed/external)
if ii <= numel(cand)
    s = cand(ii).src;
elseif ii <= nfixed
    s = 'seed';
else
    s = 'external';
end
end

% ---------------------------------------------------------------------------------------
function z = ct_newton(Dfun, z0, box, Dscale, stepmax, stallmax, stalltol, itmax, dDfun)
% 复 Newton 迭代, 带阻尼与盒约束; 收敛判据是**绝对**的 |D| < 1e-7·Dscale。
% 【v8.5】若给了 dDfun(解析 ∂D/∂ω 的句柄, 见 ct_disp_deriv), 用它代替中心差分:
%   每次迭代从 3 次 D 求值降到 **2 次**(值 1 次 + 解析导数 1 次), 且导数精确 ⇒ 迭代更少、
%   对病态点更稳(实测 CBC 'standard' 档单 k 求值次数与时间都明显下降)。
% 【v8.1】加停滞判据: 连续 stallmax 步 |D| 都没有实质下降(降不到 (1-stalltol) 倍) 就退出。
%   无根时 Newton 会在盒内游走到迭代上限(最多 3 次求值 × 80 步 × 25 起点 ≈ 6000 次求值,
%   比网格本身还贵); 真正收敛的 Newton 每步会降好几个量级。
%   默认 stallmax = 40 / stalltol = 1e-4 / itmax = 80 —— 这个取值是**实测标定**出来的(见下),
%   既保住弱根又省掉大部分空转; 可用 opt 覆盖。
%
%   【为什么是 40 而不是更小】CBC 上对两个"弱根"单元做阈值扫描(balloon_rms θ_rms=2 的
%   ky=0.75: γ=0.00606; qR factor 1.0 的 ky=0.55: γ=0.00458, 两者 γ/|ω_*| 都 ≲ 3%):
%     stallmax = 999(无判据, = v8.0 行为) → 两个弱根都找到, 无根 k 求值 ~7900
%     stallmax = 60 / 40                → **两个弱根都找到**, 无根 k 求值 ~5400
%     stallmax = 15                     → **两个弱根都丢失**(误报"无根"), 无根 k 求值 ~3400
%   因此默认取 40:相对 v8.0 把无根 k 的 Newton 空转削掉约 1/3 且结果不变;
%   若只求速度、且接受"可能漏掉 γ/|ω_*| ≲ 3% 的弱根", 可显式设 opt.stallmax = 15。
if nargin < 5 || isempty(stepmax),  stepmax  = 0.35*max(diff(box.re), diff(box.im)); end
if nargin < 6 || isempty(stallmax), stallmax = 40;   end
if nargin < 7 || isempty(stalltol), stalltol = 1e-4; end
if nargin < 8 || isempty(itmax),    itmax    = 80;   end
has_d = (nargin >= 9) && ~isempty(dDfun);
z = z0;  best = inf;  stall = 0;
for it = 1:itmax
    d0 = Dfun(z);
    if ~isfinite(d0), z = []; return; end
    a0 = abs(d0);
    if a0 < best*(1 - stalltol), best = a0;  stall = 0; else, stall = stall + 1; end
    if stall >= stallmax, break; end
    if has_d
        d1 = dDfun(z);
        if ~isfinite(d1) || abs(d1) < 1e-30        % 解析导数不可用(pole 模型等) ⇒ 退回差分
            h = 1e-6*max(1e-3, abs(z));
            d1 = (Dfun(z + h) - Dfun(z - h))/(2*h);
        end
    else
        h = 1e-6*max(1e-3, abs(z));
        d1 = (Dfun(z + h) - Dfun(z - h))/(2*h);
    end
    if ~isfinite(d1) || abs(d1) < 1e-30, break; end
    dz = -d0/d1;
    if abs(dz) > stepmax, dz = dz*stepmax/abs(dz); end
    z = z + dz;
    if abs(dz) < 1e-13*max(1, abs(z)), break; end
end
if ~isfinite(z)
    z = [];
elseif abs(Dfun(z)) > 1e-7*Dscale
    z = [];                       % 未收敛到根
end
end

% ---------------------------------------------------------------------------------------
function [roots, srcs, dg, Dscale] = search_one_box(Dwrap, dDfun, box0, opt)
%SEARCH_ONE_BOX  [v8.5] 对单个盒做完整搜索: 粗网格绕数 -> 自适应加密 -> 多起点 Newton,
%   含 v8.1 的"只有真的扩了盒才重搜"修正; 返回该盒内的根与诊断量。
box = box0;  nexpand = 0;  Dscale = 1;  nzcell = 0;  nbound = 0;  nref = 0;
nseed_used = 0;  ncand = 0;
roots = [];  srcs = {};
for iexp = 0:opt.maxexpand          % 旧循环头必须保留: 体内的 continue/break 依赖它
    % ---------------- 粗网格 ----------------
    re = linspace(box.re(1), box.re(2), opt.nre+1);
    im = linspace(box.im(1), box.im(2), opt.nim+1);
    [RE, IM] = meshgrid(re, im);
    W = RE + 1i*IM;
    Dg = zeros(size(W));
    for jj = 1:numel(W)
        Dg(jj) = Dwrap(W(jj));
        if ~isfinite(Dg(jj)), Dg(jj) = NaN; end
    end
    mags = sort(abs(Dg(isfinite(Dg(:)))));
    if isempty(mags), Dscale = 1; else, Dscale = mags(max(1,round(numel(mags)/2))); end
    if ~(Dscale > 0), Dscale = 1; end
    if isempty(opt.tol), opt.tol = 1e-6*max(diff(box.re), diff(box.im)); end

    % ---------------- (1) 粗网格绕数 ----------------
    cand = struct('ctr', {}, 'hw', {}, 'hh', {}, 'src', {});
    nzcell = 0;  nbound = 0;
    for j = 1:opt.nim
        for i = 1:opt.nre
            c = [Dg(j,i), Dg(j,i+1), Dg(j+1,i+1), Dg(j+1,i)];
            if any(~isfinite(c)), continue; end
            ph = angle(c);
            d = diff([ph, ph(1)]);
            d = mod(d + pi, 2*pi) - pi;
            winding = sum(d)/(2*pi);
            if abs(winding) >= 0.5
                nzcell = nzcell + 1;
                cand(end+1) = struct('ctr', (re(i)+re(i+1))/2 + 1i*(im(j)+im(j+1))/2, ...
                                     'hw', (re(i+1)-re(i))/2, 'hh', (im(j+1)-im(j))/2, ...
                                     'src', 'winding'); %#ok<AGROW>
                if i == 1 || i == opt.nre || j == 1 || j == opt.nim, nbound = nbound + 1; end
            end
        end
    end

    % ---------------- (4) 自动扩盒 ----------------
    % 【v8.1 修正】只有在**盒边界确实被根切到**(nbound > 0)时才值得放大重搜。
    %   旧写法用 `for iexp = 0:opt.maxexpand` + 末尾 `if isempty(roots), (不 break)`:
    %   当盒内没有根、边界也没有绕数时, 循环仍会**把同一个盒重搜 maxexpand+1 次**,
    %   实测把无根 k 的代价凭空放大 3 倍(examples: k=0.5 从 61 s 变 184 s, nD 2145→6435,
    %   而这 3 次搜索的盒完全相同)。现在改成: 没扩盒就直接跳出。
    expanded = false;
    if opt.expand > 1 && nbound > 0 && iexp < opt.maxexpand
        rc = mean(box.re);  ic = mean(box.im);
        box.re = [rc + (box.re(1)-rc)*opt.expand, rc + (box.re(2)-rc)*opt.expand];
        box.im = [ic + (box.im(1)-ic)*opt.expand, ic + (box.im(2)-ic)*opt.expand];
        nexpand = nexpand + 1;
        expanded = true;
    end
    if expanded
        continue;      % 本轮的加密/Newton 对"未放大的旧盒"没有意义, 直接用新盒重来
    end

    % ---------------- (2) 自适应加密 ----------------
    nref = min(numel(cand), opt.maxrefcell);
    for ic = 1:nref
        ctr = cand(ic).ctr;  hw = cand(ic).hw;  hh = cand(ic).hh;
        for lev = 1:opt.refine
            nsub = 2;
            xg = linspace(real(ctr)-hw, real(ctr)+hw, nsub+1);
            yg = linspace(imag(ctr)-hh, imag(ctr)+hh, nsub+1);
            best = [];  bestw = 0;
            for jj = 1:nsub
                for ii = 1:nsub
                    q = [Dwrap(xg(ii)+1i*yg(jj)), Dwrap(xg(ii+1)+1i*yg(jj)), ...
                         Dwrap(xg(ii+1)+1i*yg(jj+1)), Dwrap(xg(ii)+1i*yg(jj+1))];
                    if any(~isfinite(q)), continue; end
                    ph = angle(q);
                    dd = diff([ph, ph(1)]);
                    dd = mod(dd + pi, 2*pi) - pi;
                    wnd = sum(dd)/(2*pi);
                    if abs(wnd) >= 0.5
                        c2 = (xg(ii)+xg(ii+1))/2 + 1i*(yg(jj)+yg(jj+1))/2;
                        if isempty(best) || abs(wnd) > bestw, best = c2; bestw = abs(wnd); end
                    end
                end
            end
            if isempty(best), break; end
            ctr = best;  hw = hw/nsub;  hh = hh/nsub;
        end
        cand(ic).ctr = ctr;  cand(ic).hw = hw;  cand(ic).hh = hh;
        cand(ic).src = 'refined';
    end

    % ---------------- (3) 种子表 ----------------
    seeds = zeros(0,1);
    for ic = 1:numel(cand), seeds(end+1,1) = cand(ic).ctr; end            %#ok<AGROW>
    [~, order] = sort(abs(Dg(:)));
    for ii = 1:min(opt.nseed, numel(order))
        if isfinite(Dg(order(ii))), seeds(end+1,1) = W(order(ii)); end      %#ok<AGROW>
    end
    nfixed = numel(seeds);
    if ~isempty(opt.seeds), seeds = [seeds; opt.seeds(:)]; end
    nseed_used = numel(seeds);  ncand = numel(cand);

    % ---------------- Newton 精修 ----------------
    % 【v8.1 修正】Newton 必须用 **Dwrap**(带计数) 而不是 Dfun —— 否则 R(1).diag.nD
    %   只统计网格/加密的求值, 漏掉真正的大头(Newton 每步 3 次求值), 会低估代价 3-5 倍。
    stepmax = 0.35*max(diff(box.re), diff(box.im));
    for ii = 1:numel(seeds)
        z = ct_newton(Dwrap, seeds(ii), box, Dscale, stepmax, opt.stallmax, opt.stalltol, opt.itmax, dDfun);
        if isempty(z), continue; end
        if real(z) < box.re(1) || real(z) > box.re(2), continue; end
        if imag(z) < box.im(1) || imag(z) > box.im(2), continue; end
        if isempty(roots)
            roots = z;  srcs = {seed_src(ii, cand, nfixed)};
        elseif all(abs(roots - z) > opt.tol)
            roots = [roots; z];  srcs{end+1} = seed_src(ii, cand, nfixed); %#ok<AGROW>
        end
    end
    if ~isempty(roots), break; end
    if ~expanded, break; end       % 没扩盒 ⟹ 重搜同一个盒毫无意义(v8.1)
end
dg = struct('nexpand',nexpand, 'nzcell',nzcell, 'nboundary',nbound, 'nrefcell',nref, ...
            'nseed',nseed_used, 'ncand',ncand);
end

% ---------------------------------------------------------------------------------------
function [r, s] = merge_roots(r, s, r2, s2, tol)
%MERGE_ROOTS  把两次搜索的根取并集(按 tol 去重); 用于多尺度盒与网格收敛补充。
for i = 1:numel(r2)
    if isempty(r)
        r = r2(i);  s = s2(i);
    elseif all(abs(r - r2(i)) > tol)
        r(end+1) = r2(i);  s{end+1} = s2{i}; %#ok<AGROW>
    end
end
end

% ---------------------------------------------------------------------------------------
function d = ct_disp_deriv(M, w, disp_model)
%CT_DISP_DERIV  [v8.5] 解析 dD/dw(来自 ct_chanresp 的 JW); pole 模型返回 NaN。
[~, d] = ct_dispersion(M, w, disp_model, true);
end
