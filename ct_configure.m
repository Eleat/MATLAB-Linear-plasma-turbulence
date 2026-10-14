function cfg = ct_configure(cfg, P)
%CT_CONFIGURE  [v8.5] 配置的**单一入口**: 归一化 + 校验 + 与几何的一致性检查
%=========================================================================================
% 【为什么需要】v8.3/v8.4 反复踩同一个坑: 开关在 ct_geometry 里留了副本(G.fbw_on /
%   G.nu_eff_model / G.alpha / G.ft / G.kpar_model ...), ct_model 过去**混着读** cfg 与副本,
%   于是"改了 cfg 却没重建 G"会**静默无效**(实测 fbw_model / nu_eff_model / alpha 各一次;
%   连本仓库自己的基准脚本也踩过两次)。本函数把"配置归一化"和"配置↔几何一致性"集中处理:
%
%     ① 归一化: 把预设/别名解析成**确定的取值**, 便于打印与对比
%        · kpar_preset 'cbc_benchmark'/'none' → kpar_model/theta_rms/factor(用户显式值优先)
%        · trap_precession 'auto' → 'full'(α≠0) | 'legacy'
%        · kpar_model 别名 'shear' → 'shear_legacy'
%     ② 校验: 枚举值拼错时**立即报错**(而不是跑到一半才发现)
%     ③ 一致性: cfg 与几何指纹 G.fp 不一致时, 一次性列出**所有**差异字段并告警
%        (提示"重建 G = ct_geometry(P, cfg)"), 而不是逐个字段各自 warning。
%
% 用法:  cfg = ct_configure(cfg);              % 只归一化+校验
%        cfg = ct_configure(cfg, G);            % 再与几何 G 比对(需要 G.fp, 见 ct_geom_fp)
%        P 也可直接给: cfg = ct_configure(cfg, P);  % 内部只做前置校验, 不比指纹
%=========================================================================================

if nargin < 1 || isempty(cfg) || ~isstruct(cfg)
    error('ct_configure:badCfg', 'cfg 必须是 CT_toroidal_defaults() 生成的配置结构体');
end
% 允许把**几何结构体**当作第 2 参数(它带 fp 字段)
Gfp = [];
if nargin >= 2 && ~isempty(P) && isstruct(P) && isfield(P,'fp') && isstruct(P.fp)
    Gfp = P.fp;
end

% ---------------- ① 归一化 ----------------
% kpar_preset: 用户显式给出的 kpar_model/theta_rms/factor 优先(v8.2 政策)
if isfield(cfg,'kpar_preset') && ~isempty(cfg.kpar_preset)
    pre = lower(char(cfg.kpar_preset));
    km0 = 'balloon';
    if isfield(cfg,'kpar_model') && ~isempty(cfg.kpar_model), km0 = lower(char(cfg.kpar_model)); end
    switch pre
        case 'cbc_benchmark'
            if strcmp(km0,'balloon'), cfg.kpar_model = 'balloon_rms'; end
        case {'none','zero'}
            if strcmp(km0,'balloon'), cfg.kpar_model = 'zero'; end
        otherwise
            error('ct_configure:badPreset', ...
                  ['未识别的 cfg.kpar_preset = ''%s''; 可选: ' ...
                   'cbc_benchmark(→balloon_rms) | none(→zero) | ''''(空, 不改动)'], ...
                  char(cfg.kpar_preset));
    end
end
% 别名与大小写
if isfield(cfg,'kpar_model') && ~isempty(cfg.kpar_model)
    cfg.kpar_model = lower(char(cfg.kpar_model));
    if strcmp(cfg.kpar_model,'shear'), cfg.kpar_model = 'shear_legacy'; end
end
% trap_precession 'auto' → 依赖 α
if isfield(cfg,'trap_precession') && ~isempty(cfg.trap_precession)
    tp = lower(char(cfg.trap_precession));
    if strcmp(tp,'auto')
        al = 0;  if isfield(cfg,'alpha') && ~isempty(cfg.alpha), al = cfg.alpha; end
        if abs(al) > 0, tp = 'full'; else, tp = 'legacy'; end
    end
    cfg.trap_precession = tp;
end

% ---------------- ② 校验(枚举值) ----------------
cfg = chk_list(cfg, 'kpar_model',    {'zero','balloon','balloon_rms','qr','shear_legacy'});
cfg = chk_list(cfg, 'trap_model',    {'kinetic','bounce'});
cfg = chk_list(cfg, 'trap_precession',{'legacy','full'});
cfg = chk_list(cfg, 'passing_e',     {'kinetic','adiabatic'});
cfg = chk_list(cfg, 'disp_model',    {'kinetic','pole'});
cfg = chk_list(cfg, 'j_numerator',   {'standard','legacy'});
cfg = chk_list(cfg, 'ft_model',      {'eps','fixed'});
cfg = chk_list(cfg, 'fbw_model',     {'banana','banana_avg','off'});
cfg = chk_list(cfg, 'fbw_trap_model',{'banana','banana_avg','off','energy'});
cfg = chk_list(cfg, 'nu_eff_model',  {'zero','collisional','fixed','tglf'});
cfg = chk_list(cfg, 'cD_model',      {'thermal','fixed'});

% ---------------- ③ 与几何指纹的一致性 ----------------
% 【只告警一次的机制】ct_model 每算一个 k 都会调本函数, 若每次都 warning 会把输出刷满;
% 这里用 persistent 记住"已经报过的字段组合", 同一组合只报一次。
persistent fp_warned
if ~isempty(Gfp)
    [names, diffs] = ct_geom_diff(cfg, Gfp);
    if ~isempty(names) && ~isequal(names, fp_warned)
        fp_warned = names;
        warning('ct_configure:staleGeometry', ...
            ['以下开关与几何副本不一致(几何是按旧配置建的):\n   %s\n' ...
             '  ⇒ 开关本身的取值以 **cfg 为准**(ct_model 已按 cfg 优先), 但几何里已经算好的量\n' ...
             '     (θ 处漂移因子 g、k⊥ 因子、捕获分数 f_t、ρ_b/FBW、平行波数基准、nu_eff 等)\n' ...
             '     仍是**旧值** ⇒ 结果自相矛盾。请重建: G = ct_geometry(P, cfg)。'], ...
            strjoin(names, ', '));
        cfg.stale_geometry_fields = names;         % 便于脚本自检
        for i = 1:numel(names), cfg.stale_geometry_old.(names{i}) = diffs{i}; end
    end
end
end

% ---------------------------------------------------------------------------------------
function cfg = chk_list(cfg, fld, allowed)
%CHK_LIST  枚举值校验(拼错立即报错)并把字面值统一成小写
if ~isfield(cfg, fld) || isempty(cfg.(fld)), return; end
v = lower(char(cfg.(fld)));
if ~any(strcmp(v, allowed))
    error('ct_configure:badValue', 'cfg.%s = ''%s'' 不是合法取值; 可选: %s', ...
          fld, char(cfg.(fld)), strjoin(allowed, ' | '));
end
cfg.(fld) = v;
end
