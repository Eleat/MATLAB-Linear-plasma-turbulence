function fp = ct_geom_fp(cfg)
%CT_GEOM_FP  [v8.5] 几何指纹: 把"决定 G 里派生量"的 cfg 字段打包成一个结构体
%=========================================================================================
% 存在 G.fp 里(见 ct_geometry.m 末尾)。用途: 之后任何人拿 cfg 与 G 配用时, 由
% ct_geom_diff 一次性比对出"改了 cfg 却没重建 G"的字段 ⇒ 消除 v8.3/v8.4 反复出现的
% "静默失效"陷阱(开关在 G 里留副本, ct_model 却读 cfg)。
%
% 只收录**几何相关**字段:
%   θ      —— 漂移因子 g(θ)、k⊥ 因子、k∥ 基准都依赖它
%   shear  —— 同上(在剖面里)
%   alpha  —— s-α 压力梯度参数
%   ft_model/ft_fixed —— 捕获分数 f_t
%   no_curv / trap_geom —— 漂移是否关闭 / 捕获通道几何覆盖
%   fbw_model / fbw_factor —— 通行通道轨道宽度
%   kpar_model/kpar_factor/kpar_geom/kpar_theta_rms/kpar_sign —— 平行波数基准
%   Zeff / nu_eff_model / nu_eff_factor —— 碰撞展宽
%   polar_term / disp_model / N_perp —— 极化项与速度空间节点
%=========================================================================================
fp = struct();
all_f = {'alpha','no_curv','trap_geom','ft_model','ft_fixed','fbw_model','fbw_factor', ...
         'kpar_model','kpar_factor','kpar_geom','kpar_theta_rms','kpar_sign', ...
         'Zeff','nu_eff_model','nu_eff_factor','polar_term','disp_model','N_perp', ...
         'N_trap_x','N_trap_xi','trap_model','fbw_trap_model','trap_precession', ...
         'passing_e','j_numerator'};
for i = 1:numel(all_f)
    f = all_f{i};
    if isfield(cfg, f) && ~isempty(cfg.(f)), fp.(f) = cfg.(f); end
end
fp.shear = NaN;  fp.theta = NaN;      % 这两个在剖面里, 由 ct_geometry 填充
end
