function [names, olds] = ct_geom_diff(cfg, fp)
%CT_GEOM_DIFF  [v8.5] 比较"当前 cfg"与"几何指纹 fp"(建 G 时用的配置), 列出不一致的字段
%=========================================================================================
% 指纹由 ct_geom_fp(cfg) 生成、存在 G.fp 里(见 ct_geometry.m)。只比较**几何相关**的字段:
% 它们决定 G 的派生量(漂移因子 g(θ)、k⊥ 因子、f_t、ρ_b、平行波数基准、碰撞展宽口径等);
% 纯数值算法开关(如 root_grid/nk/kmax/plots)不在其中 —— 改它们不必重建几何。
%
% 返回: names = 不一致字段名(元胞), olds = 建 G 时的旧值(元胞, 与 names 对应)
%=========================================================================================
names = {};  olds = {};
if nargin < 2 || isempty(fp) || ~isstruct(fp), return; end
fn = fieldnames(fp);
for i = 1:numel(fn)
    f = fn{i};
    if ~isfield(cfg, f), continue; end
    a = cfg.(f);  b = fp.(f);
    same = false;
    if isnumeric(a) && isnumeric(b)
        same = isequal(size(a), size(b)) && all(abs(a(:) - b(:)) <= 1e-12*max(1, max(abs(b(:)))));
    elseif ischar(a) || isstring(a)
        same = strcmpi(char(a), char(b));
    else
        same = isequal(a, b);
    end
    if ~same
        names{end+1} = f;      %#ok<AGROW>
        olds{end+1}  = b;      %#ok<AGROW>
    end
end
end
