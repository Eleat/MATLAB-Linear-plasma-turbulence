function [d, cnt] = ct_dcount(f, w)
%CT_DCOUNT  色散函数求值计数器(供 ct_roots.m 的诊断统计使用)
%=========================================================================================
% 为什么单独一个文件: 计数需要跨调用保持状态, 用**持久变量**最省事; 但在 ct_roots.m 里用
%   "嵌套函数"会触发 MATLAB 的嵌套规则(含嵌套函数的文件里, 所有函数都必须显式 end,
%   且嵌套函数不能出现在控制流里), 易踩坑。独立文件 + persistent 既简单又安全。
%
% 用法:  ct_dcount('reset')       % 归零
%        d = ct_dcount(fh, w)     % 求值 fh(w) 并计数(返回结果)
%        n = ct_dcount()          % 读取当前计数
%=========================================================================================
persistent n
if isempty(n), n = 0; end

if nargin == 0                      % 读计数
    d = n;  cnt = n;  return;
end
if nargin == 1                      % 'reset' 或读计数
    d = n;  cnt = n;
    if (ischar(f) || isstring(f)) && strcmpi(char(f), 'reset'), n = 0; d = 0; cnt = 0; end
    return;
end
n = n + 1;                          % 求值 + 计数
d = f(w);
cnt = n;
end
