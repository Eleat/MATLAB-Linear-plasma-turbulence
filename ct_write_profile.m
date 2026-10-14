function fpath = ct_write_profile(fname, P, header)
%CT_WRITE_PROFILE  把剖面参数写成**标准剖面 txt 文件**（ct_read_profile 的唯一输入格式）
%=========================================================================================
% 本函数是"写"的一端: 把参数结构体 / 实验数据整理成规范的行 "名字 = 数值", 供
% ct_read_profile 读取, 也供自检写临时文件。**它不是另一条输入通道** ——
% 程序读剖面永远只走 ct_read_profile('xxx.txt')。
%
% 用法
%   % (1) 从参数结构体生成
%   Q = struct('rho',0.62,'R0',1.65,'B0',2.5,'a',0.6,'q',1.4,'shear',0.8, ...
%              'A_ion',2,'theta',0,'n_e',5e19,'T_e',1200,'n_i',5e19,'T_i',600, ...
%              'RL_ne',3,'RL_Te',9,'RL_ni',3,'RL_Ti',6);
%   ct_write_profile('my_profile.txt', Q);
%
%   % (2) 由已有剖面派生一个新温度(例如做 T_e/T_i 扫描)
%   P = ct_read_profile('profile_example.txt');
%   P.T_e = 1500;  P.T_i = 750;
%   ct_write_profile('scan_Te1500.txt', P);
%
%   % (3) 只写部分字段 + 自定义表头
%   ct_write_profile(tmp, struct('T_e',1000), '仅用于演示');
%
% 说明
%   * 只写 ct_read_profile 认识的字段(未知字段被忽略, 不写进文件);
%   * 值为空/NaN 的字段跳过(由读取端的默认值或派生逻辑处理);
%   * 写完后可用 ct_read_profile(fpath) 读回做往返校验(自检 [T9] 就是这么做的)。
%
% 输出
%   fpath : 实际写出的绝对路径(便于紧接其后读取/删除)
%=========================================================================================

if nargin < 1 || isempty(fname)
    error('ct_write_profile: 需要输出文件名, 例如 ct_write_profile(''my_profile.txt'', Q)');
end
if nargin < 2 || ~isstruct(P) || isempty(fieldnames(P))
    error('ct_write_profile: 第二个输入必须是含剖面字段的结构体');
end
if nargin < 3 || isempty(header)
    header = '剖面文件 (标准格式: 名字 = 数值; # 之后为注释)';
end
fname = char(fname);
if isempty(regexp(fname, '^([A-Za-z]:[\\/]|\\\\)', 'once'))
    fpath = fullfile(pwd, fname);
else
    fpath = fname;
end

% 按"几何 → 局部量 → 梯度 → 可选"的顺序写, 便于人读
order = {'rho','r','R0','B0','a','q','shear','kappa','A_ion','theta', ...
         'n_e','T_e','n_i','T_i', ...
         'RL_ne','RL_Te','RL_ni','RL_Ti', ...
         'dn_e_dr','dT_e_dr','dn_i_dr','dT_i_dr'};

[fp, nm, ext] = fileparts(fpath);   %#ok<ASGLU>
fid = fopen(fpath, 'w');
if fid < 0
    error('ct_write_profile: 无法写出 "%s" (检查目录是否存在/是否可写)', fpath);
end
fprintf(fid, '# %s\n', header);
fprintf(fid, '# 由 ct_write_profile.m 生成; 读取用 ct_read_profile(''%s'')\n', [nm ext]);
fprintf(fid, '# 必需字段: rho,R0,B0,a,q,shear,A_ion,theta,n_e,T_e + 四项 R/L(或物理梯度)\n');
fprintf(fid, '\n');
for ii = 1:numel(order)
    k = order{ii};
    if ~isfield(P, k), continue; end
    v = P.(k);
    if isempty(v) || ~isnumeric(v) || ~isfinite(v), continue; end
    fprintf(fid, '%-8s = %.10g\n', k, v);
end
fclose(fid);
end
