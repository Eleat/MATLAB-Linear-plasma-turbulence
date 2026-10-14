function P = ct_read_profile(fname)
%CT_READ_PROFILE  读取**唯一**的剖面输入：标准 txt 文件（严格校验，报错定位到行）
%=========================================================================================
% 本程序**只接受一种剖面输入**：标准文本文件。结构体 / cfg 字段两种旧方式已删除。
%
% ------------------------------- 调用方式 -------------------------------
%   P = ct_read_profile('profile_example.txt');       % 相对当前目录
%   P = ct_read_profile('D:\data\shot12345.txt');     % 或绝对路径
%   文件不在当前目录时, 程序还会在本 .m 文件所在目录里找 (便于把常用剖面放在程序包内)。
%
% ------------------------------- 文件格式 -------------------------------
%   * 每行一条:  名字 = 数值
%   * '#' 之后到行尾为注释; 空行忽略
%   * 名字必须是合法 MATLAB 变量名(字母开头, 只含字母/数字/下划线)
%   * 同一个名字只能出现一次(重复报错)
%   * 名字必须在下面的"字段表"内(拼错会报错并给出最接近的名字建议)
%   * 数值必须是可解析的十进制数(支持 1.2e3, -5, .5); 空值/非数报错
%
% ------------------------------- 字段表 -------------------------------
% 【必需】(缺任何一个都报错, 并一次性列出全部缺失项)
%   rho      r/a                                       [-]      0 < rho < 1
%   R0       大半径                                     [m]      > 0
%   B0       轴上环向磁场                                [T]      > 0
%   a        小半径                                     [m]      > 0
%   q        安全因子                                   [-]      > 0
%   shear    磁剪切  s_hat = (r/q) dq/dr                 [-]      有限
%   A_ion    主离子质量数(D=2, H=1, DT=2.5)              [-]      > 0
%   theta    极向角 (0 = 外中平面)                       [rad]    |theta| <= pi
%   n_e      局部电子密度                                [m^-3]   > 0
%   T_e      局部电子温度                                [eV]     > 0
%   梯度四项 见下(必须全部确定)
%
% 【梯度: 二选一, 同一物理量不能同时给】
%   RL_ne, RL_Te, RL_ni, RL_Ti   归一化梯度 R/L_X = -R dlnX/dr   [-]
%   或 dn_e_dr, dT_e_dr, dn_i_dr, dT_i_dr   物理梯度 dX/dr [SI]
%   规则: 给出物理梯度时自动反算 R/L_X; 两者同时给出 -> 报错(避免歧义)。
%   全部四个 R/L 都必须确定(直接给, 或由物理梯度算出), 缺一个都报错。
%   要求 RL_ne>0, RL_Te>0, RL_ni>0, RL_Ti>0 (本模型讨论的"通常剖面"区间)。
%
% 【可选】(不写则用默认值)
%   kappa    拉长比(本版按圆截面处理)                    [-]      默认 1.00
%   r        小半径 [m]; 不写则 r = rho*a                [m]      需 0 < r < a
%            若同时给了 r 与 rho, 必须满足 |r - rho*a| <= 0.05*a, 否则报错
%   n_i      局部主离子密度; 不写则 = n_e                 [m^-3]   > 0
%   T_i      局部主离子温度; 不写则 = T_e                 [eV]     > 0
%
% ------------------------------- 输出 -------------------------------
%   P : 结构体, 含上述全部字段(已补齐可选默认值), 另有派生量:
%       P.Lne,Lni,LTe,LTi       尺度长度 [m]  L_X = R0/(R/L_X)
%       P.eta_i = RL_Ti/RL_ni ,  P.eta_e = RL_Te/RL_ne
%       P.src_file              实际读取的文件名(供打印/记录)
%
% ------------------------------- 报错一览 -------------------------------
%   找不到剖面文件 / 文件为空 / 某行缺少 "=" / 字段名不合法 / 未知字段名(附建议)
%   / 值为空或不是数 / 同一字段重复定义 / 缺少必需字段(列出全部) / 数值越界
%   / 物理梯度与 R/L 同时给出 / r 与 rho*a 不一致 / 传入结构体(不再支持)
%   以上都在 ct_selftest 的 [T17] 里逐项验证。
%
% 用法:  P = ct_read_profile('my_profile.txt');
%        写文件用 ct_write_profile.m (把参数/实验数据落成标准 txt)
%=========================================================================================

% 【报错信息的写法约定 —— 请勿改成 error([...])】
%   error(['多行\n消息']) 这种"单个字符数组入参"的写法**不会**被 MATLAB 做 sprintf 处理,
%   于是 \n 与 %% 会原样打印(用户看到的是字面 "\n")。因此本文件一律写成
%       msg = sprintf([...]);
%       error('ct_read_profile:xxx', '%s', msg);
%   这样既能保证换行/百分号正确, 又给每个错误一个可 catch 的 message identifier。
if nargin < 1 || isempty(fname)
    msg = sprintf(['ct_read_profile: 必须提供剖面文件名, 例如\n' ...
                   '    P = ct_read_profile(''profile_example.txt'');\n' ...
                   '  本程序唯一的剖面输入是标准 txt 文件(格式见本文件头部说明)。']);
    error('ct_read_profile:noFileName', '%s', msg);
end
if isstruct(fname)
    msg = sprintf(['ct_read_profile: 不再支持"剖面结构体"输入。\n' ...
                   '  本程序唯一的剖面输入是标准 txt 文件: 请把字段写成 "名字 = 数值" 的行,\n' ...
                   '  保存为 txt 后用 ct_read_profile(''xxx.txt'') 读取(或用 ct_write_profile 生成)。']);
    error('ct_read_profile:structInput', '%s', msg);
end
if ~(ischar(fname) || (isstring(fname) && strlength(fname) > 0))
    error('ct_read_profile:badFileName', '文件名必须是字符型, 当前类型为 %s', class(fname));
end
fname = char(fname);

% ---------------- 1. 定位文件 ----------------
[ok, fpath] = ct_locate_profile(fname);
if ~ok
    msg = sprintf(['ct_read_profile: 找不到剖面文件 "%s"\n' ...
                   '  已尝试: ① 当前目录  %s\n' ...
                   '          ② 程序目录  %s\n' ...
                   '  请检查路径/文件名, 或把剖面文件放到上述任一目录。'], ...
                  fname, fullfile(pwd, fname), fullfile(fileparts(mfilename('fullpath')), fname));
    error('ct_read_profile:noSuchFile', '%s', msg);
end

% ---------------- 2. 逐行解析 ----------------
req = {'rho','R0','B0','a','q','shear','A_ion','theta','n_e','T_e'};
opt = struct('kappa',1.00, 'r',[], 'n_i',[], 'T_i',[], ...
             'dn_e_dr',[], 'dT_e_dr',[], 'dn_i_dr',[], 'dT_i_dr',[], ...
             'RL_ne',[], 'RL_Te',[], 'RL_ni',[], 'RL_Ti',[]);
known = [req, fieldnames(opt)'];

fid = fopen(fpath, 'r');
if fid < 0
    error('ct_read_profile:cannotOpen', '无法打开 "%s" (请检查文件是否被占用/权限)', fpath);
end
% 【必须在解析出错时也关闭文件句柄】解析循环里的 error 会把 fid 泄漏出去:
%   在 Windows 上该文件随后无法删除/覆盖(用户看到的警告"未找到文件或者权限被拒绝"
%   就是这么来的)。因此把整个循环包在 try/catch 里, 出错先 fclose 再 rethrow,
%   保证 message identifier 与文本原样保留。自检 [T17] 用"报错后文件仍可删除"检验这一点。
P = struct();  L = struct();  nline = 0;  ntot = 0;
try
    while ~feof(fid)
        ln = fgetl(fid);
        if ~ischar(ln), break; end
        nline = nline + 1;
        p = strfind(ln, '#');                                   % 去注释
        if ~isempty(p), ln = ln(1:p(1)-1); end
        ln = strtrim(ln);
        if isempty(ln), continue; end
        ntot = ntot + 1;

        eq = strfind(ln, '=');
        if isempty(eq)
            error('ct_read_profile:missingEquals', '%s 第 %d 行缺少 "=": "%s"  (格式: 名字 = 数值)', ...
                  fname, nline, ln);
        end
        key = strtrim(ln(1:eq(1)-1));
        val = strtrim(ln(eq(1)+1:end));

        if isempty(key)
            error('ct_read_profile:emptyKey', '%s 第 %d 行 "=" 左边没有字段名', fname, nline);
        end
        if ~isvarname(key)
            error('ct_read_profile:badFieldName', ...
                  '%s 第 %d 行字段名不合法: "%s"  (只允许字母/数字/下划线, 且不以数字开头)', ...
                  fname, nline, key);
        end
        if ~any(strcmp(key, known))
            sug = ct_suggest(key, known);
            msg = sprintf(['ct_read_profile: %s 第 %d 行出现未知字段 "%s"%s\n' ...
                           '  允许的字段见 ct_read_profile.m 头部字段表。'], ...
                          fname, nline, key, sug);
            error('ct_read_profile:unknownField', '%s', msg);
        end
        if isfield(P, key)
            error('ct_read_profile:duplicateField', ...
                  '%s 第 %d 行重复定义 "%s" (第 %d 行已给)', ...
                  fname, nline, key, L.(key));
        end
        if isempty(val)
            error('ct_read_profile:emptyValue', '%s 第 %d 行的值为空: "%s = "', fname, nline, key);
        end
        v = str2double(val);
        if isnan(v)
            error('ct_read_profile:notANumber', '%s 第 %d 行的值不是数: "%s = %s"', ...
                  fname, nline, key, val);
        end
        if ~isfinite(v)
            error('ct_read_profile:notFinite', '%s 第 %d 行的值不是有限数: "%s = %s"', ...
                  fname, nline, key, val);
        end
        P.(key) = v;
        L.(key) = nline;
    end
catch err
    fclose(fid);                 % 出错也要释放句柄(否则文件无法删除/覆盖)
    rethrow(err);
end
fclose(fid);

if ntot == 0
    error('ct_read_profile:noValidLine', ...
          '%s 中没有任何 "名字 = 数值" 的有效行(是否整个文件都是注释?)', fname);
end

% ---------------- 3. 补齐可选字段 ----------------
fo = fieldnames(opt);
for ii = 1:numel(fo)
    k = fo{ii};
    if ~isfield(P, k), P.(k) = opt.(k); end
end
if isempty(P.r), P.r = P.rho*P.a; end
if isempty(P.n_i), P.n_i = P.n_e; end
if isempty(P.T_i), P.T_i = P.T_e; end

% ---- [v7] 准中性提示(上一轮审查 S8) ----
%   本程序的色散关系是**单一 Z=1 主离子 + 电子**的准中性形式 D = (Γ0_i - 𝒥_i)/τ + (Γ0_e - 𝒥_e),
%   它隐含 n_i = n_e; 而 n_i 本身并不进入 D(只用 τ = T_i/T_e 与 L_ni/L_ne)。文件里给了
%   n_i ≠ n_e 时既不报错也不参与计算, 容易被误读成"已含杂质/多物种"。这里显式告警。
if ~isempty(P.n_e) && ~isempty(P.n_i) && abs(P.n_i - P.n_e) > 1e-8*max(P.n_e, realmin)
    warning('ct_read_profile:quasiNeutrality', ...
            ['n_i = %.6g 与 n_e = %.6g 相差 %.2f%%, 但本程序的色散关系隐含 n_i = n_e ' ...
             '(单一 Z=1 主离子, n_i 不进入 D)。若确有杂质/多物种, 本程序不适用。'], ...
            P.n_i, P.n_e, 100*abs(P.n_i-P.n_e)/P.n_e);
end

% ---------------- 4. 梯度: 二选一, 反算 R/L_X ----------------
pairs = {'n_e','dn_e_dr','RL_ne'; 'T_e','dT_e_dr','RL_Te'; ...
         'n_i','dn_i_dr','RL_ni'; 'T_i','dT_i_dr','RL_Ti'};
for ii = 1:size(pairs,1)
    base = pairs{ii,1};  dfld = pairs{ii,2};  rlfld = pairs{ii,3};
    has_d  = ~isempty(P.(dfld));
    has_rl = ~isempty(P.(rlfld));
    if has_d && has_rl
        msg = sprintf(['ct_read_profile: %s 中 "%s" 与 "%s" 同时给出(第 %d、%d 行)。\n' ...
                       '  同一物理量的两种写法只能选一种: 用 R/L_X 或 用 dX/dr, 以免歧义。'], ...
                      fname, rlfld, dfld, L.(rlfld), L.(dfld));
        error('ct_read_profile:ambiguousGradient', '%s', msg);
    end
    if has_d
        if P.(base) == 0
            msg = sprintf(['ct_read_profile: %s 的 %s = 0, 无法由 %s 反算 R/L_X\n' ...
                           '  请改用 R/L_X 直接给出(该写法不依赖 %s)。'], ...
                          fname, base, dfld, base);
            error('ct_read_profile:zeroBaseForGradient', '%s', msg);
        end
        P.(rlfld) = -P.R0*P.(dfld)/P.(base);      % R/L_X = -R0 (dX/dr)/X
    end
end
miss = {};
for ii = 1:size(pairs,1)
    if isempty(P.(pairs{ii,3})), miss{end+1} = pairs{ii,3}; end   %#ok<AGROW>
end
if ~isempty(miss)
    msg = sprintf(['ct_read_profile: %s 缺少梯度信息: %s\n' ...
                   '  请给出全部四项 R/L_X (RL_ne, RL_Te, RL_ni, RL_Ti), 或用对应的物理梯度\n' ...
                   '  (dn_e_dr, dT_e_dr, dn_i_dr, dT_i_dr) 代替。'], fname, strjoin(miss, ', '));
    error('ct_read_profile:missingGradient', '%s', msg);
end

% ---------------- 5. 必需字段检查 ----------------
mreq = req(~isfield(P, req));
if ~isempty(mreq)
    msg = sprintf(['ct_read_profile: %s 缺少必需字段: %s\n' ...
                   '  必需字段列表: %s'], ...
                  fname, strjoin(mreq, ', '), strjoin(req, ', '));
    error('ct_read_profile:missingRequiredField', '%s', msg);
end

% ---------------- 6. 数值范围 / 自洽性 ----------------
ct_rng(fname, L, P, 'rho',   '>', 0);      ct_rng(fname, L, P, 'rho',   '<', 1);
ct_rng(fname, L, P, 'R0',    '>', 0);      ct_rng(fname, L, P, 'B0',    '>', 0);
ct_rng(fname, L, P, 'a',     '>', 0);      ct_rng(fname, L, P, 'q',     '>', 0);
ct_rng(fname, L, P, 'A_ion', '>', 0);      ct_rng(fname, L, P, 'kappa', '>', 0);
ct_rng(fname, L, P, 'n_e',   '>', 0);      ct_rng(fname, L, P, 'T_e',   '>', 0);
ct_rng(fname, L, P, 'n_i',   '>', 0);      ct_rng(fname, L, P, 'T_i',   '>', 0);
ct_rng(fname, L, P, 'RL_ne', '>', 0);      ct_rng(fname, L, P, 'RL_ni', '>', 0);
ct_rng(fname, L, P, 'RL_Te', '>', 0);      ct_rng(fname, L, P, 'RL_Ti', '>', 0);
ct_rng(fname, L, P, 'r',     '>', 0);      ct_rng(fname, L, P, 'r',     '<', P.a);
if abs(P.theta) > pi
    msg = sprintf(['ct_read_profile: %s 的 theta = %g rad 超出 [-pi, pi] (第 %d 行)\n' ...
                   '  theta 是磁面极向角, 单位 rad。'], ...
                  fname, P.theta, ct_line(L, 'theta'));
    error('ct_read_profile:thetaOutOfRange', '%s', msg);
end
if isfield(L, 'r')
    dphi = abs(P.r - P.rho*P.a);
    if dphi > 0.05*P.a
        msg = sprintf(['ct_read_profile: %s 中 r = %.6g m 与 rho*a = %.6g m 不一致' ...
                       '(相差 %.3g m > 5%% of a)。\n' ...
                       '  r 与 rho 只需给一个: 给 rho 时自动取 r = rho*a; 给 r 时请核对 rho。'], ...
                      fname, P.r, P.rho*P.a, dphi);
        error('ct_read_profile:inconsistentR', '%s', msg);
    end
end

% ---------------- 7. 派生量 ----------------
P.Lne = P.R0/P.RL_ne;   P.Lni = P.R0/P.RL_ni;
P.LTe = P.R0/P.RL_Te;   P.LTi = P.R0/P.RL_Ti;
P.eta_i = P.RL_Ti/P.RL_ni;
P.eta_e = P.RL_Te/P.RL_ne;
P.src_file = fpath;
end

% ---------------------------------------------------------------------------------------
function [ok, fpath] = ct_locate_profile(fname)
%CT_LOCATE_PROFILE  先按给定路径找, 再在程序所在目录找; 返回规范化的绝对路径
ok = false;  fpath = '';
cand = {fname, fullfile(fileparts(mfilename('fullpath')), fname)};
if ~isempty(regexp(fname, '^([A-Za-z]:[\\/]|\\\\)', 'once'))     % 已经是绝对路径
    cand = {fname};
end
for ii = 1:numel(cand)
    if exist(cand{ii}, 'file') == 2
        d = dir(cand{ii});
        fpath = fullfile(d(1).folder, d(1).name);
        ok = true;
        return;
    end
end
end

% ---------------------------------------------------------------------------------------
function ln = ct_line(L, key)
%CT_LINE  取某字段所在行号(字段不存在时返回 0)
ln = 0;
if isfield(L, key), ln = L.(key); end
end

% ---------------------------------------------------------------------------------------
function ct_rng(fname, L, P, key, op, lim)
%CT_RNG  范围检查(报错带行号)
v = P.(key);
bad = false;
switch op
    case '>', bad = ~(v > lim);
    case '<', bad = ~(v < lim);
end
if bad
    msg = sprintf(['ct_read_profile: %s 的 %s = %g 不满足 "%s %g" (第 %d 行)\n' ...
                   '  请检查该量的物理取值范围。'], ...
                  fname, key, v, op, lim, ct_line(L, key));
    error('ct_read_profile:outOfRange', '%s', msg);
end
end

% ---------------------------------------------------------------------------------------
function sug = ct_suggest(key, known)
%CT_SUGGEST  给未知字段名一个"最接近的已知字段"建议(编辑距离 <= 3 时给出)
best = '';  bd = Inf;
for ii = 1:numel(known)
    d = ct_editdist(lower(key), lower(known{ii}));
    if d < bd, bd = d;  best = known{ii}; end
end
if bd <= 3
    sug = sprintf('  (是否想写 "%s" ?)', best);
else
    sug = '';
end
end

% ---------------------------------------------------------------------------------------
function d = ct_editdist(a, b)
%CT_EDITDIST  Levenshtein 距离(用于字段名拼写建议)
na = numel(a);  nb = numel(b);
D = zeros(na+1, nb+1);
D(:,1) = (0:na)';  D(1,:) = 0:nb;
for i = 2:na+1
    for j = 2:nb+1
        c = 1;
        if a(i-1) == b(j-1), c = 0; end
        D(i,j) = min([D(i-1,j)+1, D(i,j-1)+1, D(i-1,j-1)+c]);
    end
end
d = D(na+1, nb+1);
end
