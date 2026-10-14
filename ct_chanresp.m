function S = ct_chanresp(M, spec, omega, want_deriv)
%CT_CHANRESP  单物种的非绝热响应 𝒥_s, 以及它对三个自变量的偏导数
%=========================================================================================
% 本函数是**响应函数的唯一出处**: ct_dispersion.m(组装 D(ω))与 ct_ct_exact.m(组装 C_T)
% 都调用它, 因此"色散关系用到的 𝒥_s"与"C_T 里求导的 𝒥_s"在构造上必然一致 —— 这消除了
% 早期版本里"色散用 kinetic 分支、C_T 却用另一套闭式"的错配(见外部审查 S2)。
%
% ------------------------------- 返回的四个量 -------------------------------
%   S.J   : 𝒥_s(ω; X_s, Y_s), 使得准中性方程写作
%               D(ω) = (Γ_0(b_i) - 𝒥_i)/τ + [1 - H_e],   1 - H_e = Γ_0(b_e) - 𝒥_e
%           (与 ct_dispersion.m 的记账完全一致: 单物种 D 里 𝒥_s 恒为被减项)
%   S.JX  : ∂𝒥_s/∂X_s  在固定 Y_s 下取         (X_s ≡ ω_*s = -K_s/L_ns)
%   S.JY  : ∂𝒥_s/∂Y_s  在固定 X_s 下取         (Y_s ≡ ω_*s η_s = -K_s/L_Ts)
%   S.JW  : ∂𝒥_s/∂ω
%   S.ok  : 偏导数是否可用(单极点色散模型下为 false, 见下)
%
% ------------------------------- 分子约定 [v7 修正] -------------------------------
% **本文件 v7 起默认采用标准局域线性 GK 方程的分子**:
%
%       (ω - k∥v∥ - ω_Ds) h_s = (q_s f_0s/T_s)(ω - ω_*s^T) J_0(k⊥ρ_s)⟨φ⟩
%   ⟹  𝒥_s = ⟨ W (ω - ω_*s^T)/(ω - k∥v∥ - ω_Ds) ⟩ ,   W = Maxwell 权重 × J_0²
%
%   (Ivanov & Adkins, JPP 89, 905890504 (2023): 其式 (3) 的线性化结果就是式 (12), 推导基于
%    Abel et al. (2013) 的 GK 框架 —— 空间漂移算子 (v∥b̂ + v_d)·∇ **只作用在非绝热部分
%    h_s 上**, 不作用在绝热部分 -(q_s⟨χ⟩/T_s)f_0s 上。)
%
% 【v6 及以前的做法, 现由 cfg.j_numerator='legacy' 复现】把该算子也作用到绝热部分, 于是
%  分子里多出一项 -ω_Ds:  𝒥_s = ⟨W(ω - ω_Ds - ω_*s^T)/(ω - k∥v∥ - ω_Ds)⟩(文档① (2.6))。
%  该多出项的实测后果(见 ../复核报告_MATLAB与理论推导.md §2.1):
%    * k∥=0 且梯度为零时 1-H_i ≡ 0(机器零 2.2e-16), 即离子漂移共振通道被整体关闭;
%      标准约定在同一条件给 1-H_i = Γ_0(b_i) = 0.9152;
%    * 同一 ω、同一模型下 D 相差 O(0.6)(与 D 自身同量级);
%    * 默认剖面 k=0.2 的根由 (0.1303,0.0740) 变为 (0.0905,0.0206) c_s/a, k≥0.4 无根。
%  故默认 'standard'; 'legacy' 仅供与 v6 结果逐位对照(自检 [T23a] 专门守住这一约定)。
%
% ------------------------------- 偏导数的显式形式 -------------------------------
% 标准分子由 ct_channel.m 的两个矩 I0,I1 直接给出:
%   𝒥_s      = (ω - X_s) I0 - Y_s I1
%   ∂𝒥/∂X_s  = -I0
%   ∂𝒥/∂Y_s  = -I1
%   ∂𝒥/∂ω    = I0 - (ω - X_s) I0b + Y_s I1b
% 第三式与"商法则直接读出"的形式 X_s I0b + Y_s I1b - Ipb 恒等(用 ω = den + k∥v∥ + ω_D);
% 两种写法都不会出现额外近似, 自检 [T21c] 检验其恒等性。
% legacy 分子在此基础上再加 -I2(值)与 +I2b(ω 导数), I2 = ⟨W ω_D/den⟩。
%
% 【各分支的偏导数】
%   * 离子(kinetic)        : 直接用上面的三式。
%   * 捕获电子(kinetic)     : 同一三式, 但 ω 换成 ω + iν_eff, 再整体乘 FBW。
%   * 捕获电子(bounce)      : 见文件末尾 ct_trap_bounce; v7 起也改用标准分子,
%                             并把"FOW 因子只乘非绝热/共振部分"的记账写清:
%   * 通行电子(kinetic, k∥≠0): 𝒥_pe = (ω - X)A - Y B, A = -Z(ζ)/k∥v_te,
%                             B = (Z + ζZ')/(2k∥v_te), ζ = ω/(k∥v_te);
%                             ∂𝒥/∂X = -A、∂𝒥/∂Y = -B、
%                             ∂𝒥/∂ω = A + (ω-X)∂A/∂ω - Y ∂B/∂ω,
%                             其中 ∂A/∂ω = -Z'/(k∥v_te)²、∂B/∂ω = (Z' - ζZ - ζ²Z')/(k∥v_te)²
%                             (用 Z'' = -2(Z + ζZ'))。
%
% 【为什么单极点模型不给偏导数】disp_model='pole' 时离子与捕获电子的响应来自
%   ct_pole_integral(单极点留数+主值), 它的 ∂/∂ω 需要另一套"双重极点"的解析处理。而
%   PDF 式 (3) 的模移动项本身是**精确响应**的性质, 所以本程序规定: CT_exact 只在
%   disp_model='kinetic' 下给出; 单极点模型仍只给闭式 C_T = 3/2 - u_*。此时 ok=false,
%   调用方(ct_ct_exact.m)会据此把 CT_exact 置为 NaN 并说明原因。
%
% 用法:  S = ct_chanresp(M, 'i', omega, true)
%=========================================================================================

if nargin < 4 || isempty(want_deriv), want_deriv = false; end

S = struct('J', NaN, 'JX', NaN, 'JY', NaN, 'JW', NaN, 'ok', false, 'why', '');

% ---- 分子约定开关 (见上方说明) ----
legacy_num = false;
if isfield(M.cfg,'j_numerator') && ~isempty(M.cfg.j_numerator) ...
        && strcmpi(char(M.cfg.j_numerator), 'legacy')
    legacy_num = true;
end

model = 'kinetic';
if isfield(M.cfg, 'disp_model') && ~isempty(M.cfg.disp_model)
    model = lower(char(M.cfg.disp_model));
end
isIon = strcmpi(spec, 'i');
if isIon
    X = M.wstar_i;  eta = M.G.eta_i;
else
    X = M.wstar_e;  eta = M.G.eta_e;
end
Y = X*eta;

% =======================================================================================
% 单极点约化模型: 只给 𝒥(供 D 使用), 偏导数不支持
% ---------------------------------------------------------------------------------------
% 【v7】单极点模型是**同一个响应函数**的约化, 因此也必须跟着分子约定走:
%   standard: 𝒥_s = ⟨W(ω - ω_*s^T)/(ω - C_s u)⟩_u = -ω·I(0) + ω_*s·I(η_s)
%             (用 1/(ω-Cu) = -(1/C)/(u-u_*) 与 ct_pole_integral 的 I 定义直接得到)
%   legacy  : 𝒥_s = 1 + ω_*s·I(η_s)     (v6 形式; 它对应的分子是 ω-ω_D-ω_*^T)
%   两种形式的差异 = -(2/√π)∫√u e^{-u}·u/(u-u_*)du, 即 O(ω_D/ω) 的项; 而
%   **闭式 C_T = 3/2 - u_* 不受影响**(该差异在冻结模下与梯度无关, 自检 [T22d] 守住)。
% =======================================================================================
if strcmp(model, 'pole')
    if isIon
        Ci = -(4/3)*M.G.lam_i*M.wD0_i;
        Ii = ct_pole_integral(omega, Ci, M.G.eta_i);
        if legacy_num
            S.J = 1 + M.wstar_i*Ii;
        else
            I0i = ct_pole_integral(omega, Ci, 0);
            S.J = -omega*I0i + M.wstar_i*Ii;
        end
    else
        if M.trap_on
            [Jt, ~, ~, ~] = ct_trap_bounce(M, omega, X, Y, [], legacy_num);
        else
            Ce = -(4/3)*M.G.lam_e*M.wD0_e;
            Ie = ct_pole_integral(omega + 1i*M.nu_eff, Ce, M.G.eta_e);
            if legacy_num
                Jt = 1 + M.wstar_e*Ie;
            else
                I0e = ct_pole_integral(omega + 1i*M.nu_eff, Ce, 0);
                Jt = -(omega + 1i*M.nu_eff)*I0e + M.wstar_e*Ie;
            end
            Jt = M.FBW*Jt;
        end
        % 单极点模型下通行电子恒取绝热(与 ct_dispersion.m 原实现逐字一致)
        S.J = (1 - M.G.ft)*(M.adiab(M.b_e) - 1) + M.G.ft*Jt;
    end
    S.why = 'disp_model=''pole'': 单极点闭式不含 ∂𝒥/∂ω 的解析形式, CT_exact 需 kinetic';
    return;
end

% =======================================================================================
% 环形推广(kinetic)模型
% =======================================================================================
if isIon
    E = ct_channel(M, 'i', omega, want_deriv);
    S.J = (omega - X)*E.I0 - Y*E.I1;
    if want_deriv
        S.JX = -E.I0;
        S.JY = -E.I1;
        S.JW = E.I0 - (omega - X)*E.I0b + Y*E.I1b;
        S.ok = true;
    end
    if legacy_num                      % v6 约定: 分子多一项 -ω_D
        S.J = S.J - E.I2;
        if want_deriv, S.JW = S.JW + E.I2b; end
    end
else
    omega_e = omega + 1i*M.nu_eff;          % 捕获通道的碰撞(去捕获)展宽
    % ---- 捕获电子 ----
    if M.trap_on
        [Jt, dX, dY, dW] = ct_trap_bounce(M, omega, X, Y, want_deriv, legacy_num);
    else
        E = ct_channel(M, 'e', omega_e, want_deriv);
        Jt0 = (omega_e - X)*E.I0 - Y*E.I1;
        if legacy_num, Jt0 = Jt0 - E.I2; end
        Jt  = M.FBW*Jt0;
        if want_deriv
            dX = -M.FBW*E.I0;
            dY = -M.FBW*E.I1;
            dW =  M.FBW*(E.I0 - (omega_e - X)*E.I0b + Y*E.I1b);
            if legacy_num, dW = dW + M.FBW*E.I2b; end
        else
            dX = NaN; dY = NaN; dW = NaN;
        end
    end
    % ---- 通行电子 ----
    passingKinetic = true;
    if isfield(M.cfg, 'passing_e'), passingKinetic = ~strcmpi(M.cfg.passing_e, 'adiabatic'); end
    if ~passingKinetic || M.kpar == 0
        Jps = 0;  pX = 0;  pY = 0;  pW = 0;      % 绝热: 常数, 对梯度与 ω 都不贡献
    else
        kvt  = M.kpar*M.G.vte;
        zeta = omega/kvt;
        [Zv, Zpv] = ct_zfun(zeta);
        A = -Zv/kvt;
        B = (Zv + zeta*Zpv)/(2*kvt);
        Jps = (omega - X)*A - Y*B;
        if want_deriv
            dAdw = -Zpv/kvt^2;
            dBdw = (Zpv - zeta*Zv - zeta^2*Zpv)/kvt^2;
            pX = -A;
            pY = -B;
            pW = A + (omega - X)*dAdw - Y*dBdw;
        else
            pX = NaN; pY = NaN; pW = NaN;
        end
    end
    ft = M.G.ft;
    % 1 - H_e = (1-ft)(1 - H_ps) + ft(Γ_0(b_e) - 𝒥_t)  =>  𝒥_e = (1-ft)(Γ_0(b_e) - 1 + 𝒥_ps) + ft 𝒥_t
    S.J = (1 - ft)*(M.adiab(M.b_e) - 1 + Jps) + ft*Jt;
    if want_deriv
        S.JX = (1 - ft)*pX + ft*dX;
        S.JY = (1 - ft)*pY + ft*dY;
        S.JW = (1 - ft)*pW + ft*dW;
        S.ok = true;
    end
end
end

% ---------------------------------------------------------------------------------------
function [Jt, dX, dY, dW] = ct_trap_bounce(M, omega, X, Y, want_deriv, legacy_num)
%CT_TRAP_BOUNCE  文档① §9.4/(9.8) 的捕获区 (x, ξ0) 弹跳平均通道 —— [v7 修正]
%=========================================================================================
% 弹跳平均消去 k∥v∥ 后的**标准**捕获电子响应(标准 GK 分子 ω - ω_*e^T, 见文件头):
%
%   1 - H_e|_trap = f_t · ⟨ J_0² (ω_*e^T - ⟨ω_De⟩_b)/(ω + iν_eff - ⟨ω_De⟩_b) ⟩_t
%
% 其中 ω_*e^T = X + Y(x²-3/2), ⟨ω_De⟩_b = -ω_D0e x² B(ξ0)(文档① (9.6)), 故
%   ω + iν_eff - ⟨ω_De⟩_b = Dm = ω + iν_eff + ω_D0e x² B.
% 有限香蕉宽度因子(FOW)按"只乘非绝热(共振)响应"的记账, 即整体乘 FBW_trap:
%
%       R ≡ FBW_trap ⟨(ω_*e^T - ⟨ω_De⟩_b)/Dm⟩_t
%         = FBW_trap [ 1 - (ω+iν_eff) S0 + X S0 + Y S1 ]        (用 ⟨⟨ω_De⟩_b/Dm⟩ = (ω+iν)S0 - 1)
%
% 为与 ct_chanresp 的 (Γ_0(b_e) - 𝒥_t) 记账接口兼容, 本函数返回的是
%       𝒥_t = Γ_0(b_e) - R ,
% 于是 1 - H_e|_trap = f_t (Γ_0 - 𝒥_t) = f_t · R 恰好落在标准响应上。
%
% 【v6 及以前(现由 M.cfg.j_numerator/'legacy' 复现)的两处错误】
%   ① 驱动符号相反: 实现成 𝒥_t = FBW[1 + X S0 + Y S1], 代进 (Γ_0 - 𝒥_t) 后驱动项变成
%      -FBW⟨ω_*^T/Dm⟩, 与该分支引用的文档① (9.8)(括号内为正)相反 —— 等于把 TEM 的
%      驱动换成了阻尼。实测同一剖面 γ 差 2–8 倍(k=0.4: 0.0107 → 0.0835 c_s/a)。
%   ② 多出一项伪常数 f_t(Γ_0(b_e) - FBW_trap): 使 trapped 响应在 ω→∞ 时趋于该常数而
%      不趋于零(标准响应必须 →0)。两处都已改正; [T23b] 专门守住这两条。
%
% 【为什么不用 v6 结果"校准"FOW 因子】README §7.4b 用 FBW_trap='banana_avg' 与 TGLF 对齐
%   到 2%, 而符号错误很可能被该可调因子吸收 —— 先定标、再用同一组数据验证无法自证。
%   v7 把响应写成标准形式后, FOW 因子必须重新由轨道宽度平均给出/重新标定, 不能沿用旧值
%   当作已验证的物理量(见 ../修订说明_CORRECTIONS.md §3)。
%=========================================================================================
if nargin < 6 || isempty(legacy_num), legacy_num = false; end

Tt = M.trap;
cx = M.wD0_e*Tt.x.^2;                       % (N_x×1)   (= -⟨ω_De⟩_b / (x²B) 的符号约定见上)
Dm = (omega + 1i*M.nu_eff) + cx*Tt.B.';     % (N_x×N_xi)
gx = sum((Tt.wxi.')./Dm, 2);                % (N_x×1)   Σ_ξ w_ξ / Dm
S0 = sum(Tt.wx.*gx);
S1 = sum(Tt.wx.*(Tt.x.^2 - 1.5).*gx);
dX = NaN; dY = NaN; dW = NaN;

if legacy_num
    % ---- v6 约定(仅作逐位对照) ----
    Jt = M.FBW_trap*(1 + X*S0 + Y*S1);
    if want_deriv
        gxb = sum((Tt.wxi.')./Dm.^2, 2);
        S0b = sum(Tt.wx.*gxb);
        S1b = sum(Tt.wx.*(Tt.x.^2 - 1.5).*gxb);
        dX = M.FBW_trap*S0;
        dY = M.FBW_trap*S1;
        dW = -M.FBW_trap*(X*S0b + Y*S1b);
    end
    return;
end

% ---- 标准约定: R = FBW[1 - (ω+iν)S0 + X S0 + Y S1],  𝒥_t = Γ_0(b_e) - R ----
G0e = M.adiab(M.b_e);
Jt  = G0e - M.FBW_trap*(1 - (omega + 1i*M.nu_eff)*S0 + X*S0 + Y*S1);
if want_deriv
    gxb = sum((Tt.wxi.')./Dm.^2, 2);
    S0b = sum(Tt.wx.*gxb);                        % ∂S0/∂ω = -S0b
    S1b = sum(Tt.wx.*(Tt.x.^2 - 1.5).*gxb);       % ∂S1/∂ω = -S1b
    % ∂R/∂X = FBW·S0 ;  ∂R/∂Y = FBW·S1 ;
    % ∂R/∂ω = FBW[-(S0 + (ω+iν)∂S0/∂ω) + X ∂S0/∂ω + Y ∂S1/∂ω]
    %       = FBW[-S0 + (ω+iν-X) S0b - Y S1b]
    dX = -M.FBW_trap*S0;
    dY = -M.FBW_trap*S1;
    dW =  M.FBW_trap*(S0 - ((omega + 1i*M.nu_eff) - X)*S0b + Y*S1b);
end
end
