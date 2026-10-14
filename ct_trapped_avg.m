function T = ct_trapped_avg(eps, N_x, N_xi, shat, alpha)
%CT_TRAPPED_AVG  捕获电子区 (能量 x, 投角 ξ₀) 的**弹跳平均**求积表（[v8.4] 含 s-α 效应）
%=========================================================================================
% 实现《热扩散pinch_回旋动理学到CT系数_逐步推导.md》§9.4（弹跳平均完整推导）与 §9.5 式 (9.8)：
%
%   δn_te/n_e = f_t (eδφ/T_e) ⟨ (ω_*e^T + i ν_eff) / (ω - ⟨ω_De⟩_b + i ν_eff) ⟩_t
%
% 本函数只准备**与 k 无关**的那部分表(几何与投角结构), 逐 k 的求和放在 ct_dispersion.m。
%
% ------------------------------- 几何与变量 -------------------------------
% B(θ) = B₀(1 - ε cosθ);  陷阱参数(文档 (9.4a)):
%     κ² = (1/λ̂ - 1 + ε)/(2ε) ∈ [0,1],   λ̂ = λB₀ = v⊥²/(v²B₀)
% 在 θ=0 处取 ξ₀ ≡ v∥/v ∈ [-ξ_t, ξ_t],  ξ_t = sqrt(2ε/(1+ε)) —— 这正是文档 (9.7) 的
% 捕获份额 f_t(ε→0); 由 ξ₀ 反解:
%     λ̂ = (1-ξ₀²)/(1-ε),   κ² = (1-ε)ξ₀² / [2ε(1-ξ₀²)]
% 于是**投角分布的权重就是 ξ₀ 上的均匀权重**(速度空间测度 d³v = 2πv²dv dξ₀), 归一化后
%     Σ_b w_ξ(b) = 1,   且以 w_ξ 求平均时 ⟨1⟩_t = 1、捕获份额自动等于 ξ_t = f_t ✔
%
% ------------------------------- 弹跳平均 -------------------------------
% 文档 (9.5):  ⟨A⟩_b = (1/K)∫_0^κ A(s) ds/[sqrt(1-s²) sqrt(κ²-s²)],  s = sin(θ/2)
% 用 s = κ sinφ 化为无奇点的形式(φ ∈ [0,π/2]):
%     ⟨A⟩_b = (1/K)∫_0^{π/2} A(κ sinφ) dφ/sqrt(1-κ²sin²φ),   K = ∫dφ/sqrt(1-κ²sin²φ)
%
% ================== [v8.4] s-α 几何: α 效应如何进入 ==================
% s-α 几何里磁漂移的极向结构是(见 Kim & Han 2017 式 (5) 之后的定义)
%     g(θ) = cosθ + sinθ·(ŝθ − α sinθ) = cosθ + ŝθ sinθ − α sin²θ ,
% 其中 **α = −q²R dβ/dr**(β = 8πp/B²)是压力梯度参数(Shafranov 位移效应)。
% 关键: **θ=0 处 g = 1, α 与 ŝ 项都消失** ⇒ 本程序(单 θ 局域)在默认 θ=0 下,
% 漂移、FLR 里的 α 只能通过**捕获电子的弹跳平均**进入 —— 而弹跳平均是沿整条
% 轨道(θ ∈ [−θ_b, θ_b])做的, α 在这一步不会被 θ=0 抹掉:
%
%     ⟨ω_De⟩_b = −ω_D0e x² B(ξ₀),
%     B(ξ₀) = (2−λ̂)·⟨g⟩_b + λ̂ε·⟨cosθ·g⟩_b           [v8.4 一般形式]
%          = (2−λ̂)M₁ + λ̂εM₂c ,
%     M₁  = ⟨cosθ⟩_b + ŝ⟨θsinθ⟩_b − α⟨sin²θ⟩_b
%     M₂c = ⟨cos²θ⟩_b + ŝ⟨θsinθcosθ⟩_b − α⟨sin²θcosθ⟩_b
%
% 【v8.3 及以前(legacy)】B 只取了 ⟨cosθ⟩_b 与 ⟨cos²θ⟩_b 两项 —— 即**漏掉了 ŝθsinθ 与
%   α sin²θ 两个几何项**, α 效应完全无法进入。`cfg.trap_precession='legacy'` 可复现旧行为;
%   `'full'`(或 cfg.alpha ≠ 0 时的 'auto')用上面的完整形式。
%   α > 0 使 ⟨g⟩_b 变小(甚至**变号**)——这正是 Ma et al. 2023 所说"弱负剪切 + 大 α 使捕获
%   电子进动方向反转、从而与 ITG 共振"的机制; 在 CBC(ŝ=+0.8, α=0)上这一项为零。
%
% 数值方法: 对每个 ξ₀ 节点用 **Gauss-Legendre 在 φ∈[0,π/2] 上**求积(192 点)。
%   因为 ξ₀ 取的是 GL 节点(不含端点), κ² < 1 严格成立 ⇒ 被积函数光滑, 无端点奇性;
%   K 用同一规则算并与 AGM 椭圆积分 ct_ellip_ke 的结果比对(自检里守住)。
%
% ------------------------------- 能量节点 -------------------------------
% x = v/v_te, 权重 (2/sqrt(pi)) x² e^{-x²} dx (使 Σ w_x = 1); 上限取 6 (e^{-36} 可忽略)。
%
% ------------------------------- 输出 -------------------------------
%   T.x    (N_x×1)   能量节点;      T.wx (N_x×1)  能量权重(Σ = 1)
%   T.B    (N_xi×1)  投角因子(含 ŝ, α);  T.B_legacy (N_xi×1) 旧形式(仅 cosθ, cos²θ)
%   T.wxi  (N_xi×1)  投角权重(Σ = 1);    T.xi (N_xi×1) ξ₀ 节点
%   T.lam, T.kap, T.kap2, T.xi_t, T.eps;  T.shat, T.alpha, T.precession
%   T.g_deep / T.g_marg   深捕获(B@κ→0)与边际捕获(B@κ→1)的 ⟨g⟩_b, 供判断进动方向
%
% 用法:  T = ct_trapped_avg(eps, 32, 24)               % legacy(ŝ=α=0)
%        T = ct_trapped_avg(eps, 32, 24, 0.8, 1.5)      % s-α, ŝ=0.8, α=1.5
%=========================================================================================

if nargin < 2 || isempty(N_x),   N_x  = 32; end
if nargin < 3 || isempty(N_xi),  N_xi = 24; end
if nargin < 4 || isempty(shat),  shat = 0;  end
if nargin < 5 || isempty(alpha), alpha = 0; end
eps = max(eps, 1e-8);
precession = 'full';
if shat == 0 && alpha == 0, precession = 'legacy'; end      % 逐位复现旧行为

% ---- 能量节点(麦氏权重) ----
[xg, wg] = ct_gaussleg(N_x, 0, 6);
T.x  = xg(:);
T.wx = (2/sqrt(pi))*xg(:).^2.*exp(-xg(:).^2).*wg(:);
T.wx = T.wx/sum(T.wx);                       % 归一化(使 Σ=1, 与解析值 1 只差舍入)

% ---- 投角节点(ξ₀ 上均匀) ----
T.xi_t = sqrt(2*eps/(1 + eps));               % 捕获边界(θ=0 处)
[xg2, wq] = ct_gaussleg(N_xi, -T.xi_t, T.xi_t);
T.xi  = xg2(:);
T.wxi = wq(:)/sum(wq(:));                     % 归一化: Σ = 1

% ---- λ̂ 与 κ² ----
T.lam  = (1 - T.xi.^2)/(1 - eps);
kap2   = (1 - eps)*T.xi.^2./(2*eps*(1 - T.xi.^2));
T.kap2 = min(max(kap2, 0), 1);                % 数值上夹到 [0,1]
T.kap  = sqrt(T.kap2);

% ---- 弹跳平均的两个矩 + B ----
T.B        = zeros(N_xi,1);
T.B_legacy = zeros(N_xi,1);
T.g_deep   = zeros(N_xi,1);   % ⟨g⟩_b 本身(诊断: 是否变号)
if strcmp(precession,'legacy')
    for b = 1:N_xi
        [cb, c2b] = cos_moments(T.kap2(b));
        T.B_legacy(b) = (2 - T.lam(b))*cb + T.lam(b)*eps*c2b;
        T.B(b)        = T.B_legacy(b);
        T.g_deep(b)   = cb;
    end
else
    nq = 192;
    [phi, wphi] = ct_gaussleg(nq, 0, pi/2);
    T.mom = struct('cth',zeros(N_xi,1), 'c2',zeros(N_xi,1), 'a3',zeros(N_xi,1), ...
                   'a4',zeros(N_xi,1), 'a5',zeros(N_xi,1), 'a6',zeros(N_xi,1), ...
                   'K_num',zeros(N_xi,1), 'K_agm',zeros(N_xi,1));
    for b = 1:N_xi
        m  = T.kap2(b);
        [cb, c2b] = cos_moments(m);                       % 解析值(用于对照)
        [~, Ek]   = ct_ellip_ke(max(m,realmin));
        T.B_legacy(b) = (2 - T.lam(b))*cb + T.lam(b)*eps*c2b;
        % 数值弹跳平均: ⟨A⟩_b = Σ w_φ A(θ)/√(1−κ²sin²φ) / K
        s   = sqrt(m)*sin(phi);
        w   = wphi(:)./sqrt(max(1 - m*sin(phi(:)).^2, realmin));
        th  = 2*asin(min(max(s, -1), 1));
        K   = sum(w);
        cth = sum(w.*cos(th))/K;
        c2  = sum(w.*cos(th).^2)/K;
        a3  = sum(w.*(th.*sin(th)))/K;
        a4  = sum(w.*(sin(th).^2))/K;
        a5  = sum(w.*(th.*sin(th).*cos(th)))/K;
        a6  = sum(w.*(sin(th).^2.*cos(th)))/K;
        M1  = cth + shat*a3 - alpha*a4;
        M2c = c2  + shat*a5 - alpha*a6;
        T.B(b)      = (2 - T.lam(b))*M1 + T.lam(b)*eps*M2c;
        T.g_deep(b) = M1;
        T.mom.cth(b)=cth; T.mom.c2(b)=c2; T.mom.a3(b)=a3; T.mom.a4(b)=a4;
        T.mom.a5(b)=a5;   T.mom.a6(b)=a6;
        T.mom.K_num(b)=K; T.mom.K_agm(b)=ct_ellip_K_only(m, Ek);
    end
end

T.eps = eps;  T.shat = shat;  T.alpha = alpha;  T.precession = precession;
T.g_marg = T.g_deep(end);                     % 最靠边际捕获的节点(κ→1)
end

% ---------------------------------------------------------------------------------------
function K = ct_ellip_K_only(m, ~)
%CT_ELLIP_K_ONLY  第一类完全椭圆积分 K(m)(AGM), 仅取 K 分量
[K, ~] = ct_ellip_ke(m);
end

% ---------------------------------------------------------------------------------------
function [cb, c2b] = cos_moments(m)
%COS_MOMENTS  ⟨cosθ⟩_b 与 ⟨cos²θ⟩_b 的解析式(文档 (9.6); v8.3 及以前的实现)
%   注: m = κ²; 用 AGM 椭圆积分, 对 κ→1 也稳定。
if m <= 0
    cb = 1;  c2b = 1;  return;
end
[Kk, Ek] = ct_ellip_ke(m);
I1 = (Kk - Ek)/m;
I2 = (2*(1 + m)*I1 - Kk)/(3*m);
s2 = (Kk - Ek)/Kk;
s4 = m*m*I2/Kk;
cb  = 2*Ek/Kk - 1;
c2b = 1 - 4*s2 + 4*s4;
end
