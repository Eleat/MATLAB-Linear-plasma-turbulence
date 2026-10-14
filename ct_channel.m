function I = ct_channel(M, spec, omega, want_deriv)
%CT_CHANNEL  单通道的三个(以及 den² 族的三个)速度空间积分
%=========================================================================================
%   I0 = < W / (ω - k∥v∥ - ω_D) >
%   I1 = < W (x²-3/2) / (ω - k∥v∥ - ω_D) >
%   I2 = < W ω_D / (ω - k∥v∥ - ω_D) >
%   W  = Maxwell 权重 × J_0²(k⊥ρ_s)
%
% 【方法】ω_D 含 x∥², 因此共振分母对 x∥ 是**二次**的 —— 这正是"二次分母问题"。
%   本函数按文档② §2 的**精确**路线处理它:
%     对每个固定的 x⊥, 把 x∥ 方向的积分解析地做出(ct_xpar.m: 二次式因子化 + 部分
%     分式 + 分支正确的 f 函数), 只对 x⊥ 做一维 Gauss-Legendre 求积。
%   于是共振位置的分支结构被精确跟踪; x∥ 方向没有网格、也没有可调的路径参数
%   (旧版的"下移 σ 路径"已删除, 见 ct_model.m 文件头与自检 [T14])。
%   剩下的唯一数值参数是 x⊥ 节点数 cfg.N_perp(自检 [T16] 检验其权重归一化与 Γ0 恒等式;
%   收敛性见程序说明 §6.7)。
%
% ------------------------------- [v6] want_deriv = true -------------------------------
% 额外给出 den² 族的三个矩(want_deriv=false 时它们为 NaN, 不产生额外开销):
%   I0b = < W / den² > ,  I1b = < W (x²-3/2) / den² > ,  I2b = < W ω_D / den² >
%   Ipb = < W k∥v_t x∥ / den² >                       (供内部恒等式交叉检验)
% 它们的用途是 C_T 里的**模移动项**: 由 (1) 的商法则,
%   ∂𝒥_s/∂ω = ⟨ J0² [X_s + Y_s(u-3/2) - k∥v_ts x∥] / den² ⟩ = X_s·I0b + Y_s·I1b - Ipb
% 而用 ω = den + k∥v∥ + ω_D 可把它改写成不需要 Ipb 的等价形式
%   ∂𝒥_s/∂ω = I0 - (ω - X_s)·I0b + Y_s·I1b + I2b
% (自检 [T21c] 直接检验这两个写法恒等)。本函数两种都给出, 组装见 ct_chanresp.m。
%
% 用法:  I = ct_channel(M, spec, omega)                  % 只要 I0,I1,I2
%        I = ct_channel(M, spec, omega, true)            % 同时给 I0b,I1b,I2b,Ipb
%=========================================================================================

if nargin < 4 || isempty(want_deriv), want_deriv = false; end

if strcmpi(spec, 'i')
    kpar_vt = M.kpar*M.G.vti;
    j2  = M.j0i;
    wD0 = M.wD0_i;  geom = M.G.lam_i;
else
    kpar_vt = 0;                                             % 弹跳平均: 无平行流
    j2  = M.j0e;
    wD0 = M.wD0_e;  geom = M.G.lam_e;
end
Gs = wD0*geom;                       % ω_D = -Gs (2 x∥² + x⊥²)
yp = M.xperp;
Ay = M.wperp.*j2/sqrt(pi);           % x⊥ 权重(含 x∥ 归一化 1/sqrt(pi))

I0 = complex(0,0);  I1 = complex(0,0);  I2 = complex(0,0);
I0b = NaN;  I1b = NaN;  I2b = NaN;  Ipb = NaN;
if want_deriv
    I0b = complex(0,0);  I1b = complex(0,0);  I2b = complex(0,0);  Ipb = complex(0,0);
end
for ii = 1:numel(yp)
    S  = ct_xpar(omega, kpar_vt, Gs, yp(ii));
    y2 = yp(ii)^2;
    I0 = I0 + Ay(ii)*S.M0;
    I1 = I1 + Ay(ii)*((S.M2 + y2*S.M0) - 1.5*S.M0);
    I2 = I2 + Ay(ii)*(-Gs)*(2*S.M2 + y2*S.M0);
    if want_deriv
        I0b = I0b + Ay(ii)*S.M0b;
        I1b = I1b + Ay(ii)*((S.M2b + y2*S.M0b) - 1.5*S.M0b);
        I2b = I2b + Ay(ii)*(-Gs)*(2*S.M2b + y2*S.M0b);
        Ipb = Ipb + Ay(ii)*kpar_vt*S.M1b;
    end
end

I = struct('I0', I0, 'I1', I1, 'I2', I2, 'I0b', I0b, 'I1b', I1b, 'I2b', I2b, 'Ipb', Ipb);
end
