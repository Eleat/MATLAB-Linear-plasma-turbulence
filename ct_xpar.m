function S = ct_xpar(omega, kpar_vt, Gs, yperp)
%CT_XPAR  固定 x⊥ 时把 x∥ 方向的速度空间积分**解析地**做出来
%=========================================================================================
% 这是《DGH问题_磁漂移共振二次分母_详解.md》(文档②) §2 的精确路线:
% 共振分母对 x∥ 是二次的 —— 不要用任何"下移路径/解析延拓"的近似去绕它,
% 而是把二次式因子化、部分分式, 每个极点用分支正确的 f 函数 (ct_fz.m)。
%
% ------------------------------- 被积式 -------------------------------
%   den(x) = omega - kpar_vt x - wD(x) ,   wD = -Gs (2 x² + yperp²)
%          = 2 Gs x² - kpar_vt x + (omega + Gs yperp²)                 [Gs = ω_D0 g_s]
%          ≡ a2 x² + a1 x + a0
%
% 需要三个高斯权重的矩积分(下式的积分沿实轴, 极点不在实轴上时无需主值):
%   M0 = ∫ dx e^{-x²} / den
%   M1 = ∫ dx x e^{-x²} / den
%   M2 = ∫ dx x² e^{-x²} / den
%
% ------------------------------- 三种情形 -------------------------------
% (a) 二次 (a2 ~= 0, 两根不重合): 部分分式 + f
%       1/den = [1/(a2(r+ - r-))] ( 1/(x-r+) - 1/(x-r-) )
%       ∫ e^{-x²}/(x-z) dx = sqrt(pi) f(z)
%       M0 = pref (J+ - J-) ,  pref = 1/(a2(r+ - r-)) ,  J± = sqrt(pi) f(r±)
%       M1 = pref ( r+ J+ - r- J- )
%       M2 = ( sqrt(pi) - a1 M1 - a0 M0 ) / a2
% (b) 两根重合 (判别式 ≈ 0): 用 f 的导数 (ct_fz 的第二个输出)
% (c) 线性 (a2 = 0, a1 ~= 0): den = a1(x - z0), z0 = -a0/a1
%       M0 = sqrt(pi) f(z0)/a1 ,  M1 = sqrt(pi)(1 + z0 f(z0))/a1 ,
%       M2 = sqrt(pi)(z0 + z0² f(z0))/a1
% (d) 常数 (a2 = a1 = 0, 例如 no_curv 且 k∥ = 0): den = a0 = omega
%       M0 = sqrt(pi)/omega ,  M1 = 0 ,  M2 = sqrt(pi)/(2 omega)
%
% 输出 S: 结构体, 字段 M0, M1, M2
%
% ------------------------------- [v6] den² 族矩 M0b/M1b/M2b -------------------------------
% 为在 C_T 中计入**模频率随梯度的移动**(PDF 式 (3) 的模移动项), 需要
%     M0b = ∫ dx e^{-x²} / den² ,  M1b = ∫ dx x e^{-x²} / den² ,  M2b = ∫ dx x² e^{-x²} / den²
% 它们恰好是 M0/M1/M2 对 ω 的导数取负: ∂M*/∂ω = -M*b (因为 ∂den/∂ω = 1)。
% 这里仍然**解析地**给出(不用数值差分): 由 1/den 的部分分式再平方, 每一项都落到
% f、f' 上; 例如两根情形
%   1/den² = pref² [ 1/(x-r+)² + 1/(x-r-)² - (2/(r+-r-))(1/(x-r+) - 1/(x-r-)) ]
% 而 ∫ e^{-x²}/(x-z) dx = √π f(z)、∫ e^{-x²}/(x-z)² dx = √π f'(z), 再用
%   ∫ x e^{-x²}/(x-z)   dx = √π (1 + z f)
%   ∫ x² e^{-x²}/(x-z)  dx = √π (z + z² f)
%   ∫ x e^{-x²}/(x-z)²  dx = √π (f + z f')
%   ∫ x² e^{-x²}/(x-z)² dx = √π (1 + 2 z f + z² f')
% 即得。重根情形需要 f''、f'''(ct_fz.m 已提供)。
% 自检 [T21b] 用 ∂M*/∂ω = -M*b 的数值差分独立复核这一族公式。
%
% 用法:  S = ct_xpar(omega, kpar_vt, Gs, yperp)
%        S.M0, S.M1, S.M2, S.M0b, S.M1b, S.M2b
%
% 【相对旧做法的差别】旧版把 x∥ 沿"下移 sigma 的路径"做二维数值求积来实现所谓的
% 朗道解析延拓。对**二次**分母这是错的: 两个根恒分居实轴两侧, 且 omega 在上半平面
% 内移动时没有任何根穿越实轴, 因而不存在需要延拓的跳跃; 下移路径只会额外拾取条带
% (-sigma, 0) 内极点的残数, 其结果恰好等于文档② §2.5 明确警告的"把根代进整函数 Z"
% 的错误值 (误差 O(1))。本函数从根上消除了该问题, 并且不再需要任何速度空间网格。
%=========================================================================================

if nargin < 4 || isempty(yperp), yperp = 0; end

a2 = 2*Gs;
a1 = -kpar_vt;
a0 = omega + Gs*yperp^2;
SP = sqrt(pi);

if abs(a2) > 0
    dis  = a1^2 - 4*a2*a0;
    % 判别式的相对判据: 两根重合时部分分式退化(0/0)
    sc   = max([abs(a1)^2, abs(4*a2*a0), realmin]);
    if abs(dis) > 1e-12*sc
        sq   = sqrt(dis);
        rp   = (-a1 + sq)/(2*a2);
        rm   = (-a1 - sq)/(2*a2);
        diff = rp - rm;
        pref = 1/(a2*diff);
        [frp, dfrp] = ct_fz(rp);
        [frm, dfrm] = ct_fz(rm);
        Jp   = SP*frp;                  % ∫ e^{-x²}/(x-r+) dx
        Jm   = SP*frm;
        Gp   = SP*dfrp;                 % ∫ e^{-x²}/(x-r+)² dx
        Gm   = SP*dfrm;
        M0   = pref*(Jp - Jm);
        M1   = pref*(rp*Jp - rm*Jm);
        M2   = (SP - a1*M1 - a0*M0)/a2;
        % ---- den² 族 ----
        p2   = pref^2;
        M0b  = p2*( Gp + Gm - (2/diff)*(Jp - Jm) );
        F1p  = SP*(1 + rp*frp);   F1m = SP*(1 + rm*frm);        % ∫ x   e^{-x²}/(x-r)   dx
        G1p  = SP*(frp + rp*dfrp); G1m = SP*(frm + rm*dfrm);    % ∫ x   e^{-x²}/(x-r)²  dx
        M1b  = p2*( G1p + G1m - (2/diff)*(F1p - F1m) );
        F2p  = SP*(rp + rp^2*frp); F2m = SP*(rm + rm^2*frm);    % ∫ x²  e^{-x²}/(x-r)   dx
        G2p  = SP*(1 + 2*rp*frp + rp^2*dfrp);                   % ∫ x²  e^{-x²}/(x-r)²  dx
        G2m  = SP*(1 + 2*rm*frm + rm^2*dfrm);
        M2b  = p2*( G2p + G2m - (2/diff)*(F2p - F2m) );
    else
        % 重根: 1/den = 1/(a2 (x-r)²)  =>  1/den² = 1/(a2² (x-r)⁴)
        r  = -a1/(2*a2);
        [fr, frp, frpp, frppp] = ct_fz(r);
        M0 = (SP/a2)*frp;
        M1 = (SP/a2)*(fr + r*frp);
        M2 = (SP/a2)*(1 + 2*r*fr + r^2*frp);
        % 由 x=(x-r)+r、x²=(x-r)²+2r(x-r)+r² 展开, 配合
        %   ∫e^{-x²}/(x-r)^n dx = √π f^{(n-1)}(r)/(n-1)!
        q  = SP/a2^2;
        M0b = q*(frppp/6);
        M1b = q*(frpp/2 + r*frppp/6);
        M2b = q*(frp + r*frpp + r^2*frppp/6);
    end
elseif abs(a1) > 0
    % 线性分母: den = a1 (x - z0)
    z0 = -a0/a1;
    [fz, fzp] = ct_fz(z0);
    M0 = SP*fz/a1;
    M1 = SP*(1 + z0*fz)/a1;
    M2 = SP*(z0 + z0^2*fz)/a1;
    % den² = a1²(x-z0)²
    q  = SP/a1^2;
    M0b = q*fzp;
    M1b = q*(fz + z0*fzp);
    M2b = q*(1 + 2*z0*fz + z0^2*fzp);
else
    % 常数分母: den = a0 = omega
    M0 = SP/a0;
    M1 = 0;
    M2 = SP/(2*a0);
    M0b = SP/a0^2;
    M1b = 0;
    M2b = SP/(2*a0^2);
end

S = struct('M0', M0, 'M1', M1, 'M2', M2, 'M0b', M0b, 'M1b', M1b, 'M2b', M2b);
end
