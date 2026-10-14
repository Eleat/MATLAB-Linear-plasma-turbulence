function [f, fp, fpp, fppp] = ct_fz(z)
%CT_FZ  文档②(二次分母详解) 式 (2.1)(2.2) 的 f(z) 及其导数
%=========================================================================================
%   f(z) ≡ (1/sqrt(pi)) ∫_{-inf}^{inf} e^{-t^2}/(t-z) dt        ("沿实轴积分"的那个函数")
%
%   f(z) = Z(z)                           (Im z > 0)
%   f(z) = Z(z) - 2i sqrt(pi) e^{-z^2}    (Im z < 0)                —— 文档② (2.2)
%
% 【为什么要写这个函数】
%   二次分母的积分有**精确闭式**: 把分母因子化后部分分式, 每个极点用 f 而不是 Z
%   (文档② §2.3-2.4)。Z 是整函数、f 在实轴上有跳跃, 二者只在上半平面相同;
%   若把下半平面的根也代进整函数 Z, 误差是 O(1) 而不是小量 (文档② §2.5)。
%
% 【下半平面用反射恒等式, 不用相减】
%   由 Z 的反射关系 Z(z) = 2i sqrt(pi) e^{-z^2} - Z(-z) 可得
%
%       f(z) = -Z(-z)          (Im z < 0)          <-- 本文件采用
%
%   这与 Z(z) - 2i sqrt(pi) e^{-z^2} 恒等, 但**数值上稳定**: 当 |Im z| 较大时
%   e^{-z^2} ~ e^{(Im z)^2} 是指数大的量, 直接相减会让两个大数抵消殆尽
%   (例如 Im z = -5 时, 两项各 ~1e11, 差却只有 O(0.1), 双精度下全部有效位丢失)。
%   反射形式只调用 Im 为正的 Z, 而 ct_zfun 在该半平面是稳定的。
%
% 【校验】(ct_selftest 的 [T12])
%   * f 与"实轴直接数值求积"一致 (文档② §6 推荐的黄金基准);
%   * f(z) = -conj(f(conj(z))) 之类的基本对称性;
%   * 大 |Im z| 时反射形式仍然准确 (直接相减的形式在此处会崩掉)。
%
% 用法:  f = ct_fz(z)          % f(z)
%        [f, fp] = ct_fz(z)    % 同时给 f'(z) = Z'(z) (Im z>0) / Z'(-z) (Im z<0)
%        [f, fp, fpp, fppp] = ct_fz(z)   % [v6] 再给 f''、f'''
%
% 【v6 新增 f''、f''' 的用途】C_T 的模移动项需要 ⟨…⟩/den² 这一族积分(见 ct_xpar.m
%   的 M0b/M1b/M2b)。在**重根**退化情形下它们要落到 f''、f''' 上, 而 Z 的高阶导数
%   由 ODE 递推给出:
%        Z'  = -2(1 + z Z)
%        Z'' = -2(Z + z Z')
%        Z'''= -2(2 Z' + z Z'')
%   下半平面仍用反射关系 f(z) = -Z(-z), 于是 f'' = -Z''(-z)、f''' = +Z'''(-z)。
%=========================================================================================

if nargin < 1 || isempty(z)
    error('ct_fz: 需要输入 z');
end
if isreal(z), z = complex(z); end        % 只用 isreal(基础函数)

if imag(z) >= 0
    if nargout > 2
        [Z, Zp] = ct_zfun(z);
        f = Z;  fp = Zp;
        Zpp = -2*(Z + z*Zp);
        fpp = Zpp;
        if nargout > 3, fppp = -2*(2*Zp + z*Zpp); end
    elseif nargout > 1
        [f, fp] = ct_zfun(z);
    else
        f = ct_zfun(z);
    end
else
    % Im z < 0 : f(z) = -Z(-z) ,  f'(z) = Z'(-z) ,  f''(z) = -Z''(-z) ,  f'''(z) = Z'''(-z)
    zz = -z;
    if nargout > 2
        [Z, Zp] = ct_zfun(zz);
        f = -Z;  fp = Zp;
        Zpp = -2*(Z + zz*Zp);
        fpp = -Zpp;
        if nargout > 3, fppp = -2*(2*Zp + zz*Zpp); end
    elseif nargout > 1
        [Z, Zp] = ct_zfun(zz);
        f = -Z;  fp = Zp;
    else
        f = -ct_zfun(zz);
    end
end
end
