function [x, w] = ct_gaussleg(n, a, b)
%CT_GAUSSLEG  n 点 Gauss-Legendre 求积节点与权重(区间 [a,b])
%=========================================================================================
% 用 Jacobi 矩阵的特征分解实现(仅依赖 eig, 无需任何工具箱):
%   J = diag(beta_1..beta_{n-1}) + 次对角线,   beta_k = k/sqrt(4k^2-1)
%   节点 = eig(J),  权重 = 2*(第一行特征向量分量)^2
% 校验: sum(w) = b-a ; sum(w.*x.^2) = (b^3-a^3)/3 (见 check_CT_toroidal.js 的 [T1])
%
% 用法:  [x, w] = ct_gaussleg(n, a, b)
%=========================================================================================

k = (1:n-1)';
beta = k./sqrt(4*k.^2 - 1);
J = diag(beta, 1) + diag(beta, -1);
[V, D] = eig(J);
x = diag(D);
w = 2*(V(1,:)').^2;

% 映射到 [a,b]
c1 = 0.5*(b - a);  c2 = 0.5*(b + a);
x = c2 + c1*x;
w = c1*w;

% 升序排列, 便于与 JS 镜像实现逐点比对
[x, idx] = sort(x);
w = w(idx);
end
