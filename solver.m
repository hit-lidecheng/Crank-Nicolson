function C = solver(leftX, rightX, numX, initialT, endT, numT, icfun)
% 此方程用于求解碳浓度分布.
% C为碳浓度分布矩阵,存储着各个位置各个时间的碳浓度值.
% 欲使用此方程求解碳浓度分布,需给出以下值:
% 左边界位置leftX,右边界位置rightX,X轴格点数numX;
% 初始时间initialT,结束时间endT,时间轴格点数numT;
% 初始条件方程icfun.

% x轴格点
x = linspace(leftX, rightX, numX);
% t轴格点
t = linspace(initialT, endT, numT);
% x轴步长
stepX = (rightX - leftX) / (numX - 1);
% t轴步长
stepT = (endT - initialT) / (numT - 1);

% 扩散系数
D = 2e-11;

% 初始化C矩阵
C = zeros(numT, numX);

% 给C矩阵赋初始值
C(1,:) = icfun(x);

% 未考虑初始条件和边界条件时的系数矩阵
r = 1 / 2 * D * stepT / stepX^2;
e = ones(numX,1);
S = spdiags([e -2*e e], -1:1, numX, numX);
A = eye(numX) - r * S;
B = eye(numX) + r * S;

% ------------第一类边界条件开始------------
%b = zeros(numX,1);
% 处理左边界值,此处左边界值为0.012
%A(1,1) = 1;
%A(1,2) = 0;
%B(1,1) = 0;
%B(1,2) = 0;
%b(1) = 0.012;
% 处理右边界值,此处右边界值为0.001
%A(end,end) = 1;
%A(end,end-1) = 0;
%B(end,end) = 0;
%B(end,end-1) = 0;
%b(end) = 0.001;
% ------------第一类边界条件结束------------

% ------------第三类边界条件开始------------
b = zeros(numX,1);
% 处理左边界值,此处左边界满足第三类边界条件
beta = 3e-8;
Cp = 0.1;
A(1,1) = 2 * (stepX)^2 / D / stepT + 2 * beta * stepX / D + 2;
A(1,2) = -2;
B(1,1) = 2 * (stepX)^2 / D / stepT - 2 * beta * stepX / D - 2;
B(1,2) = 2;
b(1) = 4 * beta * stepX * Cp / D;
% 处理右边界值,此处右边界值为0.001
A(end,end) = 1;
A(end,end-1) = 0;
B(end,end) = 0;
B(end,end-1) = 0;
b(end) = 0.001;

% 边界条件为第三类边界条件时扩散方程的解析解
C0 = 0.001;
sol = (Cp - C0) * (erfc(x / sqrt(4 * D * endT)) - exp((beta * x + beta^2 * endT) / D) .* erfc(x / sqrt(4 * D * endT) + beta * sqrt(endT / D))) + C0;
% ------------第三类边界条件结束------------

% ------------C矩阵求解部分开始-------------
for n = 1:numT-1
    C(n+1,:) = (A \ (B * C(n,:)' + b))';
end
% ------------C矩阵求解部分结束-------------

% 绘图部分
figure(1)
plot(x,C(end,:),'*',x,sol);
figure(2)
plot(t,C(:,1));


end