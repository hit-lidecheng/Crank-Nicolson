function C = solver(leftX, rightX, numX, initialT, endT, numT, D, beta, Cp, C0)
%  solver函数利用Crank-Nicolson有限差分格式求解扩散方程.
%  函数返回值C为numT*numX的二维矩阵,用于存储各位置各时间的碳浓度值.
%  形参leftX,rightX,numX均为标量,分别表示左边界位置,右边界位置和位置轴格点数量.
%  形参initialT,endT,numT均为标量,分别表示初始时间,终止时间和时间轴格点数量.
%  形参D,beta,Cp均为标量,分别表示扩散系数,传递系数和碳势.
%  形参C0为1*numX的向量,表示各位置处的初始碳浓度值.


deltaX = (rightX - leftX) / (numX - 1);
deltaT = (endT - initialT) / (numT - 1);


C = zeros(numT, numX);
C(1,:) = C0;


% --------未考虑初始条件和边界条件的系数矩阵--------
r = 1 / 2 * D * deltaT / deltaX^2;
e = ones(numX,1);
S = spdiags([e -2*e e], -1:1, numX, numX);
A = eye(numX) - r * S;
B = eye(numX) + r * S;
% --------未考虑初始条件和边界条件的系数矩阵--------


% -----------------第三类边界条件-----------------
b = zeros(numX,1);

% -------------------左边界条件-------------------
A(1,1) = 2 * deltaX^2 / D / deltaT + 2 * beta * deltaX / D + 2;
A(1,2) = -2;
B(1,1) = 2 * deltaX^2 / D / deltaT - 2 * beta * deltaX / D - 2;
B(1,2) = 2;
b(1) = 4 * beta * deltaX * Cp / D;
% -------------------左边界条件-------------------

% -------------------右边界条件-------------------
A(end,end) = 1;
A(end,end-1) = 0;
B(end,end) = 0;
B(end,end-1) = 0;
b(end) = 0.0013; % 0.0013为M50NiL钢初始碳浓度
% -------------------右边界条件-------------------

% -------------------第三类边界条件-------------------


% ------------碳浓度C矩阵求解-------------
for n = 1:numT-1
    C(n+1,:) = (A \ (B * C(n,:)' + b))';
end
% ------------碳浓度C矩阵求解-------------


end
