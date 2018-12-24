% ---------------------------用户需给出以下参数的值---------------------------
leftX = 0; % leftX必须为0,因为left=0对应方程的左边界条件
rightX = 5e-5; % rightX的值需满足右边界条件,建议rightX取较大的值
numX = 101; % 划分的格点数没有特殊要求,因为Crank-Nicolson方法是稳定的,numX的值合适就好

pulses = [30,30]; % 每个强渗脉冲或扩散脉冲持续的时间,此数组应以强渗脉冲开始

D = 1e-12; % 扩散系数
beta = 4e-9; % 传递系数
Cp = 0.01; % 气氛碳势
% ---------------------------用户需给出以上参数的值---------------------------


% -------------------------强渗脉冲或扩散脉冲--------------------------------
len = length(pulses);
pulse = struct('initialT', 0,'endT', 0, 'numT', 1000);
pulse = repmat(pulse, 1, len);
for m = 1:len
    pulse(m).endT = pulses(m);
end
% -------------------------强渗脉冲或扩散脉冲--------------------------------


x = linspace(leftX, rightX, numX);

C0 = icfun(x);


% -----------------------------碳浓度C矩阵----------------------------------
for n = 1:len
    if pulse(n).endT ~= 0
        if mod(n,2) == 0
            C = solver(leftX, rightX,numX, pulse(n).initialT, pulse(n).endT, pulse(n).numT, D, 0, Cp, C0);
        else
            C = solver(leftX, rightX, numX, pulse(n).initialT, pulse(n).endT, pulse(n).numT, D, beta, Cp, C0);
        end
        C0 = C(end, :);
    end
end
% -----------------------------碳浓度C矩阵----------------------------------


sol = exactSol(Cp, 0.0013, D, beta, x, 30);


plot(x, C(end, :), '*', x, sol);
