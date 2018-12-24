% ---------------------------�û���������²�����ֵ---------------------------
leftX = 0; % leftX����Ϊ0,��Ϊleft=0��Ӧ���̵���߽�����
rightX = 5e-5; % rightX��ֵ�������ұ߽�����,����rightXȡ�ϴ��ֵ
numX = 101; % ���ֵĸ����û������Ҫ��,��ΪCrank-Nicolson�������ȶ���,numX��ֵ���ʾͺ�

pulses = [30,30]; % ÿ��ǿ���������ɢ���������ʱ��,������Ӧ��ǿ�����忪ʼ

D = 1e-12; % ��ɢϵ��
beta = 4e-9; % ����ϵ��
Cp = 0.01; % ����̼��
% ---------------------------�û���������ϲ�����ֵ---------------------------


% -------------------------ǿ���������ɢ����--------------------------------
len = length(pulses);
pulse = struct('initialT', 0,'endT', 0, 'numT', 1000);
pulse = repmat(pulse, 1, len);
for m = 1:len
    pulse(m).endT = pulses(m);
end
% -------------------------ǿ���������ɢ����--------------------------------


x = linspace(leftX, rightX, numX);

C0 = icfun(x);


% -----------------------------̼Ũ��C����----------------------------------
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
% -----------------------------̼Ũ��C����----------------------------------


sol = exactSol(Cp, 0.0013, D, beta, x, 30);


plot(x, C(end, :), '*', x, sol);
