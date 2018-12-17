function C = solver(leftX, rightX, numX, initialT, endT, numT, icfun)
% �˷����������̼Ũ�ȷֲ�.
% CΪ̼Ũ�ȷֲ�����,�洢�Ÿ���λ�ø���ʱ���̼Ũ��ֵ.
% ��ʹ�ô˷������̼Ũ�ȷֲ�,���������ֵ:
% ��߽�λ��leftX,�ұ߽�λ��rightX,X������numX;
% ��ʼʱ��initialT,����ʱ��endT,ʱ��������numT;
% ��ʼ��������icfun.

% x����
x = linspace(leftX, rightX, numX);
% t����
t = linspace(initialT, endT, numT);
% x�Ჽ��
stepX = (rightX - leftX) / (numX - 1);
% t�Ჽ��
stepT = (endT - initialT) / (numT - 1);

% ��ɢϵ��
D = 2e-11;

% ��ʼ��C����
C = zeros(numT, numX);

% ��C���󸳳�ʼֵ
C(1,:) = icfun(x);

% δ���ǳ�ʼ�����ͱ߽�����ʱ��ϵ������
r = 1 / 2 * D * stepT / stepX^2;
e = ones(numX,1);
S = spdiags([e -2*e e], -1:1, numX, numX);
A = eye(numX) - r * S;
B = eye(numX) + r * S;

% ------------��һ��߽�������ʼ------------
%b = zeros(numX,1);
% ������߽�ֵ,�˴���߽�ֵΪ0.012
%A(1,1) = 1;
%A(1,2) = 0;
%B(1,1) = 0;
%B(1,2) = 0;
%b(1) = 0.012;
% �����ұ߽�ֵ,�˴��ұ߽�ֵΪ0.001
%A(end,end) = 1;
%A(end,end-1) = 0;
%B(end,end) = 0;
%B(end,end-1) = 0;
%b(end) = 0.001;
% ------------��һ��߽���������------------

% ------------������߽�������ʼ------------
b = zeros(numX,1);
% ������߽�ֵ,�˴���߽����������߽�����
beta = 3e-8;
Cp = 0.1;
A(1,1) = 2 * (stepX)^2 / D / stepT + 2 * beta * stepX / D + 2;
A(1,2) = -2;
B(1,1) = 2 * (stepX)^2 / D / stepT - 2 * beta * stepX / D - 2;
B(1,2) = 2;
b(1) = 4 * beta * stepX * Cp / D;
% �����ұ߽�ֵ,�˴��ұ߽�ֵΪ0.001
A(end,end) = 1;
A(end,end-1) = 0;
B(end,end) = 0;
B(end,end-1) = 0;
b(end) = 0.001;

% �߽�����Ϊ������߽�����ʱ��ɢ���̵Ľ�����
C0 = 0.001;
sol = (Cp - C0) * (erfc(x / sqrt(4 * D * endT)) - exp((beta * x + beta^2 * endT) / D) .* erfc(x / sqrt(4 * D * endT) + beta * sqrt(endT / D))) + C0;
% ------------������߽���������------------

% ------------C������ⲿ�ֿ�ʼ-------------
for n = 1:numT-1
    C(n+1,:) = (A \ (B * C(n,:)' + b))';
end
% ------------C������ⲿ�ֽ���-------------

% ��ͼ����
figure(1)
plot(x,C(end,:),'*',x,sol);
figure(2)
plot(t,C(:,1));


end