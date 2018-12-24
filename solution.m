parameter = [10,20,10,20,10,30];

len = length(parameter);

pulse = struct('endT',0);
%pulse = repmat(pulse, 1, len);



for m = 1:len-1
    pulse(m).endT = parameter(m);
end










leftX = 0;
rightX = 5e-5;
numX = 51;
initialT = 0;
endT = 30;
numT = 1000;
D = 1e-12;
beta = 4e-9;
Cp = 0.01;
x = linspace(leftX, rightX, numX);
C0 = icfun(x);

C = solver(leftX, rightX, numX, initialT, endT, numT, D, beta, Cp, C0);
sol = exactSol(Cp, 0.0013, D, beta, x, endT*2);

C0 = C(end, :);
C = solver(leftX, rightX, numX, initialT, endT, numT*2, D, beta, Cp, C0);

plot(x, C(end, :), '*', x, sol);
