function C = exactSol(Cp, C0, D, beta, x, t)
C=(Cp-C0)*(erfc(x/sqrt(4*D*t))-exp((beta*x+beta^2*t)/D).*erfc(x/sqrt(4*D*t)+beta*sqrt(t/D)))+C0; 
end
