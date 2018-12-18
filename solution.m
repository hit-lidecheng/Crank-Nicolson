C = solver(0,0.005,51,0,7200,120,@icfun);
xlswrite('test.xlsx',C);