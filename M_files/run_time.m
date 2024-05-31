% measure running time
clear all;
nn = 5
for k_exp =1:nn
  t0 = cputime;
  reconstructQRnoise_v2;
  t1 = cputime;
  runtime(k_exp)= t1-t0;
end  
mean(runtime)
