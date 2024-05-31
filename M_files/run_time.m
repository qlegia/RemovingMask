% measure running time
clear all;
nn = 10
for k_exp =1:nn
  t0 = cputime;
  %reconstruct_v3_no_noise;
  %reconstructQRnoise_v2;
  reconstructSPG_grp1;
  t1 = cputime;
  runtime(k_exp)= t1-t0;
end  
mean(runtime)
