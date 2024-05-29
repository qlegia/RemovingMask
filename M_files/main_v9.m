% after running this script, you need to run ../py_files/reconstruct_map.py
% outputs ../py_files/reconstr_field_v9_gauss_noise_1e-2.png
%         ../py_files/reconstr_field_v9_gauss_noise_1e-3.png 
%         ../py_files/reconstr_field_v9_gauss_noise_1e-4.png 
%         ../py_files/error_field_v9_gauss_noise_1e-2.png
%         ../py_files/error_field_v9_gauss_noise_1e-3.png
%         ../py_files/error_field_v9_gauss_noise_1e-4.png

clear all;

% load the original random field
inst = 1;
dirname = '../linearRF/';
fname = sprintf('%sLinear_Nside2048_instance%d.mat',dirname,inst);

eval(['load ',fname]);
org_alm = alm;

% loading angular power spectrum 
% C_0 = 1
LLmax = 100;
halfL = 50;
%C_ells = ones(1,LLmax);  % C_ells(1) = C_1, C_ells(2) = C_2, etc.
%C_ells(halfL+1:LLmax) = -2*[halfL+1:LLmax]/(LLmax+1)+2;

KK = 900;
fac = 1e-1;
pow = -log10(fac);  % fac = 10^{-power}

Jmax = LLmax+KK; % L=100; K=900
halfJ = Jmax/2; 
Sig_j = fac*ones(1,Jmax); % Sig_j(1) = \Sigma_1
% with the following line we get v6_ figures,
% if it is commented out we get v6a_figures
Sig_j(halfJ+1:Jmax) = fac*(-2*[halfJ+1:Jmax]/(Jmax+1) + 2); 



%loading masked field with axial mask
%fname_mask = sprintf('%sLinear_masked_noisy_Nside2048_instance%d.mat',dirname,inst);
fname_mask = sprintf('%sMasked_noisy__1e_%d_Nside2048_instance%d.mat',dirname,pow,inst);
eval(['load ',fname_mask]);
msk_alm = alm;

Edir = '../matE_L1max100_L2max900/';

orgLmax = 100;
maskLmax = 1200; % double checked with ../py_files/rand_masked_2.py

% single sigma_k = { 10^{-15}+k*10^{-16}: k=5,10,15,20,25}
strategy = 6
switch strategy 
   case 1
      % strategy 1: run main_v8.m first
      load sig_opt_v8/L100_all_m_re.mat;
      load sig_opt_v8/L100_all_m_im.mat;
      [rec_alm, l2err] = reconstruct_image(org_alm,orgLmax,msk_alm,maskLmax,fac,opt_sigs_re, opt_sigs_im);
      l2err
   case 2
      % strategy 2: run main_v8_RSW.m first
      load sig_opt_v8rsw/L100_all_m_re.mat;
      load sig_opt_v8rsw/L100_all_m_im.mat;
      [rec_alm, l2err] = reconstruct_image(org_alm,orgLmax,msk_alm,maskLmax,fac,opt_sigs_re, opt_sigs_im);
      l2err
   case 3
      % strategy 3:
      ss2 = 1e+3; % for m>=51
      ss1_choices = [1e-16; 1e-15; 2*1e-15; 3*1e-15; 1e-14;  1e-13; 1e-12; 1; 1e+2; 1e+3];
      l2errs = zeros(1,length(ss1_choices));
      for kk=1:length(ss1_choices)
           ss1 = ss1_choices(kk);  
           opt_sigs_re = ss1*ones(LLmax+1);
           opt_sigs_re(halfL+2:end) = ss2;
           opt_sigs_im = opt_sigs_re;
           [rec_alm, l2err] = reconstruct_image(org_alm,orgLmax,msk_alm,maskLmax,fac,opt_sigs_re, opt_sigs_im);
           l2errs(kk) = l2err;
      end
   case 4
      % strategy 4:
      ss1 = 1e-15;  % when fac = 1e-4 or 1e-2 for m<=50
      ss2_choices = [1e-15; 1; 10; 1e+2; 1e+3; 1e+4; 1e+5]; % when fac = 1e-4
      %ss1 = 2*1e-15; % when fac = 1e-3
      %ss2_choices = [1e-15; 2*1e-15; 1; 10; 1e+2; 1e+3; 1e+4; 1e+5]; % when fac = 1e-3      
      l2errs = zeros(1,length(ss2_choices));
      for kk=1:length(ss2_choices)
           ss2 = ss2_choices(kk);
           opt_sigs_re = ss1*ones(LLmax+1);
           opt_sigs_re(halfL+2:end) = ss2;
           opt_sigs_im = opt_sigs_re;
           [rec_alm, l2err] = reconstruct_image(org_alm,orgLmax,msk_alm,maskLmax,fac,opt_sigs_re, opt_sigs_im);
           l2errs(kk) = l2err;
      end
   case 5
      % strategy 5:
      opt_sigs_re = 1e-15*ones(LLmax+1); % when fac=1e-4 or 1e-2
      %opt_sigs_re = 2*1e-15*ones(LLmax+1); % when fac=1e-3
      opt_sigs_im = opt_sigs_re;
      [rec_alm, l2err] = reconstruct_image(org_alm,orgLmax,msk_alm,maskLmax,fac,opt_sigs_re, opt_sigs_im);
      l2err
   case 6   
      % strategy 6:
      % single sigma_k = { 10^{-15}+k*10^{-16}: k=0,5,10,15,20,25}
      sig_choices = 1e-15 + [0:5:25]*1e-16;
      l2errs = zeros(1, length(sig_choices));
      for kk=1:length(sig_choices)
	  ss = sig_choices(kk);     
	  opt_sigs_re = ss*ones(LLmax+1);
          opt_sigs_im = opt_sigs_re;
          [rec_alm, l2err] = reconstruct_image(org_alm,orgLmax,msk_alm,maskLmax,fac,opt_sigs_re, opt_sigs_im);
          l2errs(kk) = l2err;
      end
      [l2errmin, kmin] = min(l2errs);      
end

% for strategy 6:
ss = sig_choices(kmin);
fprintf('optimal sigma= %e \n',ss);
opt_sigs_re = ss*ones(LLmax+1);
opt_sigs_im = opt_sigs_re;
[rec_alm, l2err] = reconstruct_image(org_alm,orgLmax,msk_alm,maskLmax,fac,opt_sigs_re, opt_sigs_im);
org_alm = org_alm(:); % column vector
save L100_re_const6 org_alm rec_alm maskLmax orgLmax pow 
!cp L100_re_const6.mat ../py_files

%----------------------------------------------------------------------------
% strategy 6, same sigma values for all m
%---------------------------------------------------------------------------
% fac = 1e-1
% sig = 1e-14*[0.1000    0.1500    0.2000    0.2500    0.3000    0.3500]
% l2err        24.9265   24.5892   24.5152   24.5489   24.5618   24.6012
%--------------------------------------------------------------------------
% fac = 1e-2
% sig = 1e-14*[0.1000    0.1500    0.2000    0.2500    0.3000    0.3500]
% l2err =      10.4358   10.3731   10.5226   10.7048   10.8367   11.0220
%-------------------------------------------------------------------------
% fac = 1e-3
% sig = 1e-14*[0.1000    0.1500    0.2000    0.2500    0.3000    0.3500]
% l2err =      9.6354    9.3535    9.4427    9.6428    9.7896    9.9695
%-------------------------------------------------------------------------
% fac = 1e-4
% sig = 1e-14*[0.1000    0.1500    0.2000    0.2500    0.3000    0.3500]
% l2err =      8.7794    8.7616    8.9469    9.2348    9.4465    9.6596
% =========================================================================
% results using different choices of opt_sigs, 
% fac = 1e-4
% sigmas from sig_opt_v8/     l2err = 5.3903
% sigmas from sig_opt_v8rsw/  l2err = 11.0504 
%-------------------------------------------------------------------------------------------------------------
% ss2 = 1e+3 for m>=51
% for m<=50, choose sigma = ss1
% ss1    =  1e-16     1e-15   2*1e-15  3*1e-15  1e-14     1e-13     1e-12     1         1e+2      1e+3
% l2err  =  65.5924   8.7794  8.9469   9.4465   11.0517   13.2054   15.2088   33.6583   49.9063   50.2289
%--------------------------------------------------------------------------------------------------------------
% ss1 = 1e-15  for m<=50,
% for m>=51, choose sigma=ss2
% ss2   =   1e-15;         1;      10;     1e+2;     1e+3;     1e+4;    1e+5; 
% l2err =   8.7794    8.7794    8.7794    8.7794    8.7794    8.7794    8.7794
%==============================================================================================================
% fac = 1e-3
% sigmas from sig_opt_v8/     l2err = 6.4714
% sigmas from sig_opt_v8rsw/  l2err = 11.2202
%------------------------------------------------------------------------------------------------------------
% ss2 = 1e+3
% ss1    =  1e-16     1e-15     2*1e-15   3*1e-15  1e-14     1e-13     1e-12     1         1e+2      1e+3
% l2err  = 74.2333    9.6354    9.4427    9.7896   11.2249   13.2965   15.2766   33.6785   49.9067   50.2290
%------------------------------------------------------------------------------------------------------------
% ss1 = 2*1e-15
% ss2 =   1e-15;  2*1e-15;   1;        10;       1e+2;     1e+3;     1e+4;     1e+5
% l2err = 9.4427   9.4427    9.4427    9.4427    9.4427    9.4427    9.4427    9.4427
%============================================================================================================
% fac = 1e-2
% sigmas from sig_opt_v8/     l2err = 8.0212 
% sigmas from sig_opt_v8rsw/  l2err = 12.1486
% -----------------------------------------------------------------------------------------------------------
% ss2 = 1e+3
% ss1 =    1e-16     1e-15    2*1e-15   3*1e-15   1e-14     1e-13     1e-12     1         1e+2      1e+3
% l2err = 89.9172   10.4358   10.5226   10.8367   12.1624   14.0915   15.9701   33.8767   49.9102   50.2293
%------------------------------------------------------------------------------------------------------------
% ss1 = 1e-15
% ss2 =   1e-15;   2*1e-15;   1;        10;       1e+2;     1e+3;     1e+4;     1e+5
% l2err = 10.4358   10.4358   10.4358   10.4358   10.4358   10.4358   10.4358
