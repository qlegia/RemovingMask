% solving least square problem using QR factorization
% need to run ../py_files/recons_mapQR.py to get the output figures
clear all;

tic
sigma0 = 0.0;

inst = 1;
dirname = '../mat_files/';
fname = sprintf('%sLinear_Nside2048_instance%d.mat',dirname,inst);

eval(['load ',fname]);
org_alm = alm;

% reconstructed alm
rec_alm = zeros(size(org_alm));

%loading the masked noisy mask
% noise level

fac = 1e-4;
pow = -log10(fac);  % fac = 10^{-power}

fname_mask = sprintf('%sMasked_noisy__1e_%d_Nside2048_instance%d.mat',dirname,pow,inst);
eval(['load ',fname_mask]);
msk_alm = alm;

maskLmax = 1200; %%%% cross-checked with ../py_files/rand_masked_2.py
orgLmax = 100;


% directory of all matrices E
Edir = '../mat_files/';

Lmax = 100;
err = zeros(1,Lmax+1);
rel_err = zeros(1,Lmax+1);
residual = zeros(1,Lmax+1);

% angular power spectrum of noises
halfL = 50;
C_ells = ones(1,Lmax);  % C_ells(1) = C_1, C_ells(2) = C_2, etc.
C_ells(halfL+1:Lmax) = -2*[halfL+1:Lmax]/(Lmax+1)+2;

for mm = 0:Lmax
  fname = sprintf('%sE_L1max100_L2max900_m%d.mat',Edir,mm);
  eval(['load ',fname]);

  [Jp1,Lp1] = size(E);
  J = Jp1+mm-1;
  Lmax = Lp1+mm-1;

  %construct the matrix C = Omega in the paper
  C = zeros(Lp1,Lp1);
  if (mm >= 1)
    C = diag([C_ells(mm:Lmax)]);
  else % m = 0
    C = diag([1 C_ells(1:Lmax)]);
  end
  % matrix Upsilon is just a scaled of Omega=C
  matUp = fac*C;

  % 
  ell = mm;
  i1 = getidx2(maskLmax,ell,mm);
  vec_av = msk_alm(i1:i1+maskLmax-mm);
  vec_av = vec_av(:);
  re_vec_av = real(vec_av(1:Jp1));  % re_vec_av(1) = a_{0,0}
  im_vec_av = imag(vec_av(1:Jp1));

  ell = mm;
  i2 = getidx2(orgLmax,ell,mm);
  vec_a = org_alm(i2:i2+orgLmax-mm);
  vec_a = vec_a(:);

  re_vec_a = real(vec_a);
  im_vec_a = imag(vec_a);
 
  % solving the least square problem E* vec_w = vec_av 
  % let E = QR
  [Q,R] = qr(E);
  Q1 = Q(:,1:Lmax+1);
  R1 = R(1:Lmax+1,:);
  %
  invCU = inv(C+matUp);
  re_alpha = R1\(Q1'*re_vec_av);
  %re_hata = C*invCU*re_alpha;
  re_hata = re_alpha/(1+fac);  % C*(C+matUp)^{-1} = Id/(1+fac);

  im_alpha = R1\(Q1'*im_vec_av);
  %im_hata = C*invCU*im_alpha;
  im_hata = im_alpha/(1+fac);

  re_err(mm+1)      = norm(re_hata - re_vec_a)/ sqrt( length(re_vec_av));
  re_rel_err(mm+1)  = norm(re_hata - re_vec_a)/ norm(re_vec_a);
  re_residual(mm+1) = norm(E*re_hata-re_vec_av)/sqrt( length(re_vec_av));  
  
  im_err(mm+1)      = norm(im_hata - im_vec_a)/ sqrt( length(im_vec_av));
  im_rel_err(mm+1)  = norm(im_hata - im_vec_a)/ norm(im_vec_a);
  im_residual(mm+1) = norm(E*im_hata-im_vec_av)/sqrt( length(im_vec_av));  

  % reconstructed alm (complex)
  rec_alm(i2:i2+orgLmax-mm) = re_hata+i*im_hata;
end

%plot([0:Lmax],re_err,'b',[0:Lmax],im_err,'r')
%grid on
%xlabel('$m$','interpreter','latex')
%ylabel('$\ell_2$ errors','interpreter','latex')
rec_alm = rec_alm(:); % save as a column vector
kap = 0
org_alm = org_alm(:); % column vector;
save L100_QR_noise org_alm rec_alm maskLmax orgLmax pow kap sigma0 
!cp L100_QR_noise.mat ../py_files

toc
%  fac         time (1st run)       time (2nd run)
% 0            6.242695 seconds (tic; reconstruct_v3_no_noise; toc)
% 10^{-2}      23.937409 seconds.   6.054489 seconds
% 10^{-3}      10.464860 seconds.   6.060336 seconds
% 10^{-4}      7.420839 seconds     5.980529 seconds
