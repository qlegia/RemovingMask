% solving least square problem using QR factorization
% need to run ../py_files/recons_map2A.py to get the output figures
% reconstructed_Gaussian_field_no_noise_QR.png
% error_field_no_noise_QR.png
clear all;
sigma0 = 0.0;

inst = 1;
dirname = '../linearRF/';
fname = sprintf('%sLinear_Nside2048_instance%d.mat',dirname,inst);

eval(['load ',fname]);
org_alm = alm;

% reconstructed alm
rec_alm = zeros(size(org_alm));

%loading masked field with axial mask
fname_mask = sprintf('%sLinear_masked_Nside2048_instance%d.mat',dirname,inst);
eval(['load ',fname_mask]);
msk_alm = alm;

maskLmax = 1200; %%%% cross-checked with ../py_files/rand_masked_2.py
orgLmax = 100;

% loading previously computed optimal sigma values for each m

% directory of all matrices E
Edir = '../matE_L1max100_L2max900/';

Lmax = 100;
err = zeros(1,Lmax+1);
rel_err = zeros(1,Lmax+1);
residual = zeros(1,Lmax+1);

for mm = 0:Lmax
  fname = sprintf('%sE_L1max100_L2max900_m%d.mat',Edir,mm);
  eval(['load ',fname]);

  [Jp1,Lp1] = size(E);
  J = Jp1+mm-1;
  Lmax = Lp1+mm-1;

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
  re_hata = R1\(Q1'*re_vec_av);

  im_hata = R1\(Q1'*im_vec_av);

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
save L100_re_const3 org_alm rec_alm maskLmax orgLmax sigma0 kap
!cp L100_re_const3.mat ../py_files
