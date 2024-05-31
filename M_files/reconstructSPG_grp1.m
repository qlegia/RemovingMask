% 13-March-2024
% Run the SPGL1 algorithm on the masked noisy field on the sphere
% output: L100_SPG_grp.mat
% need to run ../py_files/recons_mapSPG_grp.py to get the errors of the reconstruction
clear all;
sigma0 = 0.0;
addpath('/srv/hdrive/z9701564/old_p/CMB_Pseudo_Cell/OtherMethods/spgl1-2.1');
% load the Fourier coefficients of the original random field
inst = 1;
dirname = '../mat_files/';
fname = sprintf('%sLinear_Nside2048_instance%d.mat',dirname,inst);

eval(['load ',fname]);
org_alm = alm;


% reconstructed alm
rec_alm = zeros(size(org_alm));

% noise level
fac = 1e-4;
pow = -log10(fac);  % fac = 10^{-power}

% load the Fourier coefficients of the noise field
fname_noise = sprintf('%sNoise_1e_%d_Nside2048_instance%d.mat',dirname,pow,inst);
eval(['load ', fname_noise]);
noise_alm = alm;

% loading the Fourier coefficients of the masked noisy field
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
 
  %sigma = fac*sqrt(length(vec_av));
  % compute the tolerance
  ell = mm;
  i3 = getidx2(orgLmax, ell, mm);
  vec_eps = noise_alm(i3: i3+orgLmax-mm);
  vec_eps = vec_eps(:);
  Eeps = E*vec_eps;
  %
  % uncomment the following 2 lines for no-noise 
  %Eeps = zeros(size(Eeps)); % no-noise case
  %vec_eps = zeros(size(vec_eps)); % no-noise case

  vec_bv = vec_av(1:Jp1) + Eeps;

  re_sigma = norm(real(vec_eps));
  im_sigma = norm(real(vec_eps));


  % solving the basis pursuit denoise problem 
  % minimize ||hat_a||_2 subject to ||E*hat_a - vec_bv||_2 <= sigma
  %
  %
  opts = spgSetParms('verbosity',0);
  grp = ones(1,length(re_vec_a));
  % set up the weights
  weights = ones(1,length(re_vec_a)); 
  % let y = W*x
  % min || y ||_2 subject to || AW^{-1} y - b||_2 <= sigma
  invW = diag(1./weights);
  w_re_hata = spg_group(E*invW, real(vec_bv), grp, re_sigma, opts);
  w_im_hata = spg_group(E*invW, imag(vec_bv), grp, im_sigma, opts);
  % x = W^{-1} y
  re_hata = invW*w_re_hata;
  im_hata = invW*w_im_hata;

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
save L100_SPG_grp org_alm rec_alm maskLmax orgLmax pow kap sigma0 
!cp L100_SPG_grp.mat ../py_files

% fac    time (average after 10 runs)
% 0      2.09       
% 1e-4   2.85
% 1e-3   2.29 
% 1e-2   2.08 
