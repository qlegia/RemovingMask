function  [rec_alm, l2err] = reconstruct_image(org_alm,orgLmax,msk_alm,maskLmax,fac,opt_sigs_re, opt_sigs_im)

% reconstructed alm

rec_alm = zeros(size(org_alm));



% directory of all matrices E
Edir = '../mat_files/';

LLmax = length(opt_sigs_re);
err = zeros(1,LLmax);
% loading angular power spectrum 
% C_0 = 1
LLmax = 100;
halfL = 50;
C_ells = ones(1,LLmax);  % C_ells(1) = C_1, C_ells(2) = C_2, etc.
C_ells(halfL+1:LLmax) = -2*[halfL+1:LLmax]/(LLmax+1)+2;


for mm = 0:LLmax
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

  re_vec_a = real(vec_a(1:Lp1));
  im_vec_a = imag(vec_a(1:Lp1));
 
  sig_re = opt_sigs_re(mm+1);
  sig_im = opt_sigs_im(mm+1);
  
  matGa = eye(Jp1);
  invCU = inv(C+matUp);

  M = E'*matGa*E + sig_re*eye(Lp1); % (Lp1,Jp1)*(Jp1,Jp1)*(Jp1,Lp1)
  re_vecw = M \ (E'*matGa*re_vec_av);
  re_hata = C*invCU*re_vecw;

  M = E'*matGa*E + sig_im*eye(Lp1); % (Lp1,Jp1)*(Jp1,Jp1)*(Jp1,Lp1)
  im_vecw = M \ (E'*matGa*im_vec_av);
  im_hata = C*invCU*im_vecw;

   
  tmp          = re_hata - re_vec_a;
  re_err(mm+1) = sum(tmp.*tmp);     
  tmp          = im_hata - im_vec_a;
  im_err(mm+1) = sum(tmp.*tmp);

  % reconstructed alm (complex)
  rec_alm(i2:i2+orgLmax-mm) = re_hata+i*im_hata;
end
l2err = re_err(1) + 2*sum(re_err(2:end))+ im_err(1) + 2*sum(im_err(2:end));
l2err = sqrt(l2err);
rec_alm = rec_alm(:); % save as a column vector
org_alm = org_alm(:);
