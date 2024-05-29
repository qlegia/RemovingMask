% Plot and calculate spherical harmonic coefficients for mask
% around equator with transition from adeg to bdeg smoothness kappa

% Smmothness parameter
kappa = 3;

% Transition interval, symmetric about z = 0 (equator)
adeg = 10;
bdeg = 20;

%adeg = 5;
%bdeg = 10;
a = cos(pi/2 - pi*adeg/180);
b = cos(pi/2 - pi*bdeg/180);

fprintf('Mask: a = %.4f (%.2f degrees) to b = %.4f (%.2f degrees)', a, adeg, b, bdeg);
fprintf(', smoothness kappa = %d\n', kappa);

Lmax = 1200;
v = zeros(1,Lmax+1);
% compute the coefficients v_{ell,0} only using Gausian quadrature on [-1,1]
% and Legendre polynomials 
[Xg, Wg] = grule(3000);
% set f being the mask function given by Rob 
f = mask(Xg,a,b,kappa,0);
% set f being constant
%f = ones(size(Xg));
%f = 0.5*f+0.5; % new mask
for ell=0:Lmax
  ell	
  P_ell = legendre(ell,Xg);
  %v(ell+1) = sum(P_ell(1,:).*f.*Wg) * (2*ell+1)/2;
  v(ell+1) = sqrt(pi)*sum(P_ell(1,:).*f.*Wg) * sqrt(2*ell+1);
end

vl = v;
save v_ell0_1200_gauss vl
% P_0 = 1,            Y_{0,0} = 1/sqrt{4pi}
% P_1 = x,            Y_{1,0} = 1/sqrt{4pi} sqrt{3} cos(theta)
% P_2 = 0.5(3x^2-1)   Y_{2,0} = 1/sqrt{4pi} sqrt{5} 0.5*(3cos^2(theta)-1) 
%
% Y_{ell,0} = 1/sqrt{4pi} sqrt(2ell+1) P_{ell}(cos\theta)
%
% int_0^{2pi} int_{0}^\pi f(\cos\theta) Y_{ell,0}(\cos\theta) sin\theta d\theta d\phi  
% = 2pi /sqrt{4pi} int_{0}^pi f(\cos\theta) sqrt{2ell+1} P_{ell} (cos(theta)) \sin\theta d\theta 
%   let z = cos(theta)
%       dz = -sin(theta) dtheta 
% = sqrt{pi} int_{-1}^1 f(z) sqrt{2ell+1} P_{\ell}(z) dz 
