function f = mask(x, a, b, kappa, ipr)
% f = mask(x, a, b, kappa, ipr)
% Mask with smoothness kappa >= 0 and for b > a > 0
% f = 1 for |x| > b
% f = 0 for |x} < a
% f = p for a < |x| < b where p polynomial
% ipr = print parameter (default ipr = 0), plot if ipr > 1

if nargin < 4
    ipr = 0;
end

f = ones(size(x));
I0 = find(abs(x)<=a);
f(I0) = 0;
I = find(x > a & x < b);
x01 = (x(I)-a)/(b-a);
f01 = mask01(x01, kappa);
f(I) = f01;
I = find(x < -a & x > -b);
x01 = -(a+x(I))/(b-a);
f01 = mask01(x01, kappa);
f(I) = f01;

if ipr > 1
    
    fprintf('Mask: a = %.4f, b = %.4f, kappa = %d\n', a, b, kappa)
    figure(ipr);
    plot(x, f);
    yl = [-0.1 1.1];
    hold on
    plot(a*[1, 1], yl, 'm--', b*[1 1], yl, 'r-.');
    plot(-a*[1, 1], yl, 'm--', -b*[1 1], yl, 'r-.');
    hold off
    grid on
    tstr = sprintf('Mask: smoothness \\kappa = %d', kappa);
    title(tstr);
    
end
