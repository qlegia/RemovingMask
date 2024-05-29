function i = getidx2(LMAX,l,m)
% used to access the a_lm coefficients returned from HEALPy map2alm function
% i is the index in the 1-d array alm returned from map2alm
% the coefficients in alm are stored as follows
%  a_00, a_10, a_20,             ...., a_LMAX,0     : LMAX+1 terms
%        a_11, a_21, a_31, .....       a_LMAX,1     : LMAX   terms
%              a_22, a_32, a_42, ..... a_LMAX,2     : LMAX-1 terms
%                     ..........                      ........            
%                               a_l,m  a_LMAX,m     : LMAX-m+1 terms
%                     ......................
%                                      a_LMAX,LMAX
if (m>=0) & (m<=LMAX) & (l>=m)& (l <= LMAX)
  i = l + m*LMAX +1 - (m*(m-1)/2);
else
  error('Invalid l and m')
end

% The number of terms before a_{0,m} is
% LMAX + 1 + LMAX + LMAX-1 + LMAX-2 + .... + LMAX-m+1 
%  = m*LMAX +1 - (1+2+3+4+...+m-1)
%  = m*LMAX +1 - (m*(m-1)/2) 
% So, the ith_ position of a_{l,m} is
%    i = l + m*LMAX +1 - (m*(m-1)/2)
%
% Example 1:  LMAX = 5
%  a_00, a_10, a_20, a_30, a_40, a_50,
%        a_11, a_21, a_31, a_41, a_51,
%              a_22, a_32, a_42, a_52,
%                    a_33, a_42, a_53
%                          a_44, a_54
%                                a_55
%
% m=0, i = l+1                   l=0,..,5
% m=1, i = l+LMAX+1 = l+6        l=1,..,5
% m=2, i = l+2*LMAX + 1 - 1 = l+10, l=2,..,5
% m=3
% m=4
% m=5, i = l+5*LMAX + 1 - 10 = l+16, l=5 
