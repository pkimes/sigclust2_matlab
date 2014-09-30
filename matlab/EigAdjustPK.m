function [veig,c,lambdahat] = EigAdjustPK(data) 
%
%
% EigAdjustPK, EIGenvalue ADJUSTment procedure for d >> n case
%   Patrick Kimes' matlab function for implementing the eigenvalue
%   adjustment procedure described in:
%     'Convergence and Prediction of Principal Component Scores 
%      in High-Dimensional Settings' by S. Lee et al. AoS 2010
%
%   Started:        10/11/2012
%   Last updated:   10/11/2012
%

% Inputs:
%   data        - p x n matrix
%

%    Copyright (c) Patrick Kimes 2012

[p,n] = size(data) ;
gamma = p/n ;

paramstruct = struct('iscreenwrite',0, ...
                     'viout',1) ;
outstruct = pcaSM(data,paramstruct) ;
dstar = outstruct.veigval ;
nev = length(dstar) ;
dstar = [dstar; zeros(p-nev,1)] ;

rv = dstar / sum(dstar) ;

ssum = p ;

dhat = ssum*rv ;

epsilon = 1 ;
niter = 0 ;
while epsilon > 1e-6 ;
    
  lambdahat = rhoI(dhat,gamma) ;
  
  epsilon = abs(ssum - sum(lambdahat)) ;
  ssum = sum(lambdahat) ;
  
  dhat = ssum*rv ;  
  
  niter = niter+1 ;
  
end ;


c = sum(dstar)/ssum ;  
veig = rhoI(dhat,gamma)*c;



function [lambda] = rhoI(d,gamma)
%
%
% rhoI, RHO Inverse transformation
%   transformation to obtain consistent estimate of lambda in spike covariance
%   model with p/n -> gamma>0 asymptotics
%     'Convergence and Prediction of Principal Component Scores 
%      in High-Dimensional Settings' by S. Lee et al. AoS 2010
%

p = length(d) ;

lambda = zeros(p,1) ;
for i = 1:p ;
  di = d(i) ;
  b = gamma-1-di ;
  if di > (1+sqrt(gamma))^2 ;
    lambda(i) = (-b+sqrt(b^2-4*di)) / 2 ;
  else
    lambda(i) = 1 ;
  end ;
end ;





