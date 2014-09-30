function ClustInd = ClustIndPK(mdata,K,ClassFlags) 
% CLUSTINDPK, CLUSTerINDex
%     Index of clustering, that underlies K-means clustering
%     using:  Squared Euclidean distance
%   Patrick Kimes' modification of
%   Steve Marron's matlab function
% Inputs:
%     mdata      - d x n matrix of multivariate data
%                       (each col is a data vector)
%                       d = dimension of each data vector
%                       n = number of data vectors 
%
%     K       - number of clusters
%
%     ClassFlags - n x K logical matrix indicating elements classes
%                       There can be a '1' in at most one column for each
%                       row. OK to have an entire row or zeros.
%                           (will just leave those out of calculation)
%                       Otherwise gives a warning message,
%                       and returns an empty value of Cluster Index
% Output:
%     ClustInd   - K-means Cluster Index
%                       This is sum of within Class Sums of Squares about mean,
%                       divided by Total Sum of Squares about overall mean 
%
% Assumes path can find personal function:
%    vec2matSM.m

%    Copyright (c) J. S. Marron 2007
%    Copyright (c) P. K. Kimes 2012


d = size(mdata,1) ;
n = size(mdata,2) ;
    

%  Check Inputs
%
if (rem(K,1) ~= 0) || (K < 2) ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp('!!!   Warning from ClustIndPK:           !!!') ;
    disp('!!!   K must be integer value >=2,       !!!') ;
    disp('!!!   using K = 2                        !!!') ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;    
    K = 2 ;
end ;

if size(ClassFlags,2) ~= K ;
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  disp('!!!   Error from ClustIndPK.m:               !!!') ; 
  disp('!!!   ClassFlags dimension does not match    !!!') ;
  disp(['!!!   K = ' num2str(K) '                                 !!!']) ;
  disp('!!!   Terminating Execution                  !!!') ; 
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  ClustInd = [] ;
  return ;
end ;

if size(ClassFlags,1) ~= n ;
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  disp('!!!   Error from ClustIndPK.m:               !!!') ; 
  disp(['!!!   ClassFlags needs to have length ' num2str(n)]) ;
  disp('!!!   Terminating Execution                  !!!') ; 
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  ClustInd = [] ;
  return ;
end ;

if ~islogical(ClassFlags) ;
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  disp('!!!   Error from ClustIndPK.m:                  !!!') ; 
  disp('!!!   ClassFlags needs to be a logical matrix   !!!') ;
  disp('!!!   Terminating Execution                     !!!') ; 
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  ClustInd = [] ;
  return ;
end ;

if sum(sum(ClassFlags,2)>1) > 0 ;
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  disp('!!!   Error from ClustIndPK.m:                 !!!') ; 
  disp('!!!   data point common to multiple clusters   !!!') ;
  disp('!!!   Terminating Execution                    !!!') ; 
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  ClustInd = [] ;
  return ;
end ;


%  Compute Cluster Index
%
if d == 1 ;    %  have only d = 1

  mdataAll = mdata(:,logical(sum(ClassFlags,2))) ;
      %  need to do this, since some data points may not be in either class
  totd = sum(sum((mdataAll - mean(mdataAll,2)).^2)) ;
      %  Total sums of square distance from mean of column vectors

  classd = zeros(1,K) ;
      
  for k = 1:K ;
      classkdata = mdata(:,ClassFlags(:,k)) ;
      classd(k) = sum(sum((classkdata - mean(classkdata,2)).^2)) ;
  end ;

else     %  d > 1 

  mdataAll = mdata(:,logical(sum(ClassFlags,2))) ;
      %  need to do this, since some data points may not be in either class
  totd = sum(sum((mdataAll - vec2matSM(mean(mdataAll,2),size(mdataAll,2))).^2)) ;
      %  Total sums of square distance from mean of column vectors

  classd = zeros(1,K) ;

  for k = 1:K ;
      classkdata = mdata(:,ClassFlags(:,k)) ;
      classd(k) = sum(sum((classkdata - vec2matSM(mean(classkdata,2),size(classkdata,2))).^2)) ;
  end ;

end ;

ClustInd = sum(classd) / totd ;



