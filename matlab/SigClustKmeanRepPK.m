function [bestclass, vindex, midx] = SigClustKmeanRepPK(data,K,paramstruct) 
% SIGCLUST2MEANREPSM, statistical SIGnificance of CLUSTers,
%         computes 2-MEAN clustering REPetitions, 
%         over random restarts, for given data 
%   Steve Marron's matlab function
%     Provides fundamental tool for studying the
%     consistency of the standard 2-means clustering
%     method for splitting a given data set.
%     Does a large number of repetitions, using random 
%     starts, and outputs numerical results.
%
% Inputs:
%   data    - d x n matrix of data, each column vector is 
%                    a "d-dim data vector"
%
%   K       - number of clusters
%
%   paramstruct - a Matlab structure of input parameters
%                    Use: "help struct" and "help datatypes" to
%                         learn about these.
%                    Create one, using commands of the form:
%
%       paramstruct = struct('field1',values1,...
%                            'field2',values2,...
%                            'field3',values3) ;
%
%                          where any of the following can be used,
%                          these are optional, misspecified values
%                          revert to defaults
%
%    fields            values
%
%    nrep             Number of repetitions (random restarts).
%                     When not specified, default is 100
%
%    randstate        State of uniform random number generator
%                     When empty, or not specified, just use current seed  
%
%    randnstate       State of normal random number generator
%                     When empty, or not specified, just use current seed  
%
%    iscreenwrite     0  (default)  no screen writes
%                     1  write to screen to show progress
%
% Outputs:
%    bestclass        Best Cluster Labelling  
%                     (in sense of minimum Cluster Index, over repetitions)
%                     n x K logical matrix indicating elements classes
%
%    vindex           Vector of Cluster Index values
%                     computed over random restarts
%                     nrep x 1 vector
%
%    midx             matrix of cluster labels,
%                     computed over random restarts,
%                     nrep x n matrix,
%                     each row is vector of labelings, i.e. 1's and 2's  
%
%
% Assumes path can find personal functions:
%    vec2matSM.m

%    Copyright (c) J. S. Marron 2007
%    Copyright (c) P. K. Kimes 2012



%  First set all parameters to defaults
%
nrep = 100 ;
randstate = [] ;
randnstate = [] ;
iscreenwrite = 0 ;


%  Now update parameters as specified,
%  by parameter structure (if it is used)
%
if nargin > 1 ;   %  then paramstruct is an argument

  if isfield(paramstruct,'nrep') ;    %  then change to input value
    nrep = getfield(paramstruct,'nrep') ; 
  end ;

  if isfield(paramstruct,'randstate') ;    %  then change to input value
    randstate = getfield(paramstruct,'randstate') ; 
  end ;

  if isfield(paramstruct,'randnstate') ;    %  then change to input value
    randnstate = getfield(paramstruct,'randnstate') ; 
  end ;

  if isfield(paramstruct,'iscreenwrite') ;    %  then change to input value
    iscreenwrite = getfield(paramstruct,'iscreenwrite') ; 
  end ;

end ;    %  of resetting of input parameters


%  check inputs
%
if (rem(K,1) ~= 0) || (K < 2) ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp('!!!   Warning from SigClustKmeanRepPK:   !!!') ;
    disp('!!!   K must be integer value >=2,       !!!') ;
    disp('!!!   using K = 2                        !!!') ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;    
    K = 2 ;
end ;


%  set preliminary stuff
%
d = size(data,1) ;
         %  dimension of each data curve
n = size(data,2) ;
         %  number of data curves



%  Run nrep K-means clusterings, with random restarts
%
if ~isempty(randstate) ;
  rand('state',randstate) ;
end ;
if ~isempty(randnstate) ;
  randn('state',randnstate) ;
end ;

totd = sum(sum((data - vec2matSM(mean(data,2),size(data,2))).^2)) ;
    %  Total sum of square distance from mean of column vectors

midx = [] ;
vindex = [] ;
for irep = 1:nrep
  if iscreenwrite ~= 0 ;
    disp(['        Working on repetition ' num2str(irep) ' of ' num2str(nrep)]) ;
  end ;
  [idx,c,sumd] = kmeans(data',K,'EmptyAction','singleton') ;
  midx = [midx; idx'] ;
  vindex = [vindex; sum(sumd)/totd] ;
end ;

if iscreenwrite ~= 0 ;
  disp(['Finished K means Clustering']) ;
  disp(' ') ;
end ;



[~,imin] = min(vindex) ;

bestclass = zeros(n,K) ;
idx = sub2ind(size(bestclass), 1:n, midx(imin,:)) ;
bestclass(idx) = 1 ;
