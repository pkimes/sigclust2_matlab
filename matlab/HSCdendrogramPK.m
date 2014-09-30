function [] = HSCdendrogramPK( data, pvals, hmethod, hmetric, paramstruct )
%
%
% HSCdendrogramPK, Hierarchical SigClust Dendrogram plot
%   Patrick Kimes' matlab function for plotting annotated dendrogram
%     based on HSigClustPK p-value output.
%
%   Started:        08/03/2013
%   Last updated:   08/03/2013
%
%

% Inputs:
%
%   data          - d x n matrix
%   
%   pvals         - (n-1) vector of p-values output from HSigClustPK
%
%   hmethod       - linkage method for hierarchical clustering
%                     e.g. 'average', 'ward'
%
%   hmetric       - distance metric for hierarchical clustering
%                     e.g. 'sqeuclidean', 'euclidean'
%                     NOTE: use 'euclidean' for 'ward'
%                            Matlab automatically treats this 
%                            as 'sqeuclidean' for ward case
%
%   paramstruct   - a Matlab structure of input parameters
%
%     fields           values
%     
%     nShow            will show p-values for highest nShow joins,
%                        (default = 10)
%     pShow            will show p-values for joins with p-val < pShow
%                        in addition to top nShow, (default = 0.10)
%     labels           labels for the data if known (default = [])
%     CIused           name of the cluster index used to print in title,
%                        (default = 'unspecified')
%     dendsavestr      name to save dendrogram,
%                        (default = ['HSCdendrogram_' date])
%
%
% Outputs:
%
%   dendplot      - dendrogram colored with p-values shown
%   

%    Copyright (c) Patrick Kimes 2012,2013




%  First set all parameters to defaults
%
nShow = 10 ;
pShow = 0.1 ;
CIused = 'unspecified' ;
labels = [];
dendsavestr = ['HSCdendrogram_' date] ;

%  Now update parameters as specified,
%  by parameter structure (if it is used)
%
if nargin > 4 ;   %  then paramstruct is an argument

  if isfield(paramstruct,'nShow') ;    %  then change to input value
    nShow = paramstruct.nShow ; 
  end ;

  if isfield(paramstruct,'pShow') ;    %  then change to input value
    pShow = paramstruct.pShow ; 
  end ;

  if isfield(paramstruct,'CIused') ;    %  then change to input value
    CIused = paramstruct.CIused ; 
  end ;

  if isfield(paramstruct,'labels') ;    %  then change to input value
    labels = paramstruct.labels; 
  end ;
  if isfield(paramstruct,'dendsavestr') ;    %  then change to input value
    dendsavestr = paramstruct.dendsavestr ; 
  end ;

end ;

n = size(data,2) ;
if length(pvals) ~= n-1 ;
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  disp('!!  pvals must be vector of length n-1 !!') ;
  disp('!!  returning without any output       !!') ;
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  return ;
end ;

if length(labels) ~= n;
    labels = repmat({''}, 1, n);
end;

%  explicitly define the squared euclidean distance measure
%
if strcmp(hmetric, 'sqeuclidean') ;
  fmetric = {@(Xi,Xj) sum((ones(size(Xj,1),1)*Xi - Xj).^2,2)} ;
else
  fmetric = hmetric ;
end ;

Z = linkage(data', hmethod, fmetric) ;

% avgCI = {'2-means CI', ...
%          'linkage value',  ...
%          'modified linkage 1', ...
%          'modified linkage 2'} ;


% create annotated dendrogram plot
%
figure(1) ;
clf ;
H = dendrogram(Z, 0, 'labels', labels);
tempc = get(H,'color') ;

for i = 1:(n-1) ;
  if pvals(i) > 1 ;
    tempc(i) = {[0 0 1]} ;
  elseif pvals(i) < 0.05 ;
    tempc(i) = {[1 0 0]} ;
  else
    tempc(i) = {[0 0 0]} ;
  end ;
end ;

xdata = get(H,'XData') ;

for i = 1:(n-1) ;
  if (i>=(n-nShow) && pvals(i) < 2) || pvals(i)<pShow ;
    xcoord = mean(xdata{i}) ;
    ycoord = Z(i,3) ;
    if ycoord < Z(n-1,3)/40 ;
      ycoord = -Z(n-1,3)/40 ;
    end ;
    text(xcoord, ycoord, num2str(pvals(i), 3), ...
        'verticalalignment','top', ...
        'horizontalalignment','center', ...
        'color','blue') ;
  end ;
end ;

set(H,{'color'},tempc) ;
set(H,'LineWidth',2) ;
title(['HSigClust results, ' ...
       hmethod ' linkage, ' ...
       hmetric ' similarity, ' ...
       'CI = ' CIused ]) ;
ax = axis ;
text(ax(1) + .90*(ax(2)-ax(1)), ...
     ax(3) + .95*(ax(4)-ax(3)), ...
     ['p-val shown for last ' num2str(nShow) ' joins '], ...
     'fontsize',10) ;
text(ax(1) + .90*(ax(2)-ax(1)), ...
     ax(3) + .92*(ax(4)-ax(3)), ...
     ['and others w/ p-val < ' num2str(pShow)], ...
     'fontsize',10) ;
   
orient('landscape') ;
print('-dpdf', dendsavestr) ; 




end

