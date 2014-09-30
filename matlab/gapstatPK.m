function [gaps,SDs,best_k] = gapstatPK(data,minK,maxK,paramstruct) 
% An implementation of the gap statistic algorithm from Tibshirani, Walther, 
% and Hastie's "Estimating the number of clusters in a data set via the 
% gap statistic".
% Given a matrix `data`, where rows are observations and columns are individual
% dimensions, compute and plot the gap statistic (according to a uniform 
% reference distribution).
%
% Patrick Kimes, 2012.
%
% Inputs:
%   data    - d x n matrix of data, each column vector is 
%                    a "d-dim data vector"
%
%   minK    - minimum number of clusters to test
%             (recommended to use minK = 1)
%
%   maxK    - maximum number of clusters to test
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
%    fields           values
%
%    refdist          Method for estimating reference null distribution 
%                     1 - uniform distribution along input directions
%                         (default)
%                     2 - uniform distribution along PC directions
%
%    sim_num          Number of samples to simulate from reference
%                     (default = 100)
%
%    iprint           0  don't print errorbar plot to screen
%                     1  (default)  print error bar plot

if minK < 1 ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp('!!!   Warning from gapstatPK:          !!!') ;
    disp('!!!   minK must be integer value >=1,  !!!') ;
    disp('!!!   using minK = 1                   !!!') ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;    
    minK = 1 ;
end ;

if maxK < minK ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp('!!!   Warning from gapstatPK:            !!!') ;
    disp('!!!   maxK must be integer value >minK,  !!!') ;
    disp('!!!   using maxK = minK+5                !!!') ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;    
    maxK = minK+5 ;
end ;

% d = size(data,1) ;
n = size(data,2) ;
refdist = 1 ;
sim_num = 100 ;
iprint = 1 ;

if nargin > 1 ;   %  then paramstruct is an argument
  if isfield(paramstruct,'refdist') ;    %  then change to input value
    refdist = paramstruct.refdist ;
  end ;
  if isfield(paramstruct,'sim_num') ;    %  then change to input value
    sim_num = paramstruct.sim_num ;
  end ;
  if isfield(paramstruct,'iprint') ;    %  then change to input value
    iprint = paramstruct.iprint ;
  end ;
end ;
    

CIs = zeros(maxK-minK+1,1) ;

totd = sum(sum((data - mean(data,2)*ones(1,n)).^2)) ;

for k = minK:maxK ;
    if k == 1 ;
        CIs(k-minK+1) = totd ;
    else
        paramstruct = struct('nrep',10) ;
        labels = SigClustKmeanRepPK(data,k,paramstruct) ;
        CIs(k-minK+1) = ClustIndPK(data,k,logical(labels))*totd ;
    end ;
end ;    

CIs = log(CIs) ;

sim_CIs = zeros(sim_num,maxK-minK+1) ;

if refdist == 1 ;
    maxs = max(data,[],2) ;
    mins = min(data,[],2) ;
    Maxs = repmat(maxs,1,n) ;
    Mins = repmat(mins,1,n) ;
    for m = 1:sim_num ;
        sim_data = unifrnd(Mins,Maxs) ;
        sim_totd = sum(sum((sim_data - mean(sim_data,2)*ones(1,n)).^2)) ;
        for k = minK:maxK ;
            if k == 1 ;
                sim_CIs(m,k-minK+1) = sim_totd ;
            else
                sim_labels = SigClustKmeanRepPK(sim_data,k,paramstruct) ;
                sim_CIs(m,k-minK+1) = ClustIndPK(sim_data,k,logical(sim_labels))*totd ;
            end ;
        end ; 
    end ;
elseif refdist == 2 ;
    %%% fill in! %%%
else 
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp('!!!   Warning from gapstatPK:          !!!') ;
    disp('!!!   refdist must be 1 or 2           !!!') ;
    disp('!!!   aborting function                !!!') ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;  
    return ;
end ;

sim_CIs = log(sim_CIs) ;

ref_means = mean(sim_CIs,1) ;
SDs = std(sim_CIs,1) .* sqrt(1+1/sim_num) ;

gaps = ref_means - CIs' ;
candidates = gaps(1:end-1) > gaps(2:end)-SDs(2:end) ;
best_k = find(candidates,1,'first')+minK-1 ;

if isempty(best_k) ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp('!!!   Warning from gapstatPK:                !!!') ;
    disp('!!!   no K satisfies threshold condition     !!!') ;
    disp('!!!   selecting best K by alternate method   !!!') ;
    disp('!!!   please check code for details          !!!') ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;  
    [~,best_k] = max((gaps(1:end-1)-gaps(2:end))./SDs(2:end)) ;
    best_k = best_k+minK-1 ;
end ;
    
if iprint == 1;
    clf ;
    hold on ;
    errorbar(minK:maxK,gaps,SDs) ;
    title('Gap Statistic (Tibshirani et al., 2001)') ;
    xlabel('K') ;
    ylabel('gap') ;
    set(gca,'xtick',minK:maxK) ;
    set(gca,'xlim',[minK-.5 maxK+.5]) ;
end ;

