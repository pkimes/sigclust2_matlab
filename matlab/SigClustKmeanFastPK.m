function [BestClass, bestCI] = SigClustKmeanFastPK(data,K,paramstruct) 
% SIGCLUST2MEANFASTSM, statistical SIGnificance of CLUSTers,
%         computes 2-MEAN clustering using a FAST algorithm
%   Steve Marron's matlab function
%     Provides a fast (relative to many random restarts)
%     algorithm for the K-means clustering
%     method for splitting a given data set.
%     The main goal is to find the global minimizer, 
%     without requiring the usual large number of 
%     random restarts.  Idea is to use PCA to find 
%     suggested starting classifications for the usual
%     algorithm, and to work in appropriate orthogonal 
%     subspaces for later steps.
%     Results are not guaranteed, but are expected to be
%     of high quality for many High Dimension Low Sample
%     Size data sets, because of the theoretically
%     predicted structure of such data 
%     (e.g. relative orthogonality of clusters)
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
%    maxstep          MAXimum number of STEPs to consider.
%                     When not specified, default is 10
%
%    ioutplot         0  make no output plots
%                     1 (default) make output plots at end of each step
%                           in new coorinate system
%                     2  make output plots at end of each step, 
%                           in both old and new coordinate systems
%
%    icolor           0  fully black and white version (everywhere)
%                     Kx3 color matrix:  a color label for each class
%                         [[1 0 0]; [0 0 1]] is default for K = 2
%                         [[1 0 0]; [0 1 0]; [0 0 1]] is default for K = 3
%                         ... up to K = 7.
%                         for K > 7, equal splits between 
%                         [0 0 0] and [0 0 1] are used
%                             e.g. [[0 0 0]; [0 0 .1]; ... [0 0 1]]
%                                  for K = 11
%
%    markerstr        Character array (K x 1), of symbols,
%                         e.g. 'o', '.', '+', 'x','*'
%                         created using:  strvcat
%                             e.g. strvcat('o','+') (default for K = 2)
%                         for K > 13, 'o' is used for all classes
%
%    maxlim        Matrix of axis limits
%                        Use [] for default of all automatically chosen, 
%                            by axisSM
%                        Use 1 for symmetrically chosen, by axisSM
%                            (often preferred for centered plots, as in PCA)
%                        Otherwise, must be (size(mdir,2) + npcadir) x 2 
%                            matrix of axis limits, 
%                            with each row corresponding to a direction
%                        Note:  Use is generally not recommended,
%                        because defaults give "good visual impression
%                        of decomposition".  It is mostly intended to allow
%                        the highlighting of "visually different scales" 
%                        in data.  But this comes at the cost of reduced  
%                         detail being visible in the plots
%
%    iplotaxes        0 do not plot axes
%                     1 (default) plot axes as determined by direction vectors, 
%                           using thin line type
%
%    iplotdirvec      0 (default) do not plot direction vectors
%                     1 plot direction vectors, using thick line type
%
%    ibelowdiag       0 leave off scatterplots below diagonal
%                     1 (default) show scatterplots both above and below diagonal
%
%    titlestr         string to add to subplot titles
%                     always get the left top title indicating Step #
%                     the second top title indicating view
%                     this controls the third top title
%                     default is an empty array, [] for no additional title
%
%    titlefontsize    font size for title
%                           (only has effect when plot is made here,
%                            and when the titlecellstr is nonempty)
%                     default is empty, [], for Matlab default
%
%    labelcellstr    Vertical cell array of strings for axis labels
%                        create this using cellstr, 
%                        or {{string1; string2; ...}}
%                            Note:  These strange double brackets seems to be 
%                                needed for correct pass to subroutine
%                                It may change in later versions of Matlab
%                    Default is an empty cell array, {}, which then gives:
%                         For istep Starting Cluster plot:
%                            "Mean Diff Direction" 1st
%                            "Ortho PC _" for others
%                         For istep 2-Means Result, data dimension d = <= 4 plot:
%                            "Direction _" for each
%                         For istep 2-Means Result, data dimension d = > 4 plot:
%                            "PC _" for each
%                    For no labels, use {''; ''; ''; ''}
%                    Length of cell array must be:
%                        1 for data dimension d = 1
%                        2 for data dimension d = 2
%                        3 for data dimension d = 3
%                        4 for data dimension d = >= 4
%
%    labelfontsize    font size for axis labels
%                                    (only has effect when plot is made here,
%                                     and when a label str is nonempty)
%                           default is empty [], for Matlab default
%
%    savestr          string controlling saving of output,
%                         either a full path, or a file prefix to
%                         save in matlab's current directory
%                         Will add .ps, and save as either
%                             color postscript (for plot with Entry 3 = 1)
%                         or
%                             black&white postscript (other plots)
%                         unspecified:  results only appear on screen
%
%    iscreenwrite     0  (default)  no screen writes
%                     1  write to screen to show progress
%
% Outputs:
%    BestClass        Best Cluster Labelling  
%                     (in sense of minimum Cluster Index, over repetitions)
%                     n x K logical matrix indicating elements classes
%
%    bestCI           Best Cluster Index value
%                     scalar
%
%
% Assumes path can find personal functions:
%    scatplotSM.m
%    bwsjpiSM.m
%    kdeSM.m
%    lbinrSM.m
%    vec2matSM.m
%    pcaSM.m
%    projplot1SM.m
%    projplot2SM.m
%    bwrfphSM.m
%    bwosSM.m
%    rootfSM
%    bwrotSM.m
%    bwsnrSM.m
%    iqrSM.m
%    cquantSM.m
%    axisSM.m
%
%    ClustIndPK.m

%    Copyright (c) J. S. Marron 2007
%    Copyright (c) P. K. Kimes 2012



%  First set all parameters to defaults
%
maxstep = 10 ;
ioutplot = 1 ;
if (rem(K,1) ~= 0) || (K < 2) ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp('!!!   Warning from SigClustKmeanFastPK:  !!!') ;
    disp('!!!   K must be integer value >=2,       !!!') ;
    disp('!!!   using K = 2                        !!!') ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;    
    K = 2 ;
    icolor = [[1 0 0]; [0 0 1]] ;
elseif K <= 7 ;
    icolor = [[1 0 0];
              [0 0 1];
              [0 1 0];
              [1 1 0];
              [1 0 1];
              [0 1 1];
              [1 1 1]] ;
    icolor = icolor(1:K,:) ;
elseif K > 7 ;
    icolor = zeros(K,3) ;
    for k = 1:K ;
        icolor(k,3) = (k-1)/(K-1) ;
    end ;
else
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp('!!!   Warning from SigClustKmeanFastPK:  !!!') ;
    disp('!!!   Invalid value for K,               !!!') ;
    disp('!!!   using K = 2                        !!!') ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;  
    K = 2 ;
    icolor = [[1 0 0]; [0 0 1]] ;    
end ;
if K <= 13 ;
    markerstr = char('o','+','*','.','x','s','d','^','v','>','<','p','h''') ;
    markerstr = markerstr(1:K,:) ;
elseif K > 13 ;
    markerstr = 'o' ;
end ;
maxlim = [] ;
iplotaxes = 1 ;
iplotdirvec = 0 ;
ibelowdiag = 1 ;
titlestr = [] ;
titlefontsize = [] ;
labelcellstr = {} ;
labelfontsize = [] ;
savestr = [] ;
iscreenwrite = 0 ;
d = size(data,1) ;
n = size(data,2) ;

%  Now update parameters as specified,
%  by parameter structure (if it is used)
%
if nargin > 1 ;   %  then paramstruct is an argument

  if isfield(paramstruct,'maxstep') ;    %  then change to input value
    maxstep = paramstruct.maxstep ;
  end ;

  if isfield(paramstruct,'ioutplot') ;    %  then change to input value
    ioutplot = paramstruct.ioutplot ;
  end ;

  if isfield(paramstruct,'icolor') ;    %  then change to input value
    tempicolor = paramstruct.icolor ;
    if size(tempicolor,1) == K ;
        icolor = tempicolor ;
    else 
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      disp('!!!   Warning from SigClustKmeanFastPK:  !!!') ;
      disp('!!!   Invalid icolor   ,                 !!!') ;
      disp('!!!   using default colors               !!!') ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;        
    end ;
  end ;

  if isfield(paramstruct,'markerstr') ;    %  then change to input value
    tempmarkerstr = paramstruct.markerstr ;
    if size(tempmarkerstr,1) == K ;
        markerstr = tempmarkerstr ;
    else 
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      disp('!!!   Warning from SigClustKmeanFastPK:  !!!') ;
      disp('!!!   Invalid markerstr,                 !!!') ;
      disp('!!!   using default markers              !!!') ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;        
    end ;
  end ;

  if isfield(paramstruct,'maxlim') ;    %  then change to input value
    maxlim = paramstruct.maxlim ;
  end ;

  if isfield(paramstruct,'iplotaxes') ;    %  then change to input value
    iplotaxes = paramstruct.iplotaxes ;
  end ;

  if isfield(paramstruct,'iplotdirvec') ;    %  then change to input value
    iplotdirvec = paramstruct.iplotdirvec ;
  end ;

  if isfield(paramstruct,'ibelowdiag') ;    %  then change to input value
    ibelowdiag = paramstruct.ibelowdiag ;
  end ;

  if isfield(paramstruct,'titlestr') ;    %  then change to input value
    titlestr = paramstruct.titlestr ;
  end ;

  if isfield(paramstruct,'titlefontsize') ;    %  then change to input value
    titlefontsize = paramstruct.titlefontsize ;
  end ;

  if isfield(paramstruct,'labelcellstr') ;    %  then change to input value
    labelcellstr = paramstruct.labelcellstr ;
  end ;

  if isfield(paramstruct,'labelfontsize') ;    %  then change to input value
    labelfontsize = paramstruct.labelfontsize ;
  end ;

  if isfield(paramstruct,'savestr') ;    %  then use input value
    savestr = paramstruct.savestr ;
    if ~(ischar(savestr) || isempty(savestr)) ;    %  then invalid input, so give warning
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      disp('!!!   Warning from SigClust2meanFastSM:  !!!') ;
      disp('!!!   Invalid savestr,                   !!!') ;
      disp('!!!   using default of no save           !!!') ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      savestr = [] ;
    end ;
  end ;

  if isfield(paramstruct,'iscreenwrite') ;    %  then change to input value
    iscreenwrite = paramstruct.iscreenwrite ;
  end ;

end ;    %  of resetting of input parameters



if ioutplot ~= 0 ;  %  Will make plots, so close current windows
  close all ;
end ;

rankdata = rank(data) ;


%  Start with first K-1 PC labelling
%
if size(data,1) > 1 ;    %  d > 1, so actually do PCA

  viout = [1 0 0 0 1] ;    %  Output mpc  -  npc x n matrix of 
                           %     "principal components", 
                           %     i.e. of coefficients of projections 
                           %     of data onto eigenvectors, 
                           %     i.e. of "scores".
  if size(data,1) < (K-1) ;
      paramstruct = struct('npc',size(data,1),...
                           'viout',viout,...
                           'iscreenwrite',0) ;
  else 
      paramstruct = struct('npc',K-1,...
                           'viout',viout,...
                           'iscreenwrite',0) ;
  end ;      
  outstruct = pcaSM(data,paramstruct) ;
  veigval = outstruct.veigval ;
  vpc = outstruct.mpc ;
  %  scale PC directions by corresponding eigenvalues
  %
  vpc = vpc' * diag(veigval) ;
else 
  vpc = data ;
end ;
%  apply K-means clustering in K-1 dim PC space 
%
nidx = kmeans(vpc,K,'EmptyAction','singleton') ;
curClass = zeros(n,K) ;
idx = sub2ind(size(curClass), 1:n, nidx') ;
curClass(idx) = 1 ;
curCI = ClustIndPK(data,K,logical(curClass)) ;
if iscreenwrite == 1 ;
  disp(['    ' num2str(K) ' Means Fast, ' num2str(size(vpc,2)) '-PC Cluster Index = ' num2str(curCI)]) ;
end ;

if ioutplot == 2 ;  %  Make output plot using 
                    %    original data PCA coordinates

  figure ;
  if isempty(titlestr) ;
    titlecellstr = {{['2-Means Fast, Step ' num2str(1)] 'Starting Cluster'}} ;
  else 
    titlecellstr = {{['2-Means Fast, Step ' num2str(1)] 'Starting Cluster' titlestr}} ;
  end ;
  if ~isempty(savestr) ;
    savestrout = [savestr 'Step' num2str(1) 'StartClust'] ;
  else 
    savestrout = [] ;
  end ;
  paramstruct = struct('iMDdir',0,...
                       'icolor',icolor,...
                       'isubpopkde',1,...
                       'markerstr',markerstr,...
                       'legendcellstr',{{['K = ' num2str(K)] ['CI = ' num2str(curCI)]}},...
                       'mlegendcolor',zeros(2,3),...
                       'maxlim',maxlim,...
                       'iplotaxes',iplotaxes,...
                       'iplotdirvec',iplotdirvec,...
                       'ibelowdiag',ibelowdiag,...
                       'titlecellstr',titlecellstr,...
                       'titlefontsize',titlefontsize,...
                       'labelfontsize',labelfontsize,...
                       'savestr',savestrout,...
                       'iscreenwrite',0) ;

  %  Update to properly handle empty components
  %
  if ~(isempty(labelcellstr)) ;
    paramstruct.labelcellstr = labelcellstr{1} ;
  end ;

  SigClustLabelPlotPK(data,nidx',paramstruct) ; 

end ;



%  Compute 1st K-means clustering
%
CurMMeanStart = (data*curClass) ./ (ones(d,1)*sum(curClass)) ;
[nidx,mc] = kmeans(data',K,'EmptyAction','singleton','Start',CurMMeanStart') ;
curClass = zeros(n,K) ;
idx = sub2ind(size(curClass), 1:n, nidx') ;
curClass(idx) = 1 ;
curCI = ClustIndPK(data,K,logical(curClass)) ;
CurMMeanStart = mc' ;
if iscreenwrite == 1 ;
  disp(['    ' num2str(K) ' Means Fast, Step 1 Cluster Index = ' num2str(curCI)]) ;
end ;
MDdir = orth(CurMMeanStart - mean(CurMMeanStart,2)*ones(1,K)) ;
mdirexc = MDdir ;
    %  matrix of directions to exclude

% MDdir = CurMMeanStart(:,1) - CurMMeanStart(:,2) ;
% MDdir = MDdir / norm(MDdir) ;
% mdirexc = MDdir ;

if ioutplot ~= 0 ;  %  Make output plot using 
                    %    new PCA coordinates


  figure ;
  if isempty(titlestr) ;
    titlecellstr = {{['K-Means Fast, Step ' num2str(1)] [num2str(K) '-means Result']}} ;
  else 
    titlecellstr = {{['K-Means Fast, Step ' num2str(1)] [num2str(K) '-means Result'] titlestr}} ;
  end ;
  if ~isempty(savestr) ;
    savestrout = [savestr 'Step' num2str(1) 'Res2mean'] ;
  else 
    savestrout = [] ;
  end ;
  paramstruct = struct('iMDdir',0,...
                       'icolor',icolor,...
                       'isubpopkde',1,...
                       'markerstr',markerstr,...
                       'legendcellstr',{{['K = ' num2str(K)] ['CI = ' num2str(curCI)]}},...
                       'mlegendcolor',zeros(2,3),...
                       'maxlim',maxlim,...
                       'iplotaxes',iplotaxes,...
                       'iplotdirvec',iplotdirvec,...
                       'ibelowdiag',ibelowdiag,...
                       'titlecellstr',titlecellstr,...
                       'titlefontsize',titlefontsize,...
                       'labelfontsize',labelfontsize,...
                       'savestr',savestrout,...
                       'iscreenwrite',0) ;

  %  Update to properly handle empty components
  %
  if ~(isempty(labelcellstr)) ;
    paramstruct.labelcellstr = labelcellstr{1} ;
  end ;

  SigClustLabelPlotPK(data,nidx',paramstruct) ; 

end ;


bestCI = curCI ;
BestClass = logical(curClass) ;



for istep = 2:maxstep ;

  if size(mdirexc,2) < rankdata ;    %  then can work with non-excluded data

    %  Project data on non-excluded directions
    %
    mprojexc = mdirexc * pinv(mdirexc' * mdirexc) * mdirexc' ;
        %  projection matrix onto excluded space
    projdata = data - mprojexc * data ;


    %  PC1 labelling on nonexcluded data
    %
    viout = [1 1 0 0 1] ;    %  Output mpc  -  npc x n matrix of 
                           %     "principal components", 
                           %     i.e. of coefficients of projections 
                           %     of data onto eigenvectors, 
                           %     i.e. of "scores".
    if rankdata < (size(mdirexc,2)+K-1) ;
        paramstruct = struct('npc',(rankdata-size(mdirexc,2)),...
                             'viout',viout,...
                             'iscreenwrite',0) ;
    else 
        paramstruct = struct('npc',K-1,...
                             'viout',viout,...
                             'iscreenwrite',0) ;
    end ;      
    outstruct = pcaSM(projdata,paramstruct) ;
    veigval = outstruct.veigval ;
    vpc = outstruct.mpc ;
    vcurPCdir = outstruct.meigvec ;
    %  scale PC directions by corresponding eigenvalues
    %
    vpc = vpc' * diag(veigval) ;
    %  apply K-means clustering in K-1 dim PC space 
    %
    nidx = kmeans(vpc,K,'EmptyAction','singleton') ;
    curClass = zeros(n,K) ;
    idx = sub2ind(size(curClass), 1:n, nidx') ;
    curClass(idx) = 1 ;
    curCI = ClustIndPK(data,K,logical(curClass)) ;
    if curCI < bestCI ;
      bestCI = curCI ;
      BestClass = logical(curClass) ;
    end ;
    if iscreenwrite == 1 ;
      disp(['    ' num2str(K) ' Means Fast, Step ' num2str(istep) ' Ortho PC Cluster Index = ' num2str(curCI)]) ;
    end ;

    if ioutplot == 2 ;  %  Make output plot using 
                        %    original data PCA coordinates

      figure ;
      if isempty(titlestr) ;
        titlecellstr = {{[num2str(K) '-Means Fast, Step ' num2str(istep)] 'Starting Cluster'}} ;
      else 
        titlecellstr = {{[num2str(K) '-Means Fast, Step ' num2str(istep)] 'Starting Cluster' titlestr}} ;
      end ;
      if ~isempty(savestr) ;
        savestrout = [savestr 'Step' num2str(istep) 'StartClust'] ;
      else 
        savestrout = [] ;
      end ;
      paramstruct = struct('iMDdir',0,...
                           'icolor',icolor,...
                           'isubpopkde',1,...
                           'markerstr',markerstr,...
                           'legendcellstr',{{['K = ' num2str(K)] ['CI = ' num2str(curCI)]}},...
                           'mlegendcolor',zeros(2,3),...
                           'maxlim',maxlim,...
                           'iplotaxes',iplotaxes,...
                           'iplotdirvec',iplotdirvec,...
                           'ibelowdiag',ibelowdiag,...
                           'titlecellstr',titlecellstr,...
                           'titlefontsize',titlefontsize,...
                           'labelfontsize',labelfontsize,...
                           'savestr',savestrout,...
                           'iscreenwrite',0) ;

      %  Update to properly handle empty components
      %
      if ~(isempty(labelcellstr)) ;
        paramstruct.labelcellstr = labelcellstr{1} ;
      end ;

      SigClustLabelPlotPK(data,nidx',paramstruct) ; 

    end ;



    %  Compute next 2 means clustering
    %
    CurMMeanStart = (data*curClass) ./ (ones(d,1)*sum(curClass)) ;
    [nidx,mc] = kmeans(data',K,'EmptyAction','singleton','Start',CurMMeanStart') ;
    curClass = zeros(n,K) ;
    idx = sub2ind(size(curClass), 1:n, nidx') ;
    curClass(idx) = 1 ;
    curCI = ClustIndPK(data,K,logical(curClass)) ;
    CurMMeanStart = mc' ;
    if curCI < bestCI ;
      bestCI = curCI ;
      BestClass = logical(curClass) ;
    end ;
    if iscreenwrite == 1 ;
      disp(['    ' num2str(K) ' Means Fast, Step ' num2str(istep) ' Cluster Index = ' num2str(curCI)]) ;
    end ;
    MDdir = orth(CurMMeanStart - mean(CurMMeanStart,2)*ones(1,K)) ;
    tempresid = sum(abs(mprojexc * MDdir - MDdir)) ;
    for k = 1:size(MDdir,2) ;
        if tempresid(k) > 0 ;
                                %  this direction already in excluded subspace
          mdirexc = [mdirexc vcurPCdir(:,k)] ;
              %  matrix of directions to exclude
        else    %  have component outside current subspace
          mdirexc = [mdirexc MDdir(:,k)] ;
              %  matrix of directions to exclude
        end ;
    end ;
    
    if ioutplot ~= 0 ;  %  Make output plot using 
                        %    new PCA coordinates

      figure ;
      if isempty(titlestr) ;
        titlecellstr = {{['K-Means Fast, Step ' num2str(istep)] [num2str(K) '-means Result']}} ;
      else 
        titlecellstr = {{['K-Means Fast, Step ' num2str(istep)] [num2str(K) '-means Result'] titlestr}} ;
      end ;
      if ~isempty(savestr) ;
        savestrout = [savestr 'Step' num2str(istep) 'Res2mean'] ;
      else 
        savestrout = [] ;
      end ;
      paramstruct = struct('iMDdir',1,...
                           'icolor',icolor,...
                           'isubpopkde',1,...
                           'markerstr',markerstr,...
                           'legendcellstr',{{['K = ' num2str(K)] ['CI = ' num2str(curCI)]}},...
                           'mlegendcolor',zeros(2,3),...
                           'maxlim',maxlim,...
                           'iplotaxes',iplotaxes,...
                           'iplotdirvec',iplotdirvec,...
                           'ibelowdiag',ibelowdiag,...
                           'titlecellstr',titlecellstr,...
                           'titlefontsize',titlefontsize,...
                           'labelfontsize',labelfontsize,...
                           'savestr',savestrout,...
                           'iscreenwrite',0) ;

      %  Update to properly handle empty components
      %
      if ~(isempty(labelcellstr)) ;
        paramstruct.labelcellstr = labelcellstr{1} ;
      end ;

      SigClustLabelPlotPK(data,nidx',paramstruct) ; 

    end ;


  end ;


end ;



