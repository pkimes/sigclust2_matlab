function [pvalQ, pvalZ, cn, mCIdata] = HSigClustPK(data, paramstruct) 
%
%
% HSigClustPK, Hierarchical statistical SIGnificance of CLUSTers
%   Patrick Kimes' matlab function extending 'SigClustSM.m'.
%       - function relies on MATLAB's 'linkage' function
%         and significance is calculated using 2-means index
%           background is estimated for vdata after subtracting median for
%           each covariate (NOT using PC scores yet) - lines ~430
%
%
% Function also crashes when I use default settings - saving plots?
%   console output:
% % Error using name (line 103)
% % Cannot create output file './QQ.ps'
% % 
% % Error in print (line 209)
% %     pj = name( pj );
% % 
% % Error in qqLM (line 1450)
% %     print('-dpsc', [savestr '.ps']) ;
% % 
% % Error in HSigClustPK (line 730)
% %       qqLM(vdata,paramstruct) ;
%
%
%
%
%   Started:        09/28/2012
%   Last updated:   02/10/2014


% Inputs:
%
%   data          - d x n matrix
%   
%   paramstruct   - a Matlab structure of input parameters
%
%     fields           values
%
%     hmethod          linkage method for hierarchical clustering
%                        currently only supports:
%                                'average' (default)
%                                'ward'
%     hmetric          distance metric for hierarchical clustering
%                        currently only supports:
%                                'sqeuclidean' for 'average'
%                                'euclidean' for 'average' (default)
%                                'euclidean' for 'ward'
%                                   NOTE: Matlab automatically treats this 
%                                         as 'sqeuclidean' for ward case
%     testn            integer value for the minimum number of samples
%                        before SigClust hypothesis testing is applied
%                        (default = 10)
%     sigbackg         Background Standard Deviation, given value.
%                        When this is not given, or set to scalar 0, 
%                        estimate from Data.
%     iCovEst          Covariance Estimation Type
%                        1 - Use Ming Yuan's Soft Thresholded, Constrained MLE
%                           (default, and recommended)
%                        2 - Use Sample Covariance Estimate
%                           (recommended when diagnostics fail)
%                        3 - Use Original Background Noise Thresholded Estimate
%                           (from Liu, et al, JASA paper)
%                           (In Yuan's terminology, this is "Hard Thresholding")
%     ipvaltype        Type of p-value to compute
%                        1 - Empirical Quantile based only
%                        2 - Gaussian Fit Quantile based only
%                        3 - Both    (default) 
%                           (note:  uncomputed p-values are returned as empty,
%                                     and not shown in output plots)
%     nsim             Number of simulated Gaussian samples to use
%                        for main p-value computation
%                        (default = 1000)
%     SimRandstate     State of uniform random number generator,
%                        for Main Simulation
%                        When empty, or not specified, just use current seed  
%                        (has no effect, unless twoMtype = 1)
%     SimRandnstate    State of normal random number generator,
%                        for Main Simulation
%                        When empty, or not specified, just use current seed  
%                        (has no effect, unless twoMtype = 1)
%     datastr          String Descriptive of data set,
%                        mostly for use in graphical output,
%                        when these are requested.
%                        Often may want to end this with "data"
%                        (default is [])
%     iBGSDdiagplot    0  do not make BackGround Standard Deviation
%                           DIAGnostic PLOTs
%                      1  (default) make BackGround Standard Deviation
%                           DIAGnostic PLOTs:
%                             i.   Overlay and density plot of pixel dist'n,
%                                      with robust Gaussian fit
%                             ii.  QQ plot assessing quality of robust fit
%                                      of Gaussian Distribution
%                           This has no effect, unless sigbackg == 0
%     BGSDsavestr      string controlling saving of output for BackGround 
%                        Standard Deviation DIAGnostic PLOTs,
%                        either a full path, or a file prefix to
%                        save in matlab's current directory.
%                        Probably makes sense to end with 
%                        something like:   'AllPixel'
%                        Will add:
%                            'KDE' for plot i.
%                            'QQ' for plot ii.
%                        Will also add .ps (so don't do this again!), 
%                        and save as color postscript files.
%                        Unspecified (or empty):
%                            Results only appear on screen
%                        Only has effect when iBGSDdiagplot = 1
%     iCovEdiagplot    0  do not make COVariance Estimation
%                           DIAGnostic PLOT
%                      1  (default) make COVariance Estimation
%                           DIAGnostic PLOT
%     CovEsavestr      string controlling saving of output for 
%                        COVariance Estimation DIAGnostic PLOTs,
%                        either a full path, or a file prefix to
%                        save in matlab's current directory.
%                        Probably makes sense to end with 
%                        something like:         'EstEigVal'
%                        Will also add .ps (so don't do this again!), 
%                        and save as color postscript files.
%                        Unspecified (or empty):
%                            Results only appear on screen
%                        Only has effect when iCovEdiagplot = 1
%     ipValplot        0  do not make p Value plot
%                      1  (default) make p Value plot
%     pValsavestr      string controlling saving of output for 
%                        p Value plot,
%                        either a full path, or a file prefix to
%                        save in matlab's current directory.
%                        Probably makes sense to end with 
%                        something like:         'pVal'
%                        Will also add .ps (so don't do this again!), 
%                        and save as color postscript files.
%                        Unspecified (or empty):
%                            Results only appear on screen
%                        Only has effect when ipValplot = 1
%     legendcellstr    For p Value Plot:
%                        cell array of strings for legend (nl of them),
%                        useful for (colored) classes, create this using
%                        cellstr, or {{string1 string2 ...}}
%                        Note:  These strange double brackets seems to be needed
%                                for correct pass to subroutine
%                                It may change in later versions of Matlab
%                        CAUTION:  If are updating this field, using a command like:
%                        paramstruct = setfield(paramstruct,'legendcellstr',...
%                        Then should only use single braces in the definition of
%                        legendecellstr, i. e. {string1 string2 ...}
%                        Suggested uses:
%                         For vclass = 0:
%                          legendcellstr = {{'Testing Best 2-means Split'}}
%                         For given vclass:
%                          legendcellstr = {{'Class 1' 'vs.' 'Class 2'}}
%                             then use mlegendcolor appropriately,
%                             e.g.  mlegendcolor = [[1 0 0]; [0 0 0]; [0 0 1]]
%                        Only has effect when ipValplot = 1
%     mlegendcolor     For p Value Plot:
%                        nl x 3 color matrix, corresponding to cell legends above
%                         (not needed when legendcellstr not specified)
%                         (defaults to black when not specified)
%                        Only has effect when ipValplot = 1
%     npValplot        how many p-val plots to print (last 'n' joins)
%                        (default = 10)
%     iscreenwrite     0  (default)  no screen writes
%                              except warning messages
%                      1  write to screen to show major steps
%                      2  to show most steps, but only one write at 
%                              each step of main simulation loop
%                      3  to show all steps, not recommended since it
%                              makes multiple writes for each step of
%                              main simulation loop
%    
%
% Outputs:
%
%   pvalQ         - Simulated SigClust P-value,
%                   Based on Empirical Quantiles, 
%                   Computed by cquantSM.m
%                     average linkage: (n-1)x4
%                     ward linkage:    (n-1)x2
%
%   pvalZ         - Simulated SigClust P-value,
%                   Based on Gaussian Quantiles, 
%                   Computed by normcdf.m
%                     average linkage: (n-1)x4
%                     ward linkage:    (n-1)x2
%
%   Graphics in different Figure windows
%   When ___savestr exists,
%       Color Postscript files saved in '___savestr'.ps
%
% Assumes path can find personal functions:
%    SigClust2meanRepSM.m
%    SigClust2meanFastSM.m
%    SigClustCovEstSM.m
%    SigClustLabelPlotSM.m
%    SigClust2meanPerfSM.m
%    vec2matSM.m
%    axisSM.m
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
%    rootfSM.m
%    bwrotSM.m
%    bwsnrSM.m
%    iqrSM.m
%    cquantSM.m
%    nmfSM.m


%    Copyright (c) J. S. Marron 2007,2008
%                  Patrick Kimes 2012,2013



%  First set all parameters to defaults
%
hmethod = 'average';
hmetric = 'euclidean';
testn = 10;
sigbackg = 0;
iCovEst = 1;
ipvaltype = 3;
nsim = 1000;
SimRandstate = [];
SimRandnstate = [];
datastr = [];
iBGSDdiagplot = 1;
BGSDsavestr = [];
iCovEdiagplot = 1;
CovEsavestr = [];
ipValplot = 1;
pValsavestr = [];
legendcellstr = {};
mlegendcolor = [];
iscreenwrite = 0;
npValplot = 10;

%  Now update parameters as specified,
%  by parameter structure (if it is used)
%
if nargin > 1;   %  then paramstruct is an argument

  if isfield(paramstruct,'hmethod') ;    %  then change to input value
    hmethod = paramstruct.hmethod ; 
    if sum(strcmp(hmethod,{'average','ward'})) ~= 1 ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      disp('!!!   Warning from HSigClustPK.m:   !!!') ;
      disp('!!!   Invalid hmethod, function     !!!') ;
      disp('!!!   currently only supports       !!!') ;
      disp('!!!   "average" or "ward",          !!!') ;
      disp('!!!   using default of "average"    !!!') ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      hmethod = 'average' ;
    end ;
  end ;

  if isfield(paramstruct,'hmetric') ;    %  then change to input value
    hmetric = paramstruct.hmetric ;
    if ~(ischar(hmetric)) ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      disp('!!!   Warning from HSigClustPK.m:   !!!') ;
      disp('!!!   hmetric must be a valid       !!!') ;
      disp('!!!   string, note "sqeuclidean"    !!!') ;
      disp('!!!   can only be combined with     !!!') ;
      disp('!!!   average linkage,              !!!') ;
      disp('!!!   using default of              !!!') ;
      disp('!!!   "average" w/ "sqeuclidean"    !!!') ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      hmethod = 'average' ;
      hmetric = 'sqeuclidean' ;
    end ;
    if strcmp(hmetric,'sqeuclidean') ;
      if ~strcmp(hmethod,'average') ;
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
        disp('!!!   Warning from HSigClustPK.m:    !!!') ;
        disp('!!!   "sqeuclidean" can only be      !!!') ;
        disp('!!!   combined with average linkage, !!!') ;
        disp('!!!   using default of               !!!') ;
        disp('!!!   "average" w/ "sqeuclidean"     !!!') ;
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
        hmethod = 'average' ;
      end ;          
      hmetric = {@(Xi,Xj) sum((ones(size(Xj,1),1)*Xi - Xj).^2,2)} ;
    end ;   
  end ;
    
  if isfield(paramstruct,'testn') ;    %  then change to input value
    testn = paramstruct.testn ;
  end ;

  if isfield(paramstruct,'sigbackg') ;    %  then change to input value
    sigbackg = paramstruct.sigbackg ; 
  end ;

  if isfield(paramstruct,'iCovEst') ;    %  then change to input value
    iCovEst = paramstruct.iCovEst ; 
  end ;

  if isfield(paramstruct,'ipvaltype') ;    %  then change to input value
    ipvaltype = paramstruct.ipvaltype ; 
  end ;

  if isfield(paramstruct,'nsim') ;    %  then change to input value
    nsim = paramstruct.nsim ; 
  end ;

  if isfield(paramstruct,'SimRandstate') ;    %  then change to input value
    SimRandstate = paramstruct.SimRandstate ; 
  end ;

  if isfield(paramstruct,'SimRandnstate') ;    %  then change to input value
    SimRandnstate = paramstruct.SimRandnstate ; 
  end ;

  if isfield(paramstruct,'datastr') ;    %  then change to input value
    datastr = paramstruct.datastr ; 
  end ;

  if isfield(paramstruct,'iBGSDdiagplot') ;    %  then change to input value
    iBGSDdiagplot = paramstruct.iBGSDdiagplot ; 
  end ;

  if isfield(paramstruct,'BGSDsavestr') ;    %  then change to input value
    BGSDsavestr = paramstruct.BGSDsavestr ; 
    if ~(ischar(BGSDsavestr) || isempty(BGSDsavestr)) ;    %  then invalid input, so give warning
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      disp('!!!   Warning from HSigClustPK.m:  !!!') ;
      disp('!!!   Invalid BGSDsavestr,         !!!') ;
      disp('!!!   using default of no save     !!!') ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      BGSDsavestr = [] ;
    end ;
  end ;

  if isfield(paramstruct,'iCovEdiagplot') ;    %  then change to input value
    iCovEdiagplot = paramstruct.iCovEdiagplot ; 
  end ;

  if isfield(paramstruct,'CovEsavestr') ;    %  then change to input value
    CovEsavestr = paramstruct.CovEsavestr ; 
    if ~(ischar(CovEsavestr) || isempty(CovEsavestr)) ;    %  then invalid input, so give warning
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      disp('!!!   Warning from HSigClustPK.m:  !!!') ;
      disp('!!!   Invalid CovEsavestr,         !!!') ;
      disp('!!!   using default of no save     !!!') ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      CovEsavestr = [] ;
    end ;
  end ;

  if isfield(paramstruct,'ipValplot') ;    %  then change to input value
    ipValplot = paramstruct.ipValplot ; 
  end ;

  if isfield(paramstruct,'pValsavestr') ;    %  then change to input value
    pValsavestr = paramstruct.pValsavestr ; 
    if ~(ischar(pValsavestr) || isempty(pValsavestr)) ;    %  then invalid input, so give warning
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      disp('!!!   Warning from HSigClustPK.m:  !!!') ;
      disp('!!!   Invalid pValsavestr,         !!!') ;
      disp('!!!   using default of no save     !!!') ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      pValsavestr = [] ;
    end ;
  end ;

  if isfield(paramstruct,'npValplot') ;    %  then change to input value
    npValplot = paramstruct.npValplot ; 
  end ;
  
  if isfield(paramstruct,'legendcellstr') ;    %  then change to input value
    legendcellstr = paramstruct.legendcellstr ; 
  end ;

  if isfield(paramstruct,'mlegendcolor') ;    %  then change to input value
    mlegendcolor = paramstruct.mlegendcolor ; 
  end ;

  if isfield(paramstruct,'iscreenwrite') ;    %  then change to input value
    iscreenwrite = paramstruct.iscreenwrite ; 
  end ;

end ;    %  of resetting of input parameters
    


%  set preliminary stuff
%
d = size(data,1) ;
         %  dimension of each data curve
n = size(data,2) ;
         %  number of data curves

if iscell(hmetric) ;  %  initialize p-value matrices
  pvalQ = zeros(n-1,4) ;
  pvalZ = zeros(n-1,4) ;
else
  pvalQ = zeros(n-1,2) ;
  pvalZ = zeros(n-1,2) ;
end ;

if testn < 3 || testn > n ;
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  disp('!!!   Warning from HSigClustPK.m:    !!!') ;
  disp('!!!   testn must be value between    !!!') ;
  disp('!!!   3 and n,                       !!!') ;
  disp('!!!   using default of 10            !!!') ;
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  testn = 10 ;
end ;    



CurFigNum = 0 ;
if iscreenwrite == 2 ;
  iscreenwritemost = 1 ;
      %  for use as arguments in other calls
elseif iscreenwrite == 3 ;
  iscreenwritemost = 1 ;
      %  for use as arguments in other calls
else
  iscreenwritemost = 0 ;
      %  for use as argument in other calls
end ;

if  ~isempty(mlegendcolor) ;
  if ~(size(legendcellstr,2) == size(mlegendcolor,1)) ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp('!!!   Warning from HSigClustPK.m:       !!!') ;
    disp('!!!   legendcellstr &  mlegendcolor     !!!') ;
    disp('!!!   must have the same number         !!!') ;
    disp('!!!   of entries                        !!!') ;
    disp('!!!   Note: this could be caused by     !!!') ;
    disp('!!!   using "setfield(paramstruct..."   !!!') ;
    disp('!!!   with "lengedcellstr" defined      !!!') ;
    disp('!!!   by double braces                  !!!') ;
    disp('!!!   of entries                        !!!') ;
    disp('!!!   Resetting to no legend            !!!') ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    legendcellstr = [] ;
  end ;
end ;

if  ~isempty(mlegendcolor) ;
  if ~(size(mlegendcolor,2) == 3) ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp('!!!   Warning from SigClustSM.m:    !!!') ;
    disp('!!!   mlegendcolor                  !!!') ;
    disp('!!!   must have 3 columns           !!!') ;
    disp('!!!   Resetting to no legend        !!!') ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    legendcellstr = [] ;
  end ;
end ;



%  Compute Initial hierarchical clustering
%
Z = linkage(data', hmethod, hmetric) ;
clusters = num2cell(1:(2*n-1))' ;
cn = zeros(1,n-1) ;

if iscell(hmetric) ; %  create (n-1)x4 matrix for average linkage
                     %  with squared euclidean distance
  mCIdata = zeros(n-1,4) ;
else
  mCIdata = zeros(n-1,3) ;
end ;

for c = 1:n-1 ;
    
  %  1st cluster index, 2-means CI at each joining
  %
  clusters{n+c} = [clusters{Z(c,1:2)}] ;
  cn(c) = length(clusters{n+c}) ;
  if cn(c) < testn ;
    mCIdata(c,1) = 0 ;
  else
    tc1 = zeros(1,n) ;
    tc1(clusters{Z(c,1)}) = 1 ;
    tc1 = tc1(clusters{n+c}) ;
    tc2 = zeros(1,n) ;
    tc2(clusters{Z(c,2)}) = 1 ;      
    tc2 = tc2(clusters{n+c}) ;
    mCIdata(c,1) = ClustIndSM(data(:,clusters{n+c}), ...
                      logical(tc1),logical(tc2));
  end ;
  
  %  2nd cluster index, linkage at each joining
  %
  mCIdata(c,2) = Z(c,3) ;
  

%   if Z(c,1) <= n ; difflink1 = 0 ; nlink1 = 1 ;
%   else difflink1 = Z( Z(c,1)-n, 3) ; nlink1 = cn(Z(c,1)-n) ;
%   end ;
% 
%   if Z(c,2) <= n ; difflink2 = 0 ; nlink2 = 1 ;
%   else difflink2 = Z( Z(c,2)-n, 3) ; nlink2 = cn(Z(c,2)-n) ;
%   end ;
% 
%   %  3rd cluster index, difference of linkage and average of lower
%   %  linkage
%   mCIdata(c,3) = (nlink1/cn(c)*difflink1 + nlink2/cn(c)*difflink2) ...
%                   / Z(c,3) ;
        
  if iscell(hmetric) ;

    %  4th cluster index, 2nd modified version of linkage at each joining
    %
    if cn(c) < testn  ;
      mCIdata(c,3) = 0 ;
      mCIdata(c,4) = 0 ;
    else
      mdat = data(:,clusters{n+c}) - vec2matSM(mean(data(:,clusters{n+c}),2),cn(c)) ;
      mCIdata(c,4) = Z(c,3) / (sum(diag(mdat'*mdat))*cn(c)) ;
%       mCIdata(c,3) = sum(tc1)*sum(tc2) * mCIdata(c,4) ;
      mCIdata(c,3) = sum(tc1)*sum(tc2) * Z(c,3) ;
    end ;

  end ;    %  end of average linkage if-block    
  
end ;

clusters = clusters(n+1:end) ;
   %  cell array of sample indices corresponding to each joining



%  Estimate covariance structure of data
%
%  Start with sample eignevalues
%
msimeigval = zeros(d,n-1) ;

paramstruct = struct('iscreenwrite',iscreenwritemost, ...
                   'viout',1) ;
outstruct = pcaSM(data,paramstruct) ;
veigval = outstruct.veigval ;
nev = length(veigval) ;
veigval = [veigval; zeros(d-nev,1)] ;
   %  Pad out to have length d
varflag = 0 ;
   %  warning flag about possible bad background estimation
      
if iCovEst ~= 2 ;    %  then need to consider thresholding

    %  code for using PC scores to calculate background SD.
    %     when using PC scores, important to scale MAD*sqrt(n/d)
%   paramstruct = struct('iscreenwrite',iscreenwritemost, ...
%                        'viout',[0 0 0 0 1]) ;
%   outstruct = pcaSM(data,paramstruct) ;
%   mpc = outstruct.mpc ;
%   npix = n * nev ;
%   disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
%   disp('!!!   Message from HSigClustPK.m:                   !!!') ;
%   disp('!!!   Background variance estimated based on        !!!') ;
%   disp(['!!!   projection on ' num2str(nev) ' eigen-directions            !!!']) ;
%   disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
%   vdata = reshape(mpc',npix,1) ;    

  npix = d * n ;
  data = data - repmat(median(data,2),1,n);
  vdata = reshape(data, npix, 1);
      %  all pixel data
  asd = std(vdata);
  avar = asd^2;

  
  if sigbackg == 0 ;    %  Need to estimate background standard deviation
                        %     same background is used for all subsets 
                        %     of data


    amean = mean(vdata) ;
    amedian = median(vdata) ;
    amad = madSM(vdata) ;
        %  MAD, but on sd scale


    if iBGSDdiagplot ~= 0 ;    %  Then make BackGround Standard Deviation
                               %        DIAGnostic PLOTs

      CurFigNum = CurFigNum + 1 ;
      figure(CurFigNum) ;
      clf ;
      maxnol = 5000 ;
      paramstruct = struct('ndatovlay',maxnol, ...
                           'datovlaymax',0.65, ...
                           'datovlaymin',0.15, ...
                           'iscreenwrite',iscreenwritemost) ;
      kdeSM(vdata,paramstruct) ;
      
      if isempty(datastr) ;
        titstr = 'Distribution of All Pixel values combined (median centered)' ;
        titstr2 = 'Distribution of All Pixel values combined (better view)' ;
      else
        titstr = ['Distribution of All Pixel values combined (median centered), ' datastr] ;
        titstr2 = ['Distribution of All Pixel values combined (better view), ' datastr] ;
      end ;
      title(titstr,'FontSize',15) ;
      vax = axis ;
      
      hold on ;
      xgrid = linspace(vax(1),vax(2),401)' ;
        normden = nmfSM(xgrid,amedian,amad^2,1) ;
      plot(xgrid,normden,'r--') ;
      if length(vdata) < maxnol ;
        olstr = ['Overlay of '...
                 num2str(length(vdata)) ' data points'] ;
      else
        olstr = ['Overlay of ' num2str(maxnol) ' of '...
                 num2str(length(vdata)) ' data points'] ;
      end ;
      text(vax(1) + 0.1 * (vax(2) - vax(1)), ...
           vax(3) + 0.9 * (vax(4) - vax(3)), ...
           olstr,'FontSize',15) ;
      text(vax(1) + 0.1 * (vax(2) - vax(1)), ...
           vax(3) + 0.8 * (vax(4) - vax(3)), ...
           ['Mean = ' num2str(amean,3) ...
            ',     median = ' num2str(amedian,3)], ...
           'FontSize',15) ;
      text(vax(1) + 0.1 * (vax(2) - vax(1)), ...
           vax(3) + 0.7 * (vax(4) - vax(3)), ...
           ['s.d. = ' num2str(asd,3) ...
            ',     MAD = ' num2str(amad,3)], ...
           'FontSize',15) ;
      text(vax(1) + 0.1 * (vax(2) - vax(1)), ...
           vax(3) + 0.6 * (vax(4) - vax(3)), ...
           ['Gaussian(' num2str(amedian,3) ...
            ',' num2str(amad,3) ') density'], ...
           'FontSize',15, 'Color','r') ;
      if amad > asd ;
        text(vax(1) + 0.1 * (vax(2) - vax(1)), ...
             vax(3) + 0.55 * (vax(4) - vax(3)), ...
             'Warning: MAD > s.d., SigClust can be Anti-Conservative', ...
             'FontSize',24, 'Color','m') ;
      end ;
      hold off ;

      if ~isempty(BGSDsavestr) ;
        orient landscape ;
          savestr = [BGSDsavestr 'KDE.ps'] ;
        print('-dpsc2',savestr) ;
      end ;

        
        %   Zoomed in version of the BackGround Standard Deviation
        %   DIAGnostic PLOTs only showing range of +/- 6 MAD
        %      particularly useful when spike covariance model is true
      CurFigNum = CurFigNum + 1 ;
      figure(CurFigNum) ;
      clf ;
      paramstruct = struct('ndatovlay',maxnol, ...
                           'datovlaymax',0.65, ...
                           'datovlaymin',0.15, ...
                           'vxgrid',[-5*amad; 5*amad], ...
                           'eptflag',1, ...
                           'iscreenwrite',iscreenwritemost) ;
      kdeSM(vdata,paramstruct) ;
      title(titstr2,'FontSize',15) ;
      vax = axis ;
      axis([-6*amad, 6*amad, vax(3), vax(4)]) ;
      vax = axis ;

      hold on ;
      xgrid = linspace(vax(1),vax(2),401)' ;
      normden = nmfSM(xgrid,amedian,amad^2,1) ;
      plot(xgrid,normden,'r--') ;
      if length(vdata) < maxnol ;
      olstr = ['Overlay of '...
               num2str(length(vdata)) ' data points'] ;
      else
      olstr = ['Overlay of ' num2str(maxnol) ' of '...
               num2str(length(vdata)) ' data points'] ;
      end ;
      text(vax(1) + 0.1 * (vax(2) - vax(1)), ...
         vax(3) + 0.9 * (vax(4) - vax(3)), ...
         olstr,'FontSize',15) ;
      text(vax(1) + 0.1 * (vax(2) - vax(1)), ...
         vax(3) + 0.8 * (vax(4) - vax(3)), ...
         ['Mean = ' num2str(amean,3) ...
          ',     median = ' num2str(amedian,3)], ...
         'FontSize',15) ;
      text(vax(1) + 0.1 * (vax(2) - vax(1)), ...
         vax(3) + 0.7 * (vax(4) - vax(3)), ...
         ['s.d. = ' num2str(asd,3) ...
          ',     MAD = ' num2str(amad,3)], ...
         'FontSize',15) ;
      text(vax(1) + 0.1 * (vax(2) - vax(1)), ...
         vax(3) + 0.6 * (vax(4) - vax(3)), ...
         ['Gaussian(' num2str(amedian,3) ...
          ',' num2str(amad,3) ') density'], ...
         'FontSize',15, 'Color','r') ;
      if amad > asd ;
      text(vax(1) + 0.1 * (vax(2) - vax(1)), ...
           vax(3) + 0.55 * (vax(4) - vax(3)), ...
           'Warning: MAD > s.d., SigClust can be Anti-Conservative', ...
           'FontSize',24, 'Color','m') ;
      end ;
      hold off ;
      
      if ~isempty(BGSDsavestr) ;
        orient landscape ;
          savestr = [BGSDsavestr 'KDE_zoom.ps'] ;
        print('-dpsc2',savestr) ;
      end ;
              
        

      CurFigNum = CurFigNum + 1 ;
      figure(CurFigNum) ;
      clf ;
      savestr = [BGSDsavestr 'QQ'] ;

      if isempty(datastr) ;
        titstr = 'Robust Fit Gaussian Q-Q, All Pixel values (median centered)' ;
        titstr2 = 'Robust Fit Gaussian Q-Q, All Pixel values (better view)' ;
      else
        titstr = ['Robust Fit Gaussian Q-Q, All Pixel values (median centered), ' datastr] ;
        titstr2 = ['Robust Fit Gaussian Q-Q, All Pixel values (better view), ' datastr] ;
      end ;

      paramstruct = struct('idist',1, ...
                           'mu',amedian, ...
                           'sigma',amad, ...
                           'nsim',0, ...
                           'nsimplotval',900, ...
                           'savestr',savestr, ...
                           'titlestr',titstr, ...
                           'titlefontsize',15, ...
                           'labelfontsize',15, ...
                           'parfontsize',15, ...
                           'ishowpar',1, ...
                           'vshowq',[0.25; 0.5; 0.75], ...
                           'iscreenwrite',iscreenwritemost) ;
      qqLM(vdata,paramstruct) ;

      
        %   Zoomed in version of the BackGround Standard Deviation
        %   DIAGnostic PLOTs only showing range of +/- 6 MAD
        %      particularly useful when spike covariance model is true
      CurFigNum = CurFigNum + 1 ;
      figure(CurFigNum) ;
      clf ;
      paramstruct = struct('idist',1, ...
                           'mu',amedian, ...
                           'sigma',amad, ...
                           'nsim',0, ...
                           'nsimplotval',900, ...
                           'savestr',[savestr '_zoom'], ...
                           'top',5*amad, ...
                           'bottom',-5*amad, ...
                           'titlestr',titstr2, ...
                           'titlefontsize',15, ...
                           'labelfontsize',15, ...
                           'parfontsize',15, ...
                           'ishowpar',1, ...
                           'vshowq',[0.25; 0.5; 0.75], ...
                           'iscreenwrite',iscreenwritemost) ;
      qqLM(vdata,paramstruct) ;      

      
    end ;    %  of iBGSDdiagplot if-block

    isigbackg = amad ;

  else    %  Use input value of background standard deviation

    isigbackg = sigbackg ;

  end ;    %  of sigbackg if-block


  simbackvar = isigbackg^2 ;
      %  Background variance to use in simulations


  %  Check whether background estimate is sensible,
  %  in sense of being smaller than sample variance
  if simbackvar > avar ;    %  then have surprising background var est,
                            %      so give warning, and recommendations
    varflag = 1 ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp('!!!   Warning from SigClustSM.m:                    !!!') ;
    disp('!!!   Background Variance Estimate is Suspect       !!!') ;
    disp('!!!   Because (MAD estimated sig2BG) >              !!!') ;
    disp('!!!                  > (full matrix sample sd),     !!!') ;
    disp(['!!!   by factor of ' num2str(simbackvar / avar)]) ;
    disp('!!!   This could make SigClust Anti-Conservative    !!!') ;
    disp('!!!   Recommend Careful look at Diagnostic plots,   !!!') ;
    disp('!!!       using:   iBGSDdiagplot = 1                !!!') ;
    disp('!!!   Consider re-running SigClust test,            !!!') ;
    disp('!!!       using:   iCovEst = 2                      !!!') ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  end ;

end ;    %  of (iCovEst ~= 2) if-block



for c = 1:(n-1) ;   %  compute simulation eigenvalues for all n-1 clusters
                    %      using background estimate from full data

  paramstruct = struct('iscreenwrite',iscreenwritemost, ...
                       'viout',1) ;
  outstruct = pcaSM(data(:,clusters{c}),paramstruct) ;
  ceigval = outstruct.veigval ;
  nev = length(ceigval) ;
  ceigval = [ceigval; zeros(d-nev,1)] ; %#ok
  
  %  Compute eigenvalues to use in Gaussian Simulation
  %
  if iCovEst == 3 ;    %  Use Original Background Noise Thresholded Estimate
                       %     (from Liu, et al, JASA paper)
                       %     (In Yuan's terminology, this is "Hard Thresholding")

    msimeigval(:,c) = max(ceigval,simbackvar) ;
        %  vector of eigenvalues for simulation

  elseif iCovEst == 1 ;    %  Use Ming Yuan's Soft Thresholded, Constrained MLE

    msimeigval(:,c) = SigClustCovEstSM(ceigval,simbackvar) ;
        %  vector of eigenvalues for simulation

  else   %  no adjustment to eigenvlues, just use empiricals

    msimeigval(:,c) = ceigval ;
        %  vector of eigenvalues for simulation
        %  based on raw covariance type model

  end ;
end ;



if iCovEdiagplot == 1 ;    %  Then make COVariance Estimation DIAGnostic PLOT
                           %     plot is only created for full data and not
                           %     all n-1 subsets (clusters)

  ncut = 100 ;
  CurFigNum = CurFigNum + 1 ;
  figure(CurFigNum) ;
  clf ;

  veigvalpos = veigval(veigval > 10^(-12)) ;
  dpos = length(veigvalpos) ;

  subplot(2,2,1) ;
    plot((1:d)',veigval,'ko-') ;
      if isempty(datastr) ;
        titstr = 'Eigenvalues' ;
      else
        titstr = ['Eigenvalues, ' datastr] ;
      end ;
      title(titstr,'FontSize',12) ;
      xlabel('Component #','FontSize',12) ;
      ylabel('Eigenvalue','FontSize',12) ;
      if iCovEst == 2 ;
        vax = axisSM(veigval) ;
      else
        vax = axisSM([veigval; simbackvar]) ;
      end ;
      vax = [0 (d+1) vax(1) vax(2)] ;
      axis(vax) ;
      hold on ;
        plot([ncut + 0.5; ncut + 0.5],[vax(3); vax(4)],'g-') ;
        plot((1:d)',msimeigval(:,c),'r--','LineWidth',2) ;
        if iCovEst ~= 2 ;
          plot([0; d + 1],[simbackvar; simbackvar],'m-') ;
          text(vax(1) + 0.1 * (vax(2) - vax(1)), ...
               vax(3) + 0.9 * (vax(4) - vax(3)), ...
               ['Background variance = ' num2str(simbackvar,3)], ...
               'FontSize',12,'Color','m') ;
        end ;
        text(vax(1) + 0.1 * (vax(2) - vax(1)), ...
             vax(3) + 0.8 * (vax(4) - vax(3)), ...
             'Eigenvalues for Simulation', ...
             'FontSize',12,'Color','r') ;
        if varflag == 1 ;
          text(vax(1) + 0.1 * (vax(2) - vax(1)), ...
               vax(3) + 0.65 * (vax(4) - vax(3)), ...
               'Warning: MAD > s.d.,', ...
               'FontSize',18,'Color','m') ;
        end ;
      hold off ;
  subplot(2,2,2) ;
    plot((1:dpos)',log10(veigvalpos),'ko-') ;
      title('log_{10} Eigenvalues','FontSize',12) ;
      xlabel('Component #','FontSize',12) ;
      ylabel('log_{10}(Eigenvalue)','FontSize',12) ;
      if iCovEst == 2 ;
        vax = axisSM(log10(veigvalpos)) ;
      else
        vax = axisSM(log10([veigvalpos; simbackvar])) ;
      end ;
      vax = [0 (d+1) vax(1) vax(2)] ;
      axis(vax) ;
      hold on ;
        plot([ncut + 0.5; ncut + 0.5],[vax(3); vax(4)],'g-') ;
        plot((1:d)',log10(msimeigval(:,c)),'r--','LineWidth',2) ;
        if iCovEst ~= 2 ;
          plot([0; d + 1],log10([simbackvar; simbackvar]),'m-') ;
          text(vax(1) + 0.1 * (vax(2) - vax(1)), ...
               vax(3) + 0.9 * (vax(4) - vax(3)), ...
               ['log_{10} Background variance = ' num2str(log10(simbackvar),3)], ...
               'FontSize',12,'Color','m') ;
        end ;
        text(vax(1) + 0.1 * (vax(2) - vax(1)), ...
             vax(3) + 0.8 * (vax(4) - vax(3)), ...
             'Eigenvalues for Simulation', ...
             'FontSize',12,'Color','r') ;
        if varflag == 1 ;
          text(vax(1) + 0.1 * (vax(2) - vax(1)), ...
               vax(3) + 0.65 * (vax(4) - vax(3)), ...
               'SigClust may be Anti-Conservative', ...
               'FontSize',18,'Color','m') ;
        end ;
      hold off ;
  if length(veigval) >= ncut ;
    subplot(2,2,3) ;
      plot((1:ncut)',veigval(1:ncut),'ko-') ;
        title('Zoomed in version of above','FontSize',12) ;
        xlabel('Component #','FontSize',12) ;
        ylabel('Eigenvalue','FontSize',12) ;
        vax = axisSM(veigval(1:ncut)) ;
        vax = [0 (ncut+1) vax(1) vax(2)] ;
        axis(vax) ;
        hold on ;
          plot((1:ncut)',msimeigval(1:ncut,c),'r--','LineWidth',2) ;
          if iCovEst ~= 2 ;
            plot([0; ncut + 1],[simbackvar; simbackvar],'m-') ;
            text(vax(1) + 0.1 * (vax(2) - vax(1)), ...
                 vax(3) + 0.9 * (vax(4) - vax(3)), ...
                 ['Background variance = ' num2str(simbackvar,3)], ...
                 'FontSize',12,'Color','m') ;
          end ;
          text(vax(1) + 0.1 * (vax(2) - vax(1)), ...
               vax(3) + 0.8 * (vax(4) - vax(3)), ...
               'Eigenvalues for Simulation', ...
               'FontSize',12,'Color','r') ;
        hold off ;
    subplot(2,2,4) ;
      plot((1:min(dpos,ncut))',log10(veigvalpos(1:min(dpos,ncut))),'ko-') ;
        title('Zoomed in version of above','FontSize',12) ;
        xlabel('Component #','FontSize',12) ;
        ylabel('log_{10}(Eigenvalue)','FontSize',12) ;
        if iCovEst == 2 ;
          vax = axisSM(log10(veigvalpos(1:min(dpos,ncut)))) ;
        else
          vax = axisSM(log10([veigvalpos(1:min(dpos,ncut)); simbackvar])) ;
        end ;
        vax = [0 (ncut+1) vax(1) vax(2)] ;
        axis(vax) ;
        hold on ;
          plot((1:ncut)',log10(msimeigval(1:ncut,c)),'r--','LineWidth',2) ;
          if iCovEst ~= 2 ;
            plot([0; d + 1],log10([simbackvar; simbackvar]),'m-') ;
            text(vax(1) + 0.1 * (vax(2) - vax(1)), ...
                 vax(3) + 0.9 * (vax(4) - vax(3)), ...
                 ['log_{10} Background variance = ' num2str(log10(simbackvar),3)], ...
                 'FontSize',12,'Color','m') ;
          end ;
          text(vax(1) + 0.1 * (vax(2) - vax(1)), ...
               vax(3) + 0.8 * (vax(4) - vax(3)), ...
               'Eigenvalues for Simulation', ...
               'FontSize',12,'Color','r') ;
        hold off ;
  end ;

  if ~isempty(CovEsavestr) ;
    orient landscape ;
    savestr = [CovEsavestr '.ps'] ;
    print('-dpsc2',savestr) ;
  end ;

end ;    %  of iCovEdiagplot if-block

if  (iscreenwrite == 1)  ||  (iscreenwrite == 2)  ||  (iscreenwrite == 3)  ;
  if isempty(datastr) ;
    dispstr = 'SigClustSM:  Finished Estimation of Parameters for Gaussian simulation' ;
  else
    dispstr = ['SigClustSM:  Finished Estimation of Parameters for Gaussian simulation, for ' datastr] ;
  end ;
  disp(dispstr) ;
end ;






%  Main Simulation Loop
%
if ~isempty(SimRandstate) ;    %  Then set random number generation seed
  rand('state',SimRandstate) ; %#ok
end ;
if ~isempty(SimRandnstate) ;    %  Then set random number generation seed
  randn('state',SimRandnstate) ; %#ok
end ;

if sum(sum(msimeigval < -1e-5)) > 0 ;
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  disp('!!!   Warning from HSigClustPK.m:      !!!') ;
  disp('!!!   negative eigenvalue produced,    !!!') ;
  disp('!!!   exiting function                 !!!') ;
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  return ;
end ;
msimeigval(msimeigval<0) = 0 ;

vscale = sqrt(msimeigval) ;

if iscell(hmetric) ;  %  initialize simulation matrices
  mCIsim = repmat({zeros(nsim,n-1)},1,4) ;
else
  mCIsim = repmat({zeros(nsim,n-1)},1,3) ;
end ;


for c = 1:n-1 ;

    %  don't need to simulate data if cn(c) < testn
  if cn(c) < testn ;

    mCIsim{1}(:,c) = zeros(nsim,1) ; 
    mCIsim{2}(:,c) = zeros(nsim,1) ; 
    mCIsim{3}(:,c) = zeros(nsim,1) ; 
    if iscell(hmetric) ;
      mCIsim{4}(:,c) = zeros(nsim,1) ; 
    end ;
    
  else 
      
    for isim = 1:nsim ;
        
      if  (iscreenwrite == 2)  || (iscreenwrite == 3) ;
        disp(['      HSigClustPK:  Working on sim ' num2str(isim) ' of ' num2str(nsim)]) ;
      end ;
  
      datsim = randn(d,cn(c)) ;
      datsim = datsim .* vec2matSM(vscale(:,c),cn(c)) ;

      Zsim = linkage(datsim',hmethod,hmetric) ;
      clustersim = num2cell(1:(2*cn(c)-1))' ;

      %  1st cluster index, 2-means CI at each joining
      %
      for k = 1:cn(c)-1 ;
        clustersim{cn(c)+k} = [clustersim{Zsim(k,1:2)}] ;
      end ;
      tc1 = zeros(1,cn(c)) ;
      tc1(clustersim{Zsim(end,1)}) = 1 ;
      tc2 = zeros(1,cn(c)) ;
      tc2(clustersim{Zsim(end,2)}) = 1 ;  
      if sum(tc1+tc2) ~= cn(c) ;
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
        disp(['!!!   tc1+tc2 ~= cn(' num2str(c) ')       !!!']) ;
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
        return ;
      end ;
      mCIsim{1}(isim,c) = ClustIndSM(datsim,logical(tc1),logical(tc2)) ;

      %  2nd cluster index, linkage at each joining
      %
      mCIsim{2}(isim,c) = Zsim(cn(c)-1,3) ;
      
%       if Zsim(cn(c)-1,1) <= cn(c) ; difflink1 = 0 ;
%       else difflink1 = Zsim( Zsim(cn(c)-1,1)-cn(c), 3) ;
%       end ;
%       
%       if Zsim(cn(c)-1,2) <= cn(c) ; difflink2 = 0 ;
%       else difflink2 = Zsim( Zsim(cn(c)-1,2)-cn(c), 3) ;
%       end ;
%       
%       %  3rd cluster index, difference of linkage and average of lower
%       %  linkage
% %       mCIsim{3}(isim,c) = Zsim(cn(c)-1,3) ...
% %                             - sum(tc1)/cn(c)*difflink1 ...
% %                             - sum(tc2)/cn(c)*difflink2 ;
%       mCIsim{3}(isim,c) = (sum(tc1)/cn(c)*difflink1 ...
%                            + sum(tc2)/cn(c)*difflink2) ...
%                            / Zsim(cn(c)-1,3) ;

      if iscell(hmetric) ;
        %  4th cluster index, 2nd modified version of linkage at each joining
        %
        mdatsim = datsim - vec2matSM(mean(datsim,2),cn(c)) ;
        mCIsim{4}(isim,c) = Zsim(cn(c)-1,3) / (sum(diag(mdatsim'*mdatsim))*cn(c)) ;
%         mCIsim{3}(isim,c) = sum(tc1)*sum(tc2) * mCIsim{4}(isim,c) ;
        mCIsim{3}(isim,c) = sum(tc1)*sum(tc2) * Zsim(cn(c)-1,3) ;
      end ;

    end ;
    
  end ;  %  end of simulation if-block for cluster c with cn(c) > testn

  if  (iscreenwrite == 1)  ||  (iscreenwrite == 2)  ||  (iscreenwrite == 3)  ;
    if isempty(datastr) ;
      dispstr = 'HSigClustPK:  Finished main simulation' ;
    else
      dispstr = ['HSigClustPK:  Finished main simulation for ' datastr] ;
    end ;
    disp(dispstr) ;
    disp(['for c = ' num2str(c) ' out of (' num2str(n) '-1)']) ;
  end ;

end ;   %  end of simulation for-block 


if iscell(hmetric) ;
  nCIs = 4 ;
else
  nCIs = 2 ;
end ;

CIstr = {'2-means CI' , ...
         'Linkage' , ...
         'modified linkage 1' , ...
         'modified linkage 2'} ;  

     
for c = 1:n-1 ;     %  output results for all n-1 clustering procedures

  %  Output Final results for Cluster Index
  %
  if ipvaltype ~= 2 ;    %  then compute Quantile based p-value

    if cn(c) < testn ;
      pvalQ(c,:) = 2*ones(1,nCIs) ;
    else
%       pvalQ(c,1) = cprobSM(mCIsim{1}(:,c),mCIdata(c,1)) ;
      for j = 1:nCIs ;
        pvalQ(c,j) = 1-cprobSM(mCIsim{j}(:,c),mCIdata(c,j)) ; %  other tail
      end ;
      pvalQ(c,1) = 1-pvalQ(c,1) ;
%       pvalQ(c,3) = 1-pvalQ(c,3) ;
    end ;
    
  end ;
    
  if ipvaltype ~= 1 ;    %  then compute Gaussian based p-value

    %  vectors of simulated CI means and SDs 
    vsimmean = zeros(1,nCIs) ;
    vsimsd = zeros(1,nCIs) ;
    
    if cn(c) < testn ;
      pvalZ(c,:) = 2*ones(1,nCIs) ;
    else
      vsimmean(1) = mean(mCIsim{1}(:,c)) ;
      vsimsd(1) = std(mCIsim{1}(:,c)) ;
%       pvalZ(c,1) = normcdf((mCIdata(c,1) - vsimmean(1)) / vsimsd(1)) ;
      for j = 1:nCIs ;
        vsimmean(j) = mean(mCIsim{j}(:,c)) ;
        vsimsd(j) = std(mCIsim{j}(:,c)) ;
        pvalZ(c,j) = 1 - normcdf((mCIdata(c,j) - vsimmean(j)) / vsimsd(j)) ; %  other tail
      end ;
      pvalZ(c,1) = 1-pvalZ(c,1) ;
%       pvalZ(c,3) = 1-pvalZ(c,3) ;
    end ;
      
  end ;
  

  if ipValplot == 1 ;    %  Then make p Value plot
           
    if c >= (n-npValplot) && cn(c) >= testn ;

      CurFigNum = CurFigNum + 1 ;
      figure(CurFigNum) ;
      clf ;
                    
      for j = 1:nCIs ;
              
%         subplot(2,nCIs/2,j) ;
        subplot(2,2,j) ;
        
        vax = axisSM([mCIsim{j}(:,c); mCIdata(c,j)]) ;
        if isempty(datastr) ;
          kdetitstr = [num2str(nsim) ' Gaussian Simulated CI for: ' CIstr{j}] ;
        else
          kdetitstr = [num2str(nsim) ' Gaussian Simulated CI for: ' ...
                        CIstr{j} ', ' datastr] ;
        end ;

        kdeparamstruct = struct('vxgrid',vax, ...
                                'linecolor','k', ...
                                'dolcolor','k', ...
                                'ibigdot',1, ...
                                'titlestr',kdetitstr, ...
                                'titlefontsize',12, ...
                                'xlabelstr','Linkage value', ...
                                'ylabelstr','density', ...
                                'labelfontsize',12, ...
                                'datovlaymin',0.4, ...
                                'datovlaymax',0.6, ...
                                'iscreenwrite',iscreenwritemost) ;
        kdeSM(mCIsim{j}(:,c),kdeparamstruct) ;
        vax = axis ;
        
        hold on ;

        plot([mCIdata(c,j); mCIdata(c,j)],[vax(3); vax(4)],'Color','g') ;
        text(vax(1) + 0.5 * (vax(2) - vax(1)), ...
             vax(3) + 0.95 * (vax(4) - vax(3)), ...
             ['join # (from top) = ' num2str(n-c)],'Color','b') ;
        text(vax(1) + 0.5 * (vax(2) - vax(1)), ...
             vax(3) + 0.9 * (vax(4) - vax(3)), ...
             ['Cluster Index = ' num2str(mCIdata(c,j))],'Color','g') ;
        if ipvaltype ~= 2 ;    %  then write Quantile based p-value
          text(vax(1) + 0.5 * (vax(2) - vax(1)), ...
               vax(3) + 0.8 * (vax(4) - vax(3)), ...
               ['p-val (Quantile) = ' num2str(pvalQ(c,j))],'Color','g') ;
        end;
        if ipvaltype ~= 1 ;    %  then write Gaussian based p-value
          text(vax(1) + 0.5 * (vax(2) - vax(1)), ...
               vax(3) + 0.7 * (vax(4) - vax(3)), ...
               ['p-val (Gaussian) = ' num2str(pvalZ(c,j))],'Color','g') ;
          xgrid = linspace(vax(1),vax(2),401)' ;
          plot(xgrid,nmfSM(xgrid,vsimmean(j),vsimsd(j)^2,1),'k') ;
        end ;
        if varflag == 1 ;
          text(vax(1) + 0.5 * (vax(2) - vax(1)), ...
               vax(3) + 0.55 * (vax(4) - vax(3)), ...
               'Warning: MAD > s.d.,', ...
               'FontSize',15,'Color','m') ;
          text(vax(1) + 0.5 * (vax(2) - vax(1)), ...
               vax(3) + 0.45 * (vax(4) - vax(3)), ...
               ['by factor of ' num2str(simbackvar / avar) ','], ...
               'FontSize',15,'Color','m') ;
          text(vax(1) + 0.5 * (vax(2) - vax(1)), ...
               vax(3) + 0.35 * (vax(4) - vax(3)), ...
               'SigClust may be Anti-Conservative', ...
               'FontSize',15,'Color','m') ;
        end ;

        if ~isempty(legendcellstr) ;    %  then add legend
          nlegend = length(legendcellstr) ;
          if isempty(mlegendcolor) ;
            mlegendcolor = vec2matSM(zeros(1,3),nlegend) ;
                %  all black when unspecified
          end ;
          tx = vax(1) + 0.1 * (vax(2) - vax(1)) ;
          for ilegend = 1:nlegend ;
            ty = 0 + ((nlegend - ilegend + 1) / ...
                                 (nlegend + 1)) * (vax(4) - 0) ;
            text(tx,ty,legendcellstr(ilegend),  ...
                      'Color',mlegendcolor(ilegend,:)) ;
          end ;
        end ;
        hold off ;
        
      end ;

      if ~isempty(pValsavestr) ;
        orient landscape ;
        savestr = [pValsavestr '.ps'] ;
        print('-dpsc2',savestr,'-append') ;
      end ;
      
    end ;
    
  end ; %  of ipValplot if-block

end ; %  of cluster for-block


if (iscreenwrite == 1)  ||  (iscreenwrite == 2)  ||  (iscreenwrite == 3)  ;
  if isempty(datastr) ;
    dispstr = ['HSigClustPK:  Finished' ...
                      ' with Qp-value(' num2str(n-1) ') = ' num2str(pvalQ(n-1,:)) ...
                      '   & Zp-value(' num2str(n-1) ') = ' num2str(pvalZ(n-1,:))] ;
  else
    dispstr = ['HSigClustPK:  Finished, for ' datastr ...
                      ', with Qp-value(' num2str(n-1) ') = ' num2str(pvalQ(n-1,:)) ...
                      '   & Zp-value(' num2str(n-1) ') = ' num2str(pvalZ(n-1,:))] ;
  end ;
  disp(dispstr) ;
end ;    


% if (iscreenwrite == 1)  ||  (iscreenwrite == 2)  ||  (iscreenwrite == 3)  ;
%   if isempty(datastr) ;
%     if ipvaltype == 1 ;    %  then write Quantile based p-value
%       dispstr = ['HSigClustPK:  Finished with Qp-value(' num2str(c) ') = ' num2str(pvalQ(c))] ;
%     elseif ipvaltype == 2 ;    %  then write Gaussian based p-value
%       dispstr = ['HSigClustPK:  Finished with Zp-value(' num2str(c) ') = ' num2str(pvalZ(c))] ;
%     else    %  then write both p-values
%       dispstr = ['HSigClustPK:  Finished with Qp-value(' num2str(c) ') = ' num2str(pvalQ(c)) ...
%                                        '   & Zp-value(' num2str(c) ') = ' num2str(pvalZ(c))] ;
%     end ;
%   else
%     if ipvaltype == 1 ;    %  then write Quantile based p-value
%       dispstr = ['HSigClustPK:  Finished, for ' datastr ...
%                                       ', with Qp-value(' num2str(c) ') = ' num2str(pvalQ(c))] ;
%     elseif ipvaltype == 2 ;    %  then write Gaussian based p-value
%       dispstr = ['HSigClustPK:  Finished, for ' datastr ...
%                                       ', with Zp-value(' num2str(c) ') = ' num2str(pvalZ(c))] ;
%     else    %  then write both p-values
%       dispstr = ['HSigClustPK:  Finished, for ' datastr ...
%                                       ', with Qp-value(' num2str(c) ') = ' num2str(pvalQ(c)) ...
%                                        '   & Zp-value(' num2str(c) ') = ' num2str(pvalZ(c))] ;
%     end ;
%   end ;
%   disp(dispstr) ;
% end ;    
    


end 


