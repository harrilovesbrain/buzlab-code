function [thresh,cross,bihist,diptest,overthresh] = bz_BimodalThresh_UD(bimodaldata,varargin)
%[thresh,cross,bihist,diptest,overthresh] = bz_BimodalThresh(bimodaldata) 
%takes bimodal time series data, calculates the threshold between the modes
%(i.e. UP vs DOWN states), and returns the crossing times (i.e. UP/DOWN onset/offset times)
%
%INPUTS
%   bimodaldata     vector of bimodal data
%
%       (optional)  for optional inputs, use: ('inputoption', inputvalue)
%       'maxthresh' sets a maximum threshold
%       'Schmidt'   Schmidt trigger uses halfway points between trough and 
%                   lower/upper peaks for DOWN/UP state thresholds
%                   (default: false)
%       'maxhistbins' Maximum number of hist bins to try before giving up
%       'startbins'  minimum number of hist bins (initial hist)
%       'setthresh'     set your own threshold
%       'diptest'   (true/false) use hardigans dip test for bimodality
%                   (default: true)
%       'zInf'      include 0 and inf as the start/end (default:false)
%       'pcheck'    only calculate bimodality if p<0.05 for HDT
%       'refineDipEstimate' Once find peaks with lower bin count, re-do
%                   hist with 'maxhistbins' bin count to get finer estimate of threshold 
%
%OUTPUS
%   thresh          threshold
%   cross
%       .upints     intervals of UP states [onsets offsets]
%       .downints   intervals of DOWN states [onsets offsets]
%   bihist
%       .bins       bin centers
%       .hist       counts
%   diptest
%       .dip             hartigans dip test statistic
%       .p              p value for bimodal distribution
%
%DLevenstein Spring 2016
%-Updated Spring 2017 to include Schmidt trigger, maxthresh
%-TO DO: add hartigans dip test to determine if bimodal... also to find the
%best threshold?  Or ISO-CUT from Simons foundation?
%%

p = inputParser;
addParameter(p,'maxthresh',inf);
addParameter(p,'maxhistbins',25);
addParameter(p,'startbins',10);
addParameter(p,'Schmidt',false);
addParameter(p,'setthresh',false);
addParameter(p,'diptest',true);
addParameter(p,'zInf',false);
addParameter(p,'pcheck',false);
addParameter(p,'refineDipEstimate',false);

parse(p,varargin{:})

maxthresh = p.Results.maxthresh;
maxhistbins = p.Results.maxhistbins;
startbins = p.Results.startbins;
Schmidt = p.Results.Schmidt;
setthresh = p.Results.setthresh;
DIPTEST = p.Results.diptest;
zInf = p.Results.zInf;
pcheck = p.Results.pcheck;
refineDipEstimate = p.Results.refineDipEstimate;

%%

% Thresh input
if isnumeric(setthresh)
    thresh = setthresh;
    setthresh = true;    
end

%Run hartigansdiptest for bimodality
if DIPTEST
    %nboot = 400; %number of times for bootstrapped significance
    nboot = 1; %don't bootstrap for now; time consuming 
    [diptest.dip, diptest.p_value, ~,~]=HartigansDipSignifTest(bimodaldata,nboot);
    %[diptest.dip, diptest.p_value]=bz_hartigansdiptest(bimodaldata); % this pvalue wrong, update later

    if diptest.p_value>0.05 && pcheck 
        display('Dip test says: not bimodal')
        cross.upints = []; cross.downints = []; thresh=nan; 
        overthresh = nan(size(bimodaldata));
        [bihist.hist,bihist.bins]= hist(bimodaldata,startbins);
        return
    end
else
    diptest = [];
end

%Remove data over the threshold... this is klugey
%overmax = bimodaldata>=maxthresh;
%bimodaldata(overmax) = nan;

numpeaks = 1;
numbins = startbins; 

while numpeaks ~=2
    [bihist.hist,bihist.bins]= hist(bimodaldata,numbins);

    [PKS,LOCS] = findpeaks([0 bihist.hist 0],'NPeaks',2);
    LOCS = sort(LOCS)-1;
    numbins = numbins+1;
    numpeaks = length(LOCS);
    if numbins==maxhistbins && ~setthresh
        display('Unable to find trough')
    	cross.upints = []; cross.downints = []; thresh=nan;
        overthresh = nan(size(bimodaldata));
        return
    end
end

betweenpeaks = bihist.bins(LOCS(1):LOCS(2));
[dip,diploc] = findpeaks(-bihist.hist(LOCS(1):LOCS(2)),'NPeaks',1,'SortStr','descend');

% Optional: calculate finer gradient to find dip 
if refineDipEstimate 
    bimodaldataSub = bimodaldata(bimodaldata >= bihist.bins(LOCS(1)) & bimodaldata <= bihist.bins(LOCS(2)));
    [bihistSub.hist,bihistSub.bins]= hist(bimodaldataSub,maxhistbins); 
    [dip,diploc] = findpeaks(-bihistSub.hist,'NPeaks',1,'SortStr','descend');
end

if refineDipEstimate
    thresh = bihistSub.bins(diploc); 
elseif ~setthresh & ~refineDipEstimate
    thresh = betweenpeaks(diploc);
end
    
   
if thresh>maxthresh 
        display('Threshold over max')
    	cross.upints = []; cross.downints = []; thresh=nan;
        overthresh = nan(size(bimodaldata));
        return
end




%%

%The Schmidt trigger uses thresholds halfway between trough and peak
%on each end
if Schmidt
    threshUP = thresh + 0.5.*(betweenpeaks(end)-thresh);
    threshDOWN = thresh + 0.5.*(betweenpeaks(1)-thresh);
    
    overUP = bimodaldata>threshUP;
    overDOWN = bimodaldata>threshDOWN;
    
    crossup = find(diff(overUP)==1);
    crossdown = find(diff(overDOWN)==-1);
    
    %Delete incomplete (repeat) crossings
    allcrossings = [crossup ones(size(crossup)) ;...
        crossdown zeros(size(crossdown))];
    [~,sortorder] = sort(allcrossings(:,1));
    allcrossings = allcrossings(sortorder,:);
    updownswitch = diff(allcrossings(:,2));
    samestate = find(updownswitch==0)+1;
    allcrossings(samestate,:)=[];
    
    crossup = allcrossings(allcrossings(:,2)==1,1);
    crossdown = allcrossings(allcrossings(:,2)==0,1);
    
else
    overind = bimodaldata>thresh;
    crossup = find(diff(overind)==1);
    crossdown = find(diff(overind)==-1);
end

%%
upforup = crossup;
upfordown = crossup;
downforup = crossdown;
downfordown = crossdown;

switch zInf
    case false
        if (isempty(crossup) || isempty(crossdown)) %Only one crossing...
            cross.upints = []; cross.downints = []; overthresh = nan(size(bimodaldata)); return
        end
        if crossup(1) < crossdown(1)
            upfordown(1) = [];
        end
        if crossdown(end) > crossup(end)
            downfordown(end) = [];
        end
        if crossdown(1) < crossup(1)
            downforup(1) = [];
        end
        if crossup(end) > crossdown(end)
            upforup(end) = [];
        end
    case true
        if isempty(crossup) %Only one crossing...
            upforup = 0; upfordown = Inf;
            
        elseif isempty(crossdown)
            downfordown = 0; downforup = Inf;
            
        else
        
            switch Schmidt
                case false   %Regular threshold assumes first state 
                             %is opposite first crossing
                    if crossup(1) < crossdown(1)
                        downfordown = [0;downfordown];
                    end
                    if crossdown(end) > crossup(end)
                        upfordown = [upfordown;Inf];
                    end
                    if crossdown(1) < crossup(1)
                        upforup = [0;upforup];
                    end
                    if crossup(end) > crossdown(end)
                        downforup = [downforup;Inf];
                    end
                case true   %Sticky threshold assumes first threshold crossing  
                             %is in first state
                    if crossup(1) < crossdown(1)
                        upfordown(1) = [];
                        upforup(1) = 0;
                    end
                    if crossdown(end) > crossup(end)
                        upfordown = [upfordown;Inf];
                    end
                    if crossdown(1) < crossup(1)
                        downforup(1) = [];
                        downfordown(1) = 0;
                    end
                    if crossup(end) > crossdown(end)
                        downforup = [downforup;Inf];
                    end
            end
        
        end
        
end
    


cross.upints = [upforup downforup];
cross.downints = [downfordown upfordown];

overthresh = bz_INTtoIDX({cross.downints,cross.upints},'length',length(bimodaldata))-1;
x=find(overthresh == -1); 
overthresh(overthresh==-1) = nan;
%overthresh = logical(overthresh);



end

