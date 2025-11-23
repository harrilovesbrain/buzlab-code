

function [SlowWaves] = DetectSlowWavesMUA(basePath,MUAchan,varargin)

%% This fx uses MUA detected from .dat file to identify DOWN states

% If log10MUA is bimodal, uses dip in bimodal dist to segregate UP and DOWN

% states. If log10MUA not bimodal, will not work.

%

%INPUTS

%   basePath    -(default: cd) basePath for the recording file, in buzcode format:

%                   whateverPath/baseName/

%               folder must include file:

%                   baseName.dat

%   MUAchan   - LFP channel to use for MUA. 

%

%   (options)

%   'smoothwin' - smoothing window for extracted MU, default = .03

%   'stickyTrigger' - whether to use Schmidt triggering, default = true

%   'startbins' - starting bins for histogram for bimodality detection, default = 100

%   'maxbins' - ending bin for histogram bimodality detection, default = 100

%   'diptestLog' - test for bimodality, default = false 

%   'pcheck' - exit fx if not significant, default = false

%   'makeFig' - makes/saves smoothed MUA hist + detected dip; default true

%   'saveEvent' - save an events file

%   'refineDipEstimate' - default false, true if want to increase bin count

%               when estimating dip in narrow range 

%

%   DEPENDENCIES

%   MUAfromDat

%   bz_BimodalThresh

%   bz_BasenameFromBasepath

%

%   TO DO

%   - consider making self-contained by adding dependent fxs to bottom of

%   this



% RS 2021

%% Parms



p = inputParser;

addParameter(p,'smoothwin',.03);

addParameter(p,'stickyTrigger',true);

addParameter(p,'startbins',100);

addParameter(p,'maxbins',100);

addParameter(p,'diptestLog',true);

addParameter(p,'pcheck',false);

addParameter(p,'makeFig',true);

addParameter(p,'saveEvent',true);

addParameter(p,'refineDipEstimate',false);



parse(p,varargin{:})



smoothwin = p.Results.smoothwin;

stickyTrigger = p.Results.stickyTrigger;

startbins = p.Results.startbins;

maxbins = p.Results.maxbins;

diptestLog = p.Results.diptestLog;

pcheck = p.Results.pcheck;

makeFig = p.Results.makeFig;

saveEvent = p.Results.saveEvent;

refineDipEstimate = p.Results.refineDipEstimate;



baseName = bz_BasenameFromBasepath(basePath);



%%

% Calculate MUA

disp('Detecting MUA from dat')

[ MUA ] = MUAfromDat( basePath,'usepeaks',true,'saveMat',true,'channels',MUAchan,'SHOWFIG',false );

MUAsmooth = smooth(MUA.data,round(smoothwin.*1250),'moving' );

MUA.data = MUAsmooth;

MUA.timestamps = MUA.timestamps;

MUA.detectionParms.channels = MUAchan;

MUA.detectionParms.smoothwin = smoothwin;

%% Find UP/DOWN ints

disp('Bimodal Thresh Detection')

data = log10(MUAsmooth);

[thresh,cross,bihist,diptest,overthresh] = bz_BimodalThresh(data,'Schmidt',stickyTrigger,'startbins',startbins,'maxhistbins',maxbins,'diptest',diptestLog,'pcheck',pcheck,'refineDipEstimate',refineDipEstimate);



% Convert inds to time

UP = MUA.timestamps(cross.upints);

DOWN = MUA.timestamps(cross.downints);



%% Make struct

disp('Saving')

clear SlowWaves

SlowWaves.ints.UP = UP;

SlowWaves.ints.DOWN = DOWN;

SlowWaves.detectorinfo.detectorname = 'DetectSlowWavesMUA';

SlowWaves.detectorinfo.detectionParms.channel = MUAchan;

SlowWaves.detectorinfo.detectionParms.smoothwin = smoothwin;

SlowWaves.detectorinfo.detectionParms.schmidt = stickyTrigger;

SlowWaves.detectorinfo.detectionParms.startbins = startbins;

SlowWaves.detectorinfo.detectionParms.maxhistbins = maxbins;

SlowWaves.detectorinfo.detectionParms.diptest = diptestLog;

SlowWaves.detectorinfo.detectionParms.pcheck = pcheck;

SlowWaves.detectorinfo.detectionParms.data = MUA;



save(fullfile(basePath,[baseName '.SlowWavesMUA.events.mat']),'SlowWaves')



%% Save evt file for visual inspection with neuroscope



if saveEvent

    eventtype = ['DOWN.MUAthresh.' num2str(smoothwin) '.stickyTrigger' num2str(stickyTrigger)]; %string for you that = event type you're saving

    numeventtypes = 2; % you have 3 diff types of events, start, peak, stop

    neuroscopeheader = '.R45';

    

    basenamesave = [baseName '.' eventtype neuroscopeheader '.evt']; % you need the ROX bc neuroscope is buggy and uses this to parse files.

    lengthAll = numel(DOWN);

    events.time = zeros(1,lengthAll);

    events.time(1:2:lengthAll) = DOWN(:,1);

    events.time(2:2:lengthAll) = DOWN(:,2);

    events.description = cell(1,lengthAll);

    events.description(1:2:lengthAll) = {'start'};

    events.description(2:2:lengthAll) = {'stop'};

    SaveEvents(fullfile(basePath,basenamesave),events) %Save and look at in neuroscope

end



%% Figure



if makeFig

    

    % Plot histogram of MUA to assess bimodality

    MUAfig = figure;

    plot(bihist.bins,bihist.hist,'k','linewidth',2)

    axis tight; box off

    hold on

    plot([thresh thresh],get(gca,'ylim'),'r','linewidth',1.2)

    if diptestLog

        title(['Dist MUA(log10+Smooth): diptest ' num2str(diptest.dip) ', p = ' num2str(diptest.p_value)])

    else

        title('Dist MUA(log10+Smooth): diptest _ , pvalue _')

    end

    xlabel('MUA(log10)')

    ylabel('Count')

    set(gcf,'color','white')

    mkdir(fullfile(basePath,'DetectionFigures'))

    saveas(MUAfig,fullfile(basePath,'DetectionFigures','distMUA_DOWNdetection'),'pdf')

    

    % Plot duration UP and DOWN states

    UPdur = SlowWaves.ints.UP(:,2)-SlowWaves.ints.UP(:,1);

    DOWNdur = SlowWaves.ints.DOWN(:,2)-SlowWaves.ints.DOWN(:,1);

    [y,x] = hist(log10(UPdur),30);

    [y2,x2] = hist(log10(DOWNdur),30);

    

    % Make figure

    UPDOWNdur = figure;

    clf

    hold on

    

    yNorm = y./sum(y);

    plot(x,yNorm,'color','r','linewidth',1.7)

    y2Norm = y2./sum(y2);

    plot(x2,y2Norm,'color','k','linewidth',1.7)    

    

    axis tight

    LogScale('x',10,'nohalf',true)

    xlabel('duration(s)')

    ylabel('P(obs)')

    set(gcf,'color','white')

    saveas(UPDOWNdur,fullfile(basePath,'DetectionFigures','distDurations'),'pdf')

end





end

