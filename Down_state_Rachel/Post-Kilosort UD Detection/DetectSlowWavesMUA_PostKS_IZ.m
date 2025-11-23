

function [SlowWavesMUA_PKS] = DetectSlowWavesMUA_PostKS_IZ(basePath,MUAchan,varargin)

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


%% Load Things
load([basePath filesep baseName '.spikes.cellinfo.mat'])
load([basePath filesep baseName '.sessionInfo.mat'])
load([basePath filesep baseName '.MergePoints.events.mat'])

sr_wideband = sessionInfo.rates.wideband;
sr_lfp = sessionInfo.rates.lfp;
num_chs = length(MUAchan);
kk = 1;
session_length_seconds = MergePoints.timestamps(end,end);
session_length_samples = ceil(MergePoints.timestamps(end,end)*sr_lfp);

spike_train = zeros(session_length_samples,1);
all_spike_list = []; %list of spikes to be sorted for ISIs

for ii = 1:num_chs
    
    UD_units = find(spikes.maxWaveformCh==MUAchan(ii));
    num_UD_units = length(UD_units);
    
    if ~isempty(UD_units)
        
        for ll = 1:num_UD_units
            temp_train = zeros(session_length_samples,1);
            temp_train(round(((spikes.times{1,UD_units(ll)})*sr_lfp))) = 1;
            if length(temp_train) > length(spike_train)
                temp_train = temp_train(1:length(spike_train));
            end
            spike_train = spike_train + temp_train;
            
            all_spike_list = [all_spike_list spikes.times{1,UD_units(ll)}'];
            
        end
    end
    
end

all_spike_list = sort(all_spike_list);
ISIs = diff(all_spike_list);
log_ISIs = log10(ISIs);
% figure('units','normalized','outerposition',[0 0 1 1])
% histogram(log_ISIs)
% plot(ISIs)
% plot(spike_train)

%% Sigmoidal activation
%This is what will be saved to MUA
MUAsmooth = smoothdata(spike_train,"gaussian",round(smoothwin.*1250));

%These post-processes will not be saved.
sig_center = median(MUAsmooth((~(MUAsmooth==0))));
MUA_sigmoid = 1./((exp(1).^(-30*(MUAsmooth-(sig_center*0.375))))+1);

%Trying to fix issue with spikes zeros and negative infinity
biggest_spikes = 1;
MUA_sigmoid = MUA_sigmoid*(60/biggest_spikes);
MUA_sigmoid = MUA_sigmoid+60;
MUA_sigmoid = awgn(MUA_sigmoid,5);
%

data = log10(MUA_sigmoid);

MUA.data = MUAsmooth;

MUA.timestamps = (0:(1/1250):session_length_seconds)';

MUA.detectionParms.channels = MUAchan;

MUA.detectionParms.smoothwin = smoothwin;


%% Find UP/DOWN ints

disp('Bimodal Thresh Detection')

[thresh,cross,bihist,diptest] = bz_BimodalThresh_no_overthresh(data,'Schmidt',stickyTrigger,'startbins',startbins,'maxhistbins',maxbins,'diptest',diptestLog,'pcheck',pcheck,'refineDipEstimate',refineDipEstimate);

UP = MUA.timestamps(cross.upints);
DOWN = MUA.timestamps(cross.downints);



%% Make struct

disp('Saving')

clear SlowWavesMUA_PKS

SlowWavesMUA_PKS.ints.UP = UP;

SlowWavesMUA_PKS.ints.DOWN = DOWN;

SlowWavesMUA_PKS.detectorinfo.detectorname = 'DetectSlowWavesMUA';

SlowWavesMUA_PKS.detectorinfo.detectionParms.channel = MUAchan;

SlowWavesMUA_PKS.detectorinfo.detectionParms.smoothwin = smoothwin;

SlowWavesMUA_PKS.detectorinfo.detectionParms.schmidt = stickyTrigger;

SlowWavesMUA_PKS.detectorinfo.detectionParms.startbins = startbins;

SlowWavesMUA_PKS.detectorinfo.detectionParms.maxhistbins = maxbins;

SlowWavesMUA_PKS.detectorinfo.detectionParms.diptest = diptestLog;

SlowWavesMUA_PKS.detectorinfo.detectionParms.pcheck = pcheck;

SlowWavesMUA_PKS.detectorinfo.detectionParms.data = MUA;




save(fullfile(basePath,[baseName '.SlowWavesMUA_PKS.events.mat']),'SlowWavesMUA_PKS')



%% Save evt file for visual inspection with neuroscope



if saveEvent

    eventtype = ['PKS_DOWN.MUAthresh.' num2str(smoothwin) '.stickyTrigger' num2str(stickyTrigger)]; %string for you that = event type you're saving

    numeventtypes = 2; % you have 3 diff types of events, start, peak, stop

    neuroscopeheader = '.R46';

    

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

    mkdir(fullfile(basePath,'DetectionFigures_PKS'))

    saveas(MUAfig,fullfile(basePath,'DetectionFigures_PKS','distMUA_DOWNdetection_PKS'),'pdf')

    

    % Plot duration UP and DOWN states

    UPdur = SlowWavesMUA_PKS.ints.UP(:,2)-SlowWavesMUA_PKS.ints.UP(:,1);

    DOWNdur = SlowWavesMUA_PKS.ints.DOWN(:,2)-SlowWavesMUA_PKS.ints.DOWN(:,1);

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

    saveas(UPDOWNdur,fullfile(basePath,'DetectionFigures_PKS','distDurations_PKS'),'pdf')

end





end

