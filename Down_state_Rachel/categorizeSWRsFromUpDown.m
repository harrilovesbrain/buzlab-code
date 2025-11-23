function swrCategories = categorizeSWRsFromUpDown(basepath, label, epoch)
% categorizeSWRsFromUpDown - Categorizes ripples based on timing with up/down states
% INPUTS:
%   basepath - full path to session folder
%   label    - 'Pre' or 'Post'
%   epoch    - [start, end] of epoch to restrict analysis (e.g., MergePoints.timestamps)
% OUTPUT:
%   swrCategories - struct with categorized SWR timestamps

cd(basepath);
basename = bz_BasenameFromBasepath(basepath);

% Load data
load(fullfile(basepath, [basename '.SleepState.states.mat']));
load(fullfile(basepath, [basename '.MergePoints.events.mat']));
load(fullfile(basepath, [basename '.SlowWavesMUA_PKS.events.mat']));
load(fullfile(basepath, [basename '.ripples.events.mat']));

% Restrict to NREM
nremIntervals = SleepState.ints.NREMstate;
nremEpoch = [max(nremIntervals(:,1), epoch(1)), min(nremIntervals(:,2), epoch(2))];

% Restrict events
UPs = restrictIntervals(SlowWavesMUA_PKS.ints.UP, nremEpoch);
DOWNs = restrictIntervals(SlowWavesMUA_PKS.ints.DOWN, nremEpoch);
SWRs = restrictIntervals(ripples.timestamps, nremEpoch);

% Init
SWR_ups = []; dist_UP = []; dist_DOWN = [];

% Find UPs containing each ripple
for i = 1:size(SWRs,1)
    for j = 1:size(UPs,1)
        if SWRs(i,1) >= UPs(j,1) && SWRs(i,1) <= UPs(j,2)
            SWR_ups = [SWR_ups; SWRs(i,1)];
            dist_UP = [dist_UP; UPs(j,1)];
            break;
        end
    end
end

% Find next DOWN state after ripple
for i = 1:length(SWR_ups)
    next_down = find(DOWNs(:,1) > SWR_ups(i), 1);
    if ~isempty(next_down)
        dist_DOWN(i,1) = DOWNs(next_down,1);
    else
        dist_DOWN(i,1) = NaN;
    end
end

% Categorize
Green = SWR_ups(SWR_ups - dist_UP > 0.05 & SWR_ups - dist_UP < 0.18);
Red = SWR_ups(dist_DOWN - SWR_ups < 0.05);
Yellow = SWR_ups(SWR_ups - dist_UP > 0.18 & SWR_ups - dist_UP < 0.2);
Blue = SWR_ups(SWR_ups - dist_UP > 0.001 & SWR_ups - dist_UP < 0.05);

swrCategories.(label).Green = Green;
swrCategories.(label).Red = Red;
swrCategories.(label).Blue = Blue;
swrCategories.(label).Yellow = Yellow;

% Save output
save(fullfile(basepath, [basename '_swrCategories_' label '.mat']), 'swrCategories');
end

function out = restrictIntervals(ints, epoch)
% Helper to restrict interval list to [start, end] range
out = ints(ints(:,1) >= epoch(1) & ints(:,2) <= epoch(2), :);
end
