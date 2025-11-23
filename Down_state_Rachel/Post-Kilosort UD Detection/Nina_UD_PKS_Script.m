addpath(genpath('C:\Users\BuzsakiPC002_misi\Documents\MATLAB\Down_state_Rachel'))
%% Must be done post kilosort 
clearvars
basepath = pwd;
basename = bz_BasenameFromBasepath(basepath);
load([basepath filesep basename '.session.mat'])
smoothwin = .05; %.0.05 seems like a good smoothingwindow for these MV analyses
startbins = 40;
refineDipEstimate = true;
outlierThreshold = 0.005;
%list all channels to use for UD detection.
MUAchan = [...
    session.extracellular.spikeGroups.channels{1,1}(1:18)];
%MUAchan = [31 30 29 28 27 26 25 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 0 335 334 333 332 331 330 329 328 327 326 325 324 323 322 321 320 287 286 285 284 283 282 281 280 279 278 277 276 275 274 273 272 271 270 269 268 267 266 265 264 263 262 261 260 259 258 257 256 255 254 253 252 251 250 249 248 247 246 245 244 243 242 241 240 223 222 221 220 219 218 217 216 215 214 213 212 211 210 209 208 207 206 205 204 203 202 201 200 199 198 197 196 195 194 193 192 143 142 141 140 139 138 137 136 135 134 133 132 131 130 129 128]; %MV2
% MUAchan = [159 158 157 156 155 154 153 152 151 150 149 148 147 146 145 144 95 94 93 92 91 90 89 88 87 86 85 84 83 82 81 80 79 78 77 76 75 74 73 72 71 70 69 68 67 66 65 64 47 46 45 44 43 42 41 40 39 38 37 36 35 34 33 32 31 30 29 28 27 26 25 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 0]; %MV1
% MUAchan = [64 96]; %CL5
load([basepath filesep basename '.cell_metrics.cellinfo.mat'])
regions = {'CA1', 'DG', 'Cortex', 'Th'};
region_channels = struct();

for r = 1:length(regions)
    region_name = regions{r};
    region_ids = cell_metrics.tags.(region_name); % logical index of units
    region_channels.(region_name) = unique(cell_metrics.maxWaveformCh1(region_ids));
end
MUAchan = region_channels.DG;
[SlowWaves_PKS_DG] = DetectSlowWavesMUA_PostKS(pwd,MUAchan,'smoothwin',smoothwin,'startbins',startbins,'refineDipEstimate',refineDipEstimate);
