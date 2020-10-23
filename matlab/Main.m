%% Main script for dual mouse 
%
%
clear, clc, close all

openPath = 'B:\Social_Outputs\Matlab\social\behavior_added\'; 
openFiles = helper.getAllFiles(openPath);

meshPath = 'B:\Social_Outputs\Matlab\mesh\';
meshFiles = helper.getAllFiles(meshPath);

opaquePath = 'B:\Social_Outputs\Matlab\opaque\';
opaqueFiles = helper.getAllFiles(opaquePath);

%% constants
H = 256;
W = 256;
fs = 28.815;
imgResizeFactor = 1/2;
CM = ['Y', 'Y', 'N', 'N', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', ...
    'Y', 'Y', 'Y', 'Y', 'N', 'N', 'N', 'N', 'N', 'Y', ...
    'Y', 'Y', 'Y', 'Y', 'N', 'Y', 'N', 'Y', 'N', 'N', ...
    'Y', 'N', 'N', 'N', 'N'];

DOM = [];

%% register to single reference and get transformation

% open single frame from all files
openFrames = helper.getSingleFrames(openPath, openFiles, H, W);

% refIdx = 27;
% refIdx = [];
refIdx = 51;
plotFigs = 1;
verbose = 1;
cpIndex = [];
[tform, regMaps] = helper.getTforms(imresize(openFrames, imgResizeFactor), refIdx, plotFigs, verbose, cpIndex);

%% draw mask
mask = helper.draw_roi(mean(untile(regMaps, 128), 3),2);
% tmp = untile(regMaps, 128);
% mask = helper.draw_roi(tmp(:,:,13),2);

%% Global signal analysis
% Pipeline:
%
% 

% [socialData, corrs, GS, B] = helper.loadData(openPath, openFiles, tform, mask, ...
%     'global_signal_correlation_interbrain');

[socialData, corrs, GS, B] = helper.loadData(openPath, openFiles, tform, mask, ...
    'global_signal_correlation_interbrain_corrected');

%% make figure
helper.makeGSfig

    
%%

% statistics
t = table(CM',corrs.before,corrs.during,corrs.after,...
'VariableNames',{'CM','t0','t2','t4'});
Time = [0 2 4]';
rm = fitrm(t,'t0-t4 ~ CM','WithinDesign',Time);
anovatbl = anova(rm)

% print summary
if ~kstest(zscore(atanh(corrs.before))) && ~kstest(zscore(atanh(corrs.during))) && ~kstest(zscore(atanh(corrs.after)))
[h,p] = ttest(atanh(corrs.before), atanh(corrs.during));
disp(['together vs before; p=', num2str(p*3)])
[h,p] = ttest(atanh(corrs.after), atanh(corrs.during));
disp(['together vs after; p=', num2str(p*3)])
[h,p] = ttest(atanh(corrs.before), atanh(corrs.after));
disp(['before vs after; p=', num2str(p*3)])
end

[h,p] = ttest(atanh(corrs.during), atanh(duringCshuffle));
disp(['together vs shuffle ', num2str(p)])

%%

%% Region by region correlations

% open single frame from all files
openFrames = helper.getSingleFrames(openPath, openFiles, H, W);
meshFrames = helper.getSingleFrames(meshPath, meshFiles, H, W);
opaqueFrames = helper.getSingleFrames(opaquePath, opaqueFiles, H, W);

% refIdx = 27;
% refIdx = [];

plotFigs = 1;
verbose = 1;
%%
refIdx = 51;
[tformOpen, regMapsOpen] = helper.getTforms(imresize(openFrames, imgResizeFactor), refIdx, plotFigs, verbose);
%%
refIdx = 14;
[tformMesh, regMapsMesh] = helper.getTforms(imresize(meshFrames, imgResizeFactor), refIdx, plotFigs, verbose);
%%
refIdx = 15;
[tformOpaque, regMapsOpaque] = helper.getTforms(imresize(opaqueFrames, imgResizeFactor), refIdx, plotFigs, verbose);

%% draw mask
maskOpen = helper.draw_roi(mean(untile(regMapsOpen, 128), 3),2);
%%
maskMesh = helper.draw_roi(mean(untile(regMapsMesh, 128), 3),2);
%%
maskOpaque = helper.draw_roi(mean(untile(regMapsOpaque, 128), 3),2);

%%

scale_factor = 10.25/128;
ROIs = {'M2', 'ALM', 'wM1', 'FL', 'HL', 'aBC', 'pBC', 'V1', 'RS', 'PtA'};
[CLopen,CRopen] = get_coords(mean(untile(regMapsOpen, 128),3), scale_factor, ROIs);
[CLmesh,CRmesh] = get_coords(mean(untile(regMapsMesh, 128),3), scale_factor, ROIs);
[CLopaque,CRopaque] = get_coords(mean(untile(regMapsOpaque, 128),3), scale_factor, ROIs);
%%

[rMat.open, leftTrace.open, rightTrace.open] = ...
    helper.loadData(openPath, openFiles, tformOpen, maskOpen, ...
    'roi_signal_correlation_interbrain_corrected', ROIs, CLopen, CRopen);
%%
clc
[rMat.mesh, leftTrace.mesh, rightTrace.mesh] = ...
    helper.loadData(meshPath, meshFiles, tformMesh, maskMesh, ...
    'barrier_controls_corrected', ROIs, CLmesh, CRmesh);
%%
clc
[rMat.opaque, leftTrace.opaque, rightTrace.opaque] = ...
    helper.loadData(opaquePath, opaqueFiles, tformOpaque, maskOpaque, ...
    'barrier_controls_corrected', ROIs, CLopaque, CRopaque);

%% intrabrain




%% R by ROI
r_by_roi = [squeeze(median(median(rbao,2),3)), ...
    squeeze(median(median(rdao,2),3)), ...
    squeeze(median(median(raao,2),3))]';
figure, plot(r_by_roi)


% r_by_roi = [squeeze(median(median(rbinter,2),3)), ...
%     squeeze(median(median(rdinter,2),3)), ...
%     squeeze(median(median(rainter,2),3))]';
% figure, plot(r_by_roi)
%%
test = [];
test2 = [];
% figure
for i = 1:23
    test(i,:) = movcorr(leftTrace{i}(1:11818,6), rightTrace{i}(1:11818,6), round(10*fs));
    test2(i,:) = movcorr(leftTrace{i}(1:11818,9), rightTrace{i}(1:11818,9), round(10*fs));
%     hold on
end
    
figure, plot(xt(test,fs,2), mean(test))
hold on, plot(xt(test2,fs,2), mean(test2))

%% BEHAVIOR MODULATION
clc
% test = openFiles(1:2);
[bFrames, B] = helper.loadData(openPath, openFiles, tform, maskOpen, 'behavior_modulation_corrected');


%%

%% BARRIER CONTROLS
refIdx = [];
plotFigs = 1;
verbose = 1;
[tform, regMaps] = helper.getTforms(imresize(meshFrames, imgResizeFactor), refIdx, plotFigs, verbose);

%%

mask = helper.draw_roi(mean(untile(regMaps, 128), 3),2);

%%
clc

[rMat, leftTrace, rightTrace, ROIs, CL, CR, B] = helper.loadData(meshPath, meshFiles, tform, mask, ...
    'barrier_controls', regMaps);
