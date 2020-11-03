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
cpIndex = [13, 14, 25:28, 35, 44:45, 61, 62, 64, 67, 68];
[tform, regMaps] = helper.getTforms(imresize(openFrames, imgResizeFactor), refIdx, plotFigs, verbose, cpIndex);

%% draw mask
mask = helper.draw_roi(mean(regMaps, 3),2);
% tmp = untile(regMaps, 128);
% mask = helper.draw_roi(tmp(:,:,13),2);

%% Global signal analysis
% Pipeline:
%
% 

% [socialData, corrs, GS, B] = helper.loadData(openPath, openFiles, tform, mask, ...
%     'global_signal_correlation_interbrain');

[socialData, corrs, GS, B] = loadData(openPath, openFiles, tform, mask, ...
    'global_signal_correlation_interbrain_corrected');

%% make figure
helper.makeGSfig

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
refIdx = 9;
[tformMesh, regMapsMesh] = helper.getTforms(imresize(meshFrames, imgResizeFactor), refIdx, plotFigs, verbose);
%%
refIdx = 21;
cpindex = [];
[tformOpaque, regMapsOpaque] = helper.getTforms(imresize(opaqueFrames, imgResizeFactor), refIdx, plotFigs, verbose);

%% draw mask
maskOpen = helper.draw_roi(mean(regMapsOpen, 3),2);
%%
maskMesh = helper.draw_roi(mean(regMapsMesh, 3),2);
%%
maskOpaque = helper.draw_roi(mean(regMapsOpaque, 3),2);

%%

% open single frame from all files
openFrames = helper.getSingleFrames(openPath, openFiles, H, W);
meshFrames = helper.getSingleFrames(meshPath, meshFiles, H, W);
opaqueFrames = helper.getSingleFrames(opaquePath, opaqueFiles, H, W);

regMapsOpen = imresize(helper.applyTform(openFrames, tform.open), 1/2);
regMapsMesh = imresize(helper.applyTform(meshFrames, tform.mesh), 1/2);
regMapsOpaque = imresize(helper.applyTform(opaqueFrames, tform.opaque), 1/2);

%%
scale_factor = 10/128;
ROIs = {'M2', 'ALM', 'wM1', 'FL', 'HL', 'aBC', 'pBC', 'V1', 'RS', 'PtA'};
[CLopen,CRopen] = helper.get_coords(mean(regMapsOpen,3), scale_factor, ROIs);
[CLmesh,CRmesh] = helper.get_coords(mean(regMapsMesh,3), scale_factor, ROIs);
%%
[CLopaque,CRopaque] = helper.get_coords(mean(regMapsOpaque,3), 11/128, ROIs);

%%

figure, helper.show_coords(mean(regMapsOpaque,3) .* mask.opaque, ROIs, CLopaque, CRopaque)
%%

[rMat.open, leftTrace.open, rightTrace.open] = ...
    loadData(openPath, openFiles, tform.open, mask.open, B, ...
    'roi_signal_correlation_interbrain_corrected', CLopen, CRopen);

clc
[rMat.mesh, leftTrace.mesh, rightTrace.mesh] = ...
    loadData(meshPath, meshFiles, tform.mesh, mask.mesh, B, ...
    'barrier_controls_corrected', CLmesh, CRmesh);
%%
clc2
[rMat.opaque, leftTrace.opaque, rightTrace.opaque] = ...
    loadData(opaquePath, opaqueFiles, tform.opaque, mask.opaque, B, ...
    'barrier_controls_corrected', CLopaque, CRopaque);


%% BEHAVIOR MODULATION
clc
% test = openFiles(3);
bFrames = loadData(openPath, openFiles, tform.open, mask.open, B, 'behavior_modulation_corrected');


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

