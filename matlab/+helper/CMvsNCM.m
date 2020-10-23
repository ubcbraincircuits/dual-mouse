%% ANALYSIS FOR CAGEMATES VS NON CAGEMATES 20200910

idx = CM;
idx([15,27]) = nan;
cmidx = idx=='Y';
ncmidx = idx == 'N';

cmselfleft = catcell(4, bFrames.selfInitiatedWhiskLeftDuring(cmidx));
cmselfright = catcell(4, bFrames.selfInitiatedWhiskRightDuring(cmidx));
cmself.data = cat(4, cmselfleft, cmselfright);
cmself.order = {'cmselfleft', 'cmselfright'};
cmself.length = [size(cmselfleft, 4), size(cmselfright,4)];
clear cmselfleft cmselfright

ncmselfleft = catcell(4, bFrames.selfInitiatedWhiskLeftDuring(ncmidx));
ncmselfright = catcell(4, bFrames.selfInitiatedWhiskRightDuring(ncmidx));
ncmself.data = cat(4, ncmselfleft, ncmselfright);
ncmself.order = {'ncmselfleft', 'ncmselfright'};
ncmself.length = [size(ncmselfleft, 4), size(ncmselfright,4)];
clear ncmselfleft ncmselfright

cmon = [];
nmon = [];

for i = 1:size(cmself.data,4)
    cmon = cat(3, cmon, helper.makeMontage(cmself.data(:,:,:,i)-mean(cmself.data(:,:,1:29,i),3), fs));
end
for i = 1:size(ncmself.data,4)
    nmon = cat(3, nmon, helper.makeMontage(ncmself.data(:,:,:,i)-mean(ncmself.data(:,:,1:29,i),3), fs));
end

cmon(cmon==0) = nan;
nmon(nmon==0) = nan;

cs = squeeze(nanmean(nanmean(cmself.data,2),1));
ns = squeeze(nanmean(nanmean(ncmself.data,2),1));


cs = cs - mean(cs(1:29,:));
ns = ns - mean(ns(1:29,:));



%%


figure, 

BCMtmp = [B.whisk.left(cmidx), B.whisk.right(cmidx)];
for i = 1:length(BCMtmp)
    BCM(i,:) = BCMtmp{i}(B.times{1}(1):B.times{1}(2));
end

BNCMtmp = [B.whisk.left(ncmidx), B.whisk.right(ncmidx)];
for i = 1:length(BNCMtmp)
    BNCM(i,:) = BNCMtmp{i}(B.times{1}(1):B.times{1}(2));
end
%%
for i = 1:size(BCM,1)
    [bOn, bOff] = helper.getBehaviourEvents(BCM(i,:));
    if any(bOff-bOn) < 0, error('afasf'); end
    cmNumWhisks(i) = numel(bOn);
    cmDurWhisks{i} = bOff - bOn;
end

%%
for i = 1:size(BNCM,1)
    [bOn, bOff] = helper.getBehaviourEvents(BCM(i,:));
    if any((bOff-bOn) < 0), error('asdad'); end
    ncmNumWhisks(i) = numel(bOn);
    ncmDurWhisks{i} = bOff - bOn;
end



%%

figure
[~,I] = sort(mean(cs2(30:end,:)));

subplot(1,2,1), imagesc(cs2(:,I)'), colormap jet; colorbar, caxis([-1.5 1.5])
subplot(1,2,2)
[~,I] = sort(mean(ns2(30:end,:)));
ns3 = [ns2(:,I), nan(59, size(cs2,2)-size(ns2,2))];
imagesc(ns3'), colormap jet; colorbar, caxis([-1.5 1.5])


%%
randidx = randperm(size(cmon,3), size(nmon,3));

figure,
subplot(3,4,1), helper.uboxplot((cmNumWhisks./2)', (ncmNumWhisks./2)');
ylabel('Whisk rate (min^{-1})'), xticklabels({'CM', 'NCM'})
subplot(3,4,2), helper.uboxplot((catcell(2,cmDurWhisks)./fs)', (catcell(2,ncmDurWhisks)./fs)')
ylabel('Whisk duration ( log_{10}(s) )'), xticklabels({'CM','NCM'})
set(gca, 'YScale', 'log')

subplot(3,4,5:7),  hold on, imagesc(flipud(mean(cmon(:,:,randidx),3))), colormap jet; colorbar,
xticklabels(''), yticklabels(''), ylabel('CM n=740'), title('Mean self-initiated whisk montage')
axis([1 1280 1 128])
line([size(cmon,2)/2 size(cmon,2)/2], [1 128], 'color', [1 1 1 1], 'LineWidth', 3)
subplot(3,4,9:11), imagesc(mean(nmon,3)), colormap jet; colorbar,
xticklabels(''), yticklabels(''), ylabel('NCM n=399'), hold on
line([size(cmon,2)/2 size(cmon,2)/2], [1 128], 'color', [1 1 1 1], 'LineWidth', 3)
subplot(3,4,[8,12]), hold on
shadedErrorBar(xt(cs,fs,1)-1,mean(cs(:,randidx),2), std(cs(:,randidx),[],2)./sqrt(size(cs(:,randidx),2)), 'k', 1);
shadedErrorBar(xt(ns,fs,1)-1,mean(ns,2), std(ns,[],2)./sqrt(size(ns,2)), 'b', 1)
xlabel('time (s)')
ylabel('\DeltaF/F_0 (\sigma)')
axis([-1 1 -0.05 0.2])
title('global activation')

% 
% figure,
% subplot(1,2,1), imagesc(mean(cmon,3)), colormap jet; colorbar, caxis([-1 1])
% subplot(1,2,2), imagesc(mean(nmon,3)), colormap jet; colorbar, caxis([-1 1])
