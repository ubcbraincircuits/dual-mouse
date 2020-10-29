%% whisk triggered figure


% selfWhisk = [];
selfWhiskAvg = [];
% partnerWhisk = [];
% selfFL = [];
% partnerFL = [];
% partnerWhiskSolo = [];
% selfWhiskSolo = [];
% selfWhiskSoloAvg = [];
durmax = 0.5;
durmin = 0.1;
chooseidxanalyze = [];
for i = 1:length(CM)
    if i == 15 || i == 27, continue; end
    before_wcount(i) = size(X{i}, 4);
    before_w_rate(i) = before_wcount(i)/(B.times{i}(3)/fs);
    during_wcount(i) = size(Y{i},4);
    during_w_rate(i) = during_wcount(i)/((B.times{i}(2)-B.times{i}(1))/fs);
%     
%     before_fcount(i) = size(bFrames.selfInitiatedFLLeftBefore{i},4);
%     before_f_rate(i) = before_fcount(i)/(B.times{i}(3)/fs);
%     during_fcount(i) = size(bFrames.selfInitiatedFLLeftDuring{i},4);
%     during_f_rate(i) = during_fcount(i)/((B.times{i}(2)-B.times{i}(1))/fs);
    


%     selfWhisk = cat(4, selfWhisk,bFrames.selfInitiatedWhiskLeftDuring{i} - ...
%         mean(bFrames.selfInitiatedWhiskLeftDuring{i}(:,:,1:29,:),3));
%     selfWhisk = cat(4, selfWhisk, bFrames.selfInitiatedWhiskRightDuring{i}- ...
%         mean(bFrames.selfInitiatedWhiskRightDuring{i}(:,:,1:29,:),3));
    
    
    
    [wOn, wOff] = helper.getBehaviourEvents(B.exclusive.whiskLeftDuring{i}); % get all whisk events
    test = helper.getWhiskEvents(B.exclusive.whiskLeftDuring{i}, round(fs*2)); % get whisk events excluding ones happening in window
    [~, ~, ia] = intersect(wOn+1, test); % find index of intersection
    dur = (wOff-wOn)/fs;
    dur = dur(ia);
    chooseidx = (dur <= durmax & dur >= durmin);
    chooseidxanalyze = [chooseidxanalyze, [sum(chooseidx); sum(chooseidx)/numel(chooseidx)]];
    
    
    tmp = bFrames.selfInitiatedWhiskLeftDuring{i}(:,:,:,chooseidx);
%     tmp = bFrames.partnerInitiatedWhiskRightDuring{i}(:,:,:,chooseidx);
    selfWhiskAvg = cat(4, selfWhiskAvg, mean(tmp - ...
        mean(tmp(:,:,1:29,:),3),4));
    
    [wOn, wOff] = helper.getBehaviourEvents(B.exclusive.whiskRightDuring{i}(B.times{i}(1):B.times{i}(2))); % get all whisk events
    test = helper.getWhiskEvents(B.exclusive.whiskRightDuring{i}(B.times{i}(1):B.times{i}(2)), round(fs*2)); % get whisk events excluding ones happening in window
    [~, ~, ia] = intersect(wOn+1, test); % find index of intersection
    dur = (wOff-wOn)/fs;
    dur = dur(ia);
    chooseidx = (dur <= durmax & dur >= durmin);
    chooseidxanalyze = [chooseidxanalyze, [sum(chooseidx); sum(chooseidx)/numel(chooseidx)]];
    
    tmp = bFrames.selfInitiatedWhiskRightDuring{i}(:,:,:,chooseidx);
%     tmp = bFrames.partnerInitiatedWhiskLeftDuring{i}(:,:,:,chooseidx);
    selfWhiskAvg = cat(4, selfWhiskAvg, mean(tmp- ...
        mean(tmp(:,:,1:29,:),3),4));
    
    
%     partnerWhisk = cat(4, partnerWhisk, bFrames.partnerInitiatedWhiskLeftDuring{i}-...
%         mean(bFrames.partnerInitiatedWhiskLeftDuring{i}(:,:,1:29,:),3));
%     partnerWhisk = cat(4, partnerWhisk, bFrames.partnerInitiatedWhiskRightDuring{i}-...
%         mean(bFrames.partnerInitiatedWhiskRightDuring{i}(:,:,1:29,:),3));

%     partnerWhisk = cat(4, partnerWhisk, mean(bFrames.partnerInitiatedWhiskLeftDuring{i}-...
%         mean(bFrames.partnerInitiatedWhiskLeftDuring{i}(:,:,1:29,:),3),4));
%     partnerWhisk = cat(4, partnerWhisk, mean(bFrames.partnerInitiatedWhiskRightDuring{i}-...
%         mean(bFrames.partnerInitiatedWhiskRightDuring{i}(:,:,1:29,:),3),4));
    
%     selfFL = cat(4, selfFL, mean(bFrames.selfInitiatedFLLeftDuring{i},4));
%     selfFL = cat(4, selfFL, mean(bFrames.selfInitiatedFLRightDuring{i},4));
%     
%     partnerFL = cat(4, partnerFL, mean(bFrames.partnerInitiatedFLLeftDuring{i},4));
%     partnerFL = cat(4, partnerFL, mean(bFrames.partnerInitiatedFLRightDuring{i},4));
%     
% %     partnerWhiskSolo = cat(4, partnerWhiskSolo, mean(bFrames.partnerInitiatedWhiskRightBefore{i}, 4));
%     selfWhiskSolo = cat(4, selfWhiskSolo, bFrames.selfInitiatedWhiskLeftBefore{i} - ...
%         mean(bFrames.selfInitiatedWhiskLeftBefore{i}(:,:,1:29,:),3));


%     selfWhiskSolo = cat(4, selfWhiskSolo, bFrames.selfInitiatedWhiskLeftBefore{i} - ...
%         mean(bFrames.selfInitiatedWhiskLeftBefore{i}(:,:,1:29,:),3));
%     selfWhiskSolo = cat(4, selfWhiskSolo, bFrames.selfInitiatedWhiskLeftAfter{i} - ...
%         mean(bFrames.selfInitiatedWhiskLeftAfter{i}(:,:,1:29,:),3));
    
%     selfWhiskSoloAvg = cat(4, selfWhiskSoloAvg, mean(bFrames.selfInitiatedWhiskLeftBefore{i} - ...
%         mean(bFrames.selfInitiatedWhiskLeftBefore{i}(:,:,1:29,:),3),4));
%     selfWhiskSoloAvg = cat(4, selfWhiskSoloAvg, mean(bFrames.selfInitiatedWhiskLeftAfter{i} - ...
%         mean(bFrames.selfInitiatedWhiskLeftAfter{i}(:,:,1:29,:),3),4));
end

%%


cmbwcount =[];
cmdwcount = [];
ncmbwcount = [];
ncmdwcount = [];

cmdwdur = [];
ncmdwdur = [];

for i = 1:length(CM)
    if i == 15 || i == 27, continue; end
    if CM(i) == 'Y'
        cmbwcount = [cmbwcount, before_wcount(i)];
        
        
        
        % get all whisk events
        [wOn, wOff] = helper.getBehaviourEvents(B.exclusive.whiskLeftDuring{i}(B.times{i}(1):B.times{i}(2)));
        
        % get whisk events excluding ones happening in window
        test = helper.getWhiskEvents(B.exclusive.whiskLeftDuring{i}(B.times{i}(1):B.times{i}(2)), round(fs*2));
        
        % find index of intersection
        [~, ia, ~] = intersect(wOn+1, test);
        
        dur = (wOff-wOn)/fs;
        dur = dur(ia);
        dur = dur(logical(dur <= durmax & dur >= durmin));
        cmdwcount = [cmdwcount, sum(logical(dur <= durmax & dur >= durmin))];
        
        cmdwdur = [cmdwdur, median(dur)];
%         cmdwdur = [cmdwdur, dur];
        
        
        [wOn, wOff] = helper.getBehaviourEvents(B.exclusive.whiskRightDuring{i}(B.times{i}(1):B.times{i}(2))); 
        test = helper.getWhiskEvents(B.exclusive.whiskRightDuring{i}(B.times{i}(1):B.times{i}(2)), round(fs*2));
        [~, ia, ~] = intersect(wOn+1, test);
        
        dur = (wOff-wOn)/fs;
        dur = dur(ia);
        dur = dur(logical(dur <= durmax & dur >= durmin));
        cmdwdur = [cmdwdur, median(dur)];
%         cmdwdur = [cmdwdur, dur];
        
        cmdwcount = [cmdwcount, sum(logical(dur <= durmax & dur >= durmin))];

    else
        ncmbwcount = [ncmbwcount, before_wcount(i)];
        
        
        [wOn, wOff] = helper.getBehaviourEvents(B.exclusive.whiskLeftDuring{i}(B.times{i}(1):B.times{i}(2)));  
        test = helper.getWhiskEvents(B.exclusive.whiskLeftDuring{i}(B.times{i}(1):B.times{i}(2)), round(fs*2));
        [~, ia, ~] = intersect(wOn+1, test);
        dur = (wOff-wOn)/fs;
        dur = dur(ia);
        dur = dur(logical(dur <= durmax & dur >= durmin));
        ncmdwdur = [ncmdwdur, median(dur)];
%         ncmdwdur = [ncmdwdur, dur];
        
        ncmdwcount = [ncmdwcount, sum(logical(dur <= durmax & dur >= durmin))];

        
        [wOn, wOff] = helper.getBehaviourEvents(B.exclusive.whiskRightDuring{i}(B.times{i}(1):B.times{i}(2)));   
        test = helper.getWhiskEvents(B.exclusive.whiskRightDuring{i}(B.times{i}(1):B.times{i}(2)), round(fs*2));
        [~, ia, ~] = intersect(wOn+1, test);
        dur = (wOff-wOn)/fs;
        dur = dur(ia);
        dur = dur(logical(dur <= durmax & dur >= durmin));
        ncmdwdur = [ncmdwdur, median(dur)];
%         ncmdwdur = [ncmdwdur, dur];
        
        ncmdwcount = [ncmdwcount, sum(logical(dur <= durmax & dur >= durmin))];


    end
end

%%

idxtmp = CM;
idxtmp([15, 27]) = [];
% idxtmp = repelem(idxtmp, 2);

% 
% cmmaps = [];
% ncmmaps = [];
% cmmapssolo = [];
% ncmmapssolo = [];
cmmapsAvg = [];
% cmmapssoloAvg = [];
ncmmapsAvg = [];
% ncmmapssoloAvg = [];

counter = 1;
for i = 1:length(idxtmp)
%     N = size(bFrames.partnerInitiatedWhiskLeftDuring{i},4);
%     M = size(bFrames.partnerInitiatedWhiskRightDuring{i},4);

    N = size(bFrames.selfInitiatedWhiskLeftDuring{i},4);
    M = size(bFrames.selfInitiatedWhiskRightDuring{i},4);
%     tmp = selfWhisk(:,:,:,counter:counter+N+M-1);
    
    if idxtmp(i) == 'Y'
        cmmapsAvg = cat(4, cmmapsAvg, selfWhiskAvg(:,:,:,2*i-1:2*i) .* mask);
%         cmmaps = cat(4, cmmaps, partnerWhisk(:,:,:,2*i-1:2*i));
%         cmmapssoloAvg = cat(4, cmmapssoloAvg, selfWhiskSoloAvg(:,:,:,2*i-1:2*i));
        
%         
%         cmmaps = cat(4, cmmaps, selfWhisk(:,:,:,counter:counter+N+M-1));
%         cmmapssolo = cat(4, cmmapssolo, selfWhiskSolo(:,:,:,counter:counter+N-1));
        
%         
%         cmmaps = cat(4, cmmaps, tmp(:,:,:,[chooseidx{2*i-1}, chooseidx{2*i}]));
%         cmmapssolo = cat(4, cmmapssolo, selfWhiskSolo(:,:,:,counter:counter+N-1));
    else
        ncmmapsAvg = cat(4, ncmmapsAvg, selfWhiskAvg(:,:,:,2*i-1:2*i) .* mask);
%         ncmmaps = cat(4, ncmmaps, partnerWhisk(:,:,:,2*i-1:2*i));
%         ncmmapssoloAvg = cat(4,ncmmapssoloAvg, selfWhiskSoloAvg(:,:,:,2*i-1:2*i));
        
%         
%         ncmmaps = cat(4, ncmmaps, selfWhisk(:,:,:,counter:counter+N+M-1));
%         ncmmapssolo = cat(4, ncmmapssolo, selfWhiskSolo(:,:,:,counter:counter+N-1));


%         ncmmaps = cat(4, ncmmaps, tmp(:,:,:,[chooseidx{2*i-1}, chooseidx{2*i}]));
%         ncmmapssolo = cat(4, ncmmapssolo, selfWhiskSolo(:,:,:,counter:counter+N-1));
    end
%     counter = counter + N;
    counter = counter + N + M;
end




%% Trial average together

cmmapsAvg(cmmapsAvg == 0) = nan;
ncmmapsAvg(ncmmapsAvg==0) = nan;

figure,

domDiff = [2, 2, nan, nan, 3, 3, 1, 1, 1, 1, 1, 1, nan, nan, ...
    nan, nan, nan, nan, nan, 1, 3, 3, 1];
subplot(4,3,1),
plot(domDiff, corrs.during(1:23), 'x'), axis([0 4 -0.2 0.8]) 
xlabel('Rank Difference')
ylabel('Inter-brain correlation')

[r, p] = corrcoef(domDiff', corrs.during(1:23)', 'Rows', 'complete');
title(['r=',num2str(r(2,1)),' p=',num2str(p(2,1))])


subplot(4,3,2),
% pp = corrs.during(1:23);
% CMpp = CM(1:23);
% helper.uboxplot(pp(CMpp=='Y'), pp(CMpp=='N'));
helper.uboxplot(corrs.during(CM=='Y'), corrs.during(CM=='N'));
xticklabels({'CM', 'NCM'});
ylabel('Inter-brain correlation')


if helper.isnormal(corrs.during(CM=='Y')) && helper.isnormal(corrs.during(CM=='N'))
    [h,p] = ttest2(corrs.during(CM=='Y'), corrs.during(CM=='N'));
else
    p = ranksum(corrs.during(CM=='Y'), corrs.during(CM=='N'));
end
title(p)

subplot(4,3,4:6)

montage = [helper.makeMontage(nanmean(cmmapsAvg,4),fs); 
    helper.makeMontage(nanmean(ncmmapsAvg,4), fs)];
imagesc(montage);
colorbar; caxis([-0.1 0.4]); colormap jet;
xticks(64:128:size(montage,2))
yticks(64:128:size(montage,1))
xticklabels({'-1', '-0.78', '-0.56', '-0.33', '-0.11', '0.11', '0.33', '0.56', '0.78', '1'})
yticklabels({'CM', 'NCM'})
xlabel('Time (s)')


cmts = squeeze(nanmedian(nanmedian(cmmapsAvg,2),1));
ncmts = squeeze(nanmedian(nanmedian(ncmmapsAvg,2),1));


subplot(4,3,7), 
plot(xt(cmts,fs,1)-1, cmts, 'color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
hold on, plot(xt(cmts,fs,1)-1, mean(cmts,2), 'k', 'LineWidth', 2)
xlabel('Time (s)'), ylabel('\DeltaF/F_0 (z-score)'), axis([-1 1 -0.5 1])

subplot(4,3,8), 
plot(xt(ncmts,fs,1)-1, ncmts, 'color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
hold on, plot(xt(ncmts,fs,1)-1, mean(ncmts,2), 'k', 'LineWidth', 2)
xlabel('Time (s)'), ylabel('\DeltaF/F_0 (z-score)'), axis([-1 1 -0.5 1])

subplot(4,3,9)
p = peak2peak(cmts(30:end,:));
q = peak2peak(ncmts(30:end,:));
helper.uboxplot(p',q')
title(['p=',num2str(ranksum(p,q))]),
ylabel('mean post-whisk \DeltaF/F_0 (z-score)')
xticklabels({'CM', 'NCM'})

subplot(4,3,10)
helper.uboxplot(cmdwcount', ncmdwcount')
ylabel('# Whisking Bouts')
xticklabels({'CM', 'NCM'})
[h, p] = ttest2(cmdwcount, ncmdwcount);
title(['p=',num2str(p)])


subplot(4,3,11)
helper.uboxplot(cmdwdur', ncmdwdur')
ylabel('Whisk Duration (s)')
xticklabels({'CM', 'NCM'})
[h, p] = ttest2(cmdwdur,ncmdwdur);
title(['p=',num2str(p)])

%% Trial average solo

cmmapssoloAvg(cmmapssoloAvg == 0) = nan;
ncmmapssoloAvg(ncmmapssoloAvg==0) = nan;

cmts = squeeze(nanmedian(nanmedian(cmmapssoloAvg,2),1));
ncmts = squeeze(nanmedian(nanmedian(ncmmapssoloAvg,2),1));


[~,cmsort] = sort(mean(cmts(30:end,:)),2);
[~,ncmsort] = sort(mean(ncmts(30:end,:)),2);

figure, subplot(4,2,[1,3,5]), imagesc(cmts(:,cmsort)'),
% colorbar,
caxis([-1 1]), colormap(helper.redblue)
title('CM'), ylabel('Trial-averaged Whisk Events'), xticks([])
subplot(4,2,[2,4,6]), imagesc(ncmts(:,ncmsort)'), 
% colorbar, 
caxis([-1 1]), colormap(helper.redblue)
title('NCM'), xticks([])

subplot(4,2,7), 
plot(xt(cmts,fs,1)-1, cmts, 'color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.1)
hold on, plot(xt(cmts,fs,1)-1, mean(cmts,2), 'k', 'LineWidth', 1.5)
xlabel('time (s)'), ylabel('\DeltaF/F_0 (z-score)'), axis([-1 1 -1 1])

subplot(4,2,8), 
plot(xt(ncmts,fs,1)-1, ncmts, 'color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.1)
hold on, plot(xt(ncmts,fs,1)-1, mean(ncmts,2), 'k', 'LineWidth', 1.5)
xlabel('time (s)'), ylabel('\DeltaF/F_0 (z-score)'), axis([-1 1 -1 1])



