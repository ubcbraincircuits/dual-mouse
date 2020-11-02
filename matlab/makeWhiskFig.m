%% whisk triggered figure



mouseAvg = [];
window = round(fs*2);
durmax = 0.5;
durmin = 0;
chooseidxanalyze = [];
for i = 1:length(CM)
    if i == 15 || i == 27, continue; end
    i1 = B.times(i,1);
    i2 = B.times(i,2);
    t1 = B.times(i,3);
    t2 = B.times(i,4);
    
    before_wcount(i) = size(X{i}, 4);
    before_w_rate(i) = before_wcount(i)/(t1/fs);
    during_wcount(i) = size(Y{i},4);
    during_w_rate(i) = during_wcount(i)/((i2-i1)/fs);

    
    [wOn, wOff] = helper.getBehaviourEvents(B.exclusive.whiskLeftDuring{i}, window); % get all whisk events
    dur = (wOff-wOn)/fs;
    chooseidx = (dur <= durmax & dur >= durmin);
    chooseidxanalyze = [chooseidxanalyze, [sum(chooseidx); sum(chooseidx)/numel(chooseidx)]];
    
    
    tmp = X{i}(:,:,:,chooseidx);
    mouseAvg = cat(4, mouseAvg, mean(tmp - ...
        mean(tmp(:,:,1:29,:),3),4));
    
    [wOn, wOff] = helper.getBehaviourEvents(B.exclusive.whiskRightDuring{i}, window); % get all whisk events
    dur = (wOff-wOn)/fs;

    chooseidx = (dur <= durmax & dur >= durmin);
    chooseidxanalyze = [chooseidxanalyze, [sum(chooseidx); sum(chooseidx)/numel(chooseidx)]];
    tmp = Y{i}(:,:,:,chooseidx);
    mouseAvg = cat(4, mouseAvg, mean(tmp- ...
        mean(tmp(:,:,1:29,:),3),4));
    
    

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
        [wOn, wOff] = helper.getBehaviourEvents(B.exclusive.whiskLeftDuring{i}, window);
        dur = (wOff-wOn)/fs;
        dur = dur(logical(dur <= durmax & dur >= durmin));
        cmdwcount = [cmdwcount, sum(logical(dur <= durmax & dur >= durmin))];
        
        cmdwdur = [cmdwdur, median(dur)];
        
        
        [wOn, wOff] = helper.getBehaviourEvents(B.exclusive.whiskRightDuring{i}, window); 
        
        dur = (wOff-wOn)/fs;
        dur = dur(logical(dur <= durmax & dur >= durmin));
        cmdwdur = [cmdwdur, median(dur)];
        
        cmdwcount = [cmdwcount, sum(logical(dur <= durmax & dur >= durmin))];

    else
        ncmbwcount = [ncmbwcount, before_wcount(i)];
        
        
        [wOn, wOff] = helper.getBehaviourEvents(B.exclusive.whiskLeftDuring{i}, window);  
        dur = (wOff-wOn)/fs;
        dur = dur(logical(dur <= durmax & dur >= durmin));
        ncmdwdur = [ncmdwdur, median(dur)];      
        ncmdwcount = [ncmdwcount, sum(logical(dur <= durmax & dur >= durmin))];

        
        [wOn, wOff] = helper.getBehaviourEvents(B.exclusive.whiskRightDuring{i}, window);   
        dur = (wOff-wOn)/fs;
        dur = dur(logical(dur <= durmax & dur >= durmin));
        ncmdwdur = [ncmdwdur, median(dur)];
        ncmdwcount = [ncmdwcount, sum(logical(dur <= durmax & dur >= durmin))];


    end
end

%%

idxtmp = CM;
idxtmp([15, 27]) = [];

cmmapsAvg = [];
ncmmapsAvg = [];

for i = 1:length(idxtmp)

    if idxtmp(i) == 'Y'
        cmmapsAvg = cat(4, cmmapsAvg, mouseAvg(:,:,:,2*i-1:2*i) .* mask.open);

    else
        ncmmapsAvg = cat(4, ncmmapsAvg, mouseAvg(:,:,:,2*i-1:2*i) .* mask.open);

    end

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
helper.uboxplot(corrs.during(CM=='Y'), corrs.during(CM=='N'));
xticklabels({'CM', 'NCM'});
ylabel('Inter-brain correlation')


if helper.isnormal(corrs.during(CM=='Y')) && helper.isnormal(corrs.during(CM=='N'))
    [~,p] = ttest2(corrs.during(CM=='Y'), corrs.during(CM=='N'));
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
p = mean(cmts(30:end,:));
q = mean(ncmts(30:end,:));
helper.uboxplot(p',q')
title(['p=',num2str(ranksum(p,q))]),
ylabel('mean post-whisk \DeltaF/F_0 (z-score)')
xticklabels({'CM', 'NCM'})

subplot(4,3,10)
helper.uboxplot(cmdwcount', ncmdwcount')
ylabel('# Whisking Bouts')
xticklabels({'CM', 'NCM'})
[~, p] = ttest2(cmdwcount, ncmdwcount);
title(['p=',num2str(p)])


subplot(4,3,11)
helper.uboxplot(cmdwdur', ncmdwdur')
ylabel('Whisk Duration (s)')
xticklabels({'CM', 'NCM'})
[h, p] = ttest2(cmdwdur,ncmdwdur);
title(['p=',num2str(p)])

