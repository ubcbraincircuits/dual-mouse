num_trials = 35;
whisk_before = nan(num_trials, 1);
whisk_during = nan(num_trials, 1);
whisk_after = nan(num_trials, 1);
FL_before = nan(num_trials, 1);
FL_during = nan(num_trials, 1);
FL_after = nan(num_trials, 1);

xx = [];
yy = [];
qq = [];
zz = [];


for i = 1:length(B.whisk.left)
    if i == 15 || i ==27, continue; end
        i1 = B.times(i,1);
        i2 = B.times(i,2);
        t1 = B.times(i,3);
        t2 = B.times(i,4);
        bOn = helper.getBehaviourEvents(B.whisk.left{i}(t1-minute:t1), 0);
        whisk_before(i) = numel(bOn);
        bOn = helper.getBehaviourEvents(B.whisk.left{i}(i1:i1+minute), 0);
        whisk_during(i) = numel(bOn);
        bOn = helper.getBehaviourEvents(B.whisk.left{i}(t2:t2+minute), 0);
        whisk_after(i) = numel(bOn);
        
        
        bOn = helper.getBehaviourEvents(B.FL.left{i}(t1-minute:t1), 0);
        FL_before(i) = numel(bOn);
        bOn = helper.getBehaviourEvents(B.FL.left{i}(i1:i1+minute), 0);
        FL_during(i) = numel(bOn);
        bOn = helper.getBehaviourEvents(B.FL.left{i}(t2:t2+minute), 0);
        FL_after(i) = numel(bOn);

        
        
        
        l = ((B.FL.left{i}(i1:i1+minute)) + (B.whisk.left{i}(i1:i1+minute)))>0;
        r = ((B.FL.right{i}(i1:i1+minute)) + (B.whisk.right{i}(i1:i1+minute)))>0;
        [lags,tmp]=helper.CXCORR(single(l),single(r));
        ctr = round(numel(tmp)/2);
        tmp = circshift(tmp, ctr);
        xx = cat(1, xx, tmp);
         
        l = ((B.FL.left{i}(t1-minute:t1)) + (B.whisk.left{i}(t1-minute:t1)))>0;
        r = ((B.FL.right{i}(i1:i1+minute)) + (B.whisk.right{i}(i1:i1+minute)))>0;
        [~, tmp] = helper.CXCORR(single(l), single(r));
        tmp = circshift(tmp, round(numel(tmp)/2));
        yy = cat(1,yy, tmp);

        l = ((B.FL.left{i}(t2:t2+minute)) + (B.whisk.left{i}(t2:t2+minute)))>0;
        r = ((B.FL.right{i}(i1:i1+minute)) + (B.whisk.right{i}(i1:i1+minute)))>0;
        [~, tmp] = helper.CXCORR(single(l), single(r));
        tmp = circshift(tmp, round(numel(tmp)/2));
        zz = cat(1,zz, tmp);


        l = ((B.FL.left{i}(t1-minute:t1)) + (B.whisk.left{i}(t1-minute:t1)))>0;
        r = ((B.FL.left{i}(t2:t2+minute)) + (B.whisk.left{i}(t2:t2+minute)))>0;

        [~, tmp] = helper.CXCORR(single(l), single(r));
        tmp = circshift(tmp, ctr);
        qq = cat(1, qq, tmp);

end



figure, 
subplot(2,2,1), boxplot([whisk_before, whisk_during, whisk_after])
xticklabels({'before','during','after'})
ylabel('# whisker movements')
subplot(2,2,2), boxplot([FL_before, FL_during, FL_after])
xticklabels({'before','during','after'})
ylabel('# forelimb movements')
subplot(2,2,3:4)
hold on
plot(lags/fs - 30, mean(xx), 'k')
plot(lags/fs - 30, mean(yy), 'b')
plot(lags/fs - 30, mean(zz), 'r')
plot(lags/fs - 30, mean(qq), 'c')

ylabel('Mean cross-correlation')
xlabel('Time (s)')
legend({'S_{during} vs M_{during}', 'S_{before} vs M_{during}', 'S_{after} vs M_{during}', 'S_{before} vs S_{after}'}, 'location', 'best')


% remove nans
nanidx = isnan(whisk_before);
whisk_before = whisk_before(~nanidx);
whisk_during  = whisk_during(~nanidx);
whisk_after = whisk_after(~nanidx);
t = table(CM(~nanidx)',whisk_before,whisk_during,whisk_after,...
'VariableNames',{'CM','t0','t2','t4'});
Time = [0 2 4]';
rm = fitrm(t,'t0-t4 ~ CM','WithinDesign',Time);
anovatbl = anova(rm)

if helper.isnormal(whisk_before) && ...
helper.isnormal(whisk_during) && ...
helper.isnormal(whisk_after)
    [h1,p] = ttest(whisk_before, whisk_during);
    disp(['together vs before; p=', num2str(p*3)])
    [h2,p] = ttest(whisk_after, whisk_during);
    disp(['together vs after; p=', num2str(p*3)])
    [h3,p] = ttest(whisk_before, whisk_after);
    disp(['before vs after; p=', num2str(p*3)])
end



FL_before = FL_before(~nanidx);
FL_during  = FL_during(~nanidx);
FL_after = FL_after(~nanidx);
t = table(CM(~nanidx)',FL_before,FL_during,FL_after,...
'VariableNames',{'CM','t0','t2','t4'});
Time = [0 2 4]';
rm = fitrm(t,'t0-t4 ~ CM','WithinDesign',Time);
anovatbl = anova(rm)

if helper.isnormal(FL_before) && ...
helper.isnormal(FL_during) && ...
helper.isnormal(FL_after)
    [h1,p] = ttest(FL_before, FL_during);
    disp(['together vs before; p=', num2str(p*3)])
    [h2,p] = ttest(FL_after, FL_during);
    disp(['together vs after; p=', num2str(p*3)])
    [h3,p] = ttest(FL_before, FL_after);
    disp(['before vs after; p=', num2str(p*3)])
end
