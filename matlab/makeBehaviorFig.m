%% behaviour figure


load('behavior_data.mat')


fs = 28.815;
sum_whisks = [];
sum_fls = [];
overlap = [];
minute = round(60*fs);
for i = 1:length(B.whisk.left)
    if ~isempty(B.whisk.left{i}) && i~=27
        i1 = B.times(i,1);
        i2 = B.times(i,2);
        t1 = B.times(i,3);
        t2 = B.times(i,4);
        
        total_time = numel(1:t1) + numel(t2:length(B.whisk.left{i}));
        
        sum_whisks = [sum_whisks; 100*sum(B.whisk.left{i}([1:t1,t2:end]))/total_time];
%         sum_whisks = [sum_whisks; sum(B.whisk.left{i}(1:t1))/total_time];

        sum_fls = [sum_fls; 100*sum(B.FL.left{i}([1:t1,t2:end]))/total_time];
%         sum_fls = [sum_fls; sum(B.FL.left{i}(1:t1))/total_time];
        
        tmp = B.whisk.left{i}([1:t1,t2:end]) + B.FL.left{i}([1:t1,t2:end]);
%         tmp = B.whisk.left{i}(1:t1) + B.FL.left{i}(1:t1);
        overlap = [overlap; 100*sum(tmp>1)/total_time];
    end
end

w_separate = sum_whisks;
fl_separate = sum_fls;

figure, 
idx = 4;
right_idx = zeros(size(B.whisk.right{idx}));
right_idx(i1:i2) = 1;
brightw = B.whisk.right{idx}.*right_idx;
brightf =  B.FL.right{idx}.*right_idx;
dvec = ones(size(B.whisk.right{idx}));
dvec(i1:i2) = 0.2;
dvec(t1:i1) = linspace(1,0.2,numel(t1:i1));
dvec(i2:t2) = linspace(0.2,1,numel(i2:t2));
dvec = dvec(round(t1-90*fs):round(t2+90*fs));

subplot(4,3,7:9),  hold on, axis tight

plot(xt(dvec,fs),dvec,'k')
helper.makeEthogram(B.whisk.left{idx}(round(t1-90*fs):round(t2+90*fs)),'magenta',fs,[-1 0])
helper.makeEthogram(B.FL.left{idx}(round(t1-90*fs):round(t2+90*fs)),'cyan',fs,[-1 0])
helper.makeEthogram(brightw(round(t1-90*fs):round(t2+90*fs)),'magenta',fs,[-2 -1])
helper.makeEthogram(brightf(round(t1-90*fs):round(t2+90*fs)),'cyan',fs,[-2 -1])
% line([0 30],[-1.5 -1.5], 'Color', [0 0 0])


subplot(4,3,3), axis equal, axis off, title('Solo')
A = [mean(sum_whisks) mean(sum_fls)]; I = mean(overlap);
venn(A,I,'FaceColor',{'m','c'},'FaceAlpha',{0.6,0.6},'EdgeColor','None')
hold on, text(-4,0,[num2str(helper.roundk(A(1),1),3),"+/-",num2str(helper.roundk(std(sum_whisks),1),2)])
text(4,0,[num2str(helper.roundk(A(2),1),3),"+/-",num2str(helper.roundk(std(sum_fls),1),2)])
text(0,0,[num2str(helper.roundk(I,1),3),"+/-",num2str(helper.roundk(std(overlap),1),2)])



sum_whisks = [];
sum_fls = [];
overlap = [];
for i = 1:length(B.whisk.left)
    if ~isempty(B.whisk.left{i}) && i~=27
        i1 = B.times(i,1);
        i2 = B.times(i,2);
        t1 = B.times(i,3);
        t2 = B.times(i,4);
        
        total_time = numel(i1:i2);
        
        sum_whisks = [sum_whisks; 100*sum(B.whisk.left{i}(i1:i2))/total_time];
%         sum_whisks = [sum_whisks; 100*sum(B.whisk.right{i}(i1:i2))/total_time];
        sum_fls = [sum_fls; 100*sum(B.FL.left{i}(i1:i2))/total_time];        
%         sum_fls = [sum_fls; 100*sum(B.FL.right{i}(i1:i2))/total_time];
        tmp = B.whisk.left{i}(i1:i2) + B.FL.left{i}(i1:i2);
        overlap = [overlap; 100*sum(tmp>1)/total_time];
%         tmp = B.whisk.right{i}(i1:i2) + B.FL.right{i}(i1:i2);
%         overlap = [overlap; 100*sum(tmp>1)/total_time];
    end
end

w_together = sum_whisks;
fl_together = sum_fls;

if helper.isnormal(w_together) && helper.isnormal(w_separate)
    [~,p2c] = ttest(w_separate, w_together)
end

subplot(4,3,6), axis equal, axis off, title('Social')
A = [mean(sum_whisks) mean(sum_fls)]; I = mean(overlap);
venn(A,I,'FaceColor',{'m','c'},'FaceAlpha',{0.6,0.6},'EdgeColor','None')
hold on, text(-4,0,[num2str(helper.roundk(A(1),1),3),"+/-",num2str(helper.roundk(std(sum_whisks),1),2)])
text(4,0,[num2str(helper.roundk(A(2),1),3),"+/-",num2str(helper.roundk(std(sum_fls),1),2)])
text(0,0,[num2str(helper.roundk(I,1),3),"+/-",num2str(helper.roundk(std(overlap),1),2)])



tog_ji = [];
tog_ji_rr = [];
sep_ji = [];
for i = 1:35
    if i~=15 && i~=27
%         b_left = (B.whisk.left{i}(i1:i2) + B.FL.left{i}(i1:i2))>0;
%         b_right = (B.whisk.right{i}(i1:i2) + B.FL.right{i}(i1:i2))>0;

        b_left = (B.whisk.left{i}(i1:i1+minute) + B.FL.left{i}(i1:i1+minute))>0;
        b_right = (B.whisk.right{i}(i1:i1+minute) + B.FL.right{i}(i1:i1+minute))>0;
        
        tog_ji_rr = [tog_ji_rr; 1-pdist(double([b_left; b_right]), 'jaccard')];
        tog_ji = [tog_ji; 1-pdist(double([b_left; b_right])+1, 'jaccard')];
        
        b_before = (B.whisk.left{i}(t1-minute:t1) + B.FL.left{i}(t1-minute:t1))>0;
        b_after = (B.whisk.left{i}(t2:t2+minute) + B.FL.left{i}(t2:t2+minute))>0;
        sep_ji = [sep_ji; 1-pdist(double([b_before; b_after]), 'jaccard')];
    end
end

subplot(4,3,12)
boxplot([tog_ji_rr, sep_ji]), ylabel('Intersection Over Union'), axis([0.25 2.75 0 0.5])
xticklabels({'Social', 'Solo'});

if helper.isnormal(tog_ji_rr) && helper.isnormal(sep_ji)
    [~,p2e] = ttest(tog_ji_rr, sep_ji)
else
    [~,p2e] = signrank(tog_ji_rr, sep_ji);
end



fs = 28.815;
xx = [];
yy = [];
zz = [];
qq=[];
minute = round(60*fs);
maxlags = round(20*fs);
for i = 1:35
    if i == 15 || i == 27, continue; end
    
    i1 = B.times(i,1);
    t1 = B.times(i,3);
    t2 = B.times(i,4);


%     [lags,tmp]=helper.CXCORR(single(B.whisk.left{i}(B.times(i,1):B.times(i,1)+minute)), ...
%         single(B.whisk.right{i}(B.times(i,1):B.times(i,1)+minute)));
%     tmp = circshift(tmp, round(numel(tmp)/2));
%     xx = cat(1, xx, tmp);

    l = ((B.FL.left{i}(i1:i1+minute)) + (B.whisk.left{i}(i1:i1+minute)))>0;
    r = ((B.FL.right{i}(i1:i1+minute)) + (B.whisk.right{i}(i1:i1+minute)))>0;

    
    [lags,tmp]=helper.CXCORR(single(l),single(r));
    
    ctr = round(numel(tmp)/2);
    tmp = circshift(tmp, ctr);
    xx = cat(1, xx, tmp);
%     
    
%     [~, tmp] = helper.CXCORR(single(B.whisk.left{i}(B.times(i,1)-minute:B.times(i,1))), ...
%     single(B.whisk.right{i}(B.times(i,1):B.times(i,1)+minute)));
%     tmp = circshift(tmp, round(numel(tmp)/2));
%     yy = cat(1,yy, tmp);
%     
%     
%     [~,tmp] = helper.CXCORR(single(B.whisk.left{i}(B.times(i,4):B.times(i,4)+minute)), ...
%         single(B.whisk.right{i}(B.times(i,1):B.times(i,1)+minute)));
%     tmp = circshift(tmp, round(numel(tmp)/2));
%     zz = cat(1,zz, tmp);

    
%     [~, tmp] = helper.CXCORR(single(B.whisk.left{i}(B.times(i,3)-minute:B.times(i,3))), ...
%         single(B.whisk.left{i}(B.times(i,4):B.times(i,4)+minute)));
%     tmp = circshift(tmp, round(numel(tmp)/2));
%     qq = cat(1, qq, tmp);

    l = ((B.FL.left{i}(t1-minute:t1)) + (B.whisk.left{i}(t1-minute:t1)))>0;
    r = ((B.FL.left{i}(t2:t2+minute)) + (B.whisk.left{i}(t2:t2+minute)))>0;

    [~, tmp] = helper.CXCORR(single(l), single(r));
    tmp = circshift(tmp, ctr);
    qq = cat(1, qq, tmp);
end

if helper.isnormal(qq(:, ctr)) && helper.isnormal(xx(:,ctr))
    [~,p2d] = ttest(qq(:,ctr), xx(:,ctr))
end


lags = lags-lags(end)/2;

subplot(4,3,10:11), shadedErrorBar(lags/fs,mean(xx),std(xx)/sqrt(size(xx,1)),'k',1)
hold on
% shadedErrorBar(lags/fs,mean(yy),std(yy)/sqrt(size(yy,1)),'b',1)
% shadedErrorBar(lags/fs,mean(zz),std(zz)/sqrt(size(zz,1)),'r',1)
shadedErrorBar(lags/fs,mean(qq),std(qq)/sqrt(size(qq,1)), 'r', 1), axis tight