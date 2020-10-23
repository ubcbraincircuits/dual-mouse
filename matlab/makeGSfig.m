trial_num = 1;
trans_dur = 27.5;
together_dur = 120;
separate_dur = 90;
N = numel(openFiles);
during_period = round((trans_dur+separate_dur)*fs):round((trans_dur+separate_dur+together_dur)*fs);


figure, subplot(3, 4, 1)
tmpMaps = untile(regMaps, [128, 128]);
imagesc([tmpMaps(:,:,trial_num).*mask; tmpMaps(:,:,N+trial_num).*mask]);
colormap gray
yticks([64, 192])
% yticklabels({'Stationary Mouse', 'Moving Mouse'})
xticks('')
c = colorbar;
c.Visible = 'off';

helper.freezeColors()
minute = round(90*fs);

% % original
lgs = zscore(GS.left{trial_num});
rgs = zscore(GS.right{trial_num});

% lgs = mean(cell2mat(GS.left'),2);
% rgs = mean(cell2mat(GS.right'),2);

% lgs = LGS{trial_num}(t1-minute:t2+minute);
% rgs = RGS{trial_num}(t1-minute:t2+minute);
rgs2 = rgs+min(lgs)-max(rgs);

time = xt(rgs2, fs);
% LM = bsxfun(@times, LM, cast(mask, 'like', LM));
% RM = bsxfun(@times, RM, cast(mask, 'like', RM));

% lgs = LGS{trial_num}(i1:i2);
% rgs = RGS{trial_num}(i1:i2);

subplot(3, 4, 2:3), plot(time, lgs, 'k'), 
hold on, plot(time, rgs2, 'k')
yL = get(gca,'YLim');

patch([separate_dur, separate_dur + trans_dur, ...
    separate_dur + trans_dur, separate_dur], ...
    [yL(1), yL(1), yL(2), yL(2)], ...
    'yellow','FaceAlpha', 0.25, 'EdgeColor', 'None')

patch([separate_dur+trans_dur+together_dur, separate_dur + 2*trans_dur + together_dur, ...
    separate_dur + 2*trans_dur + together_dur, separate_dur + trans_dur+together_dur], ...
    [yL(1), yL(1), yL(2), yL(2)], ...
    'yellow','FaceAlpha', 0.25, 'EdgeColor', 'None')

patch([trans_dur + separate_dur, trans_dur + together_dur + 90, ...
    trans_dur + together_dur + 90, trans_dur + 90], ...
    [yL(1), yL(1), yL(2), yL(2)],'green','FaceAlpha', 0.25, 'EdgeColor', 'None')
set(gca,'children',flipud(get(gca,'children')))
% line([(minute + trans_time)/fs, (minute + trans_time)/fs],yL, 'Color', 'k','LineStyle','--')
% line([(length(rgs)-minute-trans_time)/fs, (length(rgs)-minute-trans_time)/fs],yL, 'Color', 'k','LineStyle','--')
% line([(minute + trans_time)/fs, (length(rgs)-minute-trans_time)/fs],[yL(1), yL(1)], 'Color', 'k','LineStyle','--')
% line([(minute + trans_time)/fs, (length(rgs)-minute-trans_time)/fs],[yL(2), yL(2)], 'Color', 'k','LineStyle','--')

line([-5 25], [min(rgs2)-0.5, min(rgs2)-0.5], 'Color', [0 0 0], 'LineWidth', 2)
line([-5 -5], [min(rgs2)-0.5 min(rgs2)+1.5], 'Color', [0 0 0], 'LineWidth', 2)
text(1, double(min(rgs2)-1.5), '30 s')
text(-40, double(min(rgs2)), {['   2 \sigma'], ['\DeltaF/F_0']})

axis([-5 time(end), min(rgs2)-0.5, max(lgs)])
axis off



subplot(3,4,5:8), hold on
% helper.annotateBehaviour(socialData, B, trial_num);
gs_left_social = GS.left{trial_num}(round((90+trans_dur)*fs):end-round((90+trans_dur)*fs)); 
gs_right_social = GS.right{trial_num}(round((90+trans_dur)*fs):end-round((90+trans_dur)*fs)); 
gs_right_social = gs_right_social + min(gs_left_social) - max(gs_right_social);

helper.makeEthogram(B.whisk.left{trial_num}(B.times{trial_num}(1):B.times{trial_num}(2)), ...
    'magenta', fs, [min(gs_left_social), max(gs_left_social)], 0.5);
helper.makeEthogram(B.FL.left{trial_num}(B.times{trial_num}(1):B.times{trial_num}(2)), ...
    'cyan', fs, [min(gs_left_social), max(gs_left_social)], 0.5);
plot(xt(gs_left_social,fs), gs_left_social, 'k'), 

helper.makeEthogram(B.whisk.right{trial_num}(B.times{trial_num}(1):B.times{trial_num}(2)), ...
    'magenta', fs, [min(gs_right_social), max(gs_right_social)], 0.5);
helper.makeEthogram(B.FL.right{trial_num}(B.times{trial_num}(1):B.times{trial_num}(2)), ...
    'cyan', fs, [min(gs_right_social), max(gs_right_social)], 0.5);
plot(xt(gs_right_social,fs), gs_right_social, 'k'), 

time = xt(gs_left_social, fs);
yL = get(gca, 'YLim');
xL = get(gca, 'XLim');
line([-2 8], [yL(1)-1, yL(1)-1], 'Color', [0 0 0], 'LineWidth', 2)
line([-2 -2], [yL(1)-1 yL(1)+1], 'Color', [0 0 0], 'LineWidth', 2)
text(1, yL(1)-1.5, '10 s')
text(-13, yL(1), {['   2 \sigma'], ['\DeltaF/F_0']})

axis([-2 time(end), yL(1)-1, yL(2)])
axis off



subplot(3,4,4)
duringCshuffle = [];
count = 1;
for i = 1:N
    for j = 1:N
        if i<j
            tmp = corrcoef(GS.left{i}(during_period), GS.right{j}(during_period));
            duringCshuffle(count) = tmp(2,1);            
            count = count +1;
        end
    end
end
test = nan(10000,4);
test(1:length(corrs.before),1) = corrs.before;
test(1:length(corrs.during),2) = corrs.during;
test(1:length(corrs.after), 3) = corrs.after;
% test(1:23,1) = corrs.before(1:23);
% test(1:23,2) = corrs.during(1:23);
% test(1:23, 3) = corrs.after(1:23);
test(1:length(duringCshuffle), 4) = duringCshuffle;
boxplot(test)
hold on
yL = get(gca, 'YLim');
line([3.5 3.5], yL, 'Color', [0 0 0], 'LineStyle', ':')
patch([1.5, 2.5, 2.5, 1.5], [yL(1), yL(1), yL(2), yL(2)], 'green', 'FaceAlpha', 0.25, 'EdgeColor', 'None')
text(1.9,yL(2), '*', 'FontSize', 14)
set(gca,'children',flipud(get(gca,'children')))

xticklabels({'Separate1','Together','Seperate2', 'Together Shuffled'})
xtickangle(30)
ylabel('PCC_{Interbrain}')
xlabel(' ')


%%% brain-behavior correlation
s90 = round(90*fs);
b_dur = [];
b_dff = [];

% crop behavior traces to 90s before and after translation
for i = 1:35
    if i == 15 || i == 27, continue; end
    i1 = B.times{i}(1);
    i2 = B.times{i}(2);
    gs_left_social = GS.left{i}(round((90+trans_dur)*fs)-1:end-round((90+trans_dur)*fs)+1);
    gs_right_social = GS.right{i}(round((90+trans_dur)*fs)-1:end-round((90+trans_dur)*fs)+1);

%     b_left = single((B.whisk.left{i}(t1-s90:t2+s90) + B.FL.left{i}(t1-s90:t2+s90))>0);
%     b_right = (B.whisk.left{i}(i1:i2) + B.FL.left{i}(i1:i2))>0;
        b_left = (B.whisk.left{i}(i1:i2) + B.FL.left{i}(i1:i2))>0;
        b_right = (B.whisk.right{i}(i1:i2) + B.FL.right{i}(i1:i2))>0;
    
    
    [b_on, b_off] = helper.getBehaviourEvents(b_left);
    assert(numel(b_on) == numel(b_off), 'assertion failed');
    for j = 1:length(b_on)
        if b_off(j) > b_on(j)
            b_dur = [b_dur; (b_off(j)-b_on(j))./fs];
%             b_dff = [b_dff; max(GS.left{i}(b_on(j):b_off(j)))];
            b_dff = [b_dff; mean(gs_left_social(b_on(j):b_off(j)))];
        else
            try
                b_dur = [b_dur; (b_off(j+1)-b_on(j))./fs];
%                 b_dff = [b_dff; max(GS.left{i}(b_on(j):b_off(j+1)))];
                b_dff = [b_dff; mean(gs_left_social(b_on(j):b_off(j+1)))];
            catch
%                 b_dur = [b_dur; (numel(GS.left{i})-b_on(j))./fs];
                b_dur = [b_dur; (numel(gs_left_social)-b_on(j))./fs];
%                 b_dff = [b_dff; max(GS.left{i}(b_on(j):end))];
                b_dff = [b_dff; mean(gs_left_social(b_on(j):end))];
            end
        end
        assert(numel(b_dur) == numel(b_dff), [num2str(i),'   ', num2str(j)])
    end
    
    [b_on, b_off] = helper.getBehaviourEvents(b_right);
    assert(numel(b_on) == numel(b_off), 'assertion failed');
    for j = 1:length(b_on)
        if b_off(j) > b_on(j)
            b_dur = [b_dur; (b_off(j)-b_on(j))./fs];
%             b_dff = [b_dff; max(GS.left{i}(b_on(j):b_off(j)))];
            b_dff = [b_dff; mean(gs_right_social(b_on(j):b_off(j)))];
        else
            try
                b_dur = [b_dur; (b_off(j+1)-b_on(j))./fs];
%                 b_dff = [b_dff; max(GS.left{i}(b_on(j):b_off(j+1)))];
                b_dff = [b_dff; mean(gs_right_social(b_on(j):b_off(j+1)))];
            catch
%                 b_dur = [b_dur; (numel(GS.left{i})-b_on(j))./fs];
                b_dur = [b_dur; (numel(gs_right_social)-b_on(j))./fs];
%                 b_dff = [b_dff; max(GS.left{i}(b_on(j):end))];
                b_dff = [b_dff; mean(gs_right_social(b_on(j):end))];
            end
        end
        assert(numel(b_dur) == numel(b_dff), [num2str(i),'   ', num2str(j)])
    end
    
    
end

b_dur = log10(b_dur);
bias = ones(size(b_dur));

% calculate coefficient
w = [b_dur, bias]\b_dff;
[r,p] = corrcoef(b_dur, b_dff);

% plot results
subplot(3,4,9)
scatter(b_dur, b_dff, 10, '.'); hold on
plot(b_dur, w(2) + w(1)*b_dur, 'k', 'LineWidth', 2)
ylabel('Mean \DeltaF/F_0 (\sigma)')
xlabel('Behavior duration ( log_{10}( seconds) )')
title(['r=',num2str(r(2,1),2),' p=',num2str(p(2,1),2)])


% subplot(3,4,9:10),
% maxlag = round(60*fs);
% cross = nan(N, 3459);
% for i = 1:N
%     [cross(i,:), lags] = xcorr(zscore(GS.left{i}), zscore(GS.right{i}), maxlag);
% end
% 
% shadedErrorBar(lags./fs, mean(cross), std(cross,[],1)./sqrt(size(cross,1)), '-k', 1)
% hold on
% 
% ylabel('Cross-correlation (AU)')
% xlabel('Time lag (s)')
% yL = get(gca, 'yLim');
% line([0 0], yL, 'color', [0 0 0], 'LineStyle', '--');
% axis([min(lags./fs), max(lags./fs), yL])

% 
% subplot(4,5,11:12)
% num_iter = 100;
% surrC = helper.getSurrogateCorr(GS, num_iter);
% histogram(surrC(:), 50, 'FaceColor', 'k'); hold on
% yL = get(gca,'YLim');
% 
% for i = 1:N
%     line([corrs.during(i) corrs.during(i)], yL, 'color', [0 1 0 0.5], 'LineWidth', 0.5)
% end
% line([median(corrs.during) median(corrs.during)], yL, 'color', [0 0 0], 'LineWidth', 3)
% text(median(corrs.during)-0.05, yL(2)+10, [num2str(invprctile(surrC(:), median(corrs.during)), 3),'%'])
% axis([-0.4 0.6 yL])
% xlabel('PCC_{Surrogate}')
% ylabel('Count')

subplot(3,4,11:12)
window = 45;
noverlap = window/2;

x = catcell(2,GS.left);
y = catcell(2, GS.right);
movingwin = [window noverlap];
params.Fs = fs;
params.trialave=1;
params.fpass = [0 1];
params.tapers = [5 9];
[C,phi,S12,S1,S2,t,f]=cohgramc(zscore(x),zscore(y),movingwin,params);
pcolor(t,f,C'), shading interp, colormap jet
xlabel('Time (s)')
ylabel('Frequency (Hz)')
c = colorbar;
c.Label.String = 'Coherence_{Interbrain}';
hold on
line([separate_dur + trans_dur, separate_dur + trans_dur], [params.fpass(1), params.fpass(2)], 'Color', [1 1 1], 'LineStyle', '--', 'LineWidth', 2)
line([separate_dur + trans_dur + together_dur, separate_dur + trans_dur + together_dur], [params.fpass(1), params.fpass(2)], 'Color', [1 1 1], 'LineStyle', '--', 'LineWidth', 2)
title([num2str(window),'s window'])

pp = gcf;
pp.Renderer = 'Painters';
