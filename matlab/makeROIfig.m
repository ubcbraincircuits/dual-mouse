%% roi_figure


R = numel(ROIs);
refIdx = 51;
rbao = rMat.open.before(1:R,1:R,:);
rdao = rMat.open.during(1:R,1:R,:);
raao = rMat.open.after(1:R,1:R,:);

rbro = rMat.open.before(R+1:end,1:R,:);
rdro = rMat.open.during(R+1:end,1:R,:);
raro = rMat.open.after(R+1:end,1:R,:);

% mesh inter
rbrm = rMat.mesh.before(R+1:end, 1:R, :);
rdrm = rMat.mesh.during(R+1:end, 1:R, :);
rarm = rMat.mesh.after(R+1:end, 1:R, :);

% opaque inter
rbrq = rMat.opaque.before(R+1:end, 1:R, :);
rdrq = rMat.opaque.during(R+1:end, 1:R, :);
rarq = rMat.opaque.after(R+1:end, 1:R, :);

figure, 
ms = helper.untile(regMapsOpen, [128 128]);
subplot(5,4,[1,5]), helper.show_coords(ms(:,:,refIdx), ROIs, CLopen, CRopen),
axis off, helper.freezeColors()


subplot(5,4,6), imagesc(mean(rbao,3)), colormap parula; 
colorbar, caxis([0.4 1])
xticks(1:size(rbao,2)); xticklabels(ROIs); xtickangle(90)
yticks(1:size(rbao,1)); yticklabels(ROIs); 

subplot(5,4,7), imagesc(mean(rdao,3)), %colormap jet; 
colorbar, caxis([0.4 1])
xticks(1:size(rbao,2)); xticklabels(ROIs); xtickangle(90)
yticks(1:size(rbao,1)); yticklabels(ROIs);

subplot(5,4,8), 
r_by_roi = [squeeze(median(median(rbao,2),3)), ...
    squeeze(median(median(rdao,2),3)), ...
    squeeze(median(median(raao,2),3))]';
plot(r_by_roi)
axis([0.75 3.25 0.6 1])
hold on
yL = get(gca, 'YLim');
patch([1.5, 2.5, 2.5, 1.5], [yL(1), yL(1), yL(2), yL(2)], 'green', 'FaceAlpha', 0.25, 'EdgeColor', 'None')
set(gca,'children',flipud(get(gca,'children')))
% xticklabels({'Separate1','Together','Seperate2', 'Together Shuffled'})
xtickangle(30)
ylabel('Avg r x ROI')
xlabel(' ')
% imagesc(mean(raintra,3)), %colormap jet; 
% colorbar, caxis([0.4 1])
% xticks(1:size(rbintra,2)); xticklabels(ROIs); xtickangle(90)
% yticks(1:size(rbintra,1)); yticklabels(ROIs);
 
subplot(5,4,2), imagesc(mean(rbro,3)), %colormap jet; 
colorbar, caxis([0 0.5])
xticks(1:size(rbao,2)); xticklabels(ROIs); xtickangle(90)
yticks(1:size(rbao,1)); yticklabels(ROIs); 

subplot(5,4,3), imagesc(mean(rdro,3)), %colormap jet; 
colorbar, caxis([0 0.5])
xticks(1:size(rbao,2)); xticklabels(ROIs); xtickangle(90)
yticks(1:size(rbao,1)); yticklabels(ROIs);

subplot(5,4,4), 
r_by_roi = [squeeze(median(median(rbro,2),3)), ...
    squeeze(median(median(rdro,2),3)), ...
    squeeze(median(median(raro,2),3))]';
plot(r_by_roi)
hold on
yL = get(gca, 'YLim');
axis([0.75 3.25 yL(1) yL(2)])
patch([1.5, 2.5, 2.5, 1.5], [yL(1), yL(1), yL(2), yL(2)], 'green', 'FaceAlpha', 0.25, 'EdgeColor', 'None')
set(gca,'children',flipud(get(gca,'children')))
% xticklabels({'Separate1','Together','Seperate2', 'Together Shuffled'})
xtickangle(30)
ylabel('Avg r x ROI')
xlabel(' ')
helper.freezeColors();





% selfleft = catcell(4, bFrames.selfInitiatedWhiskLeftDuring([1:14,16:26,28:end]));
% selfright = catcell(4, bFrames.selfInitiatedWhiskRightDuring([1:14,16:26,28:end]));
% self = cat(4, selfleft, selfright);
% clear selfleft selfright
% 
% partnerleft = catcell(4,bFrames.partnerInitiatedWhiskLeftDuring([1:14,16:26,28:end]));
% partnerright = catcell(4,bFrames.partnerInitiatedWhiskRightDuring([1:14,16:26,28:end]));
% partner = cat(4, partnerleft, partnerright);
% clear partnerleft partnerright
% 
% smon = [];
% for i = 1:size(self,4)
%     smon = cat(3, smon, helper.makeMontage(self(:,:,:,i),fs));
% end
% 
% pmon = [];
% for i = 1:size(partner,4)
%     pmon = cat(3, pmon, helper.makeMontage(partner(:,:,:,i),fs));
% end
% 
% smon(smon==0) = nan;
% pmon(pmon==0) = nan;

subplot(5,4,9:12), imagesc(mean(smon,3)), colormap jet; c = colorbar;
axis off, c.Label.String = '\DeltaF/F_0 (\sigma)';
subplot(5,4,13:16), imagesc(mean(pmon,3)), colormap jet; c=colorbar;
axis off, c.Label.String = '\DeltaF/F_0 (\sigma)';

pmo = [squeeze(median(median(rbro,2),1)), ...
    squeeze(median(median(rdro,2),1)), ...
    squeeze(median(median(raro,2),1))];

subplot(5,4,18), plot(mean(pmo), 'k', 'LineWidth', 1)
hold on, plot(pmo', 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5 0.5]),
axis([0.75 3.25 -0.4 0.8])
yL = get(gca, 'YLim');
patch([1.5, 2.5, 2.5, 1.5], [yL(1), yL(1), yL(2), yL(2)], 'green', 'FaceAlpha', 0.25, 'EdgeColor', 'None')
set(gca,'children',flipud(get(gca,'children')))
% xticklabels({'Separate1','Together','Seperate2', 'Together Shuffled'})
xtickangle(30)
ylabel('Avg r')
xlabel(' ')
title('open')


pmm = [squeeze(median(median(rbrm,2),1)), ...
    squeeze(median(median(rdrm,2),1)), ...
    squeeze(median(median(rarm,2),1))];
subplot(5,4,20),hold on, plot(mean(pmm), 'k', 'LineWidth', 1), axis([0.75 3.25 -0.4 0.8])
plot(pmm', 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5 0.5]),
yL = get(gca, 'YLim');
patch([1.5, 2.5, 2.5, 1.5], [yL(1), yL(1), yL(2), yL(2)], 'green', 'FaceAlpha', 0.25, 'EdgeColor', 'None')
set(gca,'children',flipud(get(gca,'children')))
% xticklabels({'Separate1','Together','Seperate2', 'Together Shuffled'})
xtickangle(30)
ylabel('Avg r')
xlabel(' ')
title('mesh')


pmq = [squeeze(median(median(rbrq,2),1)), ...
    squeeze(median(median(rdrq,2),1)), ...
    squeeze(median(median(rarq,2),1))];

subplot(5,4,19), hold on, plot(mean(pmq), 'k', 'LineWidth', 1), axis([0.75 3.25 -0.4 0.8])
plot(pmq', 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5 0.5]), 

title('opaque')
yL = get(gca, 'YLim');
patch([1.5, 2.5, 2.5, 1.5], [yL(1), yL(1), yL(2), yL(2)], 'green', 'FaceAlpha', 0.25, 'EdgeColor', 'None')
set(gca,'children',flipud(get(gca,'children')))
% xticklabels({'Separate1','Together','Seperate2', 'Together Shuffled'})
xtickangle(30)
ylabel('Avg r')
xlabel(' ')