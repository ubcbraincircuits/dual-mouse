f = fieldnames(bFrames);
% count = 0;
X = cell(1,35);
Y = cell(1,35);
for i = 1:35
    if i == 15 || i == 27, continue; end
%     for j = 1:length(bFrames.(f{i}))
% %         bFrames.(f{i}){j} = single(bFrames.(f{i}){j});
% %         count = count + numel(bFrames.(f{i}){j});
%     end

    X{i} = bFrames.selfInitiatedWhiskLeftBefore{i};
    Y{i} = bFrames.selfInitiatedWhiskLeftAfter{i};
end


%%




