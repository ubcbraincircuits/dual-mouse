function varargout = loadData(pathname, filenames, tform, mask, B, analysis, CL, CR)

unfiltPath = 'B:\Social_Outputs\Matlab\social\unfiltered\';
unfiltFiles = helper.getAllFiles(unfiltPath);

% Constants
N = length(filenames);
fs = 28.815;
imgResizeFactor = 1/2;

i1 = B.times(1,1);
t1 = B.times(1,3);
t2 = B.times(1,4);

minute = round(60*fs);

switch analysis       
    case 'global_signal_correlation_interbrain'
        % returns:
        %   socialData  (5D matrix of shape H x W x Frame x Trial x 2)
        %   corrs       (structure of correlation in minute long periods
        %               immediately before translation, immediately after
        %               first translation ends, or immediately after 2nd
        %               translation ends)
        %   GS          (structure that contains global signals for left 
        %               and right mouse from 90s before translation1 start 
        %               to 90s after translation2 end)
        %   B           (structure that contains binary vectors of
        %               behaviour event indexes during the together period)
        
        % pre-allocate
        socialData = zeros(128, 128, 3459, N, 2);
        C.before = zeros(length(filenames),1);
        C.during = zeros(length(filenames),1);
        C.after = zeros(length(filenames),1);
        GS.left = cell(length(filenames),1);
        GS.right = cell(length(filenames),1);

        % check minute in each trial phase
        s90 = round(90*fs);
        for i = 1:N
            tic

            % load corrected brain data
            leftGreenFile = unfiltFiles{contains(unfiltFiles, filenames{i}(1:11)) & ...
                contains(unfiltFiles, 'GREEN') & ...
                contains(unfiltFiles, 'LEFT')};
            rightGreenFile = unfiltFiles{contains(unfiltFiles, filenames{i}(1:11)) & ...
                contains(unfiltFiles, 'GREEN') & ...
                contains(unfiltFiles, 'RIGHT')};
            
            
            lDataG = flipud(helper.open_raw([unfiltPath, leftGreenFile]) );
            rDataG = helper.open_raw([unfiltPath, rightGreenFile]);
            
            [~, lDataz] = processUnfilteredData(lDataG, imgResizeFactor, tform{i}, mask);
            [~, rDataz] = processUnfilteredData(rDataG, imgResizeFactor, tform{N+i}, mask);
            

            socialData(:,:,:,i,1) = lDataz(:,:,i1:i2);
            socialData(:,:,:,i,2) = rDataz(:,:,i1:i2);


            % global signals
            GS.left{i} = squeeze(nanmean(nanmean(lDataz)));
            GS.right{i} = squeeze(nanmean(nanmean(rDataz)));

            % calculate correlations
            tmp = corrcoef(GS.left{i}(t1 - minute:t1), GS.right{i}(t1 - minute:t1));
            C.before(i) = tmp(2,1);

            tmp = corrcoef(GS.left{i}(i1:i1+minute), GS.right{i}(i1:i1+minute));
            C.during(i) = tmp(2,1);

            tmp = corrcoef(GS.left{i}(t2:t2+minute), GS.right{i}(t2:t2+minute));
            C.after(i) = tmp(2,1);
                        
            clear lData rData
            
            % crop global signals to 90s before/after translation 
            % initiation/end to remove initial dFF decay
            GS.left{i} = GS.left{i}(t1-s90:t2+s90);
            GS.right{i} = GS.right{i}(t1-s90:t2+s90);

            toc
        end
        
        varargout{1} = socialData;
        varargout{2} = C;
        varargout{3} = GS;
        
    case 'global_signal_correlation_interbrain_corrected'
        % returns:
        %   socialData  (5D matrix of shape H x W x Frame x Trial x 2)
        %   corrs       (structure of correlation in minute long periods
        %               immediately before translation, immediately after
        %               first translation ends, or immediately after 2nd
        %               translation ends)
        %   GS          (structure that contains global signals for left 
        %               and right mouse from 90s before translation1 start 
        %               to 90s after translation2 end)

        
        % pre-allocate
        socialData.green = zeros(128, 128, 3459, N, 2);
        socialData.blue = zeros(128, 128, 3459, N, 2);
        C.before = zeros(length(filenames),1);
        C.during = zeros(length(filenames),1);
        C.after = zeros(length(filenames),1);
        GS.left = cell(length(filenames),1);
        GS.right = cell(length(filenames),1);
        
        % check minute in each trial phase
        s90 = round(90*fs);
        for i = 1:N
            tic

            % load corrected brain data
            leftGreenFile = unfiltFiles{contains(unfiltFiles, filenames{i}(1:11)) & ...
                contains(unfiltFiles, 'GREEN') & ...
                contains(unfiltFiles, 'LEFT')};
            rightGreenFile = unfiltFiles{contains(unfiltFiles, filenames{i}(1:11)) & ...
                contains(unfiltFiles, 'GREEN') & ...
                contains(unfiltFiles, 'RIGHT')};
            
            leftBlueFile = unfiltFiles{contains(unfiltFiles, filenames{i}(1:11)) & ...
                contains(unfiltFiles, 'BLUE') & ...
                contains(unfiltFiles, 'LEFT')};
            rightBlueFile = unfiltFiles{contains(unfiltFiles, filenames{i}(1:11)) & ...
                contains(unfiltFiles, 'BLUE') & ...
                contains(unfiltFiles, 'RIGHT')};
            
            lDataG = flipud(helper.open_raw([unfiltPath, leftGreenFile]) );
            rDataG = helper.open_raw([unfiltPath, rightGreenFile]);
            
            lDataB = flipud(helper.open_raw([unfiltPath, leftBlueFile]) );
            rDataB = helper.open_raw([unfiltPath, rightBlueFile]);
            
            lDataG = processData(lDataG, imgResizeFactor, tform{i}, mask, 3.5);
            rDataG = processData(rDataG, imgResizeFactor, tform{N+i}, mask, 3.5);
            
            lDataB = processData(lDataB, imgResizeFactor, tform{i}, mask, 0.15);
            rDataB = processData(rDataB, imgResizeFactor, tform{N+i}, mask, 0.15);

            
            GS.left{i} = zscore(squeeze(nanmedian(nanmedian(lDataG))) - ...
                squeeze(nanmedian(nanmedian(lDataB))));
            GS.right{i} = zscore(squeeze(nanmedian(nanmedian(rDataG))) - ...
                squeeze(nanmedian(nanmedian(rDataB))));

            % calculate correlations
            tmp = corrcoef(GS.left{i}(t1 - minute:t1), GS.right{i}(t1 - minute:t1));
            C.before(i) = tmp(2,1);

            tmp = corrcoef(GS.left{i}(i1:i1+minute), GS.right{i}(i1:i1+minute));
            C.during(i) = tmp(2,1);

            tmp = corrcoef(GS.left{i}(t2:t2+minute), GS.right{i}(t2:t2+minute));
            C.after(i) = tmp(2,1);
             
            clear lData rData
            
            % crop global signals to 90s before/after translation 
            % initiation/end to remove initial dFF decay
            GS.left{i} = GS.left{i}(t1-s90:t2+s90);
            GS.right{i} = GS.right{i}(t1-s90:t2+s90);

            toc
        end
        
        varargout{1} = socialData;
        varargout{2} = C;
        varargout{3} = GS;
        
    case 'roi_signal_correlation_interbrain'               
        radius = 2;
        rMat = struct();
        leftTrace = cell(1,N);
        rightTrace = cell(1,N);
        
        for i = 1:N
            % load corrected brain data
            leftGreenFile = unfiltFiles{contains(unfiltFiles, filenames{i}(1:11)) & ...
                contains(unfiltFiles, 'GREEN') & ...
                contains(unfiltFiles, 'LEFT')};
            rightGreenFile = unfiltFiles{contains(unfiltFiles, filenames{i}(1:11)) & ...
                contains(unfiltFiles, 'GREEN') & ...
                contains(unfiltFiles, 'RIGHT')};
            
            
            lDataG = flipud(helper.open_raw([unfiltPath, leftGreenFile]));
            rDataG = helper.open_raw([unfiltPath, rightGreenFile]);
            
            
            mult = 3.5;
            filterOn = 0;
            [~, lDataz] = processData(lDataG, imgResizeFactor, tform{i}, mask, mult, filterOn);
            [~, rDataz] = processData(rDataG, imgResizeFactor, tform{N+i}, mask, mult, filterOn);
            
            
            leftTrace{i} = helper.getTimeseries(lDataz, CL, CR, radius);
            rightTrace{i} = helper.getTimeseries(rDataz, CL, CR, radius);
            
            rMat.before(:,:,i) = corrcoef([leftTrace{i}(t1-minute:t1,:), rightTrace{i}(t1-minute:t1,:)]);
            rMat.during(:,:,i) = corrcoef([leftTrace{i}(i1:i1+minute,:), rightTrace{i}(i1:i1+minute,:)]);
            rMat.after(:,:,i)  = corrcoef([leftTrace{i}(t2:t2+minute,:), rightTrace{i}(t2:t2+minute,:)]);            
            
        end
        
        varargout{1} = rMat;
        varargout{2} = leftTrace;
        varargout{3} = rightTrace;

    case 'roi_signal_correlation_interbrain_corrected'                
        radius = 2;
        rMat = struct();
        leftTrace = cell(1,N);
        rightTrace = cell(1,N);
        
        for i = 1:N           
            % load corrected brain data
            leftGreenFile = unfiltFiles{contains(unfiltFiles, filenames{i}(1:11)) & ...
                contains(unfiltFiles, 'GREEN') & ...
                contains(unfiltFiles, 'LEFT')};
            rightGreenFile = unfiltFiles{contains(unfiltFiles, filenames{i}(1:11)) & ...
                contains(unfiltFiles, 'GREEN') & ...
                contains(unfiltFiles, 'RIGHT')};
            
            leftBlueFile = unfiltFiles{contains(unfiltFiles, filenames{i}(1:11)) & ...
                contains(unfiltFiles, 'BLUE') & ...
                contains(unfiltFiles, 'LEFT')};
            rightBlueFile = unfiltFiles{contains(unfiltFiles, filenames{i}(1:11)) & ...
                contains(unfiltFiles, 'BLUE') & ...
                contains(unfiltFiles, 'RIGHT')};
            
            lDataG = flipud(helper.open_raw([unfiltPath, leftGreenFile]));
            rDataG = helper.open_raw([unfiltPath, rightGreenFile]);
            
            lDataB = flipud(helper.open_raw([unfiltPath, leftBlueFile]) );
            rDataB = helper.open_raw([unfiltPath, rightBlueFile]);
            
            mult = 3.5;
            filterOn = 1;
            lDataG = processData(lDataG, imgResizeFactor, tform{i}, mask, mult, filterOn);
            rDataG = processData(rDataG, imgResizeFactor, tform{N+i}, mask, mult, filterOn);
            
            mult = 0.15;
            lDataB = processData(lDataB, imgResizeFactor, tform{i}, mask, mult, filterOn);
            rDataB = processData(rDataB, imgResizeFactor, tform{N+i}, mask, mult, filterOn);
            
            lDataz = zscore(lDataG - lDataB, [], 3);
            rDataz = zscore(rDataG - rDataB, [], 3);
            
            leftTrace{i} = helper.getTimeseries(lDataz, CL, CR, radius);
            rightTrace{i} = helper.getTimeseries(rDataz, CL, CR, radius);
            

            rMat.before(:,:,i) = corrcoef([leftTrace{i}(t1-minute:t1,:), rightTrace{i}(t1-minute:t1,:)]);
            rMat.during(:,:,i) = corrcoef([leftTrace{i}(i1:i1+minute,:), rightTrace{i}(i1:i1+minute,:)]);
            rMat.after(:,:,i)  = corrcoef([leftTrace{i}(t2:t2+minute,:), rightTrace{i}(t2:t2+minute,:)]);            
            
        end
        
        varargout{1} = rMat;
        varargout{2} = leftTrace;
        varargout{3} = rightTrace;
        
    case 'behavior_modulation'
        % Examine effect of animal A behavior on animal B brain
        %
        % returns:
        %   bFrames     (structure containing cell arrays of length N
        %               with 4D matrices of behavior events (HxWxTxW))
      
        % constants
        window = 2; window = round(window * fs);

        for i = 1:N
            tic           
            % ignore trial with corrupted behaviour data
            if contains(filenames{i}, 'August-21_1123')
                continue
            end

            % --- Process brain data ---
            leftGreenFile = unfiltFiles{contains(unfiltFiles, filenames{i}(1:11)) & ...
                contains(unfiltFiles, 'GREEN') & ...
                contains(unfiltFiles, 'LEFT')};
            rightGreenFile = unfiltFiles{contains(unfiltFiles, filenames{i}(1:11)) & ...
                contains(unfiltFiles, 'GREEN') & ...
                contains(unfiltFiles, 'RIGHT')};
            
            lDataG = flipud(helper.open_raw([unfiltPath, leftGreenFile]) );
            rDataG = helper.open_raw([unfiltPath, rightGreenFile]);
            
            [~, lDataz] = processUnfilteredData(lDataG, imgResizeFactor, tform{i}, mask);
            [~, rDataz] = processUnfilteredData(rDataG, imgResizeFactor, tform{N+i}, mask);

            clear lDataG rDataG
            % --------------------------               
                           
            % get event initiation indices
            wEventsBeforeLeft = helper.getBehaviourEvents(B.exclusive.whiskLeftBefore{i}, window);
            wEventsDuringLeft = helper.getBehaviourEvents(B.exclusive.whiskLeftDuring{i}, window);
            wEventsDuringRight = helper.getBehaviourEvents(B.exclusive.whiskRightDuring{i}, window);
            fEventsBeforeLeft = helper.getBehaviourEvents(B.exclusive.flLeftBefore{i}, window);
            fEventsDuringLeft = helper.getBehaviourEvents(B.exclusive.flLeftDuring{i}, window);
            fEventsDuringRight = helper.getBehaviourEvents(B.exclusive.flRightDuring{i}, window);
            
            % self initiated maps
            bFrames.selfInitiatedWhiskLeftBefore{i} = helper.getWhiskFrames(lDataz, wEventsBeforeLeft, window);
            bFrames.selfInitiatedWhiskLeftDuring{i} = helper.getWhiskFrames(lDataz, wEventsDuringLeft, window);
            bFrames.selfInitiatedWhiskRightDuring{i} = helper.getWhiskFrames(rDataz, wEventsDuringRight, window);
            
            bFrames.selfInitiatedFLLeftBefore{i} = helper.getWhiskFrames(lDataz, fEventsBeforeLeft, window);
            bFrames.selfInitiatedFLLeftDuring{i} = helper.getWhiskFrames(lDataz, fEventsDuringLeft, window);
            bFrames.selfInitiatedFLRightDuring{i} = helper.getWhiskFrames(rDataz, fEventsDuringRight, window);
            
            % partner-initiated maps (Left/Right in variable name refers to
            % Left/Right maps)
            bFrames.partnerInitiatedWhiskLeftDuring{i} = helper.getWhiskFrames(lDataz, wEventsDuringRight, window);
            bFrames.partnerInitiatedWhiskRightDuring{i} = helper.getWhiskFrames(rDataz, wEventsDuringLeft, window);
            bFrames.partnerInitiatedWhiskRightBefore{i} = helper.getWhiskFrames(rDataz, wEventsBeforeLeft, window);
            
            bFrames.partnerInitiatedFLLeftDuring{i} = helper.getWhiskFrames(lDataz, fEventsDuringRight, window);
            bFrames.partnerInitiatedFLRightDuring{i} = helper.getWhiskFrames(rDataz, fEventsDuringLeft, window);
            bFrames.partnerInitiatedFLRightBefore{i} = helper.getWhiskFrames(rDataz, fEventsBeforeLeft, window);
            
            clear wEvents* fEvents*

            toc
        end
        
        varargout{1} = bFrames;
        
    case 'behavior_modulation_corrected'
        % Examine effect of animal A behavior on animal B brain
        %
        % returns:
        %   bFrames     (structure containing cell arrays of length N
        %               with 4D matrices of behavior events (HxWxTxW))

        
        % constants
        window = 2; window = round(window * fs);      

        for i = 1:N
            tic
            
            % ignore trials with corrupted behaviour data
            if i==15 || i == 27
                continue
            end
            
            % --- Process brain data ---
            leftGreenFile = unfiltFiles{contains(unfiltFiles, filenames{i}(1:11)) & ...
                contains(unfiltFiles, 'GREEN') & ...
                contains(unfiltFiles, 'LEFT')};
            rightGreenFile = unfiltFiles{contains(unfiltFiles, filenames{i}(1:11)) & ...
                contains(unfiltFiles, 'GREEN') & ...
                contains(unfiltFiles, 'RIGHT')};
            
            leftBlueFile = unfiltFiles{contains(unfiltFiles, filenames{i}(1:11)) & ...
                contains(unfiltFiles, 'BLUE') & ...
                contains(unfiltFiles, 'LEFT')};
            rightBlueFile = unfiltFiles{contains(unfiltFiles, filenames{i}(1:11)) & ...
                contains(unfiltFiles, 'BLUE') & ...
                contains(unfiltFiles, 'RIGHT')};
            
            lDataG = flipud(helper.open_raw([unfiltPath, leftGreenFile]) );
            rDataG = helper.open_raw([unfiltPath, rightGreenFile]);
            
            lDataB = flipud(helper.open_raw([unfiltPath, leftBlueFile]) );
            rDataB = helper.open_raw([unfiltPath, rightBlueFile]);
            
            mult = 3.5;
            filterOn = 1;
            lDataG = processData(lDataG, imgResizeFactor, tform{i}, mask, mult, filterOn);
            rDataG = processData(rDataG, imgResizeFactor, tform{N+i}, mask, mult, filterOn);
            
            mult = 0.15;
            lDataB = processData(lDataB, imgResizeFactor, tform{i}, mask, mult, filterOn);
            rDataB = processData(rDataB, imgResizeFactor, tform{N+i}, mask, mult, filterOn);
            
            lDataz = zscore(lDataG-lDataB, [], 3);
            rDataz = zscore(rDataG-rDataB, [], 3);

            clear lDataG rDataG lDataB rDataB
            % --------------------------               
                        
            % get event initiation indices
            wEventsBeforeLeft = helper.getBehaviourEvents(B.exclusive.whiskLeftBefore{i}, window);
            wEventsDuringLeft = helper.getBehaviourEvents(B.exclusive.whiskLeftDuring{i}, window);
            wEventsDuringRight = helper.getBehaviourEvents(B.exclusive.whiskRightDuring{i}, window);
            wEventsAfterLeft = helper.getBehaviourEvents(B.exclusive.whiskLeftAfter{i}, window);
%             fEventsBeforeLeft = helper.getWhiskEvents(B.exclusive.flLeftBefore{i}(1:bEndL), window);
%             fEventsDuringLeft = helper.getWhiskEvents(B.exclusive.flLeftDuring{i}(1:bEndL), window);
%             fEventsDuringRight = helper.getWhiskEvents(B.exclusive.flRightDuring{i}(1:bEndL), window);
            
            % self initiated maps
            bFrames.selfInitiatedWhiskLeftBefore{i} = helper.getWhiskFrames(lDataz, wEventsBeforeLeft, window);
            bFrames.selfInitiatedWhiskLeftDuring{i} = helper.getWhiskFrames(lDataz, wEventsDuringLeft, window);
            bFrames.selfInitiatedWhiskRightDuring{i} = helper.getWhiskFrames(rDataz, wEventsDuringRight, window);
            bFrames.selfInitiatedWhiskLeftAfter{i} = helper.getWhiskFrames(lDataz, wEventsAfterLeft, window);
            
%             bFrames.selfInitiatedFLLeftBefore{i} = helper.getWhiskFrames(lDataz, fEventsBeforeLeft, window);
%             bFrames.selfInitiatedFLLeftDuring{i} = helper.getWhiskFrames(lDataz, fEventsDuringLeft, window);
%             bFrames.selfInitiatedFLRightDuring{i} = helper.getWhiskFrames(rDataz, fEventsDuringRight, window);
            
            % partner-initiated maps (Left/Right in variable name refers to
            % Left/Right maps)
            bFrames.partnerInitiatedWhiskLeftDuring{i} = helper.getWhiskFrames(lDataz, wEventsDuringRight, window);
            bFrames.partnerInitiatedWhiskRightDuring{i} = helper.getWhiskFrames(rDataz, wEventsDuringLeft, window);
            bFrames.partnerInitiatedWhiskRightBefore{i} = helper.getWhiskFrames(rDataz, wEventsBeforeLeft, window);
            bFrames.partnerInitiatedWhiskRightAfter{i} = helper.getWhiskFrames(rDataz, wEventsAfterLeft, window);
            
%             bFrames.partnerInitiatedFLLeftDuring{i} = helper.getWhiskFrames(lDataz, fEventsDuringRight, window);
%             bFrames.partnerInitiatedFLRightDuring{i} = helper.getWhiskFrames(rDataz, fEventsDuringLeft, window);
%             bFrames.partnerInitiatedFLRightBefore{i} = helper.getWhiskFrames(rDataz, fEventsBeforeLeft, window);
            
            clear wEvents* fEvents*

            toc
        end
        
        varargout{1} = bFrames;
        
        
    case 'barrier_controls'
        radius = 2;
        rMat = struct();
        leftTrace = cell(1,N);
        rightTrace = cell(1,N);
        
        for i = 1:N
            load([pathname, filenames{i}],'left_dFF', 'right_dFF')
            
            lDataG = imrotate(permute(left_dFF, [2, 3, 1]), 90);
            rDataG = permute(right_dFF, [3, 2, 1]);
            
            [~, lDataz] = processData(lDataG, imgResizeFactor, tform{i}, mask);
            [~, rDataz] = processData(rDataG, imgResizeFactor, tform{N+i}, mask); 
            
            leftTrace{i} = helper.getTimeseries(lDataz, CL, CR, radius);
            rightTrace{i} = helper.getTimeseries(rDataz, CL, CR, radius);           

            rMat.before(:,:,i) = corrcoef([leftTrace{i}(t1-minute:t1,:), rightTrace{i}(t1-minute:t1,:)]);
            rMat.during(:,:,i) = corrcoef([leftTrace{i}(i1:i1+minute,:), rightTrace{i}(i1:i1+minute,:)]);
            rMat.after(:,:,i)  = corrcoef([leftTrace{i}(t2:t2+minute,:), rightTrace{i}(t2:t2+minute,:)]);            
                        
        end
        
        varargout{1} = rMat;
        varargout{2} = leftTrace;
        varargout{3} = rightTrace;     
        
    case 'barrier_controls_corrected'
        radius = 2;
        rMat = struct();
        leftTrace = cell(1,N);
        rightTrace = cell(1,N);
        
        for i = 1:N           
            load([pathname, filenames{i}],'left_dFF_green', 'right_dFF_green', 'left_dFF_blue', ...
                'right_dFF_blue')
            
            lDataG = imrotate(permute(left_dFF_green, [2, 3, 1]), 90);
            rDataG = permute(right_dFF_green, [3, 2, 1]);
            lDataB = imrotate(permute(left_dFF_blue, [2, 3, 1]), 90);
            rDataB = permute(right_dFF_blue, [3, 2, 1]);
            
            mult = 3.5;
            filterOn = 0;
            lDataG = processData(lDataG, imgResizeFactor, tform{i}, mask, mult, filterOn);
            rDataG = processData(rDataG, imgResizeFactor, tform{N+i}, mask, mult, filterOn);
            
            mult = 0.15;
            lDataB = processData(lDataB, imgResizeFactor, tform{i}, mask, mult, filterOn);
            rDataB = processData(rDataB, imgResizeFactor, tform{N+i}, mask, mult, filterOn);
            
            lDataz = zscore(lDataG - lDataB, [], 3);
            rDataz = zscore(rDataG - rDataB, [], 3);
            
            leftTrace{i} = helper.getTimeseries(lDataz, CL, CR, radius);
            rightTrace{i} = helper.getTimeseries(rDataz, CL, CR, radius);
                       
            rMat.before(:,:,i) = corrcoef([leftTrace{i}(t1-minute:t1,:), rightTrace{i}(t1-minute:t1,:)]);
            rMat.during(:,:,i) = corrcoef([leftTrace{i}(i1:i1+minute,:), rightTrace{i}(i1:i1+minute,:)]);
            rMat.after(:,:,i)  = corrcoef([leftTrace{i}(t2:t2+minute,:), rightTrace{i}(t2:t2+minute,:)]);            
                        
        end
        
        varargout{1} = rMat;
        varargout{2} = leftTrace;
        varargout{3} = rightTrace;
    
    otherwise
        error('Choose valid analysis type')
end



end



% function [data, dataz] = processData(data, imgResizeFactor, tform, mask)
% % Process brain data
% %   resize and orient properly
% %   register
% %   apply mask
% %   spatial smoothing
% %   z-score
% 
% data = imresize(data, imgResizeFactor);
% 
% if ~isempty(tform)
%     data = imwarp(data, tform, ...
%         'OutputView', imref2d([size(data,1), size(data,2)]), 'FillValues', nan);
% end            
% data = bsxfun(@times, data, cast(mask, 'like', data));
% data = imgaussfilt(data, 1);
% dataz = zscore(data,[],3);
% 
% end


function [data, dataz] = processData(data, imgResizeFactor, tform, mask, mult, filterOn)
if nargin < 6 || isempty(filterOn), filterOn = 1; end
tic

data = imresize(data, imgResizeFactor, 'bilinear');

% register left brain
if ~isempty(tform)
    data = imwarp(data, tform, ...
        'OutputView', imref2d([size(data,1), size(data,2)]), 'FillValues', nan);
end

% apply mask
data =  double(bsxfun(@times, data, cast(mask, 'like', data)));



% crop values 
if mult < 1
    data(data<-mult) = -mult;
    data(data>mult) = mult;
else    
    d = reshape(data,[],1);
    m = nanmean(d);
    s = nanstd(d);
    data(data<(m-(mult*s))) = m-(mult*s);
    data(data>m+mult*s) = m+mult*s;
end
% smooth
data(isnan(data)) = 0;
data = imgaussfilt(data, 1);

if filterOn, data = single(helper.image_filter(data)); end

if nargout > 1, dataz = zscore(data, [], 3); end
toc

end
