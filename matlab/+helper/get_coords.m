function [CL,CR, bregma] = get_coords(data,scale_factor,ROIs)
% Get coordinates for cortical regions
%
% Inputs:
%   data            (required, image or image series)
%   scale_factor    (optional, scalar - specify mm/px)
%   ROIs            (optional, cell array containing string of ROI labels)
%
% Outputs:
%   CL, CR          (Nx2 coordinates [x,y] of ROIs on left and right
%                   hemisphere)  
%   bregma          (coordinates [x,y] of bregma location)
% 
% Usage:  coords = get_coords(data,scale_factor,ROIs);

if nargin < 2 || isempty(scale_factor)
    scale_factor = 8.4/size(data,1); 
end

if nargin < 3 || isempty(ROIs)
    load ai.mat atlas
    ROIs = atlas.areatag;
end

flag = 1;
while flag == 1
    [bx,by] = get_breg(data); hold on
    
    % pre-allocate
    CL = nan(length(ROIs),2);
    CR = nan(length(ROIs),2);

    for i = 1:length(ROIs)
        switch ROIs{i}
            case 'A'
                ML=2.2932;   AP=-2.4962;
            case 'AC'
                ML=0.097951; AP=1.8536;
            case 'AL'
                ML=3.8271;   AP=-3.3393;
            case 'AM'
                ML=1.6479;   AP=-2.696;
            case 'AU'
                ML=4.5304;   AP=-2.901;
            case 'BC'
                ML=3.4569;   AP=-1.727;
            case 'FL'
                ML=2.4526;   AP=-0.5668;
            case 'HL'
                ML=1.6942;   AP=-1.1457;
            case 'L'
                ML=3.7126;   AP=-4.2615;
            case 'LI'
                ML=4.0586;   AP=-4.2293;
            case 'M1'
                ML=1.8603;   AP=0.64181;
            case 'M2'
%                 ML=0.87002;  AP=1.4205;
                ML = 1;      AP=2.5;
            case 'MO'
                ML=3.4917;   AP=0.58712;
            case 'NO'
                ML=3.8001;   AP=-0.47733;
            case 'PL'
                ML=3.5161;   AP=-5.2146;
            case 'PM'
                ML=1.6217;   AP=-3.6247;
            case 'POR'
                ML=4.2231;   AP=-4.755;
            case 'RL'
                ML=3.1712;   AP=-2.849;
            case 'RS'
                ML=0.62043 + 0.2;  AP=-2.8858;
            case 'S2'
                ML=4.3977;   AP=-1.2027;
            case 'TEA'
                ML=4.5657;   AP=-4.1622;
            case 'TR'
                ML=1.8644;   AP=-2.0204;
            case 'UN'
                ML=2.7979;   AP=-0.97112;
            case 'V1'
%                 ML=2.5168;   AP=-4.2678;
                ML=2.5618;   AP=-4;
            case 'ALM'
%                 ML=2.5;      AP=1.8;
%                 ML=1.8;      AP=2.5;
                ML=2.0;      AP=2.4;
            case 'AMA'
                ML=1.6;      AP=2.4;
            case 'OMA'
                ML=2.5;      AP=1.8;
            case 'aBC'
                AP = -1.36 + 0.575;   ML = 3.35 + 0.5;
            case 'pBC'
                AP = -1.9;            ML = 3.35 + 0.5;
            case 'wM1'
%                 AP = 1-0.3;               ML = 0.6 + 0.5;
%                 AP = 1.9;                ML = 1.4;
                  AP = 1.15;                ML = 0.9;
            case 'PtA'
                AP = -1.94;             ML = 1.2;
        end
        ML = ML/scale_factor;
        AP = AP/scale_factor;

        CL(i,:) = round([bx-ML by-AP]);
        CR(i,:) = round([bx+ML by-AP]);
        
        plot(CL(i,1),CL(i,2),'r*');
%         plot(CR(i,1),CR(i,2),'r*');
        text(CR(i,1),CR(i,2),ROIs{i},'Color','red');
    end
    
    % Prompt user to evaluate ROI placement
    answer = questdlg('Would you like try again?', ...
        'Attention', ...
        'Yes','No','Cancel','Cancel');
    
    switch answer
        case 'Yes'
            flag = 1;
            close(gcf)
        case 'No'
            flag = 0;
            bregma = round([bx, by]);
        case 'Cancel'
            error("Operation terminated by user.")
    end
    
    
end