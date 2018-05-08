% Example : PlotCorrelationMatrix(corrmat, [], 'corrmat') 
%               Max scale depends on maximum of correlation matrix value
%               Save figure as corrmat_minsc-1_maxsc1.jpg
%
%           PlotCorrelationMatrix(corrmat, [-0.5 0.5])
%               Will not save figure
%
% Note that the highly "dense" white lines will become more appropriate when the figure is saved.
%
% inputs are 	1) Correlation matrix, symmetric along diagonal line (i.e.
%                       corrmat(1,3) is equal to corrmat(3,1)
%               2) Scale limit : min and max scale limit
%                  If not specified, or specified as [], scale limit is -1*max(abs(corr_mat)) to max(abs(corr_mat))
%               3) Filename : basename. E.g. 'RWRestingCorrMat', final filename = RWRestingCorrMat_min-3_max3.jpg
%                  If not specified, figure will not be saved
%
% Major networks separated by thick white grid lines
% Thin white grid lines separate the breakdowns of major networks
%
% Ordering of major networks from left to right 
%
%       Major network     Sub-network
%     
%     1)  Default :       DefaultC
%                         DefaultB
%                         DefaultA
%         
%     2)  Control :       ContC
%                         ContB
%                         ContA
%         
%     3)  Limbic
%     
%     4)  SalVentAttn :   SalVentAttnB
%                         SalVentAttnA
%         
%     5)  DorsAttn :      DorsAttnB
%                         DorsAttnA
%         
%     6)  SomMot :        SomMotB
%                         SomMotA
%     
%     7)  Visual :        VisPeri
%                         VisCent
%
%     8)  TempPar
%
%     9)  Subcortical :   Thalamus
%                         Striatum
%                         Hippocampus
%                         Amygdala
%
% Within each sub-network, correlation matrix entries start from left
% hemisphere, then right hemisphere entries (from left to right).
    
%=========================================================================
%
%  Copyright (c) 2013 Jesisca Tandi & Thomas Yeo
%  All rights reserved.
%
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are met:
%
%    * Redistributions of source code must retain the above copyright notice,
%      this list of conditions and the following disclaimer.
%
%    * Redistributions in binary form must reproduce the above copyright notice,
%      this list of conditions and the following disclaimer in the documentation
%      and/or other materials provided with the distribution.
%
%    * Neither the names of the copyright holders nor the names of future
%      contributors may be used to endorse or promote products derived from this
%      software without specific prior written permission.
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
%ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
%WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
%ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
%(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
%LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
%ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
%SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.    
%
%=========================================================================

function curr_net = PlotCorrelationMatrix(corr_mat, scalelim, filename_prefix)

% Ensure it's symmetric
corr_mat = (corr_mat+corr_mat')/2;
% Get rearrangement indexes
[curr_net, Index] = LabelsRearrangebyNetwork;
corr_mat = corr_mat(Index,Index);

% Load colormap
load colormap_redblue.mat;

% Plot corr_mat using imagesc
figure; imagesc(corr_mat); set(gcf,'Colormap',rbmap2);
    
% Generate thin grid lines
[xline, yline, ymaj] = generateline(size(corr_mat,1));
patch(xline, yline,'w', 'edgecolor', 'w', 'Linewidth', 0.05, 'EdgeAlpha', 0.2)
patch(yline, xline,'w', 'edgecolor', 'w', 'Linewidth', 0.05, 'EdgeAlpha', 0.2)

xlim(gca,[1 size(corr_mat, 1)]);
ylim(gca,[1 size(corr_mat, 1)]);
postn = get(gca, 'Position');
postn(2) = 0.15;
set(gca, 'Position', postn);

% Set colorbar
hcol=colorbar('peer',gca,'SouthOutside');
cpos=get(hcol,'Position');
cpos(4)=cpos(4)/4; % Halve the thickness
cpos(3)=cpos(3)*0.75; % Reduce length
cpos(1)=cpos(1) + 0.1; % Move it to the center
cpos(2)=cpos(2) - 0.12; % Move it down outside the plot
set(hcol,'Position',cpos);

% Set color limit 
if ((nargin < 2) || (isempty(scalelim)))
	collim = max(max(abs(corr_mat)));
	scalelim = [-1*collim, 1*collim];
end

set(gca, 'CLim', scalelim);

axis equal;
grid off;
axis([-15 size(corr_mat, 1)+15.5 -15 size(corr_mat, 1)+15.5]);
set(gca, 'Visible', 'off');
set(gcf, 'color', 'white');

% Generate major and minor gridlines 
	minor_grid = [6 15 28 39 52 67 86 100 108 116 118 120];
	major_grid = [24 50 54 78 92 102 112 114];

	patch(xline(:,minor_grid), yline(:,minor_grid),'w', 'edgecolor', 'w', 'Linewidth', 0.06, 'EdgeAlpha', 0.4)
	patch(yline(:,minor_grid), xline(:,minor_grid),'w', 'edgecolor', 'w', 'Linewidth', 0.06, 'EdgeAlpha', 0.4)
	patch(xline(:,major_grid), ymaj(:,major_grid),'w', 'edgecolor', 'w', 'Linewidth', 1.2)
	patch(ymaj(:,major_grid), xline(:,major_grid),'w', 'edgecolor', 'w', 'Linewidth', 1.2)

    
% Generate labels on the left side
    for k = 1:numel(major_grid)
        line([-15 -0.5], [major_grid(k)+0.5 major_grid(k)+0.5], 'Color', [0.3 0.3 0.3]);
    end
    for j = 1:numel(minor_grid)
        line([-5 -0.5], [minor_grid(j)+0.5 minor_grid(j)+0.5]);
    end

    labelstoprint = {'Default' 'Control' 'Limbic' 'Sal/Vent\nAttn' 'Dors\nAttn' 'Som\nMot' 'Visual' 'Sub\ncortical'};
    labels_pos = {[-12 12.5] [-12 37.5] [-12 52.5] [-12 66.5] [-12 85.5] [-12 97.5] [-12 107.5] [-12 115.5]};
    for w = 1:numel(labelstoprint)
        q = text(labels_pos{w}(1), labels_pos{w}(2), sprintf(labelstoprint{w}));
        if sum(strcmp(labelstoprint{w}, {'Limbic' 'Visual' }))>0
            set(q, 'HorizontalAlignment', 'center', 'FontSize', 9,'VerticalAlignment', 'middle');
        elseif sum(strcmp(labelstoprint{w}, {'Sub\ncortical' }))>0
            set(q, 'Rotation', 90, 'HorizontalAlignment', 'right','VerticalAlignment', 'middle');
        else
            set(q, 'Rotation', 90, 'HorizontalAlignment', 'center','VerticalAlignment', 'middle');
        end
    end
    
    sublabelstoprint = {'A' 'B' 'C'};
    sublabels_pos = {[-2.5 19.5] [-2.5 10.5] [-2.5 3]};
    for t = 1:numel(sublabelstoprint)
        text(sublabels_pos{t}(1), sublabels_pos{t}(2), sublabelstoprint{t}, 'HorizontalAlignment', 'center','VerticalAlignment', 'middle')
    end
   
% Save figure
    if ((nargin == 3) && ~isempty(filename_prefix))
        filenamefin = [filename_prefix '_minsc' num2str(scalelim(1), '%.2f') '_maxsc' num2str(scalelim(2), '%.2f') '.jpg'];
        set(gcf, 'PaperPositionMode', 'auto');
        print(gcf, '-djpeg', '-r600', filenamefin);
    end
end






function [x,y, ymaj] = generateline(n)
x = 1.5:1:n;
x = [ x; x; repmat(nan,1,(n-1)) ];
y = [ 0.5 n+0.5 nan ].';
y = repmat(y,1,(n-1));

ymaj = [ -5 n+5.5 nan]';
% For short grid (same as the figure boundary)
% ymaj = [ 0.5 n+0.5 nan]';
ymaj = repmat(ymaj,1,(n-1));
end


function [curr_net, Index] = LabelsRearrangebyNetwork

% Load original labels
labeldir = pwd;
roi_fid = fopen(fullfile(labeldir, 'ROIslist.txt'));
roi_label = textscan(roi_fid, '%s %s %s', 'delimiter', '.'); roi_label = roi_label{1,2}; 

roi_netw = cell(numel(roi_label),1);
roi_subnet = cell(numel(roi_label),1);



for i = 1:numel(roi_label);
    netw = textscan(char(roi_label(i,1)), '%s %s %s %s', 'delimiter', '_');
    tmp = netw{1,3};
    roi_netw(i,1) = tmp;
    subnet = netw{1,4};
    if ~isempty(subnet)
        roi_subnet(i,1) = subnet;
    end
end


% Arrange new label based on template order
tmplate = {'DefaultC'; 'DefaultB';'DefaultA'; 'ContC'; 'ContB'; 'ContA'; 'Limbic'; 'SalVentAttnB'; 'SalVentAttnA'; 'DorsAttnB'; 'DorsAttnA'; 'SomMotB'; 'SomMotA'; 'VisPeri'; 'VisCent'; 'TempPar'; 'Thalamus'; 'Striatum'; 'Hippocampus'; 'Amygdala'};
tmplate2 = {'TempPole'; 'OFC'};
% initiate new label
newlabel = [];
curr_net = NaN(size(roi_netw));

for j = 1:numel(tmplate);
    ind = strmatch(tmplate(j), roi_netw);
    curr_net(ind) = j;
    % For Limbic, further separate the networks to TempPole and OFC
    if ~isempty(strmatch('Limbic', tmplate(j)));
        ind2 = [];
        for s = 1:numel(tmplate2);
            if ~isempty(strmatch(tmplate2(s), roi_subnet(ind)))
                ind2 = [ind2; ind(strmatch(tmplate2(s), roi_subnet(ind)))];
            end
            
        end
        if numel(ind) ~= numel(ind2)
            disp('Wrong Index')
        end
        ind = ind2;
    end
    
    newlabel = cat(1, newlabel,roi_label(ind));
end

% Create indexing from old labeling to new labeling
Index = zeros(numel(roi_label),1);

for k = 1:numel(roi_label);
Index(k) = strmatch(newlabel(k), roi_label, 'exact');
end

fin_label = roi_label(Index);
% Save new label
fid = fopen('ROIslist_reordered.txt', 'w');
for l = 1:numel(roi_label);
   fprintf(fid, '%s\n', char(fin_label(l)));
end
end
