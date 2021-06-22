function tmpFig = Publish_Figure(varargin)
%
% Function to create a qausi-publication ready figure for HystLab
%
% Last Modified 2019/05/07
%

%% Input checking
% Should be redundant as this function is only internally called

nPlot = nargin;

% Check we have handles
if sum(cellfun(@ishghandle, varargin)) ~= nPlot
    error('Publish_Figure:Input', 'At least one input is not a graphics object handle.');
end

PlotAxes = varargin;

%% Define some parameters

switch nPlot
    case 1
        PlotSize = 4.5;
        FigSize = [9, 7.5];
    case 2
        PlotSize = 4.5;
        FigSize = [18, 7.5];
    case 3
        PlotSize = 4;
        FigSize = [18, 7.5];
    otherwise
        error('Publish_Figure:nPlot', 'Currently more than 3 plots are not supported.');
end

FontSize1 = 8;
FontSize2 = 9;
FontSize3 = 10;
FontSize4 = 11;

FontName = 'Arial';

Msize = 4; % Marker Size
Tsize = [0.03, 0.06]; % Tick Size

Axis_Line = 1; % Axis line width
Plot_Line = 1; % Plot line width

% Make the figure
tmpFig = figure('Visible', 'off', 'Units', 'Centimeters','PaperPositionMode','auto');
oldPos = get(tmpFig, 'Position');
set(tmpFig, 'Position', [oldPos(1), oldPos(2), FigSize(1), FigSize(2)]); % make the figure bigger than needed (300x300)


%% Make the figure

% To track the axes position
PosTrack = NaN(nPlot+1,1);
PosTrack(1) = 0;

for ii = 1:nPlot
    
    newAxes=copyobj(PlotAxes{ii}, tmpFig);
    axis(newAxes, 'square');
    
    % Adjust the figure
    set(newAxes, 'FontUnits', 'Points', 'FontSize', FontSize2, 'Units', 'Centimeters', 'FontName', FontName)
    set(get(newAxes, 'XLabel'), 'FontUnits', 'Points', 'FontSize', FontSize3, 'FontName', FontName, 'Color', [0,0,0])
    set(get(newAxes, 'YLabel'), 'FontUnits', 'Points', 'FontSize', FontSize3, 'FontName', FontName, 'Color', [0,0,0]);
    set(get(newAxes, 'Title'), 'FontUnits', 'Points', 'FontSize', FontSize4, 'FontName', FontName, 'Color', [0,0,0]);
    
    % Delete the title
    % Temporary fix for HytLab until better adjustemnts are made
    set(get(newAxes, 'Title'), 'String', '');
    
    % Readjust the x-axis scale and tickmarks
    set(newAxes, 'Xlim', get(PlotAxes{ii}, 'Xlim'))
    set(newAxes, 'XTick', get(PlotAxes{ii}, 'XTick'))
    set(newAxes, 'XTickLabel', get(PlotAxes{ii}, 'XTick'), 'FontName', FontName);%, 'Color', [0,0,0])
    
        
    % Adjust size
    if ii==1
        NewPos = [1.5, 1.5, PlotSize, PlotSize];
    else
        NewPos = [PosTrack(ii) + 1.75 + PlotSize, 1.5, PlotSize, PlotSize];
    end
    
    set(newAxes, 'Position', NewPos, 'XColor', [0,0,0], 'YColor', [0,0,0], 'Box', 'off', 'TickDir', 'Out', 'TickLength', Tsize, 'LineWidth', Axis_Line);
    
    % Place a new set of axes on top to create the box
    h0 = axes('Units', 'Centimeters', 'Position', NewPos);
    set(h0, 'box', 'on', 'XTick', [], 'YTick', [], 'color', 'none', 'XColor', [0,0,0], 'YColor', [0,0,0], 'LineWidth', Axis_Line);
    
    PosTrack(ii+1) = NewPos(1);
    
    % Reset the line widths and any text handles
    C = get(newAxes, 'Children');
    for jj = 1: length(C)
        
        tmp_hand = get(C(jj));
        
        if isfield(tmp_hand, 'LineWidth')
            set(C(jj),'LineWidth',1);
        end
        
        % Set the text label font size
        if isfield(tmp_hand, 'String')
            set(C(jj), 'FontUnits', 'Points', 'FontSize', FontSize1, 'FontName', FontName);
        end
        
    end

    
end


