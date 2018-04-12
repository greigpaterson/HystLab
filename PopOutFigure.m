function PopOutFigure(AxisHandle, FigName)
%
% Creates an external figure of the selected axes
% This is design to be called from a buttondown function of an axes object
%
% Input:
%        AxisHandle - the hanlde of the axes to copy
%        FigName - the name of the figure to create or update if it already
%                  exists
%

%% Check to see if the figure exists

OpenFigHandles = findall(0,'Type','figure');

if isempty(OpenFigHandles)
    FigExist = 0;
else
    OpenFigNames = get(OpenFigHandles, 'name');
    FigInd = find(strcmp(OpenFigNames, FigName));
    
    if isempty(FigInd)
        FigExist = 0;
    else
        FigExist = FigInd;
    end
end

%% Find the Figure to plot to

if FigExist == 0
    
    % Make a new figure
    tmpFig = figure('Name', FigName, 'Visible', 'on', 'Units', 'Pixels','PaperPositionMode','auto');
    
    % Set the postion relative to the parent window
    pPos = get(get(AxisHandle, 'Parent'), 'Position');
    xPos = round(pPos(1) + pPos(3)./2 - 250);
    yPos = round(pPos(2) + pPos(4)./2 - 250);
    
    set(tmpFig, 'Position', [xPos, yPos, 500, 500]);
    
    % Set position to defaul new fig position
    %         oldPos = get(tmpFig, 'Position');
    %     set(tmpFig, 'Position', [oldPos(1), oldPos(2), 500, 500]);
    
    NewPos = [60, 60, 400, 400];
    Units = 'Pixels';
    
else
    
    % Select the existing one
    tmpFig = OpenFigHandles(FigInd);
    
    % Get some info from the exisitng axes
    C = get(tmpFig, 'Children');
    
    if length(C) ~=1
        % User has created different objects
        C = findall(C, 'type', 'axes', 'Tag', FigName);
    end
    
    NewPos = get(C, 'Position');
    Units = get(C, 'Units');
    
    clf(tmpFig); % clear the figure
    figure(tmpFig); % bring it to the front
end

%% Make the new plot

newAxes = copyobj(AxisHandle, tmpFig);

% Adjust the axes and reset the units to normalized for scaling
set(newAxes, 'ButtonDownFcn', '', 'Units', Units, 'Position', NewPos);
set(newAxes, 'Units', 'Normalized', 'Tag', FigName);

