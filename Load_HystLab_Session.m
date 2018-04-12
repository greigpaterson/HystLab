function handles = Load_HystLab_Session(handles, Session_handles)
%
% Function to load saved HystLab sessions
%

%%

fNames = fieldnames(Session_handles);
nFields = length(fNames);

for ii = 1:nFields
    handles.(fNames{ii}) = Session_handles.(fNames{ii});
end


if ~isfield(Session_handles, 'HystLab_Version')
    warning('Load_HystLab_Session:Legacy_Version', 'The loaded session is from a VERY old version of HystLab and IS NOT fully compatible with the latest release.');
    handles.Session_Fail = 1;
    return
end

Current_Version = handles.HystLab_Version; %#ok<NASGU>
Load_Version = Session_handles.HystLab_Version;
Ver_Num = str2double(strjoin(strsplit(Load_Version, '.'), ''));

%% Do some version checking
% this will alos run some updates to the loaded session to ensure thay are
% compatible with the latest version


if Ver_Num < 081
    
    % Does not include paramagnetic drift correction
    % Need to add Temp_ratio into the processing parameters
    
    % Sound a warning to the command line
    warning('Load_HystLab_Session:Legacy_Version', 'The loaded session is from an older version of HystLab and may not be fully compatible with the latest release.');
    
    try
        nData = size(handles.All_Data,1);
        
        % Make a default neasure order
        handles.All_Data_Order = ones(nData,1);
        
        % Add in the Temp_ratio parameters
        tmp_Proc = handles.All_Processing_Parameters(:,1:3);
        tmp_Proc(:,4) = NaN;
        tmp_Proc = [tmp_Proc, handles.All_Processing_Parameters(:,4:end)];
        handles.All_Processing_Parameters = tmp_Proc;
        
    catch
        
        % Give a warning if something failed
        MSG = [{'You are trying to load a pre-v0.8.1 session file.'}; {'HystLab encountered an error loading this session.'};...
            {'If data need to be recovered, please contact the authors, or use an older version of HystLab to read the session.'}];
        warndlg(MSG, 'Legacy version not supported');
        handles.Session_Fail = 1;
        return;
        
    end
    
end


if Ver_Num < 082
    
    % The matix for the fitted loop data was updated to match the format of
    % the main loop data matrix
    % Need to restrtuct the fitted data for compatibility
    
    warning('Load_HystLab_Session:Legacy_Version', 'The loaded session is from an older version of HystLab and may not be fully compatible with the latest release.');
    
    
    try
        nData = size(handles.All_Fitted_Data,1);
        
        % loop through the fitted data and update
        for ii = 1:nData
            
            old_data = handles.All_Fitted_Data{ii};
            
            new_data = [old_data(:,1), flipud(old_data(:,1)), old_data(:,2), flipud(old_data(:,3)),...
                old_data(:,4), old_data(:,5)];
            
            handles.All_Fitted_Data{ii} = new_data;
        end
        
        % Do the currently selected specimen
        
        old_data = handles.Current_Fitted_Data;
        
        new_data = [old_data(:,1), flipud(old_data(:,1)), old_data(:,2), flipud(old_data(:,3)),...
            old_data(:,4), old_data(:,5)];
        
        handles.Current_Fitted_Data = new_data;
        
    catch
        
        % Give a warning if something failed
        MSG = [{'You are trying to load a pre-v0.8.2 session file.'}; {'HystLab encountered an error loading this session.'};...
            {'If data need to be recovered, please contact the authors, or use an older version of HystLab to read the session.'}];
        warndlg(MSG, 'Legacy version not supported');
        handles.Session_Fail = 1;
        return;
        
    end
    
    % Added fit RMS to the end of the data parameters
    nData = size(handles.All_Data_Parameters,1);
    
    tmp_data = [handles.All_Data_Parameters, NaN(nData, 1)];
    
    handles.All_Data_Parameters = tmp_data;
    
end


if Ver_Num < 083
    
    % No test for high-field zero Mrh    
    try
        % add in 8 additional values for the test stats
        nData = size(handles.All_Data_Parameters,1);
        
        tmp_data = [handles.All_Data_Parameters, NaN(nData, 8)];
        
        handles.All_Data_Parameters = tmp_data;
        
    catch
        
        % Give a warning if something failed
        MSG = [{'You are trying to load a pre-v0.8.3 session file.'}; {'HystLab encountered an error loading this session.'};...
            {'If data need to be recovered, please contact the authors, or use an older version of HystLab to read the session.'}];
        warndlg(MSG, 'Legacy version not supported');
        handles.Session_Fail = 1;
        return;
        
    end
    
    
    
end


if Ver_Num < 100
% Will consider using version 1.0.0 as a cut off for backward compatibility
% All pre-v1.0.0 session will be from test and development code only
%The authors should be the only users severly affected
%     handles = [];
%     MSG = [{'You are trying to load a pre-v1.0.0 session file.'}; {'This is not supported in the current version of HystLab.'};...
%         {'If data need to be recovered, please contact the authors, or use an older version of HystLab to read the session.'}];
%     warndlg(MSG, 'Legacy version not supported');
%     handles.Session_Fail = 1;
%     return;
end


% Flag for checking failure
% All is good
handles.Session_Fail = 0;

