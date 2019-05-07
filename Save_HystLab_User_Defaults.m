function Save_HystLab_User_Defaults(Defaults)
%
% Saves user specified default paramters for HystLab
%
% Last Modified 2019/05/07
%

%%

% Get the main path
S = mfilename('fullpath');
name_len = length(mfilename());
MyPath = S(1:end-name_len);


% Do the color defaults first
Hyst_Plot = Defaults.Hyst_Plot_Color; %#ok<*NASGU>
IR_Plot = Defaults.IR_Plot_Color;
Noise_Plot = Defaults.Noise_Plot_Color;
BiPlot = Defaults.BiPlot_Color;


% Save the color file and remove the fields from Defaults
save(strcat(MyPath, 'HystLabUserColorDefaults.mat'), 'Hyst_Plot', 'IR_Plot', 'Noise_Plot', 'BiPlot');

Defaults = rmfield(Defaults, {'Hyst_Plot_Color', 'IR_Plot_Color', 'Noise_Plot_Color', 'BiPlot_Color'});


% Get the text based settings
FID = fopen( strcat(MyPath, 'HystLabUserDefaults.cfg'), 'wt');
fprintf(FID, '%s\n', 'PlotColorFile = HystLabUserColorDefaults.mat');

FN = fieldnames(Defaults);
nFields = length(fieldnames(Defaults));

for ii=1:nFields
    fprintf(FID, '%s\n', [FN{ii}, ' = ', num2str(Defaults.(FN{ii})) ] );
end

fclose(FID);

msgbox('Default values saved.', 'modal')
