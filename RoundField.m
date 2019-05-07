function OutField = RoundField(InField)
%
% Function to return a field increment
% If the field is < 1 , round the first non-zero digit down, otherwise
% round down to the 1
%
% Last Modified 2019/05/07
%

%%

if InField == 0
    error('RoundField:Zero', 'Input field step is zero. This is invalid.')
end

Step_Str = num2str(InField);
S = regexp(Step_Str, '\.', 'split');


if str2double(S{1}) == 0
    
    % Loop through each character to find the first non-zero digit
    nChar = length(S{2});
    
    for ii = 1:nChar
        
        if str2double(S{2}(ii)) ~= 0
            break;
        end
        
    end
    
    OutField = fix(InField.*10^ii) / 10^ii;
    
else
    
    % Leading number is not zero
    % Round to nearest 1 mT
    OutField = fix(InField);
    
end