function [Model_Fit] = Get_Fitted_Data(Fields, Basis, type, Params)
%
% Function to generate hysteresis, or IRM curves from a set of basis 
% function coefficients
%
% Input:
%       Fields - The field steps to evaluate the basis functions at [0 MaxField]
%       Basis - Cell array of hysteresis ([1x2] or [2x1]) or IRM [1x1] basis
%               coefficients. For hysteresis, Basis{1} is Mih and Basis{2} is Mrh.
%               Each cell is a [nBasis x 4] matrix. Column 1 is the basis
%               type (1 = sigmoid, 2 = hyperbolic, 3 = linear), Column 2 is
%               the realtive abundance of the basis, Columns 3 and 4 are
%               the basis fucntion coefficients.
%       type - String for the type of basis data to generate ('Hys' or 'IRM')
%       Params - Parameters for scaling the model fits. For type 'Hys' this
%                is [Ms, Mrs], for 'IRM' this is [Mrs]
%
% Output:
%       Mode_Fit - the sum of the basis functions evaluated at the given
%                  fields
%                  For hysteresis, [Fields(U,L), Upper, Lower, Mih, Mrh]
%                  For IRM, [Fields, IRM]
%
% Last Modified 2019/05/07
%

%% Some input checking

if nargin < 3
    error('Get_Fitted_Data:Input', 'At least 3 input arguments are required.');
end

% Check the size of Fields
nFields = size(Fields,1);

if nFields == 1
    Fields = Fields';
end


%% The many function
switch type
    
    case 'Hys'
        
        % Check we have the correct number of cells
        if sum(size(Basis) == 2) ~= 1
            error('Get_Fitted_Data:Hys_Basis', 'Missing induced or remanence basis coefficients.');
        end
        
        if sum(size(Params) == 2) ~=1
            error('Get_Fitted_Data:Hys_Params', 'Incorrect number of input parameters. Should be [Ms, Mrs].');
        end
        
%         if Params(2) > Params(1)
%             error('Get_Fitted_Data:Hys_Params', 'Mrs > Ms. Parameter order should be [Ms, Mrs].');
%         end

        % Get the bases and the number of them
        IH_Basis = Basis{1};
        RH_Basis = Basis{2};
        
        nIH = size(IH_Basis,1);
        nRH = size(RH_Basis,1);
        
        Pos_Fields = sort(Fields(Fields>=0));
        nFields = length(Pos_Fields);
        
        %% Get the IH curve
        Mih_Basis = NaN(nIH, nFields);
        
        for ii = 1:nIH
            
            Flag = IH_Basis(ii,1);
            C1 = IH_Basis(ii,3); % b_ih  or t
            C2 = IH_Basis(ii,4); % B_ih_m  or MaxField
            
            switch Flag
                case 1
                    % Sigmoidal basis
                    tmp_basis = 1 ./ (1 + exp(-C1.*(Pos_Fields - C2) ) ); % A = 0, C = 1
                case 2
                    % Hyperbolic tangent basis
                    tmp_basis = tanh(C1.*Pos_Fields) ./ tanh(C1.*C2);
                case 3
                    % Linear basis
                    tmp_basis = linspace(0, 1, nFields);
                case -3
                     % Linear basis
                    tmp_basis = linspace(1, 0, nFields);
                otherwise
                    error('Get_Fitted_Data:BasisType', 'Wrong basis flag, use either 1,2, or 3 not %s', Flag);
            end
            
            tmp_basis = tmp_basis./sum(tmp_basis);
            
            Mih_Basis(ii,:) = tmp_basis';
            
        end
        
        % Get the scaled Mrh curve
        Mih = Params(1) .* IH_Basis(:,2)'*Mih_Basis;        
        
        
        %% Get the RH curve
        
        % Check we have the correct number of cells
        if sum(size(Basis) == 1) ~= 1
            error('Get_Fitted_Data:IRM_Basis', 'Incorrect number of input basis coefficients.');
        end
        
        if sum(size(Params) == 1) ~=1
            error('Get_Fitted_Data:IRM_Params', 'Incorrect number of input parameters. Should be Mrs only.');
        end
        
        Mrh_Basis = NaN(nRH, nFields);
        
        for ii = 1:nRH
            
            Flag = RH_Basis(ii,1);
            C1 = RH_Basis(ii,3);
            C2 = RH_Basis(ii,4);
            
            switch Flag
                case 1
                    % Sigmoidal basis
                    tmp_basis = 1 + -1./ (1 + exp(-C1.*(Pos_Fields - C2) ) );  % A = 1, C = 0
                case 2
                    % Hyperbolic secant basis
                    tmp_basis = (sech(C1.*Pos_Fields)-sech(C1.*C2))./(1-sech(C1.*C2));
                case 3
                    % Linear basis
                    tmp_basis = linspace(1, 0, nFields);
                case -3
                    % Linear basis
                    tmp_basis = linspace(0, 1, nFields);
                case 4
                    % Constant basis function
                    tmp_basis = ones(1, nFields);
                otherwise
                    error('Get_Basis_Data:BasisType', 'Wrong basis flag, use either 1,2, or 3 not %s', Flag);
            end
            
            tmp_basis = tmp_basis./sum(tmp_basis);
            Mrh_Basis(ii,:) = tmp_basis';
            
        end
        
        % Get the scaled Mrh curve
        Mrh = Params(2) .*RH_Basis(:,2)'*Mrh_Basis;
        
        
        %% Get the modeled data to return
        
        if any(Fields == 0)
            % Stop the zero field step being doubled
            Mv = [[fliplr(-Mih)'; Mih(2:end)'], [fliplr(Mrh)'; Mrh(2:end)'] ];
        else
            Mv = [[fliplr(-Mih)'; Mih(1:end)'], [fliplr(Mrh)'; Mrh(1:end)'] ];
        end
        
        Top_Curve = Mv(:,1) + Mv(:,2);
        Bot_Curve = Mv(:,1) - Mv(:,2);
        
        
        % Flip to keep the correct field ordering
         Mv = flipud(Mv);
        
        if any(Fields == 0)
            % Stop the zero field step being doubled
            Full_Fields = [flipud(-Pos_Fields); Pos_Fields(2:end)];
        else
            Full_Fields = [flipud(-Pos_Fields); Pos_Fields(1:end)];
        end
        
        Model_Fit = [ flipud(Full_Fields), Full_Fields, flipud(Top_Curve), Bot_Curve, Mv ];
        
        
    case 'IRM'
        
        error('Get_Fitted_Data:IRM', 'Not yet supported.');
        
        % Get the bases and the number of them
%         IRM_Coeffs = Basis{1};
%         nIH = size(IRM_Coeffs,1);
%         
%         IRM_Basis = NaN(nIH, nFields);
%         
%         for ii = 1:nIRM
%             
%             Flag = IRM_Coeffs(ii,1);
%             C1 = IRM_Coeffs(ii,3); % b_ih  or t
%             C2 = IRM_Coeffs(ii,4); % B_ih_m  or MaxField
%             
%             switch Flag
%                 case 1
%                     % Sigmoidal basis
%                     tmp_basis = 1 ./ (1 + exp(-C1.*(Fields - C2) ) ); % A = 0, C = 1
%                 case 2
%                     % Hyperbolic tangent basis
%                     tmp_basis = tanh(C1.*Fields) ./ tanh(C1.*C2);
%                 case 3
%                     % Linear basis
%                     tmp_basis = linspace(0, 1, nFields);
%                 otherwise
%                     error('Get_Basis_Data:BasisType', 'Wrong basis flag, use either 1,2, or 3 not %s', Flag);
%             end
%             
%             tmp_basis = tmp_basis./sum(tmp_basis);
%             
%             IRM_Basis(ii,:) = tmp_basis';
%             
%         end
%         
%         IRM = IRM_Coeffs(:,2)'*IRM_Basis;
%         
%         
%         Model_Fit = [Fields, Param.*IRM];
        
    otherwise
        error('Get_Fitted_Data:Type', 'Unrecognized data type.');
end

