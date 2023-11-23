%% Script to run in parallel multiple PCA-CCA of psychosis-related diagnoses

% Copyright (C) 2023 University of Seville

% Written by Natalia García San Martín (ngarcia1@us.es)

% This file is part of Neurobiological Centiles Mapping toolkit.
%
% Neurobiological Centiles Mapping toolkit is free software: 
% you can redistribute it and/or modify it under the terms of the 
% GNU General Public License as published by the Free Software Foundation, 
% either version 3 of the License, or (at your option) any later version.
%
% Neurobiological Centiles Mapping toolkit is distributed in the hope that 
% it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Neurobiological Centiles Mapping toolkit. If not, see 
% <https://www.gnu.org/licenses/>.

clear
close all

% Change as desired
location = 'D:\OneDrive - UNIVERSIDAD DE SEVILLA\Natalia\Code\Data';
location = 'C:\Users\usuario\OneDrive - UNIVERSIDAD DE SEVILLA\Natalia\Code\Data';
directory = [location,'\PCA_CCA\'];

dx = {'SCZ_Relative','Schizoaffective_Disorder_Relative','PE_1','PE_2','PE_3'...
    'FEP_1','SCZ','SCZ_without_ASRB','SCZ_without_BSNIP','SCZ_without_LA5c','SCZ_without_MCIC','Schizoaffective_Disorder'};

name = {};

% poolobj = gcp('nocreate');
% delete(poolobj);
% parpool(length(dx))

for i=1:length(dx)
    dx_title = dx{i};

    % Create subfolders
    if ~exist(fullfile(directory,dx_title),'dir')
       mkdir([directory,dx_title])
    end
    if ~exist(fullfile([directory,dx_title],'data'),'dir')
        mkdir([directory,dx_title],'data')
    end
    
    [pval_spins, molecular_names, weights, Y_pred, molecular_maps_labels_ordered] = TEN_CCA_cent_var(name,dx_title,location,'real');
    weights_total(i,:) = weights;
    Y_pred_total(i,:) = Y_pred;
 
end

i_significative = any(weights_total~=0,2); % same for weight_total
weights_significative = weights_total(i_significative,:);

weights_table_dx = array2table(weights_significative,'RowNames',dx(i_significative)','VariableNames',molecular_names);

Y_pred_table_dx = array2table(Y_pred_total(any(Y_pred_total~=0,2),:),'RowNames',dx(any(Y_pred_total~=0,2))','VariableNames',molecular_maps_labels_ordered);


molecular_names_table = readtable([location,'\Molecular_maps\molecular_names.xlsx'],ReadVariableNames=true,ReadRowNames=true);
molecular_types_orig = molecular_names_table{:,2};
weights_table_dx =  rows2vars(weights_table_dx);
weights_table_dx.molecular_types(1:46) = molecular_types_orig;
weights_table_dx = sortrows(weights_table_dx,{'molecular_types','OriginalVariableNames'});
weights_table_dx =  rows2vars(weights_table_dx(:,1:end-1),"VariableNamesSource","OriginalVariableNames");
weights_table_dx.Properties.RowNames = weights_table_dx{:,1};
weights_table_dx = weights_table_dx(:,2:end);

weights_table_dx.Properties.VariableNames(7) = {'α_4β_2'};
weights_table_dx.Properties.VariableNames(1) = {'5-HT_{1A}'};
weights_table_dx.Properties.VariableNames(2) = {'5-HT_{1B}'};
weights_table_dx.Properties.VariableNames(3) = {'5-HT_{2A}'};
weights_table_dx.Properties.VariableNames(4) = {'5-HT_4'};
weights_table_dx.Properties.VariableNames(5) = {'5-HT_6'};
weights_table_dx.Properties.VariableNames(6) = {'5-HTT'};
weights_table_dx.Properties.VariableNames(8) = insertBefore(weights_table_dx.Properties.VariableNames(8), 3, '_');
weights_table_dx.Properties.VariableNames([9,10,13,14]) = insertBefore(weights_table_dx.Properties.VariableNames([9,10,13,14]), 2, '_');
weights_table_dx.Properties.VariableNames(12) = {'GABA'};
weights_table_dx.Properties.VariableNames(19) = insertBefore(weights_table_dx.Properties.VariableNames(19), 6, '_');
weights_table_dx.Properties.VariableNames([23,24]) = replace(weights_table_dx.Properties.VariableNames([23,24]),'_',' ');
weights_table_dx.Properties.VariableNames(27) = {'Layer I'};
weights_table_dx.Properties.VariableNames(28) = {'Layer II'};
weights_table_dx.Properties.VariableNames(29) = {'Layer III'};
weights_table_dx.Properties.VariableNames(30) = {'Layer IV'};
weights_table_dx.Properties.VariableNames(31) = {'Layer V'};
weights_table_dx.Properties.VariableNames(32) = {'Layer VI'};
weights_table_dx.Properties.VariableNames(33) = {'Gene PC1'};
weights_table_dx.Properties.VariableNames(34) = {'Myelin'};
weights_table_dx.Properties.VariableNames(35) = {'Neurotransmitter PC1'};
weights_table_dx.Properties.VariableNames(36) = {'Synapse density'};
weights_table_dx.Properties.VariableNames(37) = {'Thickness'};
weights_table_dx.Properties.VariableNames(38) = {'Developmental exp.'};
weights_table_dx.Properties.VariableNames(39) = {'Evolutionary exp.'};
weights_table_dx.Properties.VariableNames(40) = {'Scaling NIH'};
weights_table_dx.Properties.VariableNames(41) = {'Scaling PNC'};
weights_table_dx.Properties.VariableNames(42) = {'CBF'};
weights_table_dx.Properties.VariableNames(43) = {'CBV'};
weights_table_dx.Properties.VariableNames(44) = {'CMRO_2'};
weights_table_dx.Properties.VariableNames(45) = {'CMRGlu'};
weights_table_dx.Properties.VariableNames(46) = {'Glycolytic index'};

if strcmp(gca().YLabel.String,'Weight')
    save([directory,'weights_dx.mat'],"weights_table_dx")
elseif strcmp(gca().YLabel.String,'Loading')
    loadings_table_dx = weights_table_dx;
    save([directory,'loadings_dx.mat'],"loadings_table_dx")
end

writetable(Y_pred_table_dx,[directory,'Y_pred_dx.csv'],'WriteRowNames',true)
