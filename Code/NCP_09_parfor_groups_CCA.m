%% Script to run in parallel multiple PCA-CCA of psychosis-related groups 
%% classified according to their clinical profile

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

groups = {'G0','G1','G1_5','G2','G2_vs_G0','G1_vs_G0','G2_vs_G1','G2_vs_G1_5','G1_5_vs_G1','G1_5_vs_G0'};
groups = {'G2','G2_without_ASRB','G2_without_BSNIP','G2_without_LA5c','G2_without_MCIC'};
% G0 = relatives of SCZ and SAD
% G1 = PE suspected, definite and clinical
% G1_5 = FEP
% G2 = chronic SCZ and SAD

name = {};

% poolobj = gcp('nocreate');
% delete(poolobj);
% parpool(length(groups))

correl_perm_tot = [];
weights_perm_tot = [];
for i=1:length(groups)
    dx_title = groups{i};

    % Create subfolders
    if ~exist(fullfile(directory,dx_title),'dir')
       mkdir([directory,dx_title])
    end
    if ~exist(fullfile([directory,dx_title],'data'),'dir')
        mkdir([directory,dx_title],'data')
    end
    if exist(fullfile([directory,dx_title],'Copy_of_framework'),'dir')
        if exist(fullfile([directory,dx_title],'framework'),'dir')
            rmdir(fullfile([directory,dx_title],'framework'),'s')
        end
        movefile(fullfile([directory,dx_title],'Copy_of_framework\'),fullfile([directory,dx_title],'framework\'))
        pause(5)
    end
    
    [pval_spins, molecular_names, weights_real, Y_pred, correl_real, molecular_maps_labels_ordered] = NCP_10_CCA_cent_var(name,dx_title,location,'real');
    weights_total(i,:) = weights_real;
    Y_pred_total(i,:) = Y_pred;
    save('correl_real.mat',"correl_real")
    save('weights_real.mat',"weights_real")

    % Uncomment for permutation test for significant effect sizes (all except those including relatives (G0)) 
%     if contains (dx_title,'vs') && ~contains (dx_title,'G0')
%         movefile('framework\','Copy_of_framework\')
%         nperm = 1000;
%         for iperm = 1:nperm
%             [~, ~, weight_perm, ~, correl_perm, ~] = NCP_10_CCA_cent_var(name,dx_title,location,'permutation');
%             correl_perm_tot = [correl_perm_tot;correl_perm];
%             weights_perm_tot = [weights_perm_tot; weight_perm];
%             rmdir ('framework','s')
%         end
%         save('correls_perm.mat',"correl_perm_tot")
%         save('weights_perm.mat','weights_perm_tot')
%         pval_correl = sum(correl_perm_tot > correl_real)/nperm;
%         pval_weights = sum(abs(weights_perm_tot)> abs(weights_real))/nperm;
%     
%         movefile('Copy_of_framework\','framework\')
%     end
end

i_significative = any(weights_total~=0,2); % same for weight_total
weights_significative = weights_total(i_significative,:);

weights_table_groups = array2table(weights_significative,'RowNames',groups(i_significative)','VariableNames',molecular_names);

Y_pred_table_groups = array2table(Y_pred_total(any(Y_pred_total~=0,2),:),'RowNames',groups(any(Y_pred_total~=0,2))','VariableNames',molecular_maps_labels_ordered);


molecular_names_table = readtable([location,'\Molecular_maps\molecular_names.xlsx'],ReadVariableNames=true,ReadRowNames=true);
molecular_types_orig = molecular_names_table{:,2};
weights_table_groups =  rows2vars(weights_table_groups);
weights_table_groups.molecular_types(1:46) = molecular_types_orig;
weights_table_groups = sortrows(weights_table_groups,{'molecular_types','OriginalVariableNames'});
weights_table_groups =  rows2vars(weights_table_groups(:,1:end-1),"VariableNamesSource","OriginalVariableNames");
weights_table_groups.Properties.RowNames = weights_table_groups{:,1};
weights_table_groups = weights_table_groups(:,2:end);

weights_table_groups.Properties.VariableNames(7) = {'α_4β_2'};
weights_table_groups.Properties.VariableNames(1) = {'5-HT_{1A}'};
weights_table_groups.Properties.VariableNames(2) = {'5-HT_{1B}'};
weights_table_groups.Properties.VariableNames(3) = {'5-HT_{2A}'};
weights_table_groups.Properties.VariableNames(4) = {'5-HT_4'};
weights_table_groups.Properties.VariableNames(5) = {'5-HT_6'};
weights_table_groups.Properties.VariableNames(6) = {'5-HTT'};
weights_table_groups.Properties.VariableNames(8) = insertBefore(weights_table_groups.Properties.VariableNames(8), 3, '_');
weights_table_groups.Properties.VariableNames([9,10,13,14]) = insertBefore(weights_table_groups.Properties.VariableNames([9,10,13,14]), 2, '_');
weights_table_groups.Properties.VariableNames(12) = {'GABA'};
weights_table_groups.Properties.VariableNames(19) = insertBefore(weights_table_groups.Properties.VariableNames(19), 6, '_');
weights_table_groups.Properties.VariableNames([23,24]) = replace(weights_table_groups.Properties.VariableNames([23,24]),'_',' ');
weights_table_groups.Properties.VariableNames(27) = {'Layer I'};
weights_table_groups.Properties.VariableNames(28) = {'Layer II'};
weights_table_groups.Properties.VariableNames(29) = {'Layer III'};
weights_table_groups.Properties.VariableNames(30) = {'Layer IV'};
weights_table_groups.Properties.VariableNames(31) = {'Layer V'};
weights_table_groups.Properties.VariableNames(32) = {'Layer VI'};
weights_table_groups.Properties.VariableNames(33) = {'Gene PC1'};
weights_table_groups.Properties.VariableNames(34) = {'Myelin'};
weights_table_groups.Properties.VariableNames(35) = {'Neurotransmitter PC1'};
weights_table_groups.Properties.VariableNames(36) = {'Synapse density'};
weights_table_groups.Properties.VariableNames(37) = {'Thickness'};
weights_table_groups.Properties.VariableNames(38) = {'Developmental exp.'};
weights_table_groups.Properties.VariableNames(39) = {'Evolutionary exp.'};
weights_table_groups.Properties.VariableNames(40) = {'Scaling NIH'};
weights_table_groups.Properties.VariableNames(41) = {'Scaling PNC'};
weights_table_groups.Properties.VariableNames(42) = {'CBF'};
weights_table_groups.Properties.VariableNames(43) = {'CBV'};
weights_table_groups.Properties.VariableNames(44) = {'CMRO_2'};
weights_table_groups.Properties.VariableNames(45) = {'CMRGlu'};
weights_table_groups.Properties.VariableNames(46) = {'Glycolytic index'};

if strcmp(gca().YLabel.String,'Weight')
    save([directory,'weights_groups.mat'],"weights_table_groups")
elseif strcmp(gca().YLabel.String,'Loading')
    loadings_table_groups = weights_table_groups;
    save([directory,'loadings_groups.mat'],"loadings_table_groups")
end

writetable(Y_pred_table_groups,[directory,'Y_pred_groups.csv'],'WriteRowNames',true)
