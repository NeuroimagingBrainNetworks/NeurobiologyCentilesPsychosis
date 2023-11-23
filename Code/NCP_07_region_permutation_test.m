%% Script to test the brain regions that exhibit significant differences from HC (centiles; Wilcoxon rank sum test),
%% or are significantly different from 0 (effect sizes, permutation test) with FDR correction 

% Copyright (C) 2023 University of Seville

% Written by Natalia García San Martín (ngarcia1@us.es)

% This file is part of Neurobiology Centiles Psychosis toolkit.
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
location = 'D:\OneDrive - UNIVERSIDAD DE SEVILLA\Natalia\Code\Data\';
location = 'C:\Users\usuario\OneDrive - UNIVERSIDAD DE SEVILLA\Natalia\Code\Data\';
cd(location); 

% Script to calculate pvalues for each brain region

% For centiles
load("FEP.mat",'-mat');
load("PE.mat",'-mat');
load("SCZ.mat",'-mat');
load("Schizoaffective_Disorder.mat",'-mat')
load("G0.mat",'-mat');
load("G1.mat",'-mat');
% load("G1_5.mat",'-mat') % it is the same as centiles_FEP_1
load("G2.mat",'-mat');
load("G2_without_ASRB.mat",'-mat');
load("G2_without_BSNIP.mat",'-mat');
load("G2_without_LA5c.mat",'-mat');
load("G2_without_MCIC.mat",'-mat');

% Load CN
PE_CN = readtable([location,'Centiles\PE\centiles.csv']);pos_centiles=19:52;
PE_CN = PE_CN(strcmp(PE_CN.dx,'CN'),pos_centiles);
SCZ_CN = readtable([location,'Centiles\SCZ_SAD\centiles.csv']);pos_centiles=10:43;
SCZ_CN = SCZ_CN(strcmp(SCZ_CN.dx,'CN') & ismember(SCZ_CN.study,{'ABCD','ASRB','BSNIP','LA5c','MCIC','UKB'}),pos_centiles);
FEP_CN = readtable([location,'Centiles\FEP\centiles.csv']);pos_centiles=19:52;
FEP_CN = FEP_CN(strcmp(FEP_CN.dx,'CN'),pos_centiles);

cent_names = {'centiles_FEP','centiles_SCZ_and_SAD-relatives', 'centiles_PE','centiles_SCZ_and_SAD-chronic', 'centiles_SCZ_and_SAD-chronic_without_ASRB',...
    'centiles_SCZ_and_SAD-chronic_without_BSNIP','centiles_SCZ_and_SAD-chronic_without_LA5c','centiles_SCZ_and_SAD-chronic_without_MCIC','centiles_PE-suspected',...
    'centiles_PE-definite', 'centiles_PE-clinical','centiles_SCZ','centiles_SCZ-relatives','centiles_SCZ_without_ASRB', ...
    'centiles_SCZ_without_BSNIP','centiles_SCZ_without_LA5c','centiles_SCZ_without_MCIC','centiles_SAD','centiles_SAD-relatives'};

pval_cent = table();
pval_cent_global = table();
dx = who;
dx = dx(contains(dx,'centiles_'));
k = 1;
for i = 1:length(dx)
    centiles = table2array(eval(dx{i}));
    if contains(dx(i),'PE')
        group_CN = PE_CN;
    elseif contains(dx(i),'SCZ') || contains(dx(i),'Schizo')
        group_CN = SCZ_CN;
    elseif contains(dx(i),'FEP')
        group_CN = FEP_CN;
    elseif contains(dx(i),'G0')
        group_CN = SCZ_CN;
    elseif contains(dx(i),'G1') && ~contains(dx(i),'G1_5')
        group_CN = PE_CN;
    elseif contains(dx(i),'G1_5')
        group_CN = FEP_CN;
    elseif contains(dx(i),'G2')
        group_CN = SCZ_CN;
    else
        continue
    end
    
    % Wilcoxon rank sum to test the regions that exhibit significant differences from HC 
    for ir = 1:34
        [p(ir), h(ir) stats(ir)] = ranksum(centiles(:,ir),group_CN{:,ir});
    end

    pval_cent{k,:} = p; pval_cent.Properties.RowNames(k) = cent_names(i);
    pval_corr_cent(k,:) = mafdr(pval_cent{k,:},'BHFDR',true); % multipe comparisons correction (Benjamini-Hochberg False Discovery Rate)
    mean_centiles(k,:) = mean(centiles);
    mean_centiles(k,pval_corr_cent(k,:) >= 0.05) = NaN;

    % Wilcoxon rank sum to test the global volumes that exhibit significant differences from HC 
    [p_global, h_global] = ranksum(mean(centiles,2),mean(group_CN{:,:},2));
    pval_cent_global{k,:} = p_global; pval_cent_global.Properties.RowNames(k) = cent_names(i);
    volume_global_centiles = mean(mean(centiles,2));
    volume_global_CN = mean(mean(group_CN{:,:},2));

    k = k + 1;

end
pval_cent.Properties.VariableNames = SCZ_CN.Properties.VariableNames;
pval_corr_cent = array2table(pval_corr_cent,'RowNames',pval_cent.Properties.RowNames,'VariableNames',pval_cent.Properties.VariableNames);
mean_centiles = array2table(mean_centiles,'RowNames',pval_cent.Properties.RowNames,'VariableNames',pval_cent.Properties.VariableNames);

% Group selection
pval_corr_cent_global = mafdr(pval_cent_global{[1:3,7],:},'BHFDR',true); % multipe comparisons correction (Benjamini-Hochberg False Discovery Rate)
pval_corr_cent = array2table(pval_corr_cent_global,'RowNames',pval_cent_global.Properties.RowNames([1:3,7]));

% For effsizes
clear("centiles_G2_without_ASRB")
clear("centiles_G2_without_BSNIP")
clear("centiles_G2_without_LA5c")
clear("centiles_G2_without_MCIC")
load("G1_5.mat",'-mat')
load("G1_vs_G0.mat",'-mat');
load("G2_vs_G0.mat",'-mat');
load("G2_vs_G1.mat",'-mat');
load('G1_5_vs_G0','-mat')
load('G1_5_vs_G1','-mat')
load('G2_vs_G1_5','-mat')

groups = who;
groups = groups(contains(groups,'centiles_G')); groups = groups(length(groups):-1:1);

nperm = 1000;
k = 1;
for i = 1:length(groups)-1
    for j = i+1:length(groups)
        [name_group_1,name_group_2] = groups{[i,j]};
        centiles_dx_1 = eval(name_group_1);
        centiles_dx_2 = eval(name_group_2);
        effsizes_name = [replace(name_group_1,'centiles','effsizes'),'_vs_',replace(name_group_2,'centiles_','')];
        real_diff = table2array(eval(effsizes_name));
        for iperm = 1:nperm
            % Group membership is randomly reassigned across the compared
            % groups to test whether they are significantly different from 0
            [selection_dx_1,selection_dx_2,ratio_1(iperm),ratio_2(iperm)] = mix_dx(centiles_dx_1,centiles_dx_2);
            
            for ir = 1:34
                effsizes_perm(ir) = computeCohen_d(selection_dx_1(:,ir),selection_dx_2(:,ir));
            end

            perm_diff(iperm,:) = effsizes_perm;
        end
        
        % Test if it explains more variance that expected by chance
        pval_eff(k,:) = sum(abs(perm_diff) > abs(real_diff')) / nperm;
    
        pval_corr_eff(k,:) = mafdr(pval_eff(k,:),'BHFDR',true); % multipe comparisons correction (Benjamini-Hochberg False Discovery Rate)
        total_real_diff(k,:) = real_diff;        

        effsizes_names{k} = effsizes_name;
        k = k + 1;
    end
end
eff_corr = total_real_diff;
eff_corr(pval_corr_eff >= 0.05) = NaN;

pval_corr_eff = array2table(pval_corr_eff,"RowNames",effsizes_names,'VariableNames',pval_corr_cent.Properties.VariableNames);
effsizes_names = {'effsizes_SCZ_and_SAD-chronic_vs_FEP','effsizes_SCZ_and_SAD-chronic_vs_PE','effsizes_SCZ_and_SAD-chronic_vs_SCZ_and_SAD-relatives', ... 
    'effsizes_FEP_vs_PE','effsizes_FEP_vs_SCZ_and_SAD-relatives','effsizes_PE_vs_SCZ_and_SAD-relatives'};
eff_corr = array2table(eff_corr,'RowNames',effsizes_names,'VariableNames',mean_centiles.Properties.VariableNames);

p_values_corr = [pval_corr_cent;pval_corr_eff];
significant_values = [mean_centiles;eff_corr];

writetable(rows2vars(significant_values,"VariableNamingRule","preserve"),[location,'\significant_values.csv'],"WriteRowNames",true)




