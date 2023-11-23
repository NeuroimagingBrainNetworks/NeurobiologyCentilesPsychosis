%% Script to calculate the Sum of Squared Diferences between diagnoses and groups (relatives, PE, FEP and chronic)

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
location = 'D:\OneDrive - UNIVERSIDAD DE SEVILLA\Natalia\Code\Data\';
location = 'C:\Users\usuario\OneDrive - UNIVERSIDAD DE SEVILLA\Natalia\Code\Data\';
cd(location); 

load("PE.mat",'-mat');
load("FEP.mat",'-mat');
load("SCZ.mat",'-mat');
load("Schizoaffective_Disorder.mat",'-mat')
dx_names = {'SCZ-relatives','SAD-relatives','PE-suspected','PE-definite','PE-severe', 'FEP','SCZ','SAD'};

names_group={'SCZ_Relative','Schizoaffective_Disorder_Relative','PE_1','PE_2','PE_3', ...
    'FEP_1','SCZ','Schizoaffective_Disorder'};

nperm = 10000;
sse_real_total = nan(length(names_group));
sse_real_total = array2table(sse_real_total,'RowNames',names_group,'VariableNames',names_group);
PCT_pval = nan(length(names_group));
PCT_pval = array2table(PCT_pval,'RowNames',names_group,'VariableNames',names_group);
for i = 1:1:length(names_group)-1
    for j = i+1:length(names_group)
        [names_dx_1,names_dx_2] = names_group{[i,j]};
        centiles_dx_1 = eval(['centiles_',names_dx_1]);
        centiles_dx_2 = eval(['centiles_',names_dx_2]);

        % SSE
        sse_real = mse(mean(table2array(centiles_dx_1)),mean(table2array(centiles_dx_2)))*size(centiles_dx_1,2);

        for iperm = 1:nperm
            % Group membership is randomly reassigned across the compared diagnoses
            [selection_dx_1,selection_dx_2,ratio_1(iperm),ratio_2(iperm)] = mix_dx(centiles_dx_1,centiles_dx_2);
            
            selection_dx_1_m = mean(selection_dx_1);
            selection_dx_2_m = mean(selection_dx_2);

            % Permuted SSE
            sse_perm(iperm) = mse(selection_dx_1_m,selection_dx_2_m)*size(selection_dx_1,2);
    
        end 

        % Test if it explains more variance that expected by chance
        PCT_pval{i,j} = sum(sse_perm > sse_real) / nperm; 
        sse_real_total{i,j} = sse_real;

    end
end

figure;
sse_real_total_tri = triu(sse_real_total{:,:}) + triu(sse_real_total{:,:},1).';
heatmap(sse_real_total_tri,'XData',strrep(dx_names,'_',' '),'YData',strrep(dx_names,'_',' '));
title('SSE (centiles)')
figure;
PCT_pval_tri = triu(PCT_pval{:,:}) + triu(PCT_pval{:,:},1).';
heatmap(PCT_pval_tri ,'XData',strrep(dx_names,'_',' '),'YData',strrep(dx_names,'_',' '));
title('p-value (centiles)')


%% Groups %%
load("G0.mat",'-mat');
load("G1.mat",'-mat');
load("G1_5.mat",'-mat');
load("G2.mat",'-mat');
groups_names = {'SCZ_and_SAD-relatives','PE','FEP','SCZ_and_SAD-chronic'};

names_group = {'G0','G1','G1_5','G2'};
sse_real_total_groups = nan(length(names_group));
sse_real_total_groups = array2table(sse_real_total_groups,'RowNames',names_group,'VariableNames',names_group);
PCT_pval_groups = nan(length(names_group));
PCT_pval_groups = array2table(PCT_pval_groups,'RowNames',names_group,'VariableNames',names_group);
for i = 1:1:length(names_group)-1
    for j = i+1:length(names_group)
        [names_dx_1,names_dx_2] = names_group{[i,j]};
        centiles_dx_1 = eval(['centiles_',names_dx_1]);
        centiles_dx_2 = eval(['centiles_',names_dx_2]);

        % SSE
        sse_real_groups = mse(mean(table2array(centiles_dx_1)),mean(table2array(centiles_dx_2)))*size(centiles_dx_1,2);

        for iperm = 1:nperm
            % Group membership is randomly reassigned across the compared groups
            [selection_dx_1,selection_dx_2,ratio_1(iperm),ratio_2(iperm)] = mix_dx(centiles_dx_1,centiles_dx_2);
            
            selection_dx_1_m = mean(selection_dx_1);
            selection_dx_2_m = mean(selection_dx_2);

            % SSE
            sse_perm_groups(iperm) = mse(selection_dx_1_m,selection_dx_2_m)*size(selection_dx_1,2);
    
        end 
        % Test if it explains more variance that expected by chance
        PCT_pval_groups{i,j} = sum(sse_perm_groups > sse_real_groups) / nperm; 
        sse_real_total_groups{i,j} = sse_real_groups;

    end
end

figure;
sse_real_total_groups_tri = triu(sse_real_total_groups{:,:}) + triu(sse_real_total_groups{:,:},1).';
h = heatmap(sse_real_total_groups_tri,'XData',strrep(groups_names,'_',' '),'YData',strrep(groups_names,'_',' '),'FontSize',20);
h.CellLabelFormat = '%.3f'; 
title('SSE (centiles)')
figure;
PCT_pval_groups_tri = triu(PCT_pval_groups{:,:}) + triu(PCT_pval_groups{:,:},1).';
heatmap(PCT_pval_groups_tri,'XData',strrep(groups_names,'_',' '),'YData',strrep(groups_names,'_',' '));
title('p-value (centiles)')

'DONE'
