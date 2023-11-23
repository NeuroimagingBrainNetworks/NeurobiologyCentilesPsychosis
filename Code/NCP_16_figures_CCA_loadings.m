%% Script to plot the significant loadings and weights of PCA-CCA modeling

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

location = 'D:\OneDrive - UNIVERSIDAD DE SEVILLA\Natalia\Code\Data';
location = 'C:\Users\usuario\OneDrive - UNIVERSIDAD DE SEVILLA\Natalia\Code\Data';
directory = [location,'\PCA_CCA\']; % change as desired

cd(directory);

selection = 1:4; % just groups (without effect sizes)
% selection = 1:10; % groups and effect sizes
Loadings = 1; % 0 for weights, 1 for Loadings
if Loadings
    load("loadings_dx.mat")
    load("loadings_groups.mat")
    weights_table_dx = loadings_table_dx;
    loadings_table_groups = loadings_table_groups(selection,:);
    weights_table_groups = loadings_table_groups;
    text_distance = 2.5;
    y_lab = 'Loading';
else
    load("weights_dx.mat")
    load("weights_groups.mat")
    text_distance = 0.033;
    y_lab = 'Weight';
end
weights_table_groups = weights_table_groups(selection,:);

dx_names = {'SCZ-relatives','SAD-relatives','PE-suspected','PE-definite','PE-clinical', 'FEP','SCZ','SAD'};
groups_names = {'SCZ_and_SAD-relatives','PE','FEP','SCZ_and_SAD-chronic', 'SCZ_and_SAD-chronic_vs_SCZ_and_SAD-relatives', ...
    'PE_vs_SCZ_and_SAD-relatives','SCZ_and_SAD-chronic_vs_PE','SCZ_and_SAD-chronic_vs_FEP','FEP_vs_PE','FEP_vs_SCZ_and_SAD-relatives'};
groups_names = groups_names(selection);

weights_total = [weights_table_dx; weights_table_groups];
writetable(weights_total,'weights.csv','WriteRowNames',true,'WriteVariableNames',true)

molecular_names_table = readtable([location,'\Molecular_maps\molecular_names.xlsx'],ReadVariableNames=true,ReadRowNames=true);
molecular_names_table = sortrows(molecular_names_table,{'Type_numeric'});
molecular_types_orig = molecular_names_table{:,2};

%%%%%%%%%%%%%%%%%%%%%%%%%%%% DIAGNOSES (DX) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bar (weights grouped by dx)

positive_weights_dx = table2array(weights_table_dx);
negative_weights_dx = positive_weights_dx;
positive_weights_dx(positive_weights_dx < 0) = NaN;
negative_weights_dx(negative_weights_dx > 0) = NaN;
positive_weights_dx = array2table(positive_weights_dx');
negative_weights_dx = array2table(negative_weights_dx');
positive_weights_dx.Properties.VariableNames = weights_table_dx.Properties.RowNames;
negative_weights_dx.Properties.VariableNames = weights_table_dx.Properties.RowNames;
positive_weights_dx.molecular_types(1:46) = molecular_types_orig;
negative_weights_dx.molecular_types(1:46) = molecular_types_orig;

weights_table_dx =  rows2vars(weights_table_dx);
weights_table_dx.molecular_types(1:46) = molecular_types_orig;

func = @(x) sum(x,'omitnan')/numel(x);
sum_by_type_pos = varfun(func, positive_weights_dx,"GroupingVariables", 'molecular_types') ;
sum_by_type_neg = varfun(func, negative_weights_dx,"GroupingVariables", 'molecular_types') ;
figure;
bar(sum_by_type_pos{:,3:end}','grouped');
max_pos = max(max(sum_by_type_pos{:,3:end}));
hold on
bar(sum_by_type_neg{:,3:end}','grouped');
max_neg = min(min(sum_by_type_neg{:,3:end}));
ylim([-max(max_pos,abs(max_neg)) * 1.1, max(max_pos,abs(max_neg)) * 1.1])
colorsmap = [0,1,1;0.4,0.6,1;0.8,0.7,1;0.2422,0.1504,0.6603;0.8,0.2,1;0.9059,0.0980,0.4824];
set(gca,'ColorOrder',colorsmap)
xticklabels(replace(dx_names,'_',' '))
ylabel(y_lab)
legend(unique(molecular_names_table{:,1},'stable'),'Location','best')


%% Bar (weights grouped by molecular_name)
if Loadings
    weights_table_dx = loadings_table_dx;
else
    load("weights_dx.mat")
end

figure;
hold on
[~, id_order] = sort(sum(weights_table_dx{:,:},1));
b = bar(weights_table_dx{:,id_order}','stacked');
colores = viridis(size(weights_table_dx,1));
for i = 1:size(weights_table_dx,1)
    b(i).FaceColor = colores(i,:);
end
xticks(1:size(weights_table_dx,2))
xtickangle(45);
ylabel(['\bf\fontsize{24}' y_lab]); % bold
ax = gca;
ax.FontSize = 20;
legend(replace(dx_names,'_',' '),'Location','best','Box','off')
for i = 1:size(weights_table_dx,2)
    text(i,min(ylim)-text_distance-2,weights_table_dx.Properties.VariableNames(id_order(i)),'Color',colorsmap(molecular_names_table{id_order(i),2},:),'Rotation',90,'FontSize', 16)
end
set(gca,'XTickLabel',[]);
ax = gca;
ax.Position(2) = ax.Position(2) + 0.2;
ax.Position(4) = ax.Position(4) - 0.2;


% Imagesc
figure;
imagesc(weights_table_dx{:,:})
colorbar;
xticks(1:size(weights_table_dx,2))
xticklabels(weights_table_dx.Properties.VariableNames)
yticks(1:size(weights_table_dx,1))
yticklabels(replace(dx_names,'_',' '))
title(y_lab);


% Correlation dx
figure;
imagesc(corr(weights_table_dx{:,:}));
colorbar
title(replace(['Correlation (',y_lab,'_dx)'],'_',' '))
xticks(1:length(weights_table_dx.Properties.VariableNames))
xticklabels(weights_table_dx.Properties.VariableNames);
yticks(1:length(weights_table_dx.Properties.VariableNames))
yticklabels(weights_table_dx.Properties.VariableNames);

% SSD dx
weights_table_dx =  rows2vars(weights_table_dx);
weights_table_dx.Properties.RowNames = weights_table_dx{:,1}; 
weights_table_dx = removevars(weights_table_dx,'OriginalVariableNames_1');
for molecular_name = 1:width(weights_table_dx)
    for j = 1:width(weights_table_dx)
        sse_dx(molecular_name,j) = mse(weights_table_dx{:,molecular_name},weights_table_dx{:,j})*size(weights_table_dx,1);
    end
end
figure;
heatmap(sse_dx,'XData',replace(dx_names,'_',' '),'YData',replace(dx_names,'_',' '),'Title',['SSD (',y_lab,')'])

figure;
heatmap(corr(weights_table_dx{:,:}),'XData',replace(dx_names,'_',' '),'YData',replace(dx_names,'_',' '),'Title',['Correlation (',y_lab,')'])



%%%%%%%%%%%%%%%%%%%%%%%%%%%% GROUPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Bar (weights grouped by dx)
positive_weights_groups = table2array(weights_table_groups);
negative_weights_groups = positive_weights_groups;
positive_weights_groups(positive_weights_groups < 0) = NaN;
negative_weights_groups(negative_weights_groups > 0) = NaN;
positive_weights_groups = array2table(positive_weights_groups');
negative_weights_groups = array2table(negative_weights_groups');
positive_weights_groups.Properties.VariableNames = weights_table_groups.Properties.RowNames;
negative_weights_groups.Properties.VariableNames = weights_table_groups.Properties.RowNames;
positive_weights_groups.molecular_types(1:46) = molecular_types_orig;
negative_weights_groups.molecular_types(1:46) = molecular_types_orig;

weights_table_groups =  rows2vars(weights_table_groups);
weights_table_groups.molecular_types(1:46) = molecular_types_orig;

func = @(x) sum(x,'omitnan')/numel(x);
sum_by_type_pos = varfun(func, positive_weights_groups,"GroupingVariables", 'molecular_types') ;
sum_by_type_neg = varfun(func, negative_weights_groups,"GroupingVariables", 'molecular_types') ;
figure;
bar(sum_by_type_pos{:,3:end}','grouped');
max_pos = max(max(sum_by_type_pos{:,3:end}));
hold on
bar(sum_by_type_neg{:,3:end}','grouped');
max_neg = min(min(sum_by_type_neg{:,3:end}));
ylim([-max(max_pos,abs(max_neg)) * 1.1, max(max_pos,abs(max_neg)) * 1.1])
colorsmap = [0,1,1;0.4,0.6,1;0.8,0.7,1;0.2422,0.1504,0.6603;0.8,0.2,1;0.9059,0.0980,0.4824];

set(gca,'ColorOrder',colorsmap)
xticklabels(replace(groups_names,'_',' '))
xtickangle(0);
ylabel(['\bf\fontsize{24}' y_lab]); % bold
ax = gca;
ax.FontSize = 20;
legend(unique(molecular_names_table{:,1},'stable'),'Location','best','Box','off')


%% Bar (weights grouped by molecular_name)
if Loadings
    weights_table_groups = loadings_table_groups;
else
    load("weights_groups.mat")
end
weights_table_groups = weights_table_groups(selection,:);

figure;
hold on
[~, id_order] = sort(sum(weights_table_groups{:,:},1));
b = bar(weights_table_groups{:,id_order}','stacked');
colores = viridis(size(weights_table_groups,1));
for i = 1:size(weights_table_groups,1)
    b(i).FaceColor = colores(i,:);
end
xticks(1:size(weights_table_groups,2))
xtickangle(45);
ylabel(['\bf\fontsize{24}' y_lab]); % bold
ax = gca;
ax.FontSize = 20;
legend(replace(groups_names,'_',' '),'Location','best','Box','off')
for i = 1:size(weights_table_groups,2)
%     text(i,min(ylim)-1.3,weights_table_groups.Properties.VariableNames(i),'Color',colorsmap(molecular_names_table{i,2},:),'Rotation',90)
    text(i,min(ylim)-text_distance,weights_table_groups.Properties.VariableNames(id_order(i)),'Color',colorsmap(molecular_names_table{id_order(i),2},:),'Rotation',90,'FontSize', 16)
end
set(gca,'XTickLabel',[]);
ax = gca;
ax.Position(2) = ax.Position(2) + 0.2;
ax.Position(4) = ax.Position(4) - 0.2;

% Imagesc
figure;
imagesc(weights_table_groups{:,:})
colorbar;
xticks(1:size(weights_table_groups,2))
xticklabels(weights_table_groups.Properties.VariableNames)
yticks(1:size(weights_table_groups,1))
yticklabels(replace(groups_names,'_',' '))
title(y_lab)


% Correlation groups
figure;
imagesc(corr(weights_table_groups{:,:}));
colorbar
title(replace(['Correlation (',y_lab,' groups)'],'_',' '))
xticks(1:length(weights_table_groups.Properties.VariableNames))
xticklabels(replace(weights_table_groups.Properties.VariableNames,'_',' '));
yticks(1:length(weights_table_groups.Properties.VariableNames))
yticklabels(replace(weights_table_groups.Properties.VariableNames,'_',' '));


% SSE groups
weights_table_groups =  rows2vars(weights_table_groups);
weights_table_groups.Properties.RowNames = weights_table_groups{:,1}; 
weights_table_groups = removevars(weights_table_groups,'OriginalVariableNames_1');

figure;
heatmap(corr(weights_table_groups{:,:}),'XData',replace(groups_names,'_',' '),'YData',replace(groups_names,'_',' '),'Title',replace(['Corralation (',y_lab,')'],'_',' '));
