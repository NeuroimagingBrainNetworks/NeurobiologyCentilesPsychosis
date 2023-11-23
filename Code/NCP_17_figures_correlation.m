%% Script to plot (1) neurobiological similarity matrix, 
%% (2) co-vulnerability to psychosis matrix, and (3) their correlation

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
location = 'C:\Users\usuario\OneDrive - UNIVERSIDAD DE SEVILLA\Natalia\Code\Data';
cd(location)

% Load dataset
load("G0.mat",'-mat');
load("G1.mat",'-mat');
load("G1_5.mat",'-mat')
load("G2.mat",'-mat');
load("FEP.mat")
load("PE.mat")
load("SCZ.mat")
load("Schizoaffective_Disorder.mat")
total_centiles = [centiles_G0;centiles_G1;centiles_G1_5;centiles_G2];
load("weights_groups.mat")

% Load molecular maps
molecular_map_hemi_ordered = readtable([location,'\Molecular_maps\molecular_map_hemi_ordered.csv'],ReadRowNames=true);
molecular_names = weights_table_groups.Properties.VariableNames;
MOR_a_VAChT = molecular_map_hemi_ordered(:,16:19);
mGluR_5 = molecular_map_hemi_ordered(:,15);
oligo = molecular_map_hemi_ordered(:,25);
OPC = molecular_map_hemi_ordered(:,26);
thickness = molecular_map_hemi_ordered(:,35);
expansion = molecular_map_hemi_ordered(:,36:37);
glu = molecular_map_hemi_ordered(:,38:41);
scaling = molecular_map_hemi_ordered(:,42:43);
synapse = molecular_map_hemi_ordered(:,45);

molecular_map_hemi_ordered(:,15:18) = MOR_a_VAChT;
molecular_map_hemi_ordered(:,19) = mGluR_5;
molecular_map_hemi_ordered(:,25) = OPC;
molecular_map_hemi_ordered(:,26) = oligo;
molecular_map_hemi_ordered(:,36) = synapse;
molecular_map_hemi_ordered(:,37) = thickness;
molecular_map_hemi_ordered(:,38:39) = expansion;
molecular_map_hemi_ordered(:,40:41) = scaling;
molecular_map_hemi_ordered(:,42:45) = glu;

molecular_map_hemi_ordered.Properties.VariableNames = molecular_names;

%% Molecular maps correlation (neurobiological similarity)
corr_molecular = corr(molecular_map_hemi_ordered{:,:}');
for i = 1:size(corr_molecular)
    corr_molecular(i, i) = NaN;
end

figure;
imagesc(corr_molecular);
colorbar
title('Correlation of neurobiological loadigs')
% xticks(1:height(molecular_map_hemi_ordered))
% xticklabels(molecular_map_hemi_ordered.Properties.RowNames);
xlabel('Regions')
% yticks(1:height(molecular_map_hemi_ordered))
% yticklabels(molecular_map_hemi_ordered.Properties.RowNames);
ylabel('Regions')
set(gca, "FontSize",16)

% Garnet and cream colormap 
num_colors = 256;  % numbers of colors
garnet_colormap = zeros(num_colors, 3);

garnet = [128, 0, 0] / 255;
cream = [255, 180, 0] / 255; 
white = [1, 1, 1];

% Garnet color gradually assigned
for i = 1:num_colors
    t = (i - 1) / (num_colors - 1);  % From 0 to 1
    if t <= 0.5
        garnet_colormap(i, :) = (1 - 2 * t) * garnet + (2 * t) * white;
    else
        garnet_colormap(i, :) = (2 - 2 * t) * white + (2 * t - 1) * cream;
    end
end
colormap(garnet_colormap);


%% Centile correlation (co-vulnerability to psychosis)

% "Sky" colormap
sky_colormap = zeros(num_colors, 3);

% Sky colors from -1 to 1 values
azul_oscuro = [0, 0, 128] / 255;      % Dark blue (-1)
white = [1, 1, 1];                   % White (0)
azul_claro = [0, 191, 200] / 255;   % Light blue (1)

% Colors gradually assigned
for i = 1:num_colors
    t = (i - 1) / (num_colors - 1);  % From 0 to 1
    if t < 0.5
        sky_colormap(i, :) = (1 - t * 2) * azul_oscuro + t * 2 * white;
    else
        sky_colormap(i, :) = (1 - (t - 0.5) * 2) * white + (t - 0.5) * 2 * azul_claro;
    end
end

% Load CN for effect size
PE_CN = readtable([location,'\Centiles\PE\centiles.csv']);pos_centiles=19:52;
PE_CN = PE_CN(strcmp(PE_CN.dx,'CN'),pos_centiles);
SCZ_CN = readtable([location,'\Centiles\SCZ_SAD\centiles.csv']);pos_centiles=10:43;
SCZ_CN = SCZ_CN(strcmp(SCZ_CN.dx,'CN') & ismember(SCZ_CN.study,{'ABCD','ASRB','BSNIP','LA5c','MCIC','UKB'}),pos_centiles);
FEP_CN = readtable([location,'\Centiles\FEP\centiles.csv']);pos_centiles=19:52;
FEP_CN = FEP_CN(strcmp(FEP_CN.dx,'CN'),:);
CN_PS_subjects=unique(FEP_CN.participant);
for i=1:height(CN_PS_subjects)
    CN_PS_subject = CN_PS_subjects(i);
    centiles_CN_PS(i,:)=mean(table2array(FEP_CN(FEP_CN{:,'participant'}==CN_PS_subject,pos_centiles)),1);
end
FEP_CN = array2table(centiles_CN_PS,"VariableNames",SCZ_CN.Properties.VariableNames);

total_CN = [PE_CN;SCZ_CN;FEP_CN];

Y = corr(molecular_map_hemi_ordered{:,:}');

Y=Y(find(triu(ones(size(Y)),1)));

for i=1:width(total_centiles)
    eff_G2(i)=computeCohen_d(centiles_G2{:,i},total_CN{:,i});
end
for i=1:width(total_centiles)
    eff_G1_5(i)=computeCohen_d(centiles_G1_5{:,i},total_CN{:,i});
end
for i=1:width(total_centiles)
    eff_G1(i)=computeCohen_d(centiles_G1{:,i},total_CN{:,i});
end
for i=1:width(total_centiles)
    eff_G0(i)=computeCohen_d(centiles_G0{:,i},total_CN{:,i});
end
for i=1:width(total_centiles)
    eff_PE_1(i)=computeCohen_d(centiles_PE_1{:,i},total_CN{:,i});
end
for i=1:width(total_centiles)
    eff_PE_2(i)=computeCohen_d(centiles_PE_2{:,i},total_CN{:,i});
end
for i=1:width(total_centiles)
    eff_PE_3(i)=computeCohen_d(centiles_PE_3{:,i},total_CN{:,i});
end
for i=1:width(total_centiles)
    eff_FEP(i)=computeCohen_d(centiles_FEP_1{:,i},total_CN{:,i});
end
for i=1:width(total_centiles)
    eff_SAD_relative(i)=computeCohen_d(centiles_Schizoaffective_Disorder_Relative{:,i},total_CN{:,i});
end
for i=1:width(total_centiles)
    eff_SAD(i)=computeCohen_d(centiles_Schizoaffective_Disorder{:,i},total_CN{:,i});
end
for i=1:width(total_centiles)
    eff_SCZ(i)=computeCohen_d(centiles_SCZ{:,i},total_CN{:,i});
end
for i=1:width(total_centiles)
    eff_SCZ_relative(i)=computeCohen_d(centiles_SCZ_Relative{:,i},total_CN{:,i});
end

total_eff = [eff_G0;eff_G1;eff_G1_5;eff_G2];
total_eff_dx = [eff_PE_1;eff_PE_2;eff_PE_3;eff_FEP;eff_SAD_relative;eff_SAD;eff_SCZ;eff_SCZ_relative];

vulner = corr(total_eff);
vulner_dx = corr(total_eff_dx);

for i = 1:size(vulner)
    vulner(i, i) = NaN;
end
for i = 1:size(vulner_dx)
    vulner_dx(i, i) = NaN;
end

% Effect size correlation
figure;
imagesc(vulner);
colorbar
title('Correlation of regional effect sizes (groups)')
% xticks(1:height(molecular_map_hemi_ordered))
% xticklabels(molecular_map_hemi_ordered.Properties.RowNames);
xlabel('Regions')
% yticks(1:height(molecular_map_hemi_ordered))
% yticklabels(molecular_map_hemi_ordered.Properties.RowNames);
ylabel('Regions')
set(gca, "FontSize",16)
colormap(sky_colormap)

median(vulner(:),'omitnan')
std(vulner(:),'omitnan')

figure;
imagesc(vulner_dx);
colorbar
title('Correlation of regional effect sizes (dx)')
% xticks(1:height(molecular_map_hemi_ordered))
% xticklabels(molecular_map_hemi_ordered.Properties.RowNames);
xlabel('Regions')
% yticks(1:height(molecular_map_hemi_ordered))
% yticklabels(molecular_map_hemi_ordered.Properties.RowNames);
ylabel('Regions')
set(gca, "FontSize",16)
colormap(sky_colormap)


vulner = vulner(find(triu(ones(size(vulner)),1)));
vulner_dx = vulner_dx(find(triu(ones(size(vulner_dx)),1)));

[r,p] = corr(vulner,Y);

[r_dx,p_dx] = corr(vulner_dx,Y);

X_rotated = csvread([location,'\PCA_CCA\perm_sphere_10000_DK.csv'])';
molecular_maps = readtable([location,'\PCA_CCA\all_microsc_DesikanKilliany68.csv'],ReadVariableNames=true);
X_rot = X_rotated(1:34,:);
molecular_maps_labels = molecular_maps.Var2(1:34);
[molecular_maps_labels_ordered, molecular_maps_labels_order] = sort(molecular_maps_labels);
X_rot_ordered = X_rot(molecular_maps_labels_order,:)';
nperm = 10000;
for i = 1:nperm
    molecular_map_hemi_ordered_perm = molecular_map_hemi_ordered{X_rot_ordered(i,:),:};
    Y_perm = corr(molecular_map_hemi_ordered_perm');
    Y_perm = Y_perm(find(triu(ones(size(Y_perm)),1)));
    [r_perm(i),p_perm(i)]=corr(vulner,Y_perm);
    [r_perm_dx(i),p_perm_dx(i)]=corr(vulner_dx,Y_perm);
end
pval = sum(abs(r_perm) > abs(r)) / nperm;
pval_dx = sum(abs(r_perm_dx) > abs(r_dx)) / nperm;

figure;
scatter(Y,vulner,'MarkerFaceColor',[255 150 0] / 255,'MarkerEdgeColor',[255 150 0] / 255)
h = lsline;
set(h, 'LineWidth', 2,'Color','k');
ylabel('Psychosis co-vulnerability groups','FontSize',12)
xlabel('Neurobiological features','FontSize',12)
ylim([-1.3, 1.3])
set(gca, "FontSize",16)
if p < 0.05
    textString = ['r = ', sprintf('%.2f',r), '; ', '\color{red}', ' P_{spin} = ', sprintf('%.3e', p)]; 
else            
    textString = ['r = ', sprintf('%.2f',r), '; P_{spin} = ', sprintf('%.2f', p)];
end
annotation('textbox', [0.4 0.8 0.1 0.1], 'String', textString, 'FontSize', 10, 'HorizontalAlignment', 'center');


figure;
scatter(Y,vulner_dx,'MarkerFaceColor',[255 150 0] / 255,'MarkerEdgeColor',[255 150 0] / 255)
h = lsline;
set(h, 'LineWidth', 2,'Color','k');
ylabel('Psychosis co-vulnerability dx','FontSize',12)
xlabel('Neurobiological similarity','FontSize',12)
ylim([-1.2, 1.2])
set(gca, "FontSize",16)
if p_dx < 0.05
    textString = ['r = ', sprintf('%.2f',r_dx), '; ', '\color{red}', ' P_{spin} = ', sprintf('%.3e', p_dx); ];
    
else            
    textString = ['r = ', sprintf('%.2f',r_dx), '; P_{spin} = ', sprintf('%.2f', p_dx)];
end
annotation('textbox', [0.4 0.8 0.1 0.1], 'String', textString, 'FontSize', 10, 'HorizontalAlignment', 'center');

