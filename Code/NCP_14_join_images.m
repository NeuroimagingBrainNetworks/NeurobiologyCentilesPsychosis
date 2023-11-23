%% Script to merge regional brain maps of empirical centile/effect sizes 
%% with each regional brain map of their respective significant regions and 
%% predicted centiles/effect sizes into separate figures

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
location = 'D:\OneDrive - UNIVERSIDAD DE SEVILLA\Natalia\Code\Data\maps\';
location = 'C:\Users\usuario\OneDrive - UNIVERSIDAD DE SEVILLA\Natalia\Code\Data\maps\';


files = dir([location,'\empirical_maps\']);
files = files(3:4);
% files([4,6:8]) = [];

%%%% SIGNIFICANT REGIONS %%%%
for i = 1:length(files)
    images = dir(fullfile([location,'\empirical_maps\',files(i).name], '*.png'));
   
    for j = 1:length(images)

        % Load images
        image = imread([location,'\empirical_maps\',files(i).name,'\',images(j).name]);
        [rows, cols, ~] = size(image);
        middle = floor(rows / 2);
        if contains(images(j).name,'mean')
            legend_centiles = image(160+1:end,2097:end,:);
        elseif contains(images(j).name,'effsizes')
            legend_effsizes = image(160+1:end,2097:end,:);
        end
        inferior_middle = image(middle + 1:end-80, 145:2097, :);

        % Load z_values
        z_value_image = imread([location,'empirical_maps\significant_regions\',files(i).name,'\',replace(images(j).name,'mean','centiles')]);
        if ismatrix(z_value_image)
            continue
        end
        [rows, cols, ~] = size(z_value_image);
        middle = floor(rows / 2);
        z_value_inferior_middle = z_value_image(middle + 1:end-80, 145:2097, :);

        image_recontrsucted = [inferior_middle; z_value_inferior_middle];
%         image_recontrsucted = [image_recontrsucted, legend];
        figure;
        imshow(image_recontrsucted)
        title(replace(replace(images(j).name,'mean_',''),'.png',''),'Interpreter','none','Position',[975, 7, 0],'FontSize',15)
 
        print(gcf,fullfile([location,'\merged_empirical_significant\'], images(j).name),'-dpng', '-r330')
    end
end

imshow(legend_centiles)
saveas(gcf,fullfile([location,'\merged_empirical_significant\legend_centiles.png']),'png')
imshow(legend_effsizes)
saveas(gcf,fullfile([location,'\merged_empirical_significant\legend_effsizes']),'png')


%%%% PREDICTED MAPS %%%%
images_pred = dir(fullfile([location,'\predicted_maps'], '*.png'));

for i = 1:length(files)
    images = dir(fullfile([location,'\empirical_maps\',files(i).name], '*.png'));
    if any(ismember({images_pred.name},replace(replace({images.name},'effsizes_',''),'mean_','')))
        % Select the same image
        idx_pred = find(ismember({images_pred.name},replace(replace({images.name},'effsizes_',''),'mean_','')));

        for k = 1:length(idx_pred)
            index_pred = idx_pred(k);
            idx_real = find(ismember(replace(replace({images.name},'effsizes_',''),'mean_',''), images_pred(index_pred).name));

            % Load images
            image = imread([location,'\empirical_maps\',files(i).name,'\',images(idx_real).name]);
            [rows, cols, ~] = size(image);
            middle = floor(rows / 2);
            inferior_middle = image(middle + 1:end-80, 145:2097, :);
        
            % Load Y_pred
            y_pred_image = imread([location,'\predicted_maps\',images_pred(index_pred).name]);
            [rows, cols, ~] = size(y_pred_image);
            middle = floor(rows / 2);
            y_pred_inferior_middle = y_pred_image(middle + 1:end-80, 145:2097, :);
    
            image_recontrsucted = [inferior_middle; y_pred_inferior_middle];
    %         image_recontrsucted = [image_recontrsucted, legend];
            figure;
            imshow(image_recontrsucted)
            title(replace(images_pred(index_pred).name,'.png',''),'Interpreter','none','Position',[975, 7, 0],'FontSize',15)
    
            print(gcf,fullfile([location,'\merged_empirical_predicted\'], images_pred(index_pred).name),'-dpng', '-r330')
            
        end
    end
end
