%% Script to merge regional brain maps of empirical centile/effect sizes 
%% with each regional brain map of their respective significant regions and 
%% predicted centiles/effect sizes into a single figure

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
location = 'D:\OneDrive - UNIVERSIDAD DE SEVILLA\Natalia\Code\Data\maps';
location = 'C:\Users\usuario\OneDrive - UNIVERSIDAD DE SEVILLA\Natalia\Code\Data\maps';


%%%% EMPIRICAL, SIGNIFICANT, AND PREDICTED %%%%
images = dir(fullfile([location,'\merged_empirical_significant'], '*.png'));
images(contains({images.name},'legend')) = [];
   
for j = 1:length(images)

    % Load empirical and significant images
    image = imread([location,'\merged_empirical_significant\',images(j).name]);

    % Load predicted images
    try
        if contains(images(j).name, 'effsizes')
            pred_image = imread([location,'\merged_empirical_predicted\',replace(images(j).name,'effsizes_','')]);
        else
            pred_image = imread([location,'\merged_empirical_predicted\',replace(images(j).name,'mean_','')]);
        end
        
    catch
        continue
    end
    if ismatrix(pred_image)
        continue
    end
    [rows, cols, ~] = size(pred_image);
    middle = floor(rows / 2);
    pred_inferior_middle = pred_image(middle:end-200, :,:);

    image_reconstructed = [image(1:2000,:,:); pred_inferior_middle];
%         image_recontrsucted = [image_recontrsucted, legend];
    figure;
    imshow(image_reconstructed)

    
    print(gcf,fullfile([location,'\merged_empirical_significant_predicted\'], images(j).name),'-dpng', '-r330')
end

