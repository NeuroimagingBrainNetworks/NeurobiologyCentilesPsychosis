%% Function to create randomized groups by mixing patients with different 
%% diagnoses or group membership

% Copyright (C) 2023 University of Seville

% Written by Natalia García San Martín (ngarcia1@us.es)

% This file is part of Neurobiology Centiles Psychosis toolkit.
%
% Neurobiology Centiles Mapping toolkit is free software: 
% you can redistribute it and/or modify it under the terms of the 
% GNU General Public License as published by the Free Software Foundation, 
% either version 3 of the License, or (at your option) any later version.
%
% Neurobiology Centiles Mapping toolkit is distributed in the hope that 
% it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Neurobiology Centiles Mapping toolkit. If not, see 
% <https://www.gnu.org/licenses/>.

function [selection_dx_1,selection_dx_2,ratio_1,ratio_2] = mix_dx(centiles_dx_1,centiles_dx_2)
    
    if height(centiles_dx_1) <= height(centiles_dx_2)
    else
        temp = centiles_dx_1;
        centiles_dx_1 = centiles_dx_2;
        centiles_dx_2 = temp;
    end

    centiles_dx_1 = table2array(centiles_dx_1(:,sort(centiles_dx_1.Properties.VariableNames))); % IMPORTANT
    centiles_dx_2 = table2array(centiles_dx_2(:,sort(centiles_dx_2.Properties.VariableNames)));

    ratio = size(centiles_dx_1,1)/(size(centiles_dx_1,1)+size(centiles_dx_2,1));

    random_ratio = ratio*0.9 + rand(1,1) * (ratio*0.2);  %ratio_range(randi(length(ratio_range)));
    random_ratio = max([0 random_ratio]);
    random_ratio = min([1 random_ratio]);

    num_filas_matriz_1 = round(random_ratio*size(centiles_dx_1,1));
    num_filas_matriz_2 = size(centiles_dx_1,1) - num_filas_matriz_1;
  
    dx_1_index_1 = randsample(size(centiles_dx_1,1),num_filas_matriz_1);
    dx_1_index_2 = setdiff((1:size(centiles_dx_1,1))',dx_1_index_1);

    dx_2_index_1 = randsample(size(centiles_dx_2,1),num_filas_matriz_2);
    dx_2_index_2 = setdiff((1:size(centiles_dx_2,1))',dx_2_index_1);
            
    selection_dx_1 = [centiles_dx_1(dx_1_index_1,:);centiles_dx_2(dx_2_index_1,:)];
    selection_dx_2 = [centiles_dx_1(dx_1_index_2,:);centiles_dx_2(dx_2_index_2,:)];

    ratio_1 = numel(dx_1_index_1)/size(selection_dx_1,1);
    ratio_2 = numel(dx_1_index_2)/size(selection_dx_2,1);
end


