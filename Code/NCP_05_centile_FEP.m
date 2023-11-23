%% Script to obtain centiles by FEP diagnosis (first session)

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
location = 'D:\OneDrive - UNIVERSIDAD DE SEVILLA\Natalia\Code\Data\Centiles\FEP\';
location = 'C:\Users\usuario\OneDrive - UNIVERSIDAD DE SEVILLA\Natalia\Code\Data\Centiles\FEP\';
cd(location); 

bd_cam=readtable('centiles.csv');pos_centiles=19:52;
conditions_all=unique(bd_cam.dx);
conditions_PS=conditions_all(2:end,1);

regions=bd_cam.Properties.VariableNames(:,pos_centiles);

% PS session 1
ind=1;
for is=1:size(bd_cam,1)
    if bd_cam.session(is) == 1 && strcmp(bd_cam.dx(is),'PS')
        centiles_FEP_1(ind,:)=bd_cam(is,pos_centiles);
        ind=ind+1;
    end
end
total_PS_1 = bd_cam(strcmp(bd_cam.dx,'PS') & (bd_cam.session == 1),:);
total_CN_1 = bd_cam(strcmp(bd_cam.dx,'CN') & (bd_cam.session == 1),:);

save("FEP","centiles_FEP_1") 


% For mapping (step THIRTEEN)
for ir=1:34
    mean_FEP_1(ir)=mean(table2array(centiles_FEP_1(:,ir)));
end

data=mean_FEP_1';t_effsizes=array2table(data);t_effsizes.Properties.RowNames=transpose(regions);
writetable(t_effsizes,[location,'\mean_FEP.csv'],'WriteRowNames',true);

'DONE'
