%% Script to obtain centiles by PE diagnosis (suspected, definite, and clinical)

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
location = 'D:\OneDrive - UNIVERSIDAD DE SEVILLA\Natalia\Code\Data\Centiles\PE\';
location = 'C:\Users\usuario\OneDrive - UNIVERSIDAD DE SEVILLA\Natalia\Code\Data\Centiles\PE\';
cd(location); 

bd_cam=readtable('centiles.csv');pos_centiles=19:52;
conditions_all=unique(bd_cam.dx);
conditions_PE=conditions_all(2:end,1);

regions=bd_cam.Properties.VariableNames(:,pos_centiles);

sex=bd_cam.sex;
ns=size(bd_cam,1);
ind=1;indM=1;indF=1;
for is=1:size(bd_cam,1)
    if strcmp(bd_cam.dx{is},'FEP_1')
        centiles_PE_1(ind,:)=bd_cam(is,pos_centiles);
        ind=ind+1;
    end
end
total_FEP_1 = bd_cam(strcmp(bd_cam.dx,'FEP_1'),:);

ind=1;indM=1;indF=1;
for is=1:size(bd_cam,1)
    if strcmp(bd_cam.dx{is},'FEP_2')
        centiles_PE_2(ind,:)=bd_cam(is,pos_centiles);
        ind=ind+1;
    end
end
total_FEP_2 = bd_cam(strcmp(bd_cam.dx,'FEP_2'),:);

ind=1;indM=1;indF=1;
for is=1:size(bd_cam,1)
    if strcmp(bd_cam.dx{is},'FEP_3')
        centiles_PE_3(ind,:)=bd_cam(is,pos_centiles);
        ind=ind+1;
    end
end
total_FEP_3 = bd_cam(strcmp(bd_cam.dx,'FEP_3'),:);

ind=1;indM=1;indF=1;indY=1;indO=1;
for is=1:size(bd_cam,1)
    if strcmp(bd_cam.dx{is},'CN')
        centiles_CN(ind,:)=bd_cam(is,pos_centiles);
        ind=ind+1;
    end
end
total_CN = bd_cam(strcmp(bd_cam.dx,'CN'),:);

save("PE","centiles_PE_1","centiles_PE_2","centiles_PE_3")


% For mapping (step THIRTEEN)
for ir=1:34
    mean_PE_1(ir)=mean(table2array(centiles_PE_1(:,ir)));
end
for ir=1:34
    mean_PE_2(ir)=mean(table2array(centiles_PE_2(:,ir)));
end
for ir=1:34
    mean_PE_3(ir)=mean(table2array(centiles_PE_3(:,ir)));
end

data=mean_PE_1';t_mean=array2table(data);t_mean.Properties.RowNames=transpose(regions);
writetable(t_mean,[location,'\mean_PE-suspected.csv'],'WriteRowNames',true);

data=mean_PE_2';t_mean=array2table(data);t_mean.Properties.RowNames=transpose(regions);
writetable(t_mean,[location,'\mean_PE-definite.csv'],'WriteRowNames',true);

data=mean_PE_3';t_mean=array2table(data);t_mean.Properties.RowNames=transpose(regions);
writetable(t_mean,[location,'\mean_PE-clinical.csv'],'WriteRowNames',true);



'DONE'
