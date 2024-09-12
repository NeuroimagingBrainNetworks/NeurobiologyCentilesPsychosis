%% Script to obtain centiles by SCZ-SAD diagnosis (relatives and chronic)

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

clear
close all

% Change as desired
location = 'D:\OneDrive - UNIVERSIDAD DE SEVILLA\Natalia\Code\Data\Centiles\SCZ_SAD\';
location = 'C:\Users\usuario\OneDrive - UNIVERSIDAD DE SEVILLA\Natalia\Code\Data\Centiles\SCZ_SAD\';
cd(location); 

bd_cam=readtable('centiles.csv');pos_centiles=10:43;
conditions_all=unique(bd_cam.dx);
conditions=conditions_all(2:end,1);
studies=unique(bd_cam.study);

regions=bd_cam.Properties.VariableNames(:,pos_centiles);


sex=bd_cam.sex;
ind=1;indM=1;indF=1;
for is=1:size(bd_cam,1)
    if strcmp(bd_cam.dx{is},'Schizoaffective Disorder')
        centiles_Schizoaffective_Disorder(ind,:)=bd_cam(is,pos_centiles);
        ind=ind+1;
    end
end
total_Schizoaffective_Disorder = bd_cam(strcmp(bd_cam.dx,'Schizoaffective Disorder'),:);

ind=1;indM=1;indF=1;
for is=1:size(bd_cam,1)
    if strcmp(bd_cam.dx{is},'Schizoaffective Disorder - Relative')
        centiles_Schizoaffective_Disorder_Relative(ind,:)=bd_cam(is,pos_centiles);
        ind=ind+1;
    end
end
total_Schizoaffective_Disorder_Relative = bd_cam(strcmp(bd_cam.dx,'Schizoaffective Disorder - Relative'),:);

ns=size(bd_cam,1);
ind=1;indM=1;indF=1;
for is=1:size(bd_cam,1)
    if strcmp(bd_cam.dx{is},'SCZ-Relative')
        centiles_SCZ_Relative(ind,:)=bd_cam(is,pos_centiles);
        ind=ind+1;
    end
end
total_SCZ_Relative = bd_cam(strcmp(bd_cam.dx,'SCZ-Relative'),:);

ind=1;indM=1;indF=1;
for is=1:size(bd_cam,1)
    if strcmp(bd_cam.dx{is},'SCZ')
        centiles_SCZ(ind,:)=bd_cam(is,pos_centiles);
        ind=ind+1;
    end
end
total_SCZ = bd_cam(strcmp(bd_cam.dx,'SCZ'),:);

ind=1;indM=1;indF=1;indY=1;indO=1;
for is=1:size(bd_cam,1)
    if strcmp(bd_cam.dx{is},'CN' & ismember(bd_cam.study,{'ABCD','ASRB','BSNIP','LA5c','MCIC','UKB'}))
        centiles_CN(ind,:)=bd_cam(is,pos_centiles);
        ind=ind+1;
    end
end
total_CN = bd_cam(strcmp(bd_cam.dx,'CN') & ismember(bd_cam.study,{'ABCD','ASRB','BSNIP','LA5c','MCIC','UKB'}),:);

centiles_SCZ_without_ASRB = centiles_SCZ(~strcmp(total_SCZ.study,'ASRB'),:);
centiles_SCZ_without_BSNIP = centiles_SCZ(~strcmp(total_SCZ.study,'BSNIP'),:);
centiles_SCZ_without_LA5c = centiles_SCZ(~strcmp(total_SCZ.study,'LA5c'),:);
centiles_SCZ_without_MCIC = centiles_SCZ(~strcmp(total_SCZ.study,'MCIC'),:);
save("SCZ","centiles_SCZ","centiles_SCZ_Relative","centiles_SCZ_without_ASRB","centiles_SCZ_without_BSNIP","centiles_SCZ_without_LA5c","centiles_SCZ_without_MCIC")
save("Schizoaffective_Disorder","centiles_Schizoaffective_Disorder","centiles_Schizoaffective_Disorder_Relative")


% For mapping (step THIRTEEN)
for ir=1:34
    mean_Schizoaffective_Disorder(ir)=mean(table2array(centiles_Schizoaffective_Disorder(:,ir)));
end

for ir=1:34
    mean_Schizoaffective_Disorder_Relative(ir)=mean(table2array(centiles_Schizoaffective_Disorder_Relative(:,ir)));
end

for ir=1:34
    mean_SCZ_Relative(ir)=mean(table2array(centiles_SCZ_Relative(:,ir)));
end

for ir=1:34
    mean_SCZ(ir)=mean(table2array(centiles_SCZ(:,ir)));
end

for ir=1:34
    mean_SCZ_without_ASRB(ir)=mean(table2array(centiles_SCZ_without_ASRB(:,ir)));
    mean_SCZ_without_BSNIP(ir)=mean(table2array(centiles_SCZ_without_BSNIP(:,ir)));
    mean_SCZ_without_LA5c(ir)=mean(table2array(centiles_SCZ_without_LA5c(:,ir)));
    mean_SCZ_without_MCIC(ir)=mean(table2array(centiles_SCZ_without_MCIC(:,ir)));
end

data=mean_Schizoaffective_Disorder';t_effsizes=array2table(data);t_effsizes.Properties.RowNames=transpose(regions);
writetable(t_effsizes,[location,'\mean_SAD.csv'],'WriteRowNames',true);

data=mean_Schizoaffective_Disorder_Relative';t_effsizes=array2table(data);t_effsizes.Properties.RowNames=transpose(regions);
writetable(t_effsizes,[location,'\mean_SAD-relatives.csv'],'WriteRowNames',true);

data=mean_SCZ_Relative';t_effsizes=array2table(data);t_effsizes.Properties.RowNames=transpose(regions);
writetable(t_effsizes,[location,'\mean_SCZ-relatives.csv'],'WriteRowNames',true);

data=mean_SCZ';t_effsizes=array2table(data);t_effsizes.Properties.RowNames=transpose(regions);
writetable(t_effsizes,[location,'\mean_SCZ.csv'],'WriteRowNames',true);

data=mean_SCZ_without_ASRB';t_effsizes=array2table(data);t_effsizes.Properties.RowNames=transpose(regions);
writetable(t_effsizes,[location,'\mean_SCZ_without_ASRB.csv'],'WriteRowNames',true);

data=mean_SCZ_without_BSNIP';t_effsizes=array2table(data);t_effsizes.Properties.RowNames=transpose(regions);
writetable(t_effsizes,[location,'\mean_SCZ_without_BSNIP.csv'],'WriteRowNames',true);

data=mean_SCZ_without_LA5c';t_effsizes=array2table(data);t_effsizes.Properties.RowNames=transpose(regions);
writetable(t_effsizes,[location,'\mean_SCZ_without_LA5c.csv'],'WriteRowNames',true);

data=mean_SCZ_without_MCIC';t_effsizes=array2table(data);t_effsizes.Properties.RowNames=transpose(regions);
writetable(t_effsizes,[location,'\mean_SCZ_without_MCIC.csv'],'WriteRowNames',true);

'DONE'
