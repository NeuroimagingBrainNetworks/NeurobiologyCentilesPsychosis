%% Script to obtain centiles and effect sizes between groups (relatives, PE, FEP, and chronic) 

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
location = 'D:\OneDrive - UNIVERSIDAD DE SEVILLA\Natalia\Code\Data\';
location = 'C:\Users\usuario\OneDrive - UNIVERSIDAD DE SEVILLA\Natalia\Code\Data\';
cd(location); 

% Load dataset
load("PE.mat",'-mat');
load("FEP.mat",'-mat');
load("SCZ.mat",'-mat');
load("Schizoaffective_Disorder.mat",'-mat');

% Define the groups
centiles_G0 = [centiles_SCZ_Relative; centiles_Schizoaffective_Disorder_Relative];
centiles_G1 = [centiles_PE_1; centiles_PE_2; centiles_PE_3];
centiles_G1_5 = centiles_FEP_1;
centiles_G2 = [centiles_SCZ; centiles_Schizoaffective_Disorder];
centiles_G2_without_ASRB = [centiles_SCZ_without_ASRB; centiles_Schizoaffective_Disorder];
centiles_G2_without_BSNIP = [centiles_SCZ_without_BSNIP; centiles_Schizoaffective_Disorder];
centiles_G2_without_LA5c = [centiles_SCZ_without_LA5c; centiles_Schizoaffective_Disorder];
centiles_G2_without_MCIC = [centiles_SCZ_without_MCIC; centiles_Schizoaffective_Disorder];

cd([location,'Centiles\groups\'])
save('G0',"centiles_G0")
save('G1',"centiles_G1")
save('G1_5',"centiles_G1_5")
save('G2',"centiles_G2")
save('G2_without_ASRB',"centiles_G2_without_ASRB")
save('G2_without_BSNIP',"centiles_G2_without_BSNIP")
save('G2_without_LA5c',"centiles_G2_without_LA5c")
save('G2_without_MCIC',"centiles_G2_without_MCIC")

centiles_tot = [centiles_G0;centiles_G1;centiles_G1_5;centiles_G2];
regions = centiles_tot.Properties.VariableNames;

% Group vs another
for ir=1:34
    effsizes_G1_vs_G0(ir)=computeCohen_d(table2array(centiles_G1(:,ir)),table2array(centiles_G0(:,ir)));
end

for ir=1:34
    effsizes_G2_vs_G0(ir)=computeCohen_d(table2array(centiles_G2(:,ir)),table2array(centiles_G0(:,ir)));
end

for ir=1:34
    effsizes_G2_vs_G1(ir)=computeCohen_d(table2array(centiles_G2(:,ir)),table2array(centiles_G1(:,ir)));
end

for ir=1:34
    effsizes_G1_5_vs_G0(ir)=computeCohen_d(table2array(centiles_G1_5(:,ir)),table2array(centiles_G0(:,ir)));
end

for ir=1:34
    effsizes_G1_5_vs_G1(ir)=computeCohen_d(table2array(centiles_G1_5(:,ir)),table2array(centiles_G1(:,ir)));
end

for ir=1:34
    effsizes_G2_vs_G1_5(ir)=computeCohen_d(table2array(centiles_G2(:,ir)),table2array(centiles_G1_5(:,ir)));
end

effsizes_G1_vs_G0 = array2table(effsizes_G1_vs_G0');effsizes_G1_vs_G0.Properties.RowNames=transpose(regions);
effsizes_G2_vs_G0 = array2table(effsizes_G2_vs_G0');effsizes_G2_vs_G0.Properties.RowNames=transpose(regions);
effsizes_G2_vs_G1 = array2table(effsizes_G2_vs_G1');effsizes_G2_vs_G1.Properties.RowNames=transpose(regions);
effsizes_G1_5_vs_G0 = array2table(effsizes_G1_5_vs_G0');effsizes_G1_5_vs_G0.Properties.RowNames=transpose(regions);
effsizes_G1_5_vs_G1 = array2table(effsizes_G1_5_vs_G1');effsizes_G1_5_vs_G1.Properties.RowNames=transpose(regions);
effsizes_G2_vs_G1_5 = array2table(effsizes_G2_vs_G1_5');effsizes_G2_vs_G1_5.Properties.RowNames=transpose(regions);
cd([location,'Effect sizes\'])
save('G1_vs_G0',"effsizes_G1_vs_G0")
save('G2_vs_G0',"effsizes_G2_vs_G0")
save('G2_vs_G1',"effsizes_G2_vs_G1")
save('G1_5_vs_G0',"effsizes_G1_5_vs_G0")
save('G1_5_vs_G1',"effsizes_G1_5_vs_G1")
save('G2_vs_G1_5',"effsizes_G2_vs_G1_5")



% For mapping (step THIRTEEN)
for ir=1:34
    mean_G0(ir)=mean(table2array(centiles_G0(:,ir)));
end

for ir=1:34
    mean_G1(ir)=mean(table2array(centiles_G1(:,ir)));
end

for ir=1:34
    mean_G1_5(ir)=mean(table2array(centiles_G1_5(:,ir)));
end

for ir=1:34
    mean_G2(ir)=mean(table2array(centiles_G2(:,ir)));
end

for ir=1:34
    mean_G2_without_ASRB(ir)=mean(table2array(centiles_G2_without_ASRB(:,ir)));
    mean_G2_without_BSNIP(ir)=mean(table2array(centiles_G2_without_BSNIP(:,ir)));
    mean_G2_without_LA5c(ir)=mean(table2array(centiles_G2_without_LA5c(:,ir)));
    mean_G2_without_MCIC(ir)=mean(table2array(centiles_G2_without_MCIC(:,ir)));
end

data=mean_G0';t_effsizes=array2table(data);t_effsizes.Properties.RowNames=transpose(regions);
writetable(t_effsizes,[location,'\Centiles\groups\mean_SCZ_and_SAD-relatives.csv'],'WriteRowNames',true);

data=mean_G1';t_effsizes=array2table(data);t_effsizes.Properties.RowNames=transpose(regions);
writetable(t_effsizes,[location,'\Centiles\groups\mean_PE.csv'],'WriteRowNames',true);

data=mean_G1_5';t_effsizes=array2table(data);t_effsizes.Properties.RowNames=transpose(regions);
writetable(t_effsizes,[location,'\Centiles\groups\mean_FEP.csv'],'WriteRowNames',true);

data=mean_G2';t_effsizes=array2table(data);t_effsizes.Properties.RowNames=transpose(regions);
writetable(t_effsizes,[location,'\Centiles\groups\mean_SCZ_and_SAD-chronic.csv'],'WriteRowNames',true);

data=mean_G2_without_ASRB';t_effsizes=array2table(data);t_effsizes.Properties.RowNames=transpose(regions);
writetable(t_effsizes,[location,'\Centiles\groups\mean_SCZ_and_SAD-chronic_without_ASRB.csv'],'WriteRowNames',true);

data=mean_G2_without_BSNIP';t_effsizes=array2table(data);t_effsizes.Properties.RowNames=transpose(regions);
writetable(t_effsizes,[location,'\Centiles\groups\mean_SCZ_and_SAD-chronic_without_BSNIP.csv'],'WriteRowNames',true);

data=mean_G2_without_LA5c';t_effsizes=array2table(data);t_effsizes.Properties.RowNames=transpose(regions);
writetable(t_effsizes,[location,'\Centiles\groups\mean_SCZ_and_SAD-chronic_without_LA5c.csv'],'WriteRowNames',true);

data=mean_G2_without_MCIC';t_effsizes=array2table(data);t_effsizes.Properties.RowNames=transpose(regions);
writetable(t_effsizes,[location,'\Centiles\groups\mean_SCZ_and_SAD-chronic_without_MCIC.csv'],'WriteRowNames',true);

data=effsizes_G1_vs_G0;t_effsizes=data;t_effsizes.Properties.RowNames=transpose(regions);
writetable(t_effsizes,[location,'\Effect sizes\effsizes_PE_vs_SCZ_and_SAD-relatives.csv'],'WriteRowNames',true);

data=effsizes_G2_vs_G0;t_effsizes=data;t_effsizes.Properties.RowNames=transpose(regions);
writetable(t_effsizes,[location,'\Effect sizes\effsizes_SCZ_and_SAD-chronic_vs_SCZ_and_SAD-relatives.csv'],'WriteRowNames',true);

data=effsizes_G2_vs_G1;t_effsizes=data;t_effsizes.Properties.RowNames=transpose(regions);
writetable(t_effsizes,[location,'\Effect sizes\effsizes_SCZ_and_SAD-chronic_vs_PE.csv'],'WriteRowNames',true);

data=effsizes_G1_5_vs_G0;t_effsizes=data;t_effsizes.Properties.RowNames=transpose(regions);
writetable(t_effsizes,[location,'\Effect sizes\effsizes_FEP_vs_SCZ_and_SAD-relatives.csv'],'WriteRowNames',true);

data=effsizes_G1_5_vs_G1;t_effsizes=data;t_effsizes.Properties.RowNames=transpose(regions);
writetable(t_effsizes,[location,'\Effect sizes\effsizes_FEP_vs_PE.csv'],'WriteRowNames',true);

data=effsizes_G2_vs_G1_5;t_effsizes=data;t_effsizes.Properties.RowNames=transpose(regions);
writetable(t_effsizes,[location,'\Effect sizes\effsizes_SCZ_and_SAD-chronic_vs_FEP.csv'],'WriteRowNames',true);

'DONE'
