
function [pval_spins, molecular_names, x_weights, Y_pred, correl, molecular_maps_labels_ordered] = NCP_10_CCA_cent_var(name,dx_title,location,type)
    warning('off','all')
    cd([location,'\PCA_CCA\',dx_title,'\'])
    
    X_rotated = csvread([location,'\PCA_CCA\perm_sphere_10000_DK.csv'])';
    X_rot = X_rotated(1:34,:);
    
    % Load dataset
    load("FEP.mat",'-mat');
    load("PE.mat",'-mat');
    load("SCZ.mat",'-mat');
    load("Schizoaffective_Disorder.mat",'-mat')
    load("G0.mat",'-mat');
    load("G1.mat",'-mat');
    load("G1_5.mat",'-mat')
    load("G2.mat",'-mat');
    load("G2_without_ASRB.mat",'-mat');
    load("G2_without_BSNIP.mat",'-mat');
    load("G2_without_LA5c.mat",'-mat');
    load("G2_without_MCIC.mat",'-mat');
    load("G1_vs_G0.mat",'-mat');
    load("G2_vs_G0.mat",'-mat');
    load("G2_vs_G1.mat",'-mat');
    load('G1_5_vs_G0','-mat')
    load('G1_5_vs_G1','-mat')
    load('G2_vs_G1_5','-mat')
    
    try
        data = eval(['centiles_', dx_title]); 
    catch
        data = rows2vars(eval(['effsizes_',dx_title])); data = data(:,2:end);
    end
    response = table2array(data);
    
    switch type
            case 'permutation'
                centiles_dx_1 = eval(['centiles_',dx_title(1:strfind(dx_title,'_vs')-1)]);
                centiles_dx_2 = eval(['centiles_',dx_title(strfind(dx_title,'_vs')+4:end)]);
                [centiles_perm_dx_1,centiles_perm_dx_2] = mix_dx(centiles_dx_1,centiles_dx_2);
                for ir=1:34
                    effsizes_perm(ir)=computeCohen_d(centiles_perm_dx_1(:,ir),centiles_perm_dx_2(:,ir));
                end
                response = effsizes_perm;          
     end
    
    % Load molecular maps
    molecular_maps = readtable([location,'\PCA_CCA\all_microsc_DesikanKilliany68.csv'],ReadVariableNames=true);
    molecular_names_table = readtable([location,'\Molecular_maps\molecular_names.xlsx'],ReadVariableNames=true);
    molecular_names = molecular_names_table{:,1};
    
    % Important! List of regions may follow a different order in each file. 
    % We need to set the same order for all of them, let´s say alphabetical
    y_labels = data.Properties.VariableNames;
    [~, y_labels_order] = sort(y_labels);
    molecular_maps_labels = molecular_maps.Var2(1:34);
    [molecular_maps_labels_ordered, molecular_maps_labels_order] = sort(molecular_maps_labels);
    
    X_rot_ordered = X_rot(molecular_maps_labels_order,:)';
    
    % Averaging across hemispheres
    molecular_map = molecular_maps{1:end,3:end};
    molecular_map_hemi = (molecular_map(1:34,:) + molecular_map(35:68,:))/2;
    
    molecular_map_hemi_ordered = molecular_map_hemi(molecular_maps_labels_order,:);
    writetable(array2table(molecular_map_hemi_ordered,"RowNames",sort(molecular_maps_labels),"VariableNames",molecular_names),[location,'\Molecular_maps\molecular_map_hemi_ordered.csv'],"WriteRowNames",true)
    response_ordered = response(:,y_labels_order);
    
    X = molecular_map_hemi_ordered;
    Y = mean(response_ordered,1)';
    
    
    %----- Analysis
    
    % Set path for analysis
    set_path;
    
    % Project folder
    cfg.dir.project = pwd();
    
    % Machine settings
    cfg.machine.name = 'cca';
    
    cfg.machine.metric = {'correl' 'trexvarx' 'trexvary'};
    
    cfg.machine.param.name = {'VARx', 'VARy'}; % explained variance by the PCA components
   
    cfg.machine.param.VARx = 0.6:0.1:0.9; % variance of data kept in the principal components during the SVD step of CCA, PCA-CCA or RCCA
    % % note that if variance is not sufficiently large, only very few (even 0 or 1) variables might be only kept
    cfg.machine.param.VARy = 1;
    
    cfg.machine.svd.varx = 1; % variance of X kept during the SVD step of CCA, PCA-CCA or RCCA. Default is 1 for CCA and 0.99 for RCCA. 
    % Note that if variance is not sufficiently large, only very few (even 0 or 1) variables might be only kept
    cfg.machine.svd.vary = 1; % variance of Y kept during the SVD step of CCA, PCA-CCA or RCCA. Default is 1 for CCA and 0.99 for RCCA. 
    % Note that if variance is not sufficiently large, only very few (even 0 or 1) variables might be only kept, i.e., for models with 1 output variable cfg.machine.svd.vary = 1 should be used
    
    cfg.machine.alignw = 'wX';

    % Framework settings
    cfg.frwork.name = 'permutation'; % In a holdout predictive (machine learning) framework, the data is divided into training and test sets 
    % by randomly subsampling subjects. In a 'permutation' descriptive framework, the data is not splitted, focusing on in-sample statistical evaluation
    
    cfg.frwork.split.nout = length(cfg.machine.param.(cfg.machine.param.name{1})); % number of outer splits/folds

    cfg.frwork.nlevel = 1;
    
    % Deflation settings
    cfg.defl.name = 'generalized'; % ['generalized', 'pls-projection', 'pls-modeA', 'pls-regression']
    % In the case of iterative methods, once a pair of weights is obtained, the corresponding
    % associative effect is removed from the data (by a process called deflation) and new associations are sought
    
    % Environment settings
    cfg.env.comp = 'local'; %  ['local', 'cluster']
    cfg.env.save.tableHeading = {'set' 'varx' 'correl' 'pval' 'npcax'};
    
    % Number of permutations
    switch type
        case 'real'
            cfg.stat.nperm = 1000;
            cfg.stat.nboot = 1000;

        case 'permutation'
            cfg.stat.nperm = 1;
            cfg.stat.nboot = 1;
    end

    % Run all VAR
    [var_real,pval_spins,~,res] = NCP_11_run_CCA(X,Y,cfg,molecular_names_table,dx_title,X_rot_ordered,type);
    pval_corr = mafdr(pval_spins,'BHFDR',true); % multipe comparisons correction (Benjamini-Hochberg False Discovery Rate)
    
    switch type
        case 'real'
            save(['pval_spins_',dx_title,'.mat'],'pval_spins');
            save(['pval_spins_corr_',dx_title,'.mat'],'pval_corr');
    end

    % Run best VAR
    pval_spins_0_05 = pval_corr < 0.05;
    if sum(pval_spins_0_05) == 0
        VAR = 0.6;
    elseif sum(pval_spins_0_05) == 1
        VAR = cfg.machine.param.VARx(find(pval_spins_0_05));
    elseif sum(pval_spins_0_05) > 1
        VAR = cfg.machine.param.VARx(min(find(pval_spins_0_05)));
    else
        error('VAR error')
    end
    switch type
        case 'real'
            save(['VAR_',dx_title,'.mat'],'VAR');
    end
    
    if isnumeric(VAR)
        cfg.machine.param.VARx = VAR;
        cfg.frwork.split.nout = length(cfg.machine.param.(cfg.machine.param.name{1}));
        [var_real,pval_spins,weights,res,Y_pred] = NCP_11_run_CCA(X,Y,cfg,molecular_names_table,dx_title,X_rot_ordered,type);
        x_weights = weights';
        Y_pred = Y_pred';
        correl = corr(Y_pred',Y);
    else
        x_weights = zeros(1,length(molecular_names));
        Y_pred = zeros(1,length(Y));
        correl = corr(Y_pred',Y);
    end

 end
