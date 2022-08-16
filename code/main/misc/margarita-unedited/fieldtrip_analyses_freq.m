%% suplementary experiment information
clear all
clc
master_dir = "C:\Users\marga\Desktop\Research Project\scripts";
main_path = 'D:\PhD\participant_';
addpath('C:\External_Software\fieldtrip-20210807');
addpath('C:\External_Software\spm12')
ft_defaults
%cd(master_dir);
experiment_type = 'onsets-2-8';
time_period = 1; % 0 = whole period, 1 = up to 800ms
%experiment_type = 'onsets-factor';
%experiment_type = 'partitions';
partitions_interaction_only = 0;
mean_centering = 1;
tails = 0;
region_of_interest = 0;
roi_applied = 'two-tailed';
%% get ROI for each frequency band
% get_region_of_interest_electrodes(stat_a,1,experiment_type,roi_applied,1);
% get_region_of_interest_electrodes(stat_b,1,experiment_type,roi_applied,2);
% get_region_of_interest_electrodes(stat_lg,1,experiment_type,roi_applied,3);
% get_region_of_interest_electrodes(stat_mg,1,experiment_type,roi_applied,4);
% get_region_of_interest_electrodes(stat_hg,1,experiment_type,roi_applied,5);

%% Are we looking at onsets2-8 or partitions
% set up the experiment as needed
if contains(experiment_type, 'onsets-2-8')
    data_file = 't1f1';
    regressor = 'ft_statfun_depsamplesT';
    n_participants = 40;
    partition.is_partition = 0;
    partition.partition_number = 0;
    
    [thin_data,med_data,thick_data,fieldtrip_data_agg,fieldtrip_data_PGI,participant_order] = from_fieldtrip_to_spm(n_participants,main_path,data_file,partition);
    
    frequency_data.thin1 = thin_data;
    frequency_data.med1 = med_data;
    frequency_data.thick1 = thick_data;
    frequency_data.data1 = fieldtrip_data_agg;
    frequency_data.participant_order_1 = participant_order;
    frequency_data.fieldtrip_data_PGI = fieldtrip_data_PGI;

    save("D:\PhD\margarita_data\frequency_data_mi.mat", 'frequency_data', '-v7.3')
    clear frequency_data;

    if time_period == 1
        load data_PGI
        load pax_order
    else
        load PGI_whole_period
        load whole_order
    end
    cd(master_dir);
    
    n_part = numel(data_PGI);
    design_matrix =  [1:n_part 1:n_part; ones(1,n_part) 2*ones(1,n_part)]; 
    
elseif contains(experiment_type, 'onsets-factor')
    data_file = 't1f1';
    regressor = 'ft_statfun_indepsamplesregrT';
    n_participants = 40;   
    regression_type = 'discomfort';
    partition.is_partition = 0;
    partition.partition_number = 0;
    
    
    
    [thin_data,med_data,thick_data,fieldtrip_data_agg,fieldtrip_data_PGI,participant_order] = from_fieldtrip_to_spm(n_participants,main_path,data_file,partition);
    
    frequency_data.thin1 = thin_data;
    frequency_data.med1 = med_data;
    frequency_data.thick1 = thick_data;
    frequency_data.data1 = fieldtrip_data_agg;
    frequency_data.participant_order_1 = participant_order;
    frequency_data.fieldtrip_data_PGI = fieldtrip_data_PGI;

    save("D:\PhD\margarita_data\frequency_data_mi.mat", 'frequency_data', '-v7.3')
    clear frequency_data;


    if time_period == 1
        %load data_PGI
        load data_partition1
        load pax_order
    else
        load PGI_whole_period
        load whole_order
    end
    cd(master_dir);
    data_PGI = data1;
    n_part = numel(data_PGI);
    [design_matrix, data_PGI] = create_design_matrix_partitions(participant_order, data_PGI, ...
        regression_type, 1);
    
    if region_of_interest == 1
        if contains(roi_applied, 'one-tailed')
            load('C:\Users\marga\Desktop\Research Project\scripts\one_tailed_roi_28_a.mat');
            load('C:\Users\marga\Desktop\Research Project\scripts\one_tailed_roi_28_b.mat');
            load('C:\Users\marga\Desktop\Research Project\scripts\one_tailed_roi_28_lg.mat');
            load('C:\Users\marga\Desktop\Research Project\scripts\one_tailed_roi_28_mg.mat');
            load('C:\Users\marga\Desktop\Research Project\scripts\one_tailed_roi_28_hg.mat');
        elseif contains(roi_applied, 'two-tailed')
            load('C:\Users\marga\Desktop\Research Project\scripts\two_tailed_roi_28_a.mat');
            load('C:\Users\marga\Desktop\Research Project\scripts\two_tailed_roi_28_b.mat');
            load('C:\Users\marga\Desktop\Research Project\scripts\two_tailed_roi_28_lg.mat');
            load('C:\Users\marga\Desktop\Research Project\scripts\two_tailed_roi_28_mg.mat');
            load('C:\Users\marga\Desktop\Research Project\scripts\two_tailed_roi_28_hg.mat');
        end
    end
    
elseif contains(experiment_type, 'partitions')
    
    data_file = 't1f1';
    regressor = 'ft_statfun_indepsamplesregrT';
    n_participants = 40;
    partition1.is_partition = 1; % partition 1
    partition1.partition_number = '1';
    partition2.is_partition = 1; % partition 2
    partition2.partition_number = '2';
    partition3.is_partition = 1; % partition 3
    partition3.partition_number = '3';
    
     [thin_data,med_data,thick_data,fieldtrip_data_agg,fieldtrip_data_PGI,participant_order] = from_fieldtrip_to_spm(n_participants,main_path,data_file,partition1);

    frequency_data.thin1 = thin_data;
    frequency_data.med1 = med_data;
    frequency_data.thick1 = thick_data;
    frequency_data.data1 = fieldtrip_data_agg;
    frequency_data.participant_order_1 = participant_order;
    frequency_data.fieldtrip_data_PGI = fieldtrip_data_PGI;

    save("D:\PhD\margarita_data\frequency_data.mat", 'frequency_data', '-v7.3')
    clear frequency_data;
% 
% 
%     [thin_data2,med_data2,thick_data2,fieldtrip_data_agg2,fieldtrip_data_PGI2,participant_order2]  = from_fieldtrip_to_spm(n_participants,main_path,data_file,partition2);
% 
%     frequency_data.thin1 = thin_data2;
%     frequency_data.med1 = med_data2;
%     frequency_data.thick1 = thick_data2;
%     frequency_data.data1 = fieldtrip_data_agg2;
%     frequency_data.participant_order_1 = participant_order2;
%     frequency_data.fieldtrip_data_PGI = fieldtrip_data_PGI2;
% 
%     save("D:\PhD\margarita_data\frequency_data2.mat", 'frequency_data', '-v7.3')
%     clear frequency_data;
% 
%     [thin_data3,med_data3,thick_data3,fieldtrip_data_agg3,fieldtrip_data_PGI3,participant_order3]  = from_fieldtrip_to_spm(n_participants,main_path,data_file,partition3);
% 
% 
%     frequency_data.thin1 = thin_data3;
%     frequency_data.med1 = med_data3;
%     frequency_data.thick1 = thick_data3;
%     frequency_data.data1 = fieldtrip_data_agg3;
%     frequency_data.participant_order_1 = participant_order3;
%     frequency_data.fieldtrip_data_PGI = fieldtrip_data_PGI3;
% 
%     save("D:\PhD\margarita_data\frequency_data3.mat", 'frequency_data', '-v7.3')
%     clear frequency_data;
% 
%     cd("C:\Users\marga\Desktop\Research Project\participants_Austyn")
    load D:\PhD\margarita_data\frequency_data.mat
    participant_order_1 = frequency_data.participant_order_1;
    data1 = frequency_data.fieldtrip_data_PGI;
    disp("loaded p1")
    clear frequency_data;

    load D:\PhD\margarita_data\frequency_data2.mat
    participant_order_2 = frequency_data.participant_order_1;
    data2 = frequency_data.fieldtrip_data_PGI;
    disp("loaded p2")
    clear frequency_data;

    load D:\PhD\margarita_data\frequency_data3.mat
    participant_order_3 = frequency_data.participant_order_1;
    data3 = frequency_data.fieldtrip_data_PGI;
    disp("loaded p3")
    clear frequency_data;
    
    
    if partitions_interaction_only == 0
        % choose from 'headache', 'visual stress', 'discomfort'
        regression_type = 'visual stress'; 
        partition = 1;
        [design1, new_participants1] = create_design_matrix_partitions(participant_order_1, data1, ...
                regression_type, partition);
        partition = 2;
        [design2, new_participants2] = create_design_matrix_partitions(participant_order_2, data2, ...
                regression_type, partition);
        partition = 3;
        [design3, new_participants3] = create_design_matrix_partitions(participant_order_3, data3, ...
                regression_type, partition);

        data = [new_participants1, new_participants2, new_participants3];
        null_data = set_values_to_zero(data); % create null data to hack a t-test
        design_matrix = [design1, design2, design3];
        n_part = numel(new_participants1);
        if mean_centering == 1
            design_matrix = design_matrix - mean(design_matrix);
        end
    else
        design_matrix_a(1:length(data1)) = 2.72;
        design_matrix_b(1:length(data2)) = 1.65;
        design_matrix_c(1:length(data3)) = 1.00;
        design_matrix = [design_matrix_a, design_matrix_b, design_matrix_c];
        data = [data1, data2, data3];
        null_data = set_values_to_zero(data); % create null data to hack a t-test
        n_part = numel(data1);
        if mean_centering == 1
            design_matrix = design_matrix - mean(design_matrix);
        end
    end
    %plot_design_matrix(design_matrix, n_part)
end


%% TFR of aggregated data grand average
% cd("C:\Users\marga\Desktop\Research Project\participants_Austyn")
% load data_agg
% cd(master_dir);
% 
% cfg = [];
% cfg.foilim = [30 80];
% cfg.toilim = [-0.2 0.8];
% cfg.parameter = 'powspctrm';
% grand_avg = ft_freqgrandaverage(cfg, thick{:});
% grand_avg.elec = data_agg{1}.elec;
% 
% cfg = [];
% %cfg.layout = data_agg{1}.elec;
% cfg.channel = 'A23';
% cfg.baseline = [-0.1 0];
% cfg.baselinetype = 'db';
% cfg.xlim = [-0.1 0.8];
% cfg.zlim = [-18 18];
% cfg.colorbar = 'yes';
% cfg.colormap = 'jet';
% cfg.showlabels = 'yes';
% cfg.parameter = 'powspctrm';
% ft_singleplotTFR(cfg, grand_avg);
% % ft_topoplotTFR(cfg, grand_avg);

%% setup FT analysis
% we have to switch to SPM8 to use some of the functions in FT
%cd("C:\Users\marga\Documents\MATLAB");
%rmpath spm12;
%addpath spm8;
%cd(master_dir);
%data=data_PGI;
%data=old_data;

% we need to tell fieldtrip how our electrodes are structured
cfg = [];
cfg.feedback = 'no';
cfg.method = 'distance';
cfg.elec = data{1}.elec;
neighbours = ft_prepare_neighbours(cfg);

% all experiment configurations
start_latency = 0;
end_latency = 0.8;
cfg = [];
cfg.latency = [start_latency end_latency];
cfg.channel = 'eeg';
cfg.avgoverfreq = 'yes';
cfg.statistic = regressor;
cfg.method = 'montecarlo';
cfg.correctm = 'cluster';
cfg.neighbours = neighbours;
cfg.clusteralpha = 0.025;
cfg.numrandomization = 1000;
cfg.alpha = 0.05;
cfg.tail = tails; % 0 = two-tailed, 1 = one-tailed
cfg.correcttail = 'alpha';
cfg.design = design_matrix;
cfg.computeprob = 'yes';

%% run the fieldtrip analyses
if contains(experiment_type, 'onsets-2-8')
    % set values to 0 for viz and other things
    null_data = set_values_to_zero(data); % create null data to hack a t-test
    cfg.uvar = 1;
    cfg.ivar = 2;
    cfg.frequency = [8 13];
    stat_a = ft_freqstatistics(cfg, data{:}, null_data{:});
    cfg.frequency = [20 35];
    stat_b = ft_freqstatistics(cfg, data{:}, null_data{:});
    cfg.frequency = [30 45];
    stat_lg = ft_freqstatistics(cfg, data{:}, null_data{:});
    cfg.frequency = [45 60];
    stat_mg = ft_freqstatistics(cfg, data{:}, null_data{:});
    cfg.frequency = [60 80];
    stat_hg = ft_freqstatistics(cfg, data{:}, null_data{:});
elseif contains(experiment_type, 'onsets-factor')
    cfg.ivar = 1;
    if region_of_interest == 1
        cfg.frequency = [8 13];
        roi = roi_a;
        data_PGI_a = create_hacked_roi(data_PGI, roi);
        stat_a = ft_freqstatistics(cfg, data_PGI_a{:});
        cfg.frequency = [20 35];
        roi = roi_b;
        data_PGI_b = create_hacked_roi(data_PGI, roi);
        stat_b = ft_freqstatistics(cfg, data_PGI_b{:});
        cfg.frequency = [30 45];
        roi = roi_lg;
        data_PGI_lg = create_hacked_roi(data_PGI, roi);
        stat_lg = ft_freqstatistics(cfg, data_PGI_lg{:});
        cfg.frequency = [45 60];
        roi = roi_mg;
        data_PGI_mg = create_hacked_roi(data_PGI, roi);
        stat_mg = ft_freqstatistics(cfg, data_PGI_mg{:});
        cfg.frequency = [60 80];
        roi = roi_hg;
        data_PGI_hg = create_hacked_roi(data_PGI, roi);
        stat_hg = ft_freqstatistics(cfg, data_PGI_hg{:});
    else
        cfg.frequency = [8 13];
        stat_a = ft_freqstatistics(cfg, data_PGI{:});
        cfg.frequency = [20 35];
        stat_b = ft_freqstatistics(cfg, data_PGI{:});
        cfg.frequency = [30 45];
        stat_lg = ft_freqstatistics(cfg, data_PGI{:});
        cfg.frequency = [45 60];
        stat_mg = ft_freqstatistics(cfg, data_PGI{:});
        cfg.frequency = [60 80];
        stat_hg = ft_freqstatistics(cfg, data_PGI{:});
    end
elseif contains(experiment_type, 'partitions')
    cfg.ivar = 1;
    cfg.frequency = [8 13];
    stat_a = ft_freqstatistics(cfg, data{:});
    cfg.frequency = [20 35];
    stat_b = ft_freqstatistics(cfg, data{:});
    cfg.frequency = [30 45];
    stat_lg = ft_freqstatistics(cfg, data{:});
    cfg.frequency = [45 60];
    stat_mg = ft_freqstatistics(cfg, data{:});
    cfg.frequency = [60 80];
    stat_hg = ft_freqstatistics(cfg, data{:});
end

%% get peak level stats
if tails==1
    % 1 for the most positive going cluster
    [peak_stats_a_pos, all_electrode_stata_pos] = get_peak_level_stats(stat_a, 1, 'positive');
    [peak_stats_b_pos, all_electrode_statb_pos] = get_peak_level_stats(stat_b, 1, 'positive');
    [peak_stats_lg_pos, all_electrode_statlg_pos] = get_peak_level_stats(stat_lg, 1, 'positive');
    [peak_stats_mg_pos, all_electrode_statmg_pos] = get_peak_level_stats(stat_mg, 1, 'positive');
    [peak_stats_hg_pos, all_electrode_stathg_pos] = get_peak_level_stats(stat_hg, 1, 'positive');
else
    % 1 for the most positive going cluster
    [peak_stats_a_pos, all_electrode_stata_pos] = get_peak_level_stats(stat_a, 1, 'positive');
    [peak_stats_b_pos, all_electrode_statb_pos] = get_peak_level_stats(stat_b, 1, 'positive');
    [peak_stats_lg_pos, all_electrode_statlg_pos] = get_peak_level_stats(stat_lg, 1, 'positive');
    [peak_stats_mg_pos, all_electrode_statmg_pos] = get_peak_level_stats(stat_mg, 1, 'positive');
    [peak_stats_hg_pos, all_electrode_stathg_pos] = get_peak_level_stats(stat_hg, 1, 'positive');
    
    %most negative going cluster
    [peak_stats_a_neg, all_electrode_stata_neg] = get_peak_level_stats(stat_a, 1, 'negative');
    [peak_stats_b_neg, all_electrode_statb_neg] = get_peak_level_stats(stat_b, 1, 'negative');
    [peak_stats_lg_neg, all_electrode_statlg_neg] = get_peak_level_stats(stat_lg, 1, 'negative');
    [peak_stats_mg_neg, all_electrode_statmg_neg] = get_peak_level_stats(stat_mg, 1, 'negative');
    [peak_stats_hg_neg, all_electrode_stathg_neg] = get_peak_level_stats(stat_hg, 1, 'negative');
end

%% function that plots the t values through time and decides whcih electrode to plot
if numel(all_electrode_stata_pos) > 0
    peak_stats_a_pos = compute_best_electrode_from_t_values(stat_a,all_electrode_stata_pos,'positive', peak_stats_a_pos);
end
if numel(all_electrode_statb_pos) > 0
    peak_stats_b_pos = compute_best_electrode_from_t_values(stat_b,all_electrode_statb_pos,'positive', peak_stats_b_pos);
end
if numel(all_electrode_statlg_pos) > 0
    peak_stats_lg_pos = compute_best_electrode_from_t_values(stat_lg,all_electrode_statlg_pos,'positive', peak_stats_lg_pos);
end
if numel(all_electrode_statmg_pos) > 0
    peak_stats_mg_pos = compute_best_electrode_from_t_values(stat_mg,all_electrode_statmg_pos,'positive', peak_stats_mg_pos);
end
if numel(all_electrode_stathg_pos) > 0
    peak_stats_hg_pos = compute_best_electrode_from_t_values(stat_hg,all_electrode_stathg_pos,'positive', peak_stats_hg_pos);
end

if numel(all_electrode_stata_neg) > 0
    peak_stats_a_neg = compute_best_electrode_from_t_values(stat_a,all_electrode_stata_neg,'negative', peak_stats_a_neg);
end
if numel(all_electrode_statb_neg) > 0
    peak_stats_b_neg = compute_best_electrode_from_t_values(stat_b,all_electrode_statb_neg,'negative', peak_stats_b_neg);
end
if numel(all_electrode_statlg_neg) > 0
    peak_stats_lg_neg = compute_best_electrode_from_t_values(stat_lg,all_electrode_statlg_neg,'negative', peak_stats_lg_neg);
end
if numel(all_electrode_statmg_neg) > 0
    peak_stats_mg_neg = compute_best_electrode_from_t_values(stat_mg,all_electrode_statmg_neg,'negative', peak_stats_mg_neg);
end
if numel(all_electrode_stathg_neg) > 0
    peak_stats_hg_neg = compute_best_electrode_from_t_values(stat_hg,all_electrode_stathg_neg,'negative', peak_stats_hg_neg);
end

%% get topoplots through time
%data = data_PGI;
tail = 'positive';
figure(1)
create_topoplots(data,stat_a,tail);
figure(2)
create_topoplots(data,stat_b,tail);
figure(3)
create_topoplots(data,stat_lg,tail);
figure(4)
create_topoplots(data,stat_mg,tail);
figure(5)
create_topoplots(data,stat_hg,tail);

%% get cluster size through time
make_plots = 'yes';
xlim = end_latency*1000;
if tails == 1
    title = 'Most positive going cluster through time as a % of entire volume';
    figure(6)
    calculate_cluster_size(stat_a, 1, make_plots, title, xlim, 'positive');
    figure(7)
    calculate_cluster_size(stat_b, 1, make_plots, title, xlim, 'positive');
    figure(8)
    calculate_cluster_size(stat_lg, 1, make_plots, title, xlim, 'positive');
    figure(9)
    calculate_cluster_size(stat_mg, 1, make_plots, title, xlim, 'positive');
    figure(10)
    calculate_cluster_size(stat_hg, 1, make_plots, title, xlim, 'positive');
else
    %title = 'Most positive (up) and negative (down) going clusters through time as a % of entire volume';
    title = '';
    figure(1)
    %subplot(2,1,1)
%     calculate_cluster_size(stat_a, 1, make_plots, title, xlim, 'positive');
%     subplot(2,1,2)
%     calculate_cluster_size(stat_a, 1, make_plots, '', xlim, 'negative');
%     figure(2)
%     subplot(2,1,1)
%     calculate_cluster_size(stat_b, 1, make_plots, title, xlim, 'positive');
%     subplot(2,1,2)
%     calculate_cluster_size(stat_b, 1, make_plots, '', xlim, 'negative');
    figure(3)
    calculate_cluster_size(stat_lg, 1, make_plots, title, xlim, 'negative');
%     figure(4)
%     subplot(2,1,1)
%     calculate_cluster_size(stat_mg, 1, make_plots, title, xlim, 'positive');
%     subplot(2,1,2)
%     calculate_cluster_size(stat_mg, 1, make_plots, '', xlim, 'negative');
%     figure(5)
%     subplot(2,1,1)
%     calculate_cluster_size(stat_hg, 1, make_plots, title, xlim, 'positive');
%     subplot(2,1,2)
%     calculate_cluster_size(stat_hg, 1, make_plots, '', xlim, 'negative');
end

%% plot PGI and conditions, mean/intercept-onsets only
cd("C:\Users\marga\Desktop\Research Project\participants_Austyn")
load thin
load med
load thick
cd(master_dir)

pause('on')
if tails==1
    chosen_electrode={peak_stats_a_pos,peak_stats_b_pos,peak_stats_lg_pos,peak_stats_mg_pos,peak_stats_hg_pos};
    onsets_time_series(thin, med, thick, data_PGI, chosen_electrode,end_latency);
else
    %positive
    chosen_electrode={peak_stats_a_pos,peak_stats_b_pos,peak_stats_lg_pos,peak_stats_mg_pos,peak_stats_hg_pos};
    onsets_time_series(thin, med, thick, data_PGI, chosen_electrode,end_latency);
    pause;
%     %negative
%     chosen_electrode={peak_stats_a_neg,peak_stats_b_neg,peak_stats_lg_neg,peak_stats_mg_neg,peak_stats_hg_neg};
%     onsets_time_series(thin, med, thick, data_PGI, chosen_electrode,end_latency);
end

%% plot median split for factors, onsets only
pause('on')
if tails==1
    chosen_electrode={peak_stats_a_pos,peak_stats_b_pos,peak_stats_lg_pos,peak_stats_mg_pos,peak_stats_hg_pos};
    median_split_onsets_factor(master_dir, data_PGI, participant_order,chosen_electrode,regression_type,end_latency);
else
%     %positive
%     chosen_electrode={peak_stats_a_pos,peak_stats_b_pos,peak_stats_lg_pos,peak_stats_mg_pos,peak_stats_hg_pos};
%     median_split_onsets_factor(master_dir, data_PGI, participant_order,chosen_electrode,regression_type,end_latency);
%     pause;
    %negative
    chosen_electrode={peak_stats_a_neg,peak_stats_b_neg,peak_stats_lg_neg,peak_stats_mg_neg,peak_stats_hg_neg};
    median_split_onsets_factor(master_dir, data_PGI, participant_order,chosen_electrode,regression_type,end_latency)
end

%% difference of PGI times series after median split
% PGI_alpha_diff = PGI_alpha_high - PGI_alpha_low;
% PGI_beta_diff = PGI_beta_high - PGI_beta_low;
% PGI_lgamma_diff = PGI_lgamma_high - PGI_lgamma_low;
% PGI_mgamma_diff = PGI_mgamma_high - PGI_mgamma_low;
% PGI_hgamma_diff = PGI_hgamma_high - PGI_hgamma_low;
% 
% figure(1)
% subplot(2,1,1)
% plot(time,PGI_alpha_diff(19,:),'Color','#7E2F8E','LineWidth',1)
% xlim([-0.1 0.8]);title('Difference of median split PGI time series @A23, Alpha band');
% grid on
% subplot(2,1,2)
% plot(time,PGI_alpha_diff(95,:),'Color','#7E2F8E','LineWidth',1)
% xlim([-0.1 0.8]);title('Difference of median split PGI time series @D9, Alpha band');
% grid on
% 
% figure(2)
% subplot(2,1,1)
% plot(time,PGI_beta_diff(19,:),'Color','#7E2F8E','LineWidth',1)
% xlim([-0.1 0.8]);title('Difference of median split PGI time series @A23, Beta band');
% grid on
% subplot(2,1,2)
% plot(time,PGI_beta_diff(16,:),'Color','#7E2F8E','LineWidth',1)
% xlim([-0.1 0.8]);title('Difference of median split PGI time series @A20, Beta band');
% grid on
% 
% figure(3)
% subplot(2,1,1)
% plot(time,PGI_lgamma_diff(19,:),'Color','#7E2F8E','LineWidth',1)
% xlim([-0.1 0.8]);title('Difference of median split PGI time series @A23, Low gamma band');
% grid on
% subplot(2,1,2)
% plot(time,PGI_lgamma_diff(65,:),'Color','#7E2F8E','LineWidth',1)
% xlim([-0.1 0.8]);title('Difference of median split PGI time series @C11, Low gamma band');
% grid on
% 
% figure(4)
% subplot(2,1,1)
% plot(time,PGI_mgamma_diff(19,:),'Color','#7E2F8E','LineWidth',1)
% xlim([-0.1 0.8]);title('Difference of median split PGI time series @A23, Medium gamma band');
% grid on
% subplot(2,1,2)
% plot(time,PGI_mgamma_diff(37,:),'Color','#7E2F8E','LineWidth',1)
% xlim([-0.1 0.8]);title('Difference of median split PGI time series @B15, Medium gamma band');
% grid on
% 
% figure(5)
% subplot(2,1,1)
% plot(time,PGI_hgamma_diff(19,:),'Color','#7E2F8E','LineWidth',1)
% xlim([-0.1 0.8]);title('Difference of median split PGI time series @A23, High gamma band');
% grid on
% subplot(2,1,2)
% plot(time,PGI_hgamma_diff(12,:),'Color','#7E2F8E','LineWidth',1)
% xlim([-0.1 0.8]);title('Difference of median split PGI time series @A16, High gamma band');
% grid on

%% plot median split for factors, partitions
%calculate median split for PGI and conditions
[PGI_p1,PGI_p2,PGI_p3,thin_p1,thin_p2,thin_p3,...
    med_p1,med_p2,med_p3,thick_p1,thick_p2,thick_p3] = median_split_partitions_factor(master_dir,...
    data1,data2,data3,participant_order_1,regression_type,end_latency);

%% plotting for median split partitions
pause('on')
if tails == 1
    chosen_electrode={peak_stats_a_pos,peak_stats_b_pos,peak_stats_lg_pos,peak_stats_mg_pos,peak_stats_hg_pos};
    time_series_partitions_median_split(data,PGI_p1,PGI_p2,PGI_p3,thin_p1,thin_p2,thin_p3,...
        med_p1,med_p2,med_p3,thick_p1,thick_p2,thick_p3,chosen_electrode,regression_type,end_latency);
else
    %positive
    chosen_electrode={peak_stats_a_pos,peak_stats_b_pos,peak_stats_lg_pos,peak_stats_mg_pos,peak_stats_hg_pos};
    time_series_partitions_median_split(data,PGI_p1,PGI_p2,PGI_p3,thin_p1,thin_p2,thin_p3,...
        med_p1,med_p2,med_p3,thick_p1,thick_p2,thick_p3,chosen_electrode,regression_type,end_latency);
%     pause; 
%     %negative
%     chosen_electrode={peak_stats_a_neg,peak_stats_b_neg,peak_stats_lg_neg,peak_stats_mg_neg,peak_stats_hg_neg};
%     time_series_partitions_median_split(data,PGI_p1,PGI_p2,PGI_p3,thin_p1,thin_p2,thin_p3,...
%         med_p1,med_p2,med_p3,thick_p1,thick_p2,thick_p3,chosen_electrode,regression_type,end_latency);
end

%% time series plotting for 3 partitions, PGI & conditions, partition interaction only
cd("C:\Users\marga\Desktop\Research Project\participants_Austyn")
load thin_partition1
load med_partition1
load thick_partition1
load thin_partition2
load med_partition2
load thick_partition2
load thin_partition3
load med_partition3
load thick_partition3
cd(master_dir);
pause('on')
if tails == 1
    chosen_electrode={peak_stats_a_pos,peak_stats_b_pos,peak_stats_lg_pos,peak_stats_mg_pos,peak_stats_hg_pos};
    time_series_partitions(data1,data2,data3,thin1,thin2,thin3,med1,med2,med3,thick1,thick2,thick3,chosen_electrode,end_latency);
else
    %positive
    chosen_electrode={peak_stats_a_pos,peak_stats_b_pos,peak_stats_lg_pos,peak_stats_mg_pos,peak_stats_hg_pos};
    time_series_partitions(data1,data2,data3,thin1,thin2,thin3,med1,med2,med3,thick1,thick2,thick3,chosen_electrode,end_latency);
%     pause; 
%     %negative
%     chosen_electrode={peak_stats_a_neg,peak_stats_b_neg,peak_stats_lg_neg,peak_stats_mg_neg,peak_stats_hg_neg};
%     time_series_partitions(data1,data2,data3,thin1,thin2,thin3,med1,med2,med3,thick1,thick2,thick3,chosen_electrode,end_latency);
end

%% plottint time series for onsets 2-8
function onsets_time_series(thin, med, thick, data_PGI,chosen_electrode,end_latency)
    
    PGI_grand_avgs = calculate_freq_avg(data_PGI,end_latency);
    thin_grand_avgs = calculate_avg_conditions(thin,end_latency);
    med_grand_avgs = calculate_avg_conditions(med,end_latency);
    thick_grand_avgs = calculate_avg_conditions(thick,end_latency);

    PGI_alpha = squeeze(mean(PGI_grand_avgs{1}.powspctrm,2,'omitnan'));
    PGI_beta = squeeze(mean(PGI_grand_avgs{2}.powspctrm,2,'omitnan'));
    PGI_lgamma = squeeze(mean(PGI_grand_avgs{3}.powspctrm,2,'omitnan'));
    PGI_mgamma = squeeze(mean(PGI_grand_avgs{4}.powspctrm,2,'omitnan'));
    PGI_hgamma = squeeze(mean(PGI_grand_avgs{5}.powspctrm,2,'omitnan'));

    thin_alpha = squeeze(mean(thin_grand_avgs{1}.powspctrm,2,'omitnan'));
    thin_beta = squeeze(mean(thin_grand_avgs{2}.powspctrm,2,'omitnan'));
    thin_lgamma = squeeze(mean(thin_grand_avgs{3}.powspctrm,2,'omitnan'));
    thin_mgamma = squeeze(mean(thin_grand_avgs{4}.powspctrm,2,'omitnan'));
    thin_hgamma = squeeze(mean(thin_grand_avgs{5}.powspctrm,2,'omitnan'));

    med_alpha = squeeze(mean(med_grand_avgs{1}.powspctrm,2,'omitnan'));
    med_beta = squeeze(mean(med_grand_avgs{2}.powspctrm,2,'omitnan'));
    med_lgamma = squeeze(mean(med_grand_avgs{3}.powspctrm,2,'omitnan'));
    med_mgamma = squeeze(mean(med_grand_avgs{4}.powspctrm,2,'omitnan'));
    med_hgamma = squeeze(mean(med_grand_avgs{5}.powspctrm,2,'omitnan'));

    thick_alpha = squeeze(mean(thick_grand_avgs{1}.powspctrm,2,'omitnan'));
    thick_beta = squeeze(mean(thick_grand_avgs{2}.powspctrm,2,'omitnan'));
    thick_lgamma = squeeze(mean(thick_grand_avgs{3}.powspctrm,2,'omitnan'));
    thick_mgamma = squeeze(mean(thick_grand_avgs{4}.powspctrm,2,'omitnan'));
    thick_hgamma = squeeze(mean(thick_grand_avgs{5}.powspctrm,2,'omitnan'));

    %plot results for all frequency bands (different electrodes)
    fig = figure(1);
    time = thick_grand_avgs{1}.time;
    elec_label = chosen_electrode{1}.electrode;
    elec_id = find(strcmp(data_PGI{1}.label,elec_label));
    subplot(2,1,1)
    plot(time,PGI_alpha(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-1 4]); grid on
    legend('PGI');
    subplot(2,1,2)
    hold on
    plot(time,thin_alpha(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med_alpha(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick_alpha(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-1 4]); grid on
    legend('Thin','Medium','Thick'); 
    hold off

    %use same labels
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,'Power (in dB)'); 
    xlabel(han,'Time (in sec)');
    mtitle = strcat('Time series of power @',elec_label,' ,Onsets 2-8, Alpha band');
    title(han,mtitle);

    fig = figure(2);
    time = thick_grand_avgs{1}.time;
    elec_label = chosen_electrode{2}.electrode;
    elec_id = find(strcmp(data_PGI{1}.label,elec_label));
    subplot(2,1,1)
    plot(time,PGI_beta(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-1 4]); grid on
    legend('PGI');
    subplot(2,1,2)
    hold on
    plot(time,thin_beta(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med_beta(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick_beta(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-1 4]); grid on
    legend('Thin','Medium','Thick');
    hold off
    
    %use same labels
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,'Power (in dB)'); 
    xlabel(han,'Time (in sec)');
    mtitle = strcat('Time series of power @',elec_label,' ,Onsets 2-8, Beta band');
    title(han,mtitle);

    fig = figure(3);
    time = thick_grand_avgs{1}.time;
    elec_label = chosen_electrode{3}.electrode;
    elec_id = find(strcmp(data_PGI{1}.label,elec_label));
    subplot(2,1,1)
    plot(time,PGI_lgamma(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-1 4]); grid on
    legend('PGI');
    subplot(2,1,2)
    hold on
    plot(time,thin_lgamma(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med_lgamma(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick_lgamma(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-1 4]); grid on
    legend('Thin','Medium','Thick'); 
    hold off
    
    %use same labels
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,'Power (in dB)'); 
    xlabel(han,'Time (in sec)');
    mtitle = strcat('Time series of power @',elec_label,' ,Onsets 2-8, Low gamma band');
    title(han,mtitle);

    fig = figure(4);
    time = thick_grand_avgs{1}.time;
    elec_label = chosen_electrode{4}.electrode;
    elec_id = find(strcmp(data_PGI{1}.label,elec_label));
    subplot(2,1,1)
    plot(time,PGI_mgamma(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-1 4]); grid on
    legend('PGI');
    subplot(2,1,2)
    hold on
    plot(time,thin_mgamma(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med_mgamma(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick_mgamma(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-1 4]); grid on
    legend('Thin','Medium','Thick');
    hold off
    
    %use same labels
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,'Power (in dB)'); 
    xlabel(han,'Time (in sec)');
    mtitle = strcat('Time series of power @',elec_label,' ,Onsets 2-8, Medium gamma band');
    title(han,mtitle);

    fig = figure(5);
    time = thick_grand_avgs{1}.time;
    elec_label = chosen_electrode{5}.electrode;
    elec_id = find(strcmp(data_PGI{1}.label,elec_label));
    subplot(2,1,1)
    plot(time,PGI_hgamma(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-1 4]); grid on
    legend('PGI');
    subplot(2,1,2)
    hold on
    plot(time,thin_hgamma(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med_hgamma(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick_hgamma(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-1 4]); grid on
    legend('Thin','Medium','Thick'); 
    hold off
    
    %use same labels
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,'Power (in dB)'); 
    xlabel(han,'Time (in sec)');
    mtitle = strcat('Time series of power @',elec_label,' ,Onsets 2-8, High gamma band');
    title(han,mtitle);
end

%% plotting time series for partitions, median split
function time_series_partitions_median_split(data,data1,data2,data3,thin1,thin2,thin3,med1,med2,med3,thick1,thick2,thick3,chosen_electrode,regression_type,end_latency)
    %plot results for all frequency bands (different electrodes)
    time = data{1}.time(1:451);
    fig = figure(1);
    elec_label = chosen_electrode{1}.electrode;
    elec_id = find(strcmp(data{1}.label,elec_label));
    subplot(5,2,1)
    hold on
    plot(time,data1{1,1}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,data2{1,1}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,data3{1,1}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-3 3]); 
    grid on
    legend('PGI partition 1','PGI partition 2','PGI partition 3'); title(strcat('Low ',{' '},regression_type))
    hold off
    subplot(5,2,2)
    hold on
    plot(time,data1{2,1}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,data2{2,1}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,data3{2,1}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-3 3]); 
    grid on
    legend('PGI partition 1','PGI partition 2','PGI partition 3'); title(strcat('High ',{' '},regression_type))
    hold off
    subplot(5,2,3)
    hold on
    plot(time,med1{1,1}(elec_id,:),'Color','#A2142F','LineWidth',0.8)
    plot(time,med2{1,1}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,med3{1,1}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-3 3]); 
    grid on
    legend('Medium P1','Medium P2','Medium P3'); title(strcat('Low ',{' '},regression_type))
    hold off
    subplot(5,2,4)
    hold on
    plot(time,med1{2,1}(elec_id,:),'Color','#A2142F','LineWidth',0.8)
    plot(time,med2{2,1}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,med3{2,1}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-3 3]); 
    grid on
    legend('Medium P1','Medium P2','Medium P3'); title(strcat('High ',{' '},regression_type))
    hold off
    subplot(5,2,5)
    hold on
    plot(time,thin1{1,1}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med1{1,1}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick1{1,1}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-4 2]); 
    grid on
    legend('Thin','Medium','Thick'); title(strcat('Partition 1,','Low ',{' '},regression_type));
    hold off
    subplot(5,2,6)
    hold on
    plot(time,thin1{2,1}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med1{2,1}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick1{2,1}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-4 2]); 
    grid on
    legend('Thin','Medium','Thick'); title(strcat('Partition 1,','High ',{' '},regression_type));
    hold off
    subplot(5,2,7)
    hold on
    plot(time,thin2{1,1}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med2{1,1}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick2{1,1}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-4 2]); 
    grid on
    legend('Thin','Medium','Thick'); title(strcat('Partition 2,','Low ',{' '},regression_type));
    hold off
    subplot(5,2,8)
    hold on
    plot(time,thin2{2,1}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med2{2,1}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick2{2,1}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-4 2]); 
    grid on
    legend('Thin','Medium','Thick'); title(strcat('Partition 2,','High ',{' '},regression_type));
    hold off
    subplot(5,2,9)
    hold on
    plot(time,thin3{1,1}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med3{1,1}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick3{1,1}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-4 2]); 
    grid on
    legend('Thin','Medium','Thick'); title(strcat('Partition 3,','Low ',{' '},regression_type));
    hold off
    subplot(5,2,10)
    hold on
    plot(time,thin3{2,1}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med3{2,1}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick3{2,1}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-4 2]); 
    grid on
    legend('Thin','Medium','Thick'); title(strcat('Partition 3,','High ',{' '},regression_type));
    hold off

    %use same labels
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,'Power (in dB)'); 
    xlabel(han,'Time (in sec)');
    mtitle = strcat('Time series of power @',elec_label,' ,Partitions, Alpha band');
    title(han,mtitle);

    fig = figure(2);
    elec_label = chosen_electrode{2}.electrode;
    elec_id = find(strcmp(data{1}.label,elec_label));
    subplot(5,2,1)
    hold on
    plot(time,data1{1,2}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,data2{1,2}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,data3{1,2}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-1 2]); 
    grid on
    legend('PGI partition 1','PGI partition 2','PGI partition 3'); title(strcat('Low ',{' '},regression_type))
    hold off
    subplot(5,2,2)
    hold on
    plot(time,data1{2,2}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,data2{2,2}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,data3{2,2}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-1 2]); 
    grid on
    legend('PGI partition 1','PGI partition 2','PGI partition 3'); title(strcat('High ',{' '},regression_type))
    hold off
    subplot(5,2,3)
    hold on
    plot(time,med1{1,2}(elec_id,:),'Color','#A2142F','LineWidth',0.8)
    plot(time,med2{1,2}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,med3{1,2}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-3 3]); 
    grid on
    legend('Medium P1','Medium P2','Medium P3'); title(strcat('Low ',{' '},regression_type))
    hold off
    subplot(5,2,4)
    hold on
    plot(time,med1{2,2}(elec_id,:),'Color','#A2142F','LineWidth',0.8)
    plot(time,med2{2,2}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,med3{2,2}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-3 3]); 
    grid on
    legend('Medium P1','Medium P2','Medium P3'); title(strcat('High ',{' '},regression_type))
    hold off
    subplot(5,2,5)
    hold on
    plot(time,thin1{1,2}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med1{1,2}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick1{1,2}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); 
    grid on
    legend('Thin','Medium','Thick'); title(strcat('Partition 1,','Low ',{' '},regression_type));
    hold off
    subplot(5,2,6)
    hold on
    plot(time,thin1{2,2}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med1{2,2}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick1{2,2}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); 
    grid on
    legend('Thin','Medium','Thick'); title(strcat('Partition 1,','High ',{' '},regression_type));
    hold off
    subplot(5,2,7)
    hold on
    plot(time,thin2{1,2}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med2{1,2}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick2{1,2}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); 
    grid on
    legend('Thin','Medium','Thick'); title(strcat('Partition 2,','Low ',{' '},regression_type));
    hold off
    subplot(5,2,8)
    hold on
    plot(time,thin2{2,2}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med2{2,2}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick2{2,2}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); 
    grid on
    legend('Thin','Medium','Thick'); title(strcat('Partition 2,','High ',{' '},regression_type));
    hold off
    subplot(5,2,9)
    hold on
    plot(time,thin3{1,2}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med3{1,2}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick3{1,2}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); 
    grid on
    legend('Thin','Medium','Thick'); title(strcat('Partition 3,','Low ',{' '},regression_type));
    hold off
    subplot(5,2,10)
    hold on
    plot(time,thin3{2,2}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med3{2,2}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick3{2,2}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); 
    grid on
    legend('Thin','Medium','Thick'); title(strcat('Partition 3,','High ',{' '},regression_type));
    hold off

    %use same labels
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,'Power (in dB)'); 
    xlabel(han,'Time (in sec)');
    mtitle = strcat('Time series of power @',elec_label,' ,Partitions, Beta band');
    title(han,mtitle);

    fig = figure(3);
    elec_label = chosen_electrode{3}.electrode;
    elec_id = find(strcmp(data{1}.label,elec_label));
    subplot(5,2,1)
    hold on
    plot(time,data1{1,3}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,data2{1,3}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,data3{1,3}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); 
    grid on
    legend('PGI partition 1','PGI partition 2','PGI partition 3'); title(strcat('Low ',{' '},regression_type))
    hold off
    subplot(5,2,2)
    hold on
    plot(time,data1{2,3}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,data2{2,3}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,data3{2,3}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); 
    grid on
    legend('PGI partition 1','PGI partition 2','PGI partition 3'); title(strcat('High ',{' '},regression_type))
    hold off
    subplot(5,2,3)
    hold on
    plot(time,med1{1,3}(elec_id,:),'Color','#A2142F','LineWidth',0.8)
    plot(time,med2{1,3}(elec_id,:),'Color','r','LineWidth',0.8)
    plot(time,med3{1,3}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); 
    grid on
    legend('Medium P1','Medium P2','Medium P3'); title(strcat('Low ',{' '},regression_type))
    hold off
    subplot(5,2,4)
    hold on
    plot(time,med1{2,3}(elec_id,:),'Color','#A2142F','LineWidth',0.8)
    plot(time,med2{2,3}(elec_id,:),'Color','r','LineWidth',0.8)
    plot(time,med3{2,3}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); 
    grid on
    legend('Medium P1','Medium P2','Medium P3'); title(strcat('High ',{' '},regression_type))
    hold off
    subplot(5,2,5)
    hold on
    plot(time,thin1{1,3}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med1{1,3}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick1{1,3}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-1 3]); 
    grid on
    legend('Thin','Medium','Thick'); title(strcat('Partition 1,','Low ',{' '},regression_type));
    hold off
    subplot(5,2,6)
    hold on
    plot(time,thin1{2,3}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med1{2,3}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick1{2,3}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-1 3]); 
    grid on
    legend('Thin','Medium','Thick'); title(strcat('Partition 1,','High ',{' '},regression_type));
    hold off
    subplot(5,2,7)
    hold on
    plot(time,thin2{1,3}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med2{1,3}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick2{1,3}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-1 3]); 
    grid on
    legend('Thin','Medium','Thick'); title(strcat('Partition 2,','Low ',{' '},regression_type));
    hold off
    subplot(5,2,8)
    hold on
    plot(time,thin2{2,3}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med2{2,3}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick2{2,3}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-1 3]); 
    grid on
    legend('Thin','Medium','Thick'); title(strcat('Partition 2,','High ',{' '},regression_type));
    hold off
    subplot(5,2,9)
    hold on
    plot(time,thin3{1,3}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med3{1,3}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick3{1,3}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-1 3]); 
    grid on
    legend('Thin','Medium','Thick'); title(strcat('Partition 3,','Low ',{' '},regression_type));
    hold off
    subplot(5,2,10)
    hold on
    plot(time,thin3{2,3}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med3{2,3}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick3{2,3}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-1 3]); 
    grid on
    legend('Thin','Medium','Thick'); title(strcat('Partition 3,','High ',{' '},regression_type));
    hold off

    %use same labels
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,'Power (in dB)'); 
    xlabel(han,'Time (in sec)');
    mtitle = strcat('Time series of power @',elec_label,' ,Partitions, Low gamma band');
    title(han,mtitle);

    fig = figure(4);
    elec_label = chosen_electrode{4}.electrode;
    elec_id = find(strcmp(data{1}.label,elec_label));
    subplot(5,2,1)
    hold on
    plot(time,data1{1,4}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,data2{1,4}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,data3{1,4}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); 
    grid on
    legend('PGI partition 1','PGI partition 2','PGI partition 3'); title(strcat('Low ',{' '},regression_type))
    hold off
    subplot(5,2,2)
    hold on
    plot(time,data1{2,4}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,data2{2,4}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,data3{2,4}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); 
    grid on
    legend('PGI partition 1','PGI partition 2','PGI partition 3'); title(strcat('High ',{' '},regression_type))
    hold off
    subplot(5,2,3)
    hold on
    plot(time,med1{1,4}(elec_id,:),'Color','#A2142F','LineWidth',0.8)
    plot(time,med2{1,4}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,med3{1,4}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-3 3]); 
    grid on
    legend('Medium P1','Medium P2','Medium P3'); title(strcat('Low ',{' '},regression_type))
    hold off
    subplot(5,2,4)
    hold on
    plot(time,med1{2,4}(elec_id,:),'Color','#A2142F','LineWidth',0.8)
    plot(time,med2{2,4}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,med3{2,4}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-3 3]); 
    grid on
    legend('Medium P1','Medium P2','Medium P3'); title(strcat('High ',{' '},regression_type))
    hold off
    subplot(5,2,5)
    hold on
    plot(time,thin1{1,4}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med1{1,4}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick1{1,4}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); 
    grid on
    legend('Thin','Medium','Thick'); title(strcat('Partition 1,','Low ',{' '},regression_type));
    hold off
    subplot(5,2,6)
    hold on
    plot(time,thin1{2,4}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med1{2,4}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick1{2,4}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); 
    grid on
    legend('Thin','Medium','Thick'); title(strcat('Partition 1,','High ',{' '},regression_type));
    hold off
    subplot(5,2,7)
    hold on
    plot(time,thin2{1,4}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med2{1,4}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick2{1,4}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); 
    grid on
    legend('Thin','Medium','Thick'); title(strcat('Partition 2,','Low ',{' '},regression_type));
    hold off
    subplot(5,2,8)
    hold on
    plot(time,thin2{2,4}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med2{2,4}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick2{2,4}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); 
    grid on
    legend('Thin','Medium','Thick'); title(strcat('Partition 2,','High ',{' '},regression_type));
    hold off
    subplot(5,2,9)
    hold on
    plot(time,thin3{1,4}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med3{1,4}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick3{1,4}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); 
    grid on
    legend('Thin','Medium','Thick'); title(strcat('Partition 3,','Low ',{' '},regression_type));
    hold off
    subplot(5,2,10)
    hold on
    plot(time,thin3{2,4}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med3{2,4}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick3{2,4}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); 
    grid on
    legend('Thin','Medium','Thick'); title(strcat('Partition 3,','High ',{' '},regression_type));
    hold off

    %use same labels
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,'Power (in dB)'); 
    xlabel(han,'Time (in sec)');
    mtitle = strcat('Time series of power @',elec_label,' ,Partitions, Medium gamma band');
    title(han,mtitle);

    fig = figure(5);
    elec_label = chosen_electrode{5}.electrode;
    elec_id = find(strcmp(data{1}.label,elec_label));
    subplot(5,2,1)
    hold on
    plot(time,data1{1,5}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,data2{1,5}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,data3{1,5}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 4]); 
    grid on
    legend('PGI partition 1','PGI partition 2','PGI partition 3'); title(strcat('Low ',{' '},regression_type))
    hold off
    subplot(5,2,2)
    hold on
    plot(time,data1{2,5}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,data2{2,5}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,data3{2,5}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 4]); 
    grid on
    legend('PGI partition 1','PGI partition 2','PGI partition 3'); title(strcat('High ',{' '},regression_type))
    hold off
    subplot(5,2,3)
    hold on
    plot(time,med1{1,5}(elec_id,:),'Color','#A2142F','LineWidth',0.8)
    plot(time,med2{1,5}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,med3{1,5}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-3 3]); 
    grid on
    legend('Medium P1','Medium P2','Medium P3'); title(strcat('Low ',{' '},regression_type))
    hold off
    subplot(5,2,4)
    hold on
    plot(time,med1{2,5}(elec_id,:),'Color','#A2142F','LineWidth',0.8)
    plot(time,med2{2,5}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,med3{2,5}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-3 3]); 
    grid on
    legend('Medium P1','Medium P2','Medium P3'); title(strcat('High ',{' '},regression_type))
    hold off
    subplot(5,2,5)
    hold on
    plot(time,thin1{1,5}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med1{1,5}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick1{1,5}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-3 4]); 
    grid on
    legend('Thin','Medium','Thick'); title(strcat('Partition 1,','Low ',{' '},regression_type));
    hold off
    subplot(5,2,6)
    hold on
    plot(time,thin1{2,5}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med1{2,5}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick1{2,5}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-3 4]); 
    grid on
    legend('Thin','Medium','Thick'); title(strcat('Partition 1,','High ',{' '},regression_type));
    hold off
    subplot(5,2,7)
    hold on
    plot(time,thin2{1,5}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med2{1,5}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick2{1,5}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-3 4]); 
    grid on
    legend('Thin','Medium','Thick'); title(strcat('Partition 2,','Low ',{' '},regression_type));
    hold off
    subplot(5,2,8)
    hold on
    plot(time,thin2{2,5}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med2{2,5}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick2{2,5}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-3 4]); 
    grid on
    legend('Thin','Medium','Thick'); title(strcat('Partition 2,','High ',{' '},regression_type));
    hold off
    subplot(5,2,9)
    hold on
    plot(time,thin3{1,5}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med3{1,5}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick3{1,5}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-3 4]); 
    grid on
    legend('Thin','Medium','Thick'); title(strcat('Partition 3,','Low ',{' '},regression_type));
    hold off
    subplot(5,2,10)
    hold on
    plot(time,thin3{2,5}(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med3{2,5}(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick3{2,5}(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-3 4]); 
    grid on
    legend('Thin','Medium','Thick'); title(strcat('Partition 3,','High ',{' '},regression_type));
    hold off

    %use same labels
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,'Power (in dB)'); 
    xlabel(han,'Time (in sec)');
    mtitle = strcat('Time series of power @',elec_label,' ,Partitions, High gamma band');
    title(han,mtitle);
end

%% plotting time series for partitions
function time_series_partitions(data1,data2,data3,thin1,thin2,thin3,med1,med2,med3,thick1,thick2,thick3,chosen_electrode,end_latency)
    PGI_grand_avgs1 = calculate_freq_avg(data1,end_latency);
    PGI_grand_avgs2 = calculate_freq_avg(data2,end_latency);
    PGI_grand_avgs3 = calculate_freq_avg(data3,end_latency);

    thin_grand_avgs1 = calculate_avg_conditions(thin1,end_latency);
    thin_grand_avgs2 = calculate_avg_conditions(thin2,end_latency);
    thin_grand_avgs3 = calculate_avg_conditions(thin3,end_latency);

    med_grand_avgs1 = calculate_avg_conditions(med1,end_latency);
    med_grand_avgs2 = calculate_avg_conditions(med2,end_latency);
    med_grand_avgs3 = calculate_avg_conditions(med3,end_latency);

    thick_grand_avgs1 = calculate_avg_conditions(thick1,end_latency);
    thick_grand_avgs2 = calculate_avg_conditions(thick2,end_latency);
    thick_grand_avgs3 = calculate_avg_conditions(thick3,end_latency);

    PGI_alpha1 = squeeze(mean(PGI_grand_avgs1{1}.powspctrm,2,'omitnan'));
    PGI_beta1 = squeeze(mean(PGI_grand_avgs1{2}.powspctrm,2,'omitnan'));
    PGI_lgamma1 = squeeze(mean(PGI_grand_avgs1{3}.powspctrm,2,'omitnan'));
    PGI_mgamma1 = squeeze(mean(PGI_grand_avgs1{4}.powspctrm,2,'omitnan'));
    PGI_hgamma1 = squeeze(mean(PGI_grand_avgs1{5}.powspctrm,2,'omitnan'));
    PGI_alpha2 = squeeze(mean(PGI_grand_avgs2{1}.powspctrm,2,'omitnan'));
    PGI_beta2 = squeeze(mean(PGI_grand_avgs2{2}.powspctrm,2,'omitnan'));
    PGI_lgamma2 = squeeze(mean(PGI_grand_avgs2{3}.powspctrm,2,'omitnan'));
    PGI_mgamma2 = squeeze(mean(PGI_grand_avgs2{4}.powspctrm,2,'omitnan'));
    PGI_hgamma2 = squeeze(mean(PGI_grand_avgs2{5}.powspctrm,2,'omitnan'));
    PGI_alpha3 = squeeze(mean(PGI_grand_avgs3{1}.powspctrm,2,'omitnan'));
    PGI_beta3 = squeeze(mean(PGI_grand_avgs3{2}.powspctrm,2,'omitnan'));
    PGI_lgamma3 = squeeze(mean(PGI_grand_avgs3{3}.powspctrm,2,'omitnan'));
    PGI_mgamma3 = squeeze(mean(PGI_grand_avgs3{4}.powspctrm,2,'omitnan'));
    PGI_hgamma3 = squeeze(mean(PGI_grand_avgs3{5}.powspctrm,2,'omitnan'));

    thin_alpha1 = squeeze(mean(thin_grand_avgs1{1}.powspctrm,2,'omitnan'));
    thin_beta1 = squeeze(mean(thin_grand_avgs1{2}.powspctrm,2,'omitnan'));
    thin_lgamma1 = squeeze(mean(thin_grand_avgs1{3}.powspctrm,2,'omitnan'));
    thin_mgamma1 = squeeze(mean(thin_grand_avgs1{4}.powspctrm,2,'omitnan'));
    thin_hgamma1 = squeeze(mean(thin_grand_avgs1{5}.powspctrm,2,'omitnan'));
    thin_alpha2 = squeeze(mean(thin_grand_avgs2{1}.powspctrm,2,'omitnan'));
    thin_beta2 = squeeze(mean(thin_grand_avgs2{2}.powspctrm,2,'omitnan'));
    thin_lgamma2 = squeeze(mean(thin_grand_avgs2{3}.powspctrm,2,'omitnan'));
    thin_mgamma2 = squeeze(mean(thin_grand_avgs2{4}.powspctrm,2,'omitnan'));
    thin_hgamma2 = squeeze(mean(thin_grand_avgs2{5}.powspctrm,2,'omitnan'));
    thin_alpha3 = squeeze(mean(thin_grand_avgs3{1}.powspctrm,2,'omitnan'));
    thin_beta3 = squeeze(mean(thin_grand_avgs3{2}.powspctrm,2,'omitnan'));
    thin_lgamma3 = squeeze(mean(thin_grand_avgs3{3}.powspctrm,2,'omitnan'));
    thin_mgamma3 = squeeze(mean(thin_grand_avgs3{4}.powspctrm,2,'omitnan'));
    thin_hgamma3 = squeeze(mean(thin_grand_avgs3{5}.powspctrm,2,'omitnan'));

    med_alpha1 = squeeze(mean(med_grand_avgs1{1}.powspctrm,2,'omitnan'));
    med_beta1 = squeeze(mean(med_grand_avgs1{2}.powspctrm,2,'omitnan'));
    med_lgamma1 = squeeze(mean(med_grand_avgs1{3}.powspctrm,2,'omitnan'));
    med_mgamma1 = squeeze(mean(med_grand_avgs1{4}.powspctrm,2,'omitnan'));
    med_hgamma1 = squeeze(mean(med_grand_avgs1{5}.powspctrm,2,'omitnan'));
    med_alpha2 = squeeze(mean(med_grand_avgs2{1}.powspctrm,2,'omitnan'));
    med_beta2 = squeeze(mean(med_grand_avgs2{2}.powspctrm,2,'omitnan'));
    med_lgamma2 = squeeze(mean(med_grand_avgs2{3}.powspctrm,2,'omitnan'));
    med_mgamma2 = squeeze(mean(med_grand_avgs2{4}.powspctrm,2,'omitnan'));
    med_hgamma2 = squeeze(mean(med_grand_avgs2{5}.powspctrm,2,'omitnan'));
    med_alpha3 = squeeze(mean(med_grand_avgs3{1}.powspctrm,2,'omitnan'));
    med_beta3 = squeeze(mean(med_grand_avgs3{2}.powspctrm,2,'omitnan'));
    med_lgamma3 = squeeze(mean(med_grand_avgs3{3}.powspctrm,2,'omitnan'));
    med_mgamma3 = squeeze(mean(med_grand_avgs3{4}.powspctrm,2,'omitnan'));
    med_hgamma3 = squeeze(mean(med_grand_avgs3{5}.powspctrm,2,'omitnan'));

    thick_alpha1 = squeeze(mean(thick_grand_avgs1{1}.powspctrm,2,'omitnan'));
    thick_beta1 = squeeze(mean(thick_grand_avgs1{2}.powspctrm,2,'omitnan'));
    thick_lgamma1 = squeeze(mean(thick_grand_avgs1{3}.powspctrm,2,'omitnan'));
    thick_mgamma1 = squeeze(mean(thick_grand_avgs1{4}.powspctrm,2,'omitnan'));
    thick_hgamma1 = squeeze(mean(thick_grand_avgs1{5}.powspctrm,2,'omitnan'));
    thick_alpha2 = squeeze(mean(thick_grand_avgs2{1}.powspctrm,2,'omitnan'));
    thick_beta2 = squeeze(mean(thick_grand_avgs2{2}.powspctrm,2,'omitnan'));
    thick_lgamma2 = squeeze(mean(thick_grand_avgs2{3}.powspctrm,2,'omitnan'));
    thick_mgamma2 = squeeze(mean(thick_grand_avgs2{4}.powspctrm,2,'omitnan'));
    thick_hgamma2 = squeeze(mean(thick_grand_avgs2{5}.powspctrm,2,'omitnan'));
    thick_alpha3 = squeeze(mean(thick_grand_avgs3{1}.powspctrm,2,'omitnan'));
    thick_beta3 = squeeze(mean(thick_grand_avgs3{2}.powspctrm,2,'omitnan'));
    thick_lgamma3 = squeeze(mean(thick_grand_avgs3{3}.powspctrm,2,'omitnan'));
    thick_mgamma3 = squeeze(mean(thick_grand_avgs3{4}.powspctrm,2,'omitnan'));
    thick_hgamma3 = squeeze(mean(thick_grand_avgs3{5}.powspctrm,2,'omitnan'));

    %plot results for all frequency bands (different electrodes)
    fig = figure(1);
    time = thick_grand_avgs1{1}.time;
    elec_label = chosen_electrode{1}.electrode;
    elec_id = find(strcmp(data1{1}.label,elec_label));
    subplot(5,1,1)
    hold on
    plot(time,PGI_alpha1(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,PGI_alpha2(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,PGI_alpha3(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); grid on
    legend('PGI partition 1','PGI partition 2','PGI partition 3');
    hold off
    subplot(5,1,2)
    hold on
    plot(time,med_alpha1(elec_id,:),'Color','#A2142F','LineWidth',0.8)
    plot(time,med_alpha2(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,med_alpha3(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); 
    grid on
    legend('Medium P1','Medium P2','Medium P3'); 
    hold off
    subplot(5,1,3)
    hold on
    plot(time,thin_alpha1(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med_alpha1(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick_alpha1(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); grid on
    legend('Thin','Medium','Thick'); title('Partition 1');
    hold off
    subplot(5,1,4)
    hold on
    plot(time,thin_alpha2(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med_alpha2(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick_alpha2(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); grid on
    legend('Thin','Medium','Thick'); title('Partition 2');
    hold off
    subplot(5,1,5)
    hold on
    plot(time,thin_alpha3(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med_alpha3(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick_alpha3(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); grid on
    legend('Thin','Medium','Thick'); title('Partition 3');
    hold off

    %use same labels
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,'Power (in dB)'); 
    xlabel(han,'Time (in sec)');
    mtitle = strcat('Time series of power @',elec_label,' ,Partitions, Alpha band');
    title(han,mtitle);

    fig = figure(2);
    time = thick_grand_avgs1{1}.time;
    elec_label = chosen_electrode{2}.electrode;
    elec_id = find(strcmp(data1{1}.label,elec_label));
    subplot(5,1,1)
    hold on
    plot(time,PGI_beta1(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,PGI_beta2(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,PGI_beta3(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); grid on
    legend('PGI partition 1','PGI partition 2','PGI partition 3');
    hold off
    subplot(5,1,2)
    hold on
    plot(time,med_beta1(elec_id,:),'Color','#A2142F','LineWidth',0.8)
    plot(time,med_beta2(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,med_beta3(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); 
    grid on
    legend('Medium P1','Medium P2','Medium P3'); 
    hold off
    subplot(5,1,3)
    hold on
    plot(time,thin_beta1(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med_beta1(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick_beta1(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); grid on
    legend('Thin','Medium','Thick'); title('Partition 1');
    hold off
    subplot(5,1,4)
    hold on
    plot(time,thin_beta2(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med_beta2(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick_beta2(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); grid on
    legend('Thin','Medium','Thick'); title('Partition 2');
    hold off
    subplot(5,1,5)
    hold on
    plot(time,thin_beta3(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med_beta3(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick_beta3(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); grid on
    legend('Thin','Medium','Thick'); title('Partition 3');
    hold off

    %use same labels
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,'Power (in dB)'); 
    xlabel(han,'Time (in sec)');
    mtitle = strcat('Time series of power @',elec_label,' ,Partitions, Beta band');
    title(han,mtitle);

    fig = figure(3);
    time = thick_grand_avgs1{1}.time;
    elec_label = chosen_electrode{3}.electrode;
    elec_id = find(strcmp(data1{1}.label,elec_label));
    subplot(5,1,1)
    hold on
    plot(time,PGI_lgamma1(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,PGI_lgamma2(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,PGI_lgamma3(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); grid on
    legend('PGI partition 1','PGI partition 2','PGI partition 3');
    hold off
    subplot(5,1,2)
    hold on
    plot(time,med_lgamma1(elec_id,:),'Color','#A2142F','LineWidth',0.8)
    plot(time,med_lgamma2(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,med_lgamma3(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); 
    grid on
    legend('Medium P1','Medium P2','Medium P3'); 
    hold off
    subplot(5,1,3)
    hold on
    plot(time,thin_lgamma1(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med_lgamma1(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick_lgamma1(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); grid on
    legend('Thin','Medium','Thick'); title('Partition 1');
    hold off
    subplot(5,1,4)
    hold on
    plot(time,thin_lgamma2(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med_lgamma2(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick_lgamma2(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); grid on
    legend('Thin','Medium','Thick'); title('Partition 2');
    hold off
    subplot(5,1,5)
    hold on
    plot(time,thin_lgamma3(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med_lgamma3(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick_lgamma3(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); grid on
    legend('Thin','Medium','Thick'); title('Partition 3');
    hold off

    %use same labels
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,'Power (in dB)'); 
    xlabel(han,'Time (in sec)');
    mtitle = strcat('Time series of power @',elec_label,' ,Partitions, Low gamma band');
    title(han,mtitle);

    fig = figure(4);
    time = thick_grand_avgs1{1}.time;
    elec_label = chosen_electrode{4}.electrode;
    elec_id = find(strcmp(data1{1}.label,elec_label));
    subplot(5,1,1)
    hold on
    plot(time,PGI_mgamma1(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,PGI_mgamma2(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,PGI_mgamma3(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); grid on
    legend('PGI partition 1','PGI partition 2','PGI partition 3');
    hold off
    subplot(5,1,2)
    hold on
    plot(time,med_mgamma1(elec_id,:),'Color','#A2142F','LineWidth',0.8)
    plot(time,med_mgamma2(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,med_mgamma3(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); 
    grid on
    legend('Medium P1','Medium P2','Medium P3'); 
    hold off
    subplot(5,1,3)
    hold on
    plot(time,thin_mgamma1(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med_mgamma1(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick_mgamma1(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); grid on
    legend('Thin','Medium','Thick'); title('Partition 1');
    hold off
    subplot(5,1,4)
    hold on
    plot(time,thin_mgamma2(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med_mgamma2(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick_mgamma2(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); grid on
    legend('Thin','Medium','Thick'); title('Partition 2');
    hold off
    subplot(5,1,5)
    hold on
    plot(time,thin_mgamma3(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med_mgamma3(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick_mgamma3(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); grid on
    legend('Thin','Medium','Thick'); title('Partition 3');
    hold off

    %use same labels
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,'Power (in dB)'); 
    xlabel(han,'Time (in sec)');
    mtitle = strcat('Time series of power @',elec_label,' ,Partitions, Medium gamma band');
    title(han,mtitle);

    fig = figure(5);
    time = thick_grand_avgs1{1}.time;
    elec_label = chosen_electrode{5}.electrode;
    elec_id = find(strcmp(data1{1}.label,elec_label));
    subplot(5,1,1)
    hold on
    plot(time,PGI_hgamma1(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,PGI_hgamma2(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,PGI_hgamma3(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); grid on
    legend('PGI partition 1','PGI partition 2','PGI partition 3');
    hold off
    subplot(5,1,2)
    hold on
    plot(time,med_hgamma1(elec_id,:),'Color','#A2142F','LineWidth',0.8)
    plot(time,med_hgamma2(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,med_hgamma3(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); 
    grid on
    legend('Medium P1','Medium P2','Medium P3'); 
    hold off
    subplot(5,1,3)
    hold on
    plot(time,thin_hgamma1(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med_hgamma1(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick_hgamma1(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); grid on
    legend('Thin','Medium','Thick'); title('Partition 1');
    hold off
    subplot(5,1,4)
    hold on
    plot(time,thin_hgamma2(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med_hgamma2(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick_hgamma2(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); grid on
    legend('Thin','Medium','Thick'); title('Partition 2');
    hold off
    subplot(5,1,5)
    hold on
    plot(time,thin_hgamma3(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med_hgamma3(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick_hgamma3(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    xlim([-0.1 end_latency]); ylim([-2 2]); grid on
    legend('Thin','Medium','Thick'); title('Partition 3');
    hold off

    %use same labels
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,'Power (in dB)'); 
    xlabel(han,'Time (in sec)');
    mtitle = strcat('Time series of power @',elec_label,' ,Partitions, High gamma band');
    title(han,mtitle);
end

%% find median split for factors, partitions
function [PGI_series1,PGI_series2,PGI_series3,thin_series1,thin_series2,thin_series3,...
    med_series1,med_series2,med_series3,thick_series1,thick_series2,thick_series3] = median_split_partitions_factor(master_dir,...
    data1,data2,data3,participant_order,regression_type,end_latency)
    cd("C:\Users\marga\Desktop\Research Project\participants_Austyn")
    load thin_partition1
    load med_partition1
    load thick_partition1
    load thin_partition2
    load med_partition2
    load thick_partition2
    load thin_partition3
    load med_partition3
    load thick_partition3
    cd(master_dir);
    
    partition=1;
    [PGI_data_high1, PGI_data_low1] = get_median_split(data1, participant_order, regression_type, partition);
    PGI_grand_avgs_high1 = calculate_freq_avg(PGI_data_high1,end_latency);
    PGI_grand_avgs_low1 = calculate_freq_avg(PGI_data_low1,end_latency);
    [thin_data_high1, thin_data_low1] = get_median_split(thin1,participant_order,regression_type,partition);
    [med_data_high1, med_data_low1] = get_median_split(med1,participant_order,regression_type,partition);
    [thick_data_high1, thick_data_low1] = get_median_split(thick1,participant_order,regression_type,partition);
    thin_grand_avgs_high1 = calculate_avg_conditions(thin_data_high1,end_latency);
    thin_grand_avgs_low1 = calculate_avg_conditions(thin_data_low1,end_latency);
    med_grand_avgs_high1 = calculate_avg_conditions(med_data_high1,end_latency);
    med_grand_avgs_low1 = calculate_avg_conditions(med_data_low1,end_latency);
    thick_grand_avgs_high1 = calculate_avg_conditions(thick_data_high1,end_latency);
    thick_grand_avgs_low1 = calculate_avg_conditions(thick_data_low1,end_latency);
    
    partition=2;
    [PGI_data_high2, PGI_data_low2] = get_median_split(data2, participant_order, regression_type, partition);
    PGI_grand_avgs_high2 = calculate_freq_avg(PGI_data_high2,end_latency);
    PGI_grand_avgs_low2 = calculate_freq_avg(PGI_data_low2,end_latency);
    [thin_data_high2, thin_data_low2] = get_median_split(thin2,participant_order,regression_type,partition);
    [med_data_high2, med_data_low2] = get_median_split(med2,participant_order,regression_type,partition);
    [thick_data_high2, thick_data_low2] = get_median_split(thick2,participant_order,regression_type,partition);
    thin_grand_avgs_high2 = calculate_avg_conditions(thin_data_high2,end_latency);
    thin_grand_avgs_low2 = calculate_avg_conditions(thin_data_low2,end_latency);
    med_grand_avgs_high2 = calculate_avg_conditions(med_data_high2,end_latency);
    med_grand_avgs_low2 = calculate_avg_conditions(med_data_low2,end_latency);
    thick_grand_avgs_high2 = calculate_avg_conditions(thick_data_high2,end_latency);
    thick_grand_avgs_low2 = calculate_avg_conditions(thick_data_low2,end_latency);
    
    partition=3;
    [PGI_data_high3, PGI_data_low3] = get_median_split(data3, participant_order, regression_type, partition);
    PGI_grand_avgs_high3 = calculate_freq_avg(PGI_data_high3,end_latency);
    PGI_grand_avgs_low3 = calculate_freq_avg(PGI_data_low3,end_latency);
    [thin_data_high3, thin_data_low3] = get_median_split(thin3,participant_order,regression_type,partition);
    [med_data_high3, med_data_low3] = get_median_split(med3,participant_order,regression_type,partition);
    [thick_data_high3, thick_data_low3] = get_median_split(thick3,participant_order,regression_type,partition);
    thin_grand_avgs_high3 = calculate_avg_conditions(thin_data_high3,end_latency);
    thin_grand_avgs_low3 = calculate_avg_conditions(thin_data_low3,end_latency);
    med_grand_avgs_high3 = calculate_avg_conditions(med_data_high3,end_latency);
    med_grand_avgs_low3 = calculate_avg_conditions(med_data_low3,end_latency);
    thick_grand_avgs_high3 = calculate_avg_conditions(thick_data_high3,end_latency);
    thick_grand_avgs_low3 = calculate_avg_conditions(thick_data_low3,end_latency);
    
    %partition 1
    PGI_alpha_high1 = squeeze(mean(PGI_grand_avgs_high1{1}.powspctrm,2,'omitnan'));
    PGI_alpha_low1 = squeeze(mean(PGI_grand_avgs_low1{1}.powspctrm,2,'omitnan'));
    PGI_beta_high1 = squeeze(mean(PGI_grand_avgs_high1{2}.powspctrm,2,'omitnan'));
    PGI_beta_low1 = squeeze(mean(PGI_grand_avgs_low1{2}.powspctrm,2,'omitnan'));
    PGI_lgamma_high1 = squeeze(mean(PGI_grand_avgs_high1{3}.powspctrm,2,'omitnan'));
    PGI_lgamma_low1 = squeeze(mean(PGI_grand_avgs_low1{3}.powspctrm,2,'omitnan'));
    PGI_mgamma_high1 = squeeze(mean(PGI_grand_avgs_high1{4}.powspctrm,2,'omitnan'));
    PGI_mgamma_low1 = squeeze(mean(PGI_grand_avgs_low1{4}.powspctrm,2,'omitnan'));
    PGI_hgamma_high1 = squeeze(mean(PGI_grand_avgs_high1{5}.powspctrm,2,'omitnan'));
    PGI_hgamma_low1 = squeeze(mean(PGI_grand_avgs_low1{5}.powspctrm,2,'omitnan'));

    thin_alpha_high1 = squeeze(mean(thin_grand_avgs_high1{1}.powspctrm,2,'omitnan'));
    thin_alpha_low1 = squeeze(mean(thin_grand_avgs_low1{1}.powspctrm,2,'omitnan'));
    thin_beta_high1 = squeeze(mean(thin_grand_avgs_high1{2}.powspctrm,2,'omitnan'));
    thin_beta_low1 = squeeze(mean(thin_grand_avgs_low1{2}.powspctrm,2,'omitnan'));
    thin_lgamma_high1 = squeeze(mean(thin_grand_avgs_high1{3}.powspctrm,2,'omitnan'));
    thin_lgamma_low1 = squeeze(mean(thin_grand_avgs_low1{3}.powspctrm,2,'omitnan'));
    thin_mgamma_high1 = squeeze(mean(thin_grand_avgs_high1{4}.powspctrm,2,'omitnan'));
    thin_mgamma_low1 = squeeze(mean(thin_grand_avgs_low1{4}.powspctrm,2,'omitnan'));
    thin_hgamma_high1 = squeeze(mean(thin_grand_avgs_high1{5}.powspctrm,2,'omitnan'));
    thin_hgamma_low1 = squeeze(mean(thin_grand_avgs_low1{5}.powspctrm,2,'omitnan'));

    med_alpha_high1 = squeeze(mean(med_grand_avgs_high1{1}.powspctrm,2,'omitnan'));
    med_alpha_low1 = squeeze(mean(med_grand_avgs_low1{1}.powspctrm,2,'omitnan'));
    med_beta_high1 = squeeze(mean(med_grand_avgs_high1{2}.powspctrm,2,'omitnan'));
    med_beta_low1 = squeeze(mean(med_grand_avgs_low1{2}.powspctrm,2,'omitnan'));
    med_lgamma_high1 = squeeze(mean(med_grand_avgs_high1{3}.powspctrm,2,'omitnan'));
    med_lgamma_low1 = squeeze(mean(med_grand_avgs_low1{3}.powspctrm,2,'omitnan'));
    med_mgamma_high1 = squeeze(mean(med_grand_avgs_high1{4}.powspctrm,2,'omitnan'));
    med_mgamma_low1 = squeeze(mean(med_grand_avgs_low1{4}.powspctrm,2,'omitnan'));
    med_hgamma_high1 = squeeze(mean(med_grand_avgs_high1{5}.powspctrm,2,'omitnan'));
    med_hgamma_low1 = squeeze(mean(med_grand_avgs_low1{5}.powspctrm,2,'omitnan'));

    thick_alpha_high1 = squeeze(mean(thick_grand_avgs_high1{1}.powspctrm,2,'omitnan'));
    thick_alpha_low1 = squeeze(mean(thick_grand_avgs_low1{1}.powspctrm,2,'omitnan'));
    thick_beta_high1 = squeeze(mean(thick_grand_avgs_high1{2}.powspctrm,2,'omitnan'));
    thick_beta_low1 = squeeze(mean(thick_grand_avgs_low1{2}.powspctrm,2,'omitnan'));
    thick_lgamma_high1 = squeeze(mean(thick_grand_avgs_high1{3}.powspctrm,2,'omitnan'));
    thick_lgamma_low1 = squeeze(mean(thick_grand_avgs_low1{3}.powspctrm,2,'omitnan'));
    thick_mgamma_high1 = squeeze(mean(thick_grand_avgs_high1{4}.powspctrm,2,'omitnan'));
    thick_mgamma_low1 = squeeze(mean(thick_grand_avgs_low1{4}.powspctrm,2,'omitnan'));
    thick_hgamma_high1 = squeeze(mean(thick_grand_avgs_high1{5}.powspctrm,2,'omitnan'));
    thick_hgamma_low1 = squeeze(mean(thick_grand_avgs_low1{5}.powspctrm,2,'omitnan'));

    %partition 2
    PGI_alpha_high2 = squeeze(mean(PGI_grand_avgs_high2{1}.powspctrm,2,'omitnan'));
    PGI_alpha_low2 = squeeze(mean(PGI_grand_avgs_low2{1}.powspctrm,2,'omitnan'));
    PGI_beta_high2 = squeeze(mean(PGI_grand_avgs_high2{2}.powspctrm,2,'omitnan'));
    PGI_beta_low2 = squeeze(mean(PGI_grand_avgs_low2{2}.powspctrm,2,'omitnan'));
    PGI_lgamma_high2 = squeeze(mean(PGI_grand_avgs_high2{3}.powspctrm,2,'omitnan'));
    PGI_lgamma_low2 = squeeze(mean(PGI_grand_avgs_low2{3}.powspctrm,2,'omitnan'));
    PGI_mgamma_high2 = squeeze(mean(PGI_grand_avgs_high2{4}.powspctrm,2,'omitnan'));
    PGI_mgamma_low2 = squeeze(mean(PGI_grand_avgs_low2{4}.powspctrm,2,'omitnan'));
    PGI_hgamma_high2 = squeeze(mean(PGI_grand_avgs_high2{5}.powspctrm,2,'omitnan'));
    PGI_hgamma_low2 = squeeze(mean(PGI_grand_avgs_low2{5}.powspctrm,2,'omitnan'));

    thin_alpha_high2 = squeeze(mean(thin_grand_avgs_high2{1}.powspctrm,2,'omitnan'));
    thin_alpha_low2 = squeeze(mean(thin_grand_avgs_low2{1}.powspctrm,2,'omitnan'));
    thin_beta_high2 = squeeze(mean(thin_grand_avgs_high2{2}.powspctrm,2,'omitnan'));
    thin_beta_low2 = squeeze(mean(thin_grand_avgs_low2{2}.powspctrm,2,'omitnan'));
    thin_lgamma_high2 = squeeze(mean(thin_grand_avgs_high2{3}.powspctrm,2,'omitnan'));
    thin_lgamma_low2 = squeeze(mean(thin_grand_avgs_low2{3}.powspctrm,2,'omitnan'));
    thin_mgamma_high2 = squeeze(mean(thin_grand_avgs_high2{4}.powspctrm,2,'omitnan'));
    thin_mgamma_low2 = squeeze(mean(thin_grand_avgs_low2{4}.powspctrm,2,'omitnan'));
    thin_hgamma_high2 = squeeze(mean(thin_grand_avgs_high2{5}.powspctrm,2,'omitnan'));
    thin_hgamma_low2 = squeeze(mean(thin_grand_avgs_low2{5}.powspctrm,2,'omitnan'));

    med_alpha_high2 = squeeze(mean(med_grand_avgs_high2{1}.powspctrm,2,'omitnan'));
    med_alpha_low2 = squeeze(mean(med_grand_avgs_low2{1}.powspctrm,2,'omitnan'));
    med_beta_high2 = squeeze(mean(med_grand_avgs_high2{2}.powspctrm,2,'omitnan'));
    med_beta_low2 = squeeze(mean(med_grand_avgs_low2{2}.powspctrm,2,'omitnan'));
    med_lgamma_high2 = squeeze(mean(med_grand_avgs_high2{3}.powspctrm,2,'omitnan'));
    med_lgamma_low2 = squeeze(mean(med_grand_avgs_low2{3}.powspctrm,2,'omitnan'));
    med_mgamma_high2 = squeeze(mean(med_grand_avgs_high2{4}.powspctrm,2,'omitnan'));
    med_mgamma_low2 = squeeze(mean(med_grand_avgs_low2{4}.powspctrm,2,'omitnan'));
    med_hgamma_high2 = squeeze(mean(med_grand_avgs_high2{5}.powspctrm,2,'omitnan'));
    med_hgamma_low2 = squeeze(mean(med_grand_avgs_low2{5}.powspctrm,2,'omitnan'));

    thick_alpha_high2 = squeeze(mean(thick_grand_avgs_high2{1}.powspctrm,2,'omitnan'));
    thick_alpha_low2 = squeeze(mean(thick_grand_avgs_low2{1}.powspctrm,2,'omitnan'));
    thick_beta_high2 = squeeze(mean(thick_grand_avgs_high2{2}.powspctrm,2,'omitnan'));
    thick_beta_low2 = squeeze(mean(thick_grand_avgs_low2{2}.powspctrm,2,'omitnan'));
    thick_lgamma_high2 = squeeze(mean(thick_grand_avgs_high2{3}.powspctrm,2,'omitnan'));
    thick_lgamma_low2 = squeeze(mean(thick_grand_avgs_low2{3}.powspctrm,2,'omitnan'));
    thick_mgamma_high2 = squeeze(mean(thick_grand_avgs_high2{4}.powspctrm,2,'omitnan'));
    thick_mgamma_low2 = squeeze(mean(thick_grand_avgs_low2{4}.powspctrm,2,'omitnan'));
    thick_hgamma_high2 = squeeze(mean(thick_grand_avgs_high2{5}.powspctrm,2,'omitnan'));
    thick_hgamma_low2 = squeeze(mean(thick_grand_avgs_low2{5}.powspctrm,2,'omitnan'));
    
    %partition 3
    PGI_alpha_high3 = squeeze(mean(PGI_grand_avgs_high3{1}.powspctrm,2,'omitnan'));
    PGI_alpha_low3 = squeeze(mean(PGI_grand_avgs_low3{1}.powspctrm,2,'omitnan'));
    PGI_beta_high3 = squeeze(mean(PGI_grand_avgs_high3{2}.powspctrm,2,'omitnan'));
    PGI_beta_low3 = squeeze(mean(PGI_grand_avgs_low3{2}.powspctrm,2,'omitnan'));
    PGI_lgamma_high3 = squeeze(mean(PGI_grand_avgs_high3{3}.powspctrm,2,'omitnan'));
    PGI_lgamma_low3 = squeeze(mean(PGI_grand_avgs_low3{3}.powspctrm,2,'omitnan'));
    PGI_mgamma_high3 = squeeze(mean(PGI_grand_avgs_high3{4}.powspctrm,2,'omitnan'));
    PGI_mgamma_low3 = squeeze(mean(PGI_grand_avgs_low3{4}.powspctrm,2,'omitnan'));
    PGI_hgamma_high3 = squeeze(mean(PGI_grand_avgs_high3{5}.powspctrm,2,'omitnan'));
    PGI_hgamma_low3 = squeeze(mean(PGI_grand_avgs_low3{5}.powspctrm,2,'omitnan'));

    thin_alpha_high3 = squeeze(mean(thin_grand_avgs_high3{1}.powspctrm,2,'omitnan'));
    thin_alpha_low3 = squeeze(mean(thin_grand_avgs_low3{1}.powspctrm,2,'omitnan'));
    thin_beta_high3 = squeeze(mean(thin_grand_avgs_high3{2}.powspctrm,2,'omitnan'));
    thin_beta_low3 = squeeze(mean(thin_grand_avgs_low3{2}.powspctrm,2,'omitnan'));
    thin_lgamma_high3 = squeeze(mean(thin_grand_avgs_high3{3}.powspctrm,2,'omitnan'));
    thin_lgamma_low3 = squeeze(mean(thin_grand_avgs_low3{3}.powspctrm,2,'omitnan'));
    thin_mgamma_high3 = squeeze(mean(thin_grand_avgs_high3{4}.powspctrm,2,'omitnan'));
    thin_mgamma_low3 = squeeze(mean(thin_grand_avgs_low3{4}.powspctrm,2,'omitnan'));
    thin_hgamma_high3 = squeeze(mean(thin_grand_avgs_high3{5}.powspctrm,2,'omitnan'));
    thin_hgamma_low3 = squeeze(mean(thin_grand_avgs_low3{5}.powspctrm,2,'omitnan'));

    med_alpha_high3 = squeeze(mean(med_grand_avgs_high3{1}.powspctrm,2,'omitnan'));
    med_alpha_low3 = squeeze(mean(med_grand_avgs_low3{1}.powspctrm,2,'omitnan'));
    med_beta_high3 = squeeze(mean(med_grand_avgs_high3{2}.powspctrm,2,'omitnan'));
    med_beta_low3 = squeeze(mean(med_grand_avgs_low3{2}.powspctrm,2,'omitnan'));
    med_lgamma_high3 = squeeze(mean(med_grand_avgs_high3{3}.powspctrm,2,'omitnan'));
    med_lgamma_low3 = squeeze(mean(med_grand_avgs_low3{3}.powspctrm,2,'omitnan'));
    med_mgamma_high3 = squeeze(mean(med_grand_avgs_high3{4}.powspctrm,2,'omitnan'));
    med_mgamma_low3 = squeeze(mean(med_grand_avgs_low3{4}.powspctrm,2,'omitnan'));
    med_hgamma_high3 = squeeze(mean(med_grand_avgs_high3{5}.powspctrm,2,'omitnan'));
    med_hgamma_low3 = squeeze(mean(med_grand_avgs_low3{5}.powspctrm,2,'omitnan'));

    thick_alpha_high3 = squeeze(mean(thick_grand_avgs_high3{1}.powspctrm,2,'omitnan'));
    thick_alpha_low3 = squeeze(mean(thick_grand_avgs_low3{1}.powspctrm,2,'omitnan'));
    thick_beta_high3 = squeeze(mean(thick_grand_avgs_high3{2}.powspctrm,2,'omitnan'));
    thick_beta_low3 = squeeze(mean(thick_grand_avgs_low3{2}.powspctrm,2,'omitnan'));
    thick_lgamma_high3 = squeeze(mean(thick_grand_avgs_high3{3}.powspctrm,2,'omitnan'));
    thick_lgamma_low3 = squeeze(mean(thick_grand_avgs_low3{3}.powspctrm,2,'omitnan'));
    thick_mgamma_high3 = squeeze(mean(thick_grand_avgs_high3{4}.powspctrm,2,'omitnan'));
    thick_mgamma_low3 = squeeze(mean(thick_grand_avgs_low3{4}.powspctrm,2,'omitnan'));
    thick_hgamma_high3 = squeeze(mean(thick_grand_avgs_high3{5}.powspctrm,2,'omitnan'));
    thick_hgamma_low3 = squeeze(mean(thick_grand_avgs_low3{5}.powspctrm,2,'omitnan'));
    
    PGI_series1 = [{PGI_alpha_low1},{PGI_beta_low1},{PGI_lgamma_low1},{PGI_mgamma_low1},{PGI_hgamma_low1};
        {PGI_alpha_high1},{PGI_beta_high1},{PGI_lgamma_high1},{PGI_mgamma_high1},{PGI_hgamma_high1}];
    PGI_series2 = [{PGI_alpha_low2},{PGI_beta_low2},{PGI_lgamma_low2},{PGI_mgamma_low2},{PGI_hgamma_low2};
        {PGI_alpha_high2},{PGI_beta_high2},{PGI_lgamma_high2},{PGI_mgamma_high2},{PGI_hgamma_high2}];
    PGI_series3 = [{PGI_alpha_low3},{PGI_beta_low3},{PGI_lgamma_low3},{PGI_mgamma_low3},{PGI_hgamma_low3};
        {PGI_alpha_high3},{PGI_beta_high3},{PGI_lgamma_high3},{PGI_mgamma_high3},{PGI_hgamma_high3}];
    
    thin_series1 = [{thin_alpha_low1},{thin_beta_low1},{thin_lgamma_low1},{thin_mgamma_low1},{thin_hgamma_low1};
        {thin_alpha_high1},{thin_beta_high1},{thin_lgamma_high1},{thin_mgamma_high1},{thin_hgamma_high1}];
    thin_series2 = [{thin_alpha_low2},{thin_beta_low2},{thin_lgamma_low2},{thin_mgamma_low2},{thin_hgamma_low2};
        {thin_alpha_high2},{thin_beta_high2},{thin_lgamma_high2},{thin_mgamma_high2},{thin_hgamma_high2}];
    thin_series3 = [{thin_alpha_low3},{thin_beta_low3},{thin_lgamma_low3},{thin_mgamma_low3},{thin_hgamma_low3};
        {thin_alpha_high3},{thin_beta_high3},{thin_lgamma_high3},{thin_mgamma_high3},{thin_hgamma_high3}];
    
    med_series1 = [{med_alpha_low1},{med_beta_low1},{med_lgamma_low1},{med_mgamma_low1},{med_hgamma_low1};
        {med_alpha_high1},{med_beta_high1},{med_lgamma_high1},{med_mgamma_high1},{med_hgamma_high1}];
    med_series2 = [{med_alpha_low2},{med_beta_low2},{med_lgamma_low2},{med_mgamma_low2},{med_hgamma_low2};
        {med_alpha_high2},{med_beta_high2},{med_lgamma_high2},{med_mgamma_high2},{med_hgamma_high2}];
    med_series3 = [{med_alpha_low3},{med_beta_low3},{med_lgamma_low3},{med_mgamma_low3},{med_hgamma_low3};
        {med_alpha_high3},{med_beta_high3},{med_lgamma_high3},{med_mgamma_high3},{med_hgamma_high3}];
    
    thick_series1 = [{thick_alpha_low1},{thick_beta_low1},{thick_lgamma_low1},{thick_mgamma_low1},{thick_hgamma_low1};
        {thick_alpha_high1},{thick_beta_high1},{thick_lgamma_high1},{thick_mgamma_high1},{thick_hgamma_high1}];
    thick_series2 = [{thick_alpha_low2},{thick_beta_low2},{thick_lgamma_low2},{thick_mgamma_low2},{thick_hgamma_low2};
        {thick_alpha_high2},{thick_beta_high2},{thick_lgamma_high2},{thick_mgamma_high2},{thick_hgamma_high2}];
    thick_series3 = [{thick_alpha_low3},{thick_beta_low3},{thick_lgamma_low3},{thick_mgamma_low3},{thick_hgamma_low3};
        {thick_alpha_high3},{thick_beta_high3},{thick_lgamma_high3},{thick_mgamma_high3},{thick_hgamma_high3}];
end

%% plot median split for factor stats, onsets only
function median_split_onsets_factor(master_dir, data_PGI, participant_order,chosen_electrode,regression_type,end_latency)
    cd("C:\Users\marga\Desktop\Research Project\participants_Austyn")
    load data_conditions
    cd(master_dir);
    
    [PGI_data_high, PGI_data_low] = get_median_split(data_PGI, participant_order, regression_type, 1);
    PGI_grand_avgs_high = calculate_freq_avg(PGI_data_high,end_latency);
    PGI_grand_avgs_low = calculate_freq_avg(PGI_data_low,end_latency);
    for i = 1:32
        thin_data{i} = data_conditions{i}(1);
        med_data{i} = data_conditions{i}(2);
        thick_data{i} = data_conditions{i}(3);
    end
    [thin_data_high, thin_data_low] = get_median_split(thin_data,participant_order,regression_type,1);
    [med_data_high, med_data_low] = get_median_split(med_data,participant_order,regression_type,1);
    [thick_data_high, thick_data_low] = get_median_split(thick_data,participant_order,regression_type,1);
    thin_grand_avgs_high = calculate_avg_conditions(thin_data_high,end_latency);
    thin_grand_avgs_low = calculate_avg_conditions(thin_data_low,end_latency);
    med_grand_avgs_high = calculate_avg_conditions(med_data_high,end_latency);
    med_grand_avgs_low = calculate_avg_conditions(med_data_low,end_latency);
    thick_grand_avgs_high = calculate_avg_conditions(thick_data_high,end_latency);
    thick_grand_avgs_low = calculate_avg_conditions(thick_data_low,end_latency);

    %% plot time series of power of median split (PGI & conditions)
    PGI_alpha_high = squeeze(mean(PGI_grand_avgs_high{1}.powspctrm,2,'omitnan'));
    PGI_alpha_low = squeeze(mean(PGI_grand_avgs_low{1}.powspctrm,2,'omitnan'));
    PGI_beta_high = squeeze(mean(PGI_grand_avgs_high{2}.powspctrm,2,'omitnan'));
    PGI_beta_low = squeeze(mean(PGI_grand_avgs_low{2}.powspctrm,2,'omitnan'));
    PGI_lgamma_high = squeeze(mean(PGI_grand_avgs_high{3}.powspctrm,2,'omitnan'));
    PGI_lgamma_low = squeeze(mean(PGI_grand_avgs_low{3}.powspctrm,2,'omitnan'));
    PGI_mgamma_high = squeeze(mean(PGI_grand_avgs_high{4}.powspctrm,2,'omitnan'));
    PGI_mgamma_low = squeeze(mean(PGI_grand_avgs_low{4}.powspctrm,2,'omitnan'));
    PGI_hgamma_high = squeeze(mean(PGI_grand_avgs_high{5}.powspctrm,2,'omitnan'));
    PGI_hgamma_low = squeeze(mean(PGI_grand_avgs_low{5}.powspctrm,2,'omitnan'));

    thin_alpha_high = squeeze(mean(thin_grand_avgs_high{1}.powspctrm,2,'omitnan'));
    thin_alpha_low = squeeze(mean(thin_grand_avgs_low{1}.powspctrm,2,'omitnan'));
    thin_beta_high = squeeze(mean(thin_grand_avgs_high{2}.powspctrm,2,'omitnan'));
    thin_beta_low = squeeze(mean(thin_grand_avgs_low{2}.powspctrm,2,'omitnan'));
    thin_lgamma_high = squeeze(mean(thin_grand_avgs_high{3}.powspctrm,2,'omitnan'));
    thin_lgamma_low = squeeze(mean(thin_grand_avgs_low{3}.powspctrm,2,'omitnan'));
    thin_mgamma_high = squeeze(mean(thin_grand_avgs_high{4}.powspctrm,2,'omitnan'));
    thin_mgamma_low = squeeze(mean(thin_grand_avgs_low{4}.powspctrm,2,'omitnan'));
    thin_hgamma_high = squeeze(mean(thin_grand_avgs_high{5}.powspctrm,2,'omitnan'));
    thin_hgamma_low = squeeze(mean(thin_grand_avgs_low{5}.powspctrm,2,'omitnan'));

    med_alpha_high = squeeze(mean(med_grand_avgs_high{1}.powspctrm,2,'omitnan'));
    med_alpha_low = squeeze(mean(med_grand_avgs_low{1}.powspctrm,2,'omitnan'));
    med_beta_high = squeeze(mean(med_grand_avgs_high{2}.powspctrm,2,'omitnan'));
    med_beta_low = squeeze(mean(med_grand_avgs_low{2}.powspctrm,2,'omitnan'));
    med_lgamma_high = squeeze(mean(med_grand_avgs_high{3}.powspctrm,2,'omitnan'));
    med_lgamma_low = squeeze(mean(med_grand_avgs_low{3}.powspctrm,2,'omitnan'));
    med_mgamma_high = squeeze(mean(med_grand_avgs_high{4}.powspctrm,2,'omitnan'));
    med_mgamma_low = squeeze(mean(med_grand_avgs_low{4}.powspctrm,2,'omitnan'));
    med_hgamma_high = squeeze(mean(med_grand_avgs_high{5}.powspctrm,2,'omitnan'));
    med_hgamma_low = squeeze(mean(med_grand_avgs_low{5}.powspctrm,2,'omitnan'));

    thick_alpha_high = squeeze(mean(thick_grand_avgs_high{1}.powspctrm,2,'omitnan'));
    thick_alpha_low = squeeze(mean(thick_grand_avgs_low{1}.powspctrm,2,'omitnan'));
    thick_beta_high = squeeze(mean(thick_grand_avgs_high{2}.powspctrm,2,'omitnan'));
    thick_beta_low = squeeze(mean(thick_grand_avgs_low{2}.powspctrm,2,'omitnan'));
    thick_lgamma_high = squeeze(mean(thick_grand_avgs_high{3}.powspctrm,2,'omitnan'));
    thick_lgamma_low = squeeze(mean(thick_grand_avgs_low{3}.powspctrm,2,'omitnan'));
    thick_mgamma_high = squeeze(mean(thick_grand_avgs_high{4}.powspctrm,2,'omitnan'));
    thick_mgamma_low = squeeze(mean(thick_grand_avgs_low{4}.powspctrm,2,'omitnan'));
    thick_hgamma_high = squeeze(mean(thick_grand_avgs_high{5}.powspctrm,2,'omitnan'));
    thick_hgamma_low = squeeze(mean(thick_grand_avgs_low{5}.powspctrm,2,'omitnan'));

    %plot results for all frequency bands (chosen electrode)
    time = thin_grand_avgs_high{1}.time;
    figure(1)
    subplot(2,1,1)
    elec_label = chosen_electrode{1}.electrode;
    elec_id = find(strcmp(data_PGI{1}.label,elec_label));
    mtitle = strcat('High ',{' '},regression_type,': Alpha band, Electrode ',{' '});
    mtitle = strcat(mtitle,elec_label);
    title(mtitle);
    hold on
    plot(time,thin_alpha_high(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med_alpha_high(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick_alpha_high(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    ylim([-2.5 1.5]); xlim([-0.1 end_latency]);
    grid on
    yyaxis right
    plot(time,PGI_alpha_high(elec_id,:),'Color','#7E2F8E','LineWidth',1)
    ylim([-6 1]); xlim([-0.1 end_latency]);
    legend('Thin','Medium','Thick','PGI');
    hold off
    subplot(2,1,2)
    mtitle = strcat('Low ',{' '},regression_type,': Alpha band, Electrode ',{' '});
    mtitle = strcat(mtitle,elec_label);
    title(mtitle);
    hold on
    plot(time,thin_alpha_low(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med_alpha_low(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick_alpha_low(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    ylim([-2.5 1.5]); xlim([-0.1 end_latency]);
    grid on
    yyaxis right
    plot(time,PGI_alpha_low(elec_id,:),'Color','#7E2F8E','LineWidth',1)
    ylim([-6 1]); xlim([-0.1 end_latency]);
    legend('Thin','Medium','Thick','PGI');
    hold off

    figure(2)
    subplot(2,1,1)
    elec_label = chosen_electrode{2}.electrode;
    elec_id = find(strcmp(data_PGI{1}.label,elec_label));
    mtitle = strcat('High ',{' '},regression_type,': Beta band, Electrode ',{' '});
    mtitle = strcat(mtitle,elec_label);
    title(mtitle);
    hold on
    plot(time,thin_beta_high(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med_beta_high(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick_beta_high(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    ylim([-1 2.5]); xlim([-0.1 end_latency]);
    grid on
    yyaxis right
    plot(time,PGI_beta_high(elec_id,:),'Color','#7E2F8E','LineWidth',1)
    ylim([-3.5 2]); xlim([-0.1 end_latency]);
    legend('Thin','Medium','Thick','PGI');
    hold off
    subplot(2,1,2)
    mtitle = strcat('Low ',{' '},regression_type,': Beta band, Electrode ',{' '});
    mtitle = strcat(mtitle,elec_label);
    title(mtitle);
    hold on
    plot(time,thin_beta_low(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med_beta_low(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick_beta_low(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    ylim([-1 2.5]); xlim([-0.1 end_latency]);
    grid on
    yyaxis right
    plot(time,PGI_beta_low(elec_id,:),'Color','#7E2F8E','LineWidth',1)
    ylim([-3.5 2]); xlim([-0.1 end_latency]);
    legend('Thin','Medium','Thick','PGI');
    hold off

    figure(3)
    subplot(2,1,1)
    elec_label = chosen_electrode{3}.electrode;
    elec_id = find(strcmp(data_PGI{1}.label,elec_label));
    mtitle = strcat('High ',{' '},regression_type,': Low gamma band, Electrode ',{' '});
    mtitle = strcat(mtitle,elec_label);
    title(mtitle);
    hold on
    plot(time,thin_lgamma_high(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med_lgamma_high(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick_lgamma_high(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    ylim([-1 2.5]); xlim([-0.1 end_latency]);
    grid on
    yyaxis right
    plot(time,PGI_lgamma_high(elec_id,:),'Color','#7E2F8E','LineWidth',1)
    ylim([-3.5 2]); xlim([-0.1 end_latency]);
    legend('Thin','Medium','Thick','PGI');
    hold off
    subplot(2,1,2)
    mtitle = strcat('Low ',{' '},regression_type,': Low gamma band, Electrode ',{' '});
    mtitle = strcat(mtitle,elec_label);
    title(mtitle);
    hold on
    plot(time,thin_lgamma_low(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med_lgamma_low(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick_lgamma_low(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    ylim([-1 2.5]); xlim([-0.1 end_latency]);
    grid on
    yyaxis right
    plot(time,PGI_lgamma_low(elec_id,:),'Color','#7E2F8E','LineWidth',1)
    ylim([-3.5 2]); xlim([-0.1 end_latency]);
    legend('Thin','Medium','Thick','PGI');
    hold off

    figure(4)
    subplot(2,1,1)
    elec_label = chosen_electrode{4}.electrode;
    elec_id = find(strcmp(data_PGI{1}.label,elec_label));
    mtitle = strcat('High ',{' '},regression_type,': Medium gamma band, Electrode ',{' '});
    mtitle = strcat(mtitle,elec_label);
    title(mtitle);
    hold on
    plot(time,thin_mgamma_high(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med_mgamma_high(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick_mgamma_high(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    ylim([-1 2.5]); xlim([-0.1 end_latency]);
    grid on
    yyaxis right
    plot(time,PGI_mgamma_high(elec_id,:),'Color','#7E2F8E','LineWidth',1)
    ylim([-3.5 2]); xlim([-0.1 end_latency]);
    legend('Thin','Medium','Thick','PGI');
    hold off
    subplot(2,1,2)
    mtitle = strcat('Low ',{' '},regression_type,': Medium gamma band, Electrode ',{' '});
    mtitle = strcat(mtitle,elec_label);
    title(mtitle);
    hold on
    plot(time,thin_mgamma_low(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med_mgamma_low(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick_mgamma_low(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    ylim([-1 2.5]); xlim([-0.1 end_latency]);
    grid on
    yyaxis right
    plot(time,PGI_mgamma_low(elec_id,:),'Color','#7E2F8E','LineWidth',1)
    ylim([-3.5 2]); xlim([-0.1 end_latency]);
    legend('Thin','Medium','Thick','PGI');
    hold off

    figure(5)
    subplot(2,1,1)
    elec_label = chosen_electrode{5}.electrode;
    elec_id = find(strcmp(data_PGI{1}.label,elec_label));
    mtitle = strcat('High ',{' '},regression_type,': High gamma band, Electrode ',{' '});
    mtitle = strcat(mtitle,elec_label);
    title(mtitle);
    hold on
    plot(time,thin_hgamma_high(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med_hgamma_high(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick_hgamma_high(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    ylim([-1 2.5]); xlim([-0.1 end_latency]);
    grid on
    yyaxis right
    plot(time,PGI_hgamma_high(elec_id,:),'Color','#7E2F8E','LineWidth',1)
    ylim([-3.5 2]); xlim([-0.1 end_latency]);
    legend('Thin','Medium','Thick','PGI');
    hold off
    subplot(2,1,2)
    mtitle = strcat('Low ',{' '},regression_type,': High gamma band, Electrode ',{' '});
    mtitle = strcat(mtitle,elec_label);
    title(mtitle);
    hold on
    plot(time,thin_hgamma_low(elec_id,:),'Color','#0072BD','LineWidth',0.8)
    plot(time,med_hgamma_low(elec_id,:),'Color','#D95319','LineWidth',0.8)
    plot(time,thick_hgamma_low(elec_id,:),'Color','#EDB120','LineWidth',0.8)
    ylim([-1 2.5]); xlim([-0.1 end_latency]);
    grid on
    yyaxis right
    plot(time,PGI_hgamma_low(elec_id,:),'Color','#7E2F8E','LineWidth',1)
    ylim([-3.5 2]); xlim([-0.1 end_latency]);
    legend('Thin','Medium','Thick','PGI');
    hold off
end

%% topoplots over time, highlighted electrodes
function create_topoplots(data_PGI,stat,tail)
    %fixes the interpolation issue of NaNs
    stat.stat(isnan(stat.stat))=0;

    start_latency = 0; end_latency = 0.8; 
    timestep = 0.05; % in seconds (50ms)
    sampling_rate = 512;
    sample_count  = length(stat.time);

    j = [start_latency:timestep:end_latency]; % Temporal endpoints (in seconds) of the ERP average computed in each subplot
    [i1,i2] = match_str(data_PGI{1}.label, stat.label);
    
    if strcmp(tail,'positive') && ~isempty(stat.posclusters)
        pos_cluster_pvals = [stat.posclusters(:).prob];
        pos_clust = find(pos_cluster_pvals < 0.05);
%         for i = 1:length(pos_clust)
%             if stat.posclusters(i).clusterstat < 2000 %ask about size
%                 pos_clust(i)=0;
%             end
%         end
%         pos_clust = pos_clust(pos_clust~=0);
        pos = ismember(stat.posclusterslabelmat, pos_clust);
        if isempty(pos_clust)
            clust = zeros(118,1,401); %401 before
        else
            clust = pos;
        end
    elseif strcmp(tail,'negative') && ~isempty(stat.negclusters)
        neg_cluster_pvals = [stat.negclusters(:).prob];
        neg_clust = find(neg_cluster_pvals < 0.05);
        neg = ismember(stat.negclusterslabelmat, neg_clust);
        if isempty(neg_clust)
            clust = zeros(118,1,401); %401 before
        else
            clust = neg;
        end
    end

    max_iter = numel(j)-1;
    for k = 1:max_iter
         subplot(2,8,k);
         cfg = [];
         cfg.parameter = 'stat';
         cfg.xlim = [j(k) j(k+1)];
         cfg.zlim = [-6 6];
         %cfg.colorbar = 'yes';
         clust_int = zeros(118,1);
         t_idx = find(stat.time>=j(k) & stat.time<=j(k+1));
         clust_int(i1) = all(clust(i2,1,t_idx),3);
         cfg.marker ='on';
         cfg.markersize = 1;
         cfg.highlight = 'on';
         cfg.highlightchannel = find(clust_int);
         cfg.highlightcolor = {'r',[0 0 1]};
         cfg.highlightsize = 2;
         cfg.comment = 'no';
         %cfg.commentpos = 'title';
         ft_topoplotTFR(cfg, stat);
    end
end

%% plot the design matrix
function plot_design_matrix(design_matrix, n_participants)
    plot(design_matrix(1:n_participants), 'color', 'r', 'LineWidth', 1);
    hold on;
    plot(design_matrix(n_participants+1:n_participants*2), 'color', 'g', 'LineWidth', 1);
    plot(design_matrix((n_participants*2)+1:n_participants*3), 'color', 'b', 'LineWidth', 1);
    xlabel('Participants');
    ylabel('Interaction');
    legend({'Partition 1', 'Partition 2', 'Partition 3'},'Location','northwest')
    set(gcf,'Position',[100 100 1000 1000])  
end

%% Create a hacked ROI in the data
function new_data = create_hacked_roi(data, roi)
    roi_clusters = roi.clusters;
    roi_freq = roi.freq;
    roi_time = roi.time;

    start_latency = NaN;
    end_latency = NaN;
    new_data = {};
    for idx_i = 1:numel(data)
        each_participant = data{idx_i};
        participant_data = each_participant.powspctrm;
        participant_freq = each_participant.freq;
        participant_time = each_participant.time;
        [electrodes, freq, time] = size(participant_data);
        new_participant_data = NaN(electrodes, freq, time);
        time_x = 1:1:time;
        time_x = time_x * 2;

        for roi_idx = 1:numel(roi_time)
            t = roi_time(roi_idx)*1000;
            [~,idx]=min(abs(time_x-t));
            %idx = idx + 100; % for baselining period
            clusters_at_t = roi_clusters(:,roi_idx);

            if isnan(start_latency) && sum(clusters_at_t)>=1
                start_latency =t;
            elseif roi_idx == numel(roi_time)
                end_latency = t; 
            end

            for electrode_idx=1:numel(clusters_at_t)
                if clusters_at_t(electrode_idx) == 1
                    new_participant_data(electrode_idx,:,idx) = participant_data(electrode_idx,:,idx);
                end
            end

        end

        %[~,start_idx]=min(abs(time_x-start_latency));
        %[~,end_idx]=min(abs(time_x-end_latency));
        %new_participant_data = new_participant_data(:,start_idx:end_idx);
        each_participant.powspctrm = new_participant_data;  

        %time = each_participant.time;
        %[~,start_idx]=min(abs(time-start_latency/1000));
        %[~,end_idx]=min(abs(time-end_latency/1000));

        new_data{idx_i} = each_participant;
    end
end

%% based on onsets 2-8 mean intercept effect, save the electrodes 
% which are at the peak level - used for ROI analysis
% saves the neighbourhood of electrodes for use later on
function get_region_of_interest_electrodes(stat, desired_cluster, experiment_type, roi_applied,freq_band)
    roi_clusters = squeeze(stat.posclusterslabelmat);
    roi_clusters(roi_clusters>desired_cluster) = 0;
    time = stat.time;
    if freq_band == 1
        roi_a.clusters = roi_clusters;
        roi_a.freq = 8:13;
        roi_a.time = time;
    elseif freq_band == 2
        roi_b.clusters = roi_clusters;
        roi_b.freq = 20:35;
        roi_b.time = time;
    elseif freq_band == 3
        roi_lg.clusters = roi_clusters;
        roi_lg.freq = 30:45;
        roi_lg.time = time;
    elseif freq_band == 4
        roi_mg.clusters = roi_clusters;
        roi_mg.freq = 45:60;
        roi_mg.time = time;
    elseif freq_band == 5
        roi_hg.clusters = roi_clusters;
        roi_hg.freq = 60:80;
        roi_hg.time = time;
    end
    
    if freq_band == 1
        if contains(experiment_type, 'onsets-2-8')  && contains(roi_applied, 'one-tailed')
            save 'C:\Users\marga\Desktop\Research Project\scripts\one_tailed_roi_28_a.mat' roi_a; 
        elseif contains(experiment_type, 'onsets-2-8')  && contains(roi_applied, 'two-tailed')
            save 'C:\Users\marga\Desktop\Research Project\scripts\two_tailed_roi_28_a.mat' roi_a; 
        end
    elseif freq_band == 2
        if contains(experiment_type, 'onsets-2-8')  && contains(roi_applied, 'one-tailed')
            save 'C:\Users\marga\Desktop\Research Project\scripts\one_tailed_roi_28_b.mat' roi_b; 
        elseif contains(experiment_type, 'onsets-2-8')  && contains(roi_applied, 'two-tailed')
            save 'C:\Users\marga\Desktop\Research Project\scripts\two_tailed_roi_28_b.mat' roi_b; 
        end
    elseif freq_band == 3
        if contains(experiment_type, 'onsets-2-8')  && contains(roi_applied, 'one-tailed')
            save 'C:\Users\marga\Desktop\Research Project\scripts\one_tailed_roi_28_lg.mat' roi_lg; 
        elseif contains(experiment_type, 'onsets-2-8')  && contains(roi_applied, 'two-tailed')
            save 'C:\Users\marga\Desktop\Research Project\scripts\two_tailed_roi_28_lg.mat' roi_lg; 
        end
    elseif freq_band == 4
        if contains(experiment_type, 'onsets-2-8')  && contains(roi_applied, 'one-tailed')
            save 'C:\Users\marga\Desktop\Research Project\scripts\one_tailed_roi_28_mg.mat' roi_mg; 
        elseif contains(experiment_type, 'onsets-2-8')  && contains(roi_applied, 'two-tailed')
            save 'C:\Users\marga\Desktop\Research Project\scripts\two_tailed_roi_28_mg.mat' roi_mg; 
        end
    elseif freq_band == 5
        if contains(experiment_type, 'onsets-2-8')  && contains(roi_applied, 'one-tailed')
            save 'C:\Users\marga\Desktop\Research Project\scripts\one_tailed_roi_28_hg.mat' roi_hg; 
        elseif contains(experiment_type, 'onsets-2-8')  && contains(roi_applied, 'two-tailed')
            save 'C:\Users\marga\Desktop\Research Project\scripts\two_tailed_roi_28_hg.mat' roi_hg; 
        end
    end
end

%% calculate grand average for each frequency band
function grand_avgs = calculate_freq_avg(data,end_latency)
    freq_bins = {[8 13],[20 35],[30 45],[45 60],[60 80]};
    grand_avgs = cell(1,5);
    for f = 1:5
        cfg = [];
        cfg.channel   = 'all';
        cfg.foilim = freq_bins{f};
        cfg.toilim = [-0.1 end_latency];
        cfg.latency   = 'all';
        cfg.parameter = 'powspctrm';
        grand_avg = ft_freqgrandaverage(cfg, data{:});
        grand_avg.elec = data{1}.elec;
        grand_avgs{f} = grand_avg;
    end
end

%% grand average for the 3 conditions separately
function grand_avgs_condition = calculate_avg_conditions(data_condition,end_latency) 
    grand_avgs_condition = cell(1,5);
    freq_bins = {[8 13],[20 35],[30 45],[45 60],[60 80]};
    for f = 1:5
        cfg = [];
        cfg.channel   = 'all';
        cfg.foilim = freq_bins{f};
        cfg.toilim = [-0.1 end_latency];
        cfg.latency   = 'all';
        cfg.parameter = 'powspctrm';
        grand_avg_thin = ft_freqgrandaverage(cfg, data_condition{:});
        grand_avgs_condition{1,f} = grand_avg_thin;
    end
end

%% calculate median splits
function [data_high, data_low] = get_median_split(data, participant_order, regression_type, partition)
    function split_data = get_participants(data, all_ids, current_ids)
        cnt = 1;
        split_data = {};
        for i=1:numel(data)
            participant = data{i};
            id = all_ids(i);
            
            if ismember(id, current_ids)
                split_data{cnt} = participant;
                cnt = cnt+1;
            end     
        end      
    end    

    scores = return_scores(regression_type);
    
    if partition == 1
        ratings = scores.one;
    elseif partition == 2
        ratings = scores.two;
    elseif partition == 3
        ratings = scores.three;
    end
    
    sorted(:,1) = ratings(:,2);
    sorted(:,2) = ratings(:,1);
    sorted = flipud(sortrows(sorted));
    
    n = ceil(numel(sorted(:,1))/2);
    high = sorted(1:n,:);
    high_ids = high(:,2);
    low = sorted(n+1:end,:);
    low_ids = low(:,2);
    
    data_high = get_participants(data, participant_order, high_ids);
    data_low = get_participants(data, participant_order, low_ids);
end

%% used to calculate the cluster size through time
function calculate_cluster_size(stat, desired_cluster, make_plots, ptitle, xlim_t, type)
    if contains(type, 'positive')
        cluster_labelling_mtx = stat.posclusterslabelmat;
    elseif contains(type, 'negative')
        cluster_labelling_mtx = stat.negclusterslabelmat;
    end
    
    time_mtx = stat.time;
    
    [electrodes, time] = size(cluster_labelling_mtx);
    
    cluster_size_through_time = zeros(2,time);
    
    for i = 1:time
        t = time_mtx(i);
        electrodes_in_time = cluster_labelling_mtx(:,i);
        clusters_in_time = find(electrodes_in_time==desired_cluster);
        cluster_size_through_time(1,i) = t;
        cluster_size_through_time(2,i) = numel(clusters_in_time)/electrodes;
    end
    
    if make_plots == 'yes'
        area(cluster_size_through_time(1,:)*1000, cluster_size_through_time(2,:)*100);
        ylim([0 100]);
        grid on;
        xlabel('Time (ms)');
        ylabel('Percentage of cluster');
        xlim([0,xlim_t])
        title(ptitle);
    end
    
end

%% compute best electrode
function peak_stat_info = compute_best_electrode_from_t_values(stat, electrode_stats,tail, peak_stat_info)
    significant_masks = stat.mask;
    electrodes = electrode_stats.electrodes;
    
    if strcmp(tail, 'positive')
        cluster_labels = stat.posclusterslabelmat;
    else
        cluster_labels = stat.negclusterslabelmat;
    end
    
    time = stat.time;
    most_significant_electrode = electrode_stats.electrodes{1};
    previous_best = 0;
    
    for i=1:numel(electrodes)
        electrode_idx = find(strcmp(stat.label,electrodes(i)));
        significant_time_series = significant_masks(electrode_idx,:);
        cluster_time_series = cluster_labels(electrode_idx,:);
        raw_time_series = zeros(1,numel(cluster_time_series));
        significance_count = 0;
        for j=1:numel(significant_time_series)
            if cluster_time_series(j)==1
                is_significant = significant_time_series(j);
                significance_count = significance_count + is_significant;
                raw_time_series(j) = significance_count;
            else
                raw_time_series(j) = significance_count;
            end
        end
        if previous_best < significance_count && strcmp(tail, 'positive')
            e = electrodes(i);
            most_significant_electrode = e;
            previous_best = significance_count;
        elseif previous_best > significance_count && strcmp(tail, 'negative')
            e = electrodes(i);
            most_significant_electrode = e;
            previous_best = significance_count;
        end
%         plot(time, raw_time_series);
%         hold on;
    end
    
   [col, ~] = size(electrode_stats);
    
   for i = 1:col
      if strcmp(electrode_stats.electrodes{i}, most_significant_electrode)
          peak_stat_info.electrode = most_significant_electrode;
          peak_stat_info.time = electrode_stats.time(i);
          peak_stat_info.t_value = electrode_stats.t_value(i);
      end
   end
%     xlabel('Time (ms)');
%     ylabel('Cumulative Frequency (If T-value is significant');
%     legend(electrodes,'Location','northwest')
%     title(strcat('Sustained Significant Electrode:',most_significant_electrode)); 
end

%% used to generate 'peak' level stats
function [peak_level_electrode, all_electrode_stats] = get_peak_level_stats(stat, cluster, type)
    peak_level_stats = containers.Map;

    time = stat.time;
    electrodes = stat.label;
    t_values = stat.stat;
    if contains(type, 'positive')
        cluster_matrix_locations = stat.posclusterslabelmat;
    elseif contains(type, 'negative')
        cluster_matrix_locations = stat.negclusterslabelmat;
    end
    
    [rows, columns] = size(cluster_matrix_locations);
    
    for col = 1:columns
        for row = 1:rows
            
            t_value = t_values(row, col);
            t = time(col);
            electrode = electrodes{row};

            which_cluster = cluster_matrix_locations(row, col);
            if which_cluster == cluster
                keys = peak_level_stats.keys;
                if any(strcmp(keys,electrode))
                    previous_values = peak_level_stats(electrode);
                    previous_t = previous_values{2};
                    if t_value > previous_t && contains(type, 'positive')
                        peak_level_stats(electrode) = {t, t_value};
                    elseif t_value < previous_t && contains(type, 'negative')
                        peak_level_stats(electrode) = {t, t_value};
                    end
                else
                    peak_level_stats(electrode) = {0, 0};
                end
            end
        end 
    end
    peak_level_electrode = get_most_significant_electrode(peak_level_stats, type);
    all_electrode_stats = get_all_electrode_stats(peak_level_stats, type);
end

function mtx = get_all_electrode_stats(peak_level_stats, type)
    keys = peak_level_stats.keys;
    electrodes = cell(numel(keys), 1);
    time = zeros(numel(keys), 1);
    t_value = zeros(numel(keys), 1);
    
    for i = 1:numel(keys)
        electrode = keys{i};
        stats = peak_level_stats(electrode);
        tv = stats{2};
        t = stats{1};
        
        electrodes{i} = electrode;
        time(i) = t;
        t_value(i) = tv;
    end
    
    mtx = table(electrodes, time, t_value);
    if strcmp(type, 'positive')
        mtx = sortrows(mtx, {'t_value'}, {'descend'});
    elseif strcmp(type, 'negative')
        mtx = sortrows(mtx, {'t_value'}, {'ascend'});
    end
end

function most_significant = get_most_significant_electrode(peak_level_stats, type)
    keys = peak_level_stats.keys;
    
    most_significant.electrode = '';
    most_significant.time = 0;
    most_significant.t_value = 0;
    for i = 1:numel(keys)
        electrode = keys{i};
        stats = peak_level_stats(electrode);
        t_value = stats{2};
        time = stats{1};
        if t_value > most_significant.t_value && contains(type, 'positive')
            most_significant.electrode = electrode;
            most_significant.time = time;
            most_significant.t_value = t_value;
        elseif t_value < most_significant.t_value && contains(type, 'negative')
            most_significant.electrode = electrode;
            most_significant.time = time;
            most_significant.t_value = t_value; 
        end
    end
end

%% set all my values to 0 to hack a T test
function data = set_values_to_zero(data)
    for idx = 1:numel(data)
        participant_data = data{1,idx};
        spectrum = participant_data.powspctrm;
        spectrum(:) = 0;
        participant_data.powspctrm = spectrum;
        data{1,idx} = participant_data;
    end
end
%% function to create the design matrix for all experiments
function [design, new_participants] = create_design_matrix_partitions(participant_order, participants, ...
    regression_type, partition)
    
    scores = return_scores(regression_type);
    
    if partition == 1
        ratings = scores.one;
    elseif partition == 2
        ratings = scores.two;
    elseif partition == 3
        ratings = scores.three;
    end
    
    new_participants = {};
    cnt = 1;
    for j = 1:numel(participant_order)
        participant = participant_order(j);
        score = ratings(find(ratings(:,1)==participant),2);
        
        if numel(score) > 0
            design(1,cnt)=score;
            new_participants{cnt} = participants{j};
            cnt = cnt + 1;
        end
    end
    
end

function scores = return_scores(regression_type)
    if contains(regression_type, 'headache')
        scores.one = [
        1,-0.227;2,-0.052;3,-0.721;4,0.531;5,1.72;6,-1.176;7,-0.797;
        8,0.199;9,-0.692;10,-0.358;11,-0.585;12,2.041;13,-0.266;14,1.31;
        16,0.002;17,-0.15;20,-0.416;21,1.144;22,-0.565;23,-0.728;24,-0.435;
        25,-0.699;26,-1.35;28,0.363;29,0.042;30,-0.47;31,-0.704;32,-0.732;
        33,-0.881;34,-0.796;37,1.227;38,-0.337;39,1.077;40,-0.167
        ];
    
        scores.two = [
        1,-0.227;2,-0.052;3,-0.721;4,0.531;5,1.72;6,-1.176;7,-0.797;
        8,0.199;9,-0.692;10,-0.358;11,-0.585;12,2.041;13,-0.266;14,1.31;
        16,0.002;17,-0.15;20,-0.416;21,1.144;22,-0.565;23,-0.728;24,-0.435;
        25,-0.699;26,-1.35;28,0.363;29,0.042;30,-0.47;31,-0.704;32,-0.732;
        33,-0.881;34,-0.796;37,1.227;38,-0.337;39,1.077;40,-0.167
        ];
    
        scores.three = [
        1,-0.227;2,-0.052;3,-0.721;4,0.531;5,1.72;6,-1.176;7,-0.797;
        8,0.199;9,-0.692;10,-0.358;11,-0.585;12,2.041;13,-0.266;14,1.31;
        16,0.002;17,-0.15;20,-0.416;21,1.144;22,-0.565;23,-0.728;24,-0.435;
        25,-0.699;26,-1.35;28,0.363;29,0.042;30,-0.47;31,-0.704;32,-0.732;
        33,-0.881;34,-0.796;37,1.227;38,-0.337;39,1.077;40,-0.167
        ];
        
        min_n = min(scores.one);
        scores.one(:,2) = scores.one(:,2) - min_n(2);
        scores.two(:,2) = scores.two(:,2) - min_n(2);
        scores.three(:,2) = scores.three(:,2) - min_n(2);
        
        scores.one(:,2) = scores.one(:,2) * 2.72;
        scores.two(:,2) = scores.two(:,2) * 1.65;
        scores.three(:,2) = scores.three(:,2) * 1.0;
        
        [n_participants, ~] = size(scores.one);
        
        for k=1:n_participants
            p1 = scores.one(k,1);
            p2 = scores.two(k,1);
            p3 = scores.three(k,1);
            
            if p1 == p2 && p2 == p3                
                to_remove = scores.three(k,2);
                scores.one(k,2) = scores.one(k,2) - to_remove;
                scores.two(k,2) = scores.two(k,2) - to_remove;
                scores.three(k,2) = scores.three(k,2) - to_remove;
            end
        end 
        
%         scores.one = [
%         1,0.007750386;2,0.67380154;3,-0.137934691;4,0.626158247;5,4.169589335;6,-0.324759417;7,-0.420847527;
%         8,0.378004967;9,-0.189543248;10,0.665549906;11,0.314848605;12,5.080009397;13,0.973407875;14,4.164128033;
%         16,1.687829269;17,0.283028132;20,-0.540358266;21,1.792972384;22,-0.201718237;23,-2.03537896;24,-0.196166288;
%         25,-1.056735762;26,-1.018688066;28,0.72589105;29,0.763670478;30,-0.1676338;31,-1.131270977;32,-0.334399666;
%         33,-0.294824707;37,2.106357136;38,0.588400524;39,0.552057051;40,-1.531878069];
%         
%         scores.two = [
%         1,-0.252823403;2,-0.001060599;3,-0.307922469;4,-0.019147587;5,1.319999324;6,-0.37845027;7,-0.414806971;
%         8,-0.112952699;9,-0.32734547;10,-0.0042185;11,-0.136784091;12,1.664007262;13,0.112132493;14,1.317928033;
%         16,0.382139997;17,-0.148776409;20,-0.459934006;21,0.421851639;22,-0.332005732;23,-1.024886369;24,-0.329900185;
%         25,-0.655136114;26,-0.640715183;28,0.018557252;29,0.032869376;30,-0.319092125;31,-0.683246065;32,-0.382087721;
%         33,-0.3672123;37,0.540262025;38,-0.033412518;39,-0.047133675;40,-0.834656025];
%         
%         scores.three = [
%         1,-0.411071967;2,-0.411071967;3,-0.411071967;4,-0.411071967;5,-0.411071967;6,-0.411071967;7,-0.411071967;
%         8,-0.411071967;9,-0.411071967;10,-0.411071967;11,-0.411071967;12,-0.411071967;13,-0.411071967;14,-0.411071967;
%         16,-0.411071967;17,-0.411071967;20,-0.411071967;21,-0.411071967;22,-0.411071967;23,-0.411071967;24,-0.411071967;
%         25,-0.411071967;26,-0.411071967;28,-0.411071967;29,-0.411071967;30,-0.411071967;31,-0.411071967;32,-0.411071967;
%         33,-0.411071967;37,-0.411071967;38,-0.411071967;39,-0.411071967;40,-0.411071967];
        
    elseif contains(regression_type, 'visual stress')
        scores.one = [
        1,0.323;2,-0.109;3,-0.51;4,1.134;5,-0.639;6,-1.215;7,-0.33;8,0.752;
        9,-0.39;10,-0.722;11,-0.769;12,-1.063;13,-0.899;14,-1.467;16,-1.199;
        17,0.154;20,0.587;21,1.001;22,-0.117;23,1.721;24,0.141;25,0.622;
        26,-0.748;28,0.674;29,-0.024;30,0.036;31,0.7;32,-0.3;33,-0.65;
        34,0.026;37,0.799;38,-0.588;39,2.333;40,2.267;
        ];
    
        scores.two = [
        1,0.323;2,-0.109;3,-0.51;4,1.134;5,-0.639;6,-1.215;7,-0.33;8,0.752;
        9,-0.39;10,-0.722;11,-0.769;12,-1.063;13,-0.899;14,-1.467;16,-1.199;
        17,0.154;20,0.587;21,1.001;22,-0.117;23,1.721;24,0.141;25,0.622;
        26,-0.748;28,0.674;29,-0.024;30,0.036;31,0.7;32,-0.3;33,-0.65;
        34,0.026;37,0.799;38,-0.588;39,2.333;40,2.267;
        ];
    
        scores.three = [
        1,0.323;2,-0.109;3,-0.51;4,1.134;5,-0.639;6,-1.215;7,-0.33;8,0.752;
        9,-0.39;10,-0.722;11,-0.769;12,-1.063;13,-0.899;14,-1.467;16,-1.199;
        17,0.154;20,0.587;21,1.001;22,-0.117;23,1.721;24,0.141;25,0.622;
        26,-0.748;28,0.674;29,-0.024;30,0.036;31,0.7;32,-0.3;33,-0.65;
        34,0.026;37,0.799;38,-0.588;39,2.333;40,2.267;
        ];
    
        min_n = min(scores.one);
        scores.one(:,2) = scores.one(:,2) - min_n(2);
        scores.two(:,2) = scores.two(:,2) - min_n(2);
        scores.three(:,2) = scores.three(:,2) - min_n(2);
        
        scores.one(:,2) = scores.one(:,2) * 2.72;
        scores.two(:,2) = scores.two(:,2) * 1.65;
        scores.three(:,2) = scores.three(:,2) * 1.0;
        
        [n_participants, ~] = size(scores.one);
        
        for k=1:n_participants
            p1 = scores.one(k,1);
            p2 = scores.two(k,1);
            p3 = scores.three(k,1);
            
            if p1 == p2 && p2 == p3                
                to_remove = scores.three(k,2);
                scores.one(k,2) = scores.one(k,2) - to_remove;
                scores.two(k,2) = scores.two(k,2) - to_remove;
                scores.three(k,2) = scores.three(k,2) - to_remove;
            end
        end     
        
        
    elseif contains(regression_type, 'discomfort')
        scores.one = [
        1,-0.264;2,0.446;3,-0.498;4,1.777;5,-0.556;6,0.872;7,-0.685;8,0.928;
        9,-0.806;10,-0.875;11,0.391;12,-0.761;13,-0.69;14,1.608;16,1.14;
        17,1.536;20,0.122;21,0.084;22,0.617;23,-1.48;24,2.284;25,-0.809;26,-0.557;
        28,-0.933;29,0.379;30,-0.631;31,2.147;32,-1.499;33,1.22;34,-0.797;
        37,-0.613;38,-1.026;39,-0.877;40,0.444;
        ];
    
        scores.two = [
        1,-0.264;2,0.446;3,-0.498;4,1.777;5,-0.556;6,0.872;7,-0.685;8,0.928;
        9,-0.806;10,-0.875;11,0.391;12,-0.761;13,-0.69;14,1.608;16,1.14;
        17,1.536;20,0.122;21,0.084;22,0.617;23,-1.48;24,2.284;25,-0.809;26,-0.557;
        28,-0.933;29,0.379;30,-0.631;31,2.147;32,-1.499;33,1.22;34,-0.797;
        37,-0.613;38,-1.026;39,-0.877;40,0.444;
        ];
    
        scores.three = [
        1,-0.264;2,0.446;3,-0.498;4,1.777;5,-0.556;6,0.872;7,-0.685;8,0.928;
        9,-0.806;10,-0.875;11,0.391;12,-0.761;13,-0.69;14,1.608;16,1.14;
        17,1.536;20,0.122;21,0.084;22,0.617;23,-1.48;24,2.284;25,-0.809;26,-0.557;
        28,-0.933;29,0.379;30,-0.631;31,2.147;32,-1.499;33,1.22;34,-0.797;
        37,-0.613;38,-1.026;39,-0.877;40,0.444;
        ];
    
        min_n = min(scores.one);
        scores.one(:,2) = scores.one(:,2) - min_n(2);
        scores.two(:,2) = scores.two(:,2) - min_n(2);
        scores.three(:,2) = scores.three(:,2) - min_n(2);
        
        scores.one(:,2) = scores.one(:,2) * 2.72;
        scores.two(:,2) = scores.two(:,2) * 1.65;
        scores.three(:,2) = scores.three(:,2) * 1.0;
        
        [n_participants, ~] = size(scores.one);
        
        for k=1:n_participants
            p1 = scores.one(k,1);
            p2 = scores.two(k,1);
            p3 = scores.three(k,1);
            
            if p1 == p2 && p2 == p3                
                to_remove = scores.three(k,2);
                scores.one(k,2) = scores.one(k,2) - to_remove;
                scores.two(k,2) = scores.two(k,2) - to_remove;
                scores.three(k,2) = scores.three(k,2) - to_remove;
            end
        end     
%         scores.one = [
%         1,-0.108255668;2,1.230945421;3,0.29906353;4,2.419686619;5,-1.438369996;6,3.366738286;7,-0.069516957;
%         8,1.404643457;9,-0.317741562;10,-0.485463307;11,1.877163926;12,-1.77964955;13,-0.133164122;14,3.031899157;
%         16,2.993830705;17,3.027880878;20,0.540954307;21,-0.847472812;22,1.888775558;23,-2.627430661;24,4.521337032;
%         25,-0.881617779;26,0.770653733;28,-1.867665377;29,1.002567227;30,-0.40947651;31,4.162733455;32,-1.531847366;
%         33,3.44408667;37,-1.991999175;38,-0.834492896;39,-3.19529217;40,-0.015083387];
%         
%         scores.two = [
%         1,-0.32035253;2,0.185725439;3,-0.166412743;4,0.63503181;5,-0.823033363;6,0.992878779;7,-0.305679307;
%         8,0.251396289;9,-0.399563034;10,-0.462892958;11,0.430015513;12,-0.951940453;13,-0.329776322;14,0.866337876;
%         16,0.851927607;17,0.864797823;20,-0.075009595;21,-0.599726118;22,0.434406293;23,-1.272386282;24,1.429201861;
%         25,-0.612569055;26,0.011833712;28,-0.985194577;29,0.099441499;30,-0.434173526;31,1.293718563;32,-0.858316972;
%         33,1.022098328;37,-1.03216192;38,-0.594831504;39,-1.486896968;40,-0.285154745];
%         
%         scores.three = [
%         1,-0.449133931;2,-0.449133931;3,-0.449133931;4,-0.449133931;5,-0.449133931;6,-0.449133931;7,-0.449133931;
%         8,-0.449133931;9,-0.449133931;10,-0.449133931;11,-0.449133931;12,-0.449133931;13,-0.449133931;14,-0.449133931;
%         16,-0.449133931;17,-0.449133931;20,-0.449133931;21,-0.449133931;22,-0.449133931;23,-0.449133931;24,-0.449133931;
%         25,-0.449133931;26,-0.449133931;28,-0.449133931;29,-0.449133931;30,-0.449133931;31,-0.449133931;32,-0.449133931;
%         33,-0.449133931;37,-0.449133931;38,-0.449133931;39,-0.449133931;40,-0.449133931];
    
    end
end

%% return the SPM data in a fieldtrip format
 function [thin_data,med_data,thick_data,fieldtrip_data_agg,fieldtrip_data_PGI,participant_order] = from_fieldtrip_to_spm(n_participants, main_path, filename, partition)
    fieldtrip_data_agg = {};
    fieldtrip_data_PGI = {};
    thin_data = {}; med_data = {}; thick_data = {};
    participant_order = [];
    cnt = 1;
    for participant = 24:26
        disp(participant)
        
        participant_main_path = strcat(main_path, int2str(participant));
        if exist(participant_main_path, 'dir')
            if participant < 10
               p = strcat('0', int2str(participant));
            else
               p = int2str(participant);
            end

            % data structure
            data_structure = strcat('spmeeg_P', p);  
            data_structure = strcat(data_structure, '_075_80hz_rejected_tempesta.mat');
            data_structure = strcat(filename, data_structure); 
            cd(participant_main_path); 
            cd("SPM_ARCHIVE")

            if isfile(data_structure)
                load(data_structure);
            else
                continue;
            end

            med = {}:
            thin = {};
            thick = {};

            spm_eeg = meeg(D);
            fieldtrip_raw = spm_eeg.ftraw;
            
            n_trials = size(D.trials);
            n_trials = n_trials(2);
            mt=1; tht=1; tt=1;
            for index_i = 1:n_trials
                label = D.trials(index_i).label;
                if partition.is_partition
                    if contains(label, partition.partition_number) && contains(label, 'medium')
                        med(mt) = fieldtrip_raw.trial(index_i);
                        mt=mt+1;
                    elseif contains(label, partition.partition_number) && contains(label, 'thin')
                        thin(tht) = fieldtrip_raw.trial(index_i);
                        tht=tht+1;
                    elseif contains(label, partition.partition_number) && contains(label, 'thick')
                        thick(tt) = fieldtrip_raw.trial(index_i);
                        tt=tt+1;
                    end
                else
                     if contains(label, 'medium')
                        med(mt) = fieldtrip_raw.trial(index_i);
                        mt=mt+1;
                    elseif contains(label, 'thin')
                        thin(tht) = fieldtrip_raw.trial(index_i);
                        tht=tht+1;
                    elseif contains(label, 'thick')
                        thick(tt) = fieldtrip_raw.trial(index_i);
                        tt=tt+1;
                     end
                end
            end
            
            if isempty(med) || isempty(thin) || isempty(thick)
                continue;
            end
            
            % update the fieldtrip structure with fields of information
            raw_med = [];
            raw_med.label = fieldtrip_raw.label;
            raw_med.elec = fieldtrip_raw.elec;
            raw_med.trial = med;
            raw_med.time = fieldtrip_raw.time(1:length(raw_med.trial));
            raw_med.dimord = 'chan_time';
            %raw_med = remove_electrodes(raw_med);
            
            raw_thin = [];
            raw_thin.label = fieldtrip_raw.label;
            raw_thin.elec = fieldtrip_raw.elec;
            raw_thin.trial = thin;
            raw_thin.time = fieldtrip_raw.time(1:length(raw_thin.trial));
            raw_thin.dimord = 'chan_time';
            %raw_thin = remove_electrodes(raw_thin);
            
            raw_thick = [];
            raw_thick.label = fieldtrip_raw.label;
            raw_thick.elec = fieldtrip_raw.elec;
            raw_thick.trial = thick;
            raw_thick.time = fieldtrip_raw.time(1:length(raw_thick.trial));
            raw_thick.dimord = 'chan_time';
            %raw_thick = remove_electrodes(raw_thick);
            
%             %apply a notch filter at 50Hz
%             cfg = [];
%             cfg.channel = 'all';
%             cfg.bsfilter = 'yes';
%             cfg.bsfreq = [49 51];
%             TFRwave_med = ft_preprocessing(cfg,raw_med);
%             TFRwave_thin = ft_preprocessing(cfg,raw_thin);
%             TFRwave_thick = ft_preprocessing(cfg,raw_thick);
            
            %wavelet decomposition
            cfg = [];
            cfg.channel = 'all';
            cfg.method = 'wavelet';
            cfg.width = 5;
            cfg.output = 'pow';
            cfg.pad = 'nextpow2';
            cfg.foi = 5:80;
            cfg.toi = -0.5:0.002:1.2;
            cfg.keeptrials = 'yes';
            TFRwave_med = ft_freqanalysis(cfg, raw_med);
            TFRwave_thin = ft_freqanalysis(cfg, raw_thin);
            TFRwave_thick = ft_freqanalysis(cfg, raw_thick);
            
            TFRwave_med.info = 'medium'; 
            TFRwave_thin.info = 'thin'; 
            TFRwave_thick.info = 'thick'; 
            
            %crop the epoch
            TFRwave_med.time = TFRwave_med.time(200:end); %200 before, 46 after
            TFRwave_med.powspctrm = TFRwave_med.powspctrm(:,:,:,200:end);
            TFRwave_thin.time = TFRwave_thin.time(200:end); 
            TFRwave_thin.powspctrm = TFRwave_thin.powspctrm(:,:,:,200:end);
            TFRwave_thick.time = TFRwave_thick.time(200:end);
            TFRwave_thick.powspctrm = TFRwave_thick.powspctrm(:,:,:,200:end);
           
            
            %average across trials
            cfg = [];
            avg_TFRwave_med = ft_freqdescriptives(cfg, TFRwave_med);
            avg_TFRwave_thin = ft_freqdescriptives(cfg, TFRwave_thin);
            avg_TFRwave_thick = ft_freqdescriptives(cfg, TFRwave_thick);
            
            %baseline rescale
            cfg = [];
            cfg.baselinetype = 'db'; 
            cfg.baseline = [-0.1 0];
            avg_TFRwave_med = ft_freqbaseline(cfg,avg_TFRwave_med);
            avg_TFRwave_thin = ft_freqbaseline(cfg,avg_TFRwave_thin);
            avg_TFRwave_thick = ft_freqbaseline(cfg,avg_TFRwave_thick);
            
            TFRwave_aggregated = avg_TFRwave_med;
            TFRwave_aggregated.info = 'Aggregated conditions';
            TFRwave_aggregated.powspctrm = (avg_TFRwave_med.powspctrm+avg_TFRwave_thin.powspctrm+avg_TFRwave_thick.powspctrm)/3;
            TFRwave_aggregated.elec = fieldtrip_raw.elec;
            
            TFRwave_PGI = avg_TFRwave_med;
            TFRwave_PGI.info = 'PGI'; 
            TFRwave_PGI.powspctrm = avg_TFRwave_med.powspctrm-(avg_TFRwave_thin.powspctrm+avg_TFRwave_thick.powspctrm)/2; %%!!!!%%%
            TFRwave_PGI.elec = fieldtrip_raw.elec;
            
            % update object with fieldtrip data
            fieldtrip_data_agg{cnt} = TFRwave_aggregated;
            fieldtrip_data_PGI{cnt} = TFRwave_PGI;
            thin_data{cnt} = avg_TFRwave_thin;
            med_data{cnt} = avg_TFRwave_med;
            thick_data{cnt} = avg_TFRwave_thick;
            participant_order(cnt) = participant;
            cnt = cnt + 1;
        end
    end
end

function fieldtrip_raw = remove_electrodes(fieldtrip_raw)
    %to_remove = {'A11', 'A12', 'A13', 'A14', 'A24', 'A25', 'A26','A27', 'B8', 'B9','EXG1', 'EXG2', 'EXG3', 'EXG4', 'EXG5', 'EXG6', 'HEOG', 'VEOG'};
    to_remove = {'EXG1', 'EXG2', 'EXG3', 'EXG4', 'EXG5', 'EXG6', 'HEOG', 'VEOG'};
    electorode_information = fieldtrip_raw.elec;
    
    new_chanpos = [];
    new_chantype = {};
    new_chanunit = {};
    new_elecpos = [];
    new_labels = {};
    
    cnt = 1;
    for idx = 1:numel(electorode_information.label)
        electrode =  electorode_information.label{idx};
        if ~ismember(electrode, to_remove)
            new_chanpos(cnt,:) = electorode_information.chanpos(idx,:);
            new_chantype{cnt} = electorode_information.chantype{idx};
            new_chanunit{cnt} = electorode_information.chanunit{idx};
            new_elecpos(cnt,:) = electorode_information.elecpos(idx,:);
            new_labels{cnt} = electorode_information.label{idx};
            cnt = cnt + 1;
        end  
    end
    
    
    fieldtrip_raw.elec.chanpos = new_chanpos;
    fieldtrip_raw.elec.chantype = new_chantype;
    fieldtrip_raw.elec.chanunit = new_chanunit;
    fieldtrip_raw.elec.elecpos = new_elecpos;
    fieldtrip_raw.elec.label = new_labels;
   
    cnt = 1;
    for idx = 1:numel(fieldtrip_raw.label)
        elec = fieldtrip_raw.label{idx};
        if ~ismember(elec, to_remove)
            new_elec{cnt} = elec;
            for i=1:length(fieldtrip_raw.trial)
                new_trial{1,i}(cnt, :) = fieldtrip_raw.trial{1,i}(idx,:);
            end
            cnt = cnt+1;
        end
    end
    
    fieldtrip_raw.label = new_elec;
    fieldtrip_raw.trial = new_trial;
end
