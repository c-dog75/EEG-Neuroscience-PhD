n_participants = 40;
addpath('/Users/cihandogan/Documents/Research/spm12')
main_path = '/Users/cihandogan/Documents/Research/preprocessing/after_spm_script/participant_';
addpath('/Users/cihandogan/Documents/Research/fieldtrip-20230118');
ft_defaults;
clear matlabbatch
ft_defaults;

spm eeg;

for participant =  [27, 36,]%1:n_participants
    %% get the correct file path
    clear matlabbatch
    disp(strcat('Procesisng participant..',int2str(participant)));
    participant_main_path = strcat(main_path, int2str(participant));    
    
    
    
    if exist(participant_main_path, 'dir')
        if participant < 10
            p = strcat('P0', int2str(participant));
        else
            p = 'P';
            p = strcat(p,int2str(participant));
        end
        
        cd(participant_main_path);
        
        p = strcat(p, '-Deci_ready_for_ft.set');
        participant_file_path = strcat(participant_main_path, '/');
        participant_file_path = strcat(participant_file_path, p);
        
        %% convert to SPM ready file

        if isfile(participant_file_path)
            D = spm_eeg_convert(participant_file_path);
        else
            continue;
        end


        
        %% clean up the labels
        label = D.chanlabels;
        for i = 1:numel(label)
            lbl = label{i};
            lbl(lbl=='"')=[];
            k = strfind(lbl, '(');
            if ~isempty(k)
                lbl = lbl(1:(k-1));
            end
            label{i} = lbl;
        end
        
        D = chanlabels(D, ':', label);
        
        if isfield(D, 'origchantypes')
            D=D.rmfield('origchantypes');
        end
        
        
        S = [];
        S.D =D;
        S.task = 'defaulttype';
        S.save = 1;
        D = spm_eeg_prep(S);
        
        S = [];
        S.D = D;
        S.task = 'defaulteegsens';
        S.save = 1;
        D = spm_eeg_prep(S);
        
        S = [];
        S.task = 'project3D';
        S.modality = 'EEG';
        S.updatehistory = 0;
        S.D = D;
        
        D = spm_eeg_prep(S);
        
        D = timeonset(D, -0.5); % time onset flag is wrong
        
        save(D);
        disp('COMPLETED PARTICIPANT');
        disp(participant);
        
    end
end