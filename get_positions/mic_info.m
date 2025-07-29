%% Rousettus Spatial Navigation in Big Flight Room
% use to get the mic location and mic vector,also fill in the 32 matrix
% with NaN
% Xiaoyan Yin Oct-20,2022
% Shivam modified on Mar 9

function savepath = mic_info(date_data,trl_num, mic_pos_folder, save_folder, mic_info_file, n_mics, n_markers)

%SETTING DEFAULT PARAMETERS

    if ~exist('landmark_folder','var')
         % parameter does not exist, so default it to something
          mic_pos_folder = fullfile(cd,'../../Data/Preprocessing/Vicon/Mic_Position/');
    end

    if ~exist('save_folder','var')
         % parameter does not exist, so default it to something
          save_folder = fullfile(cd,'../../Data/Preprocessing/Vicon/Mic_Position/');
    end

    if ~exist('mic_info_file','var')
         % parameter does not exist, so default it to something
        %get all the files in the mic_info folder

        files = dir(fullfile(cd,'mic_info_files', '*.mat'));
        mic_info_dates = zeros(size(files));

        for k = 1:numel(files)
            mic_info_dates(k) = str2double((files(k).name(end-11:end-4)));
        end

        diff = mic_info_dates - str2double(date_data); diff(diff>0) = Inf;
        [~, i] = min(abs(diff));

        mic_info_file = strjoin(['mic_info_files/mic_info_',string(mic_info_dates(i)),'.mat'],'');

    end
    
    if ~exist('n_mics','var')
        n_mics = 32;
        disp("number of mics: ")
        disp(n_mics)
    end

    if ~exist('n_markers','var')
        n_markers = 3;
        disp("number of markers: ")
        disp(n_markers)
    end

    n_channels = 32;
    gain = 10;

    load([mic_pos_folder 'mic_pos_' date_data '_' trl_num '.mat'],'mic_pos') %mic_pos, extract_mics.m
    
    tip_indx = n_markers*(1:n_mics);
    vec_indx = [tip_indx; tip_indx-1];
    
    tips1 = mic_pos(tip_indx,:);
    mic_vec1 = mic_pos(vec_indx(1,:),:) - mic_pos(vec_indx(2,:),:);
    
    mic_vec = mic_vec1;
    mic_loc = tips1;
    mic_gain = gain.*ones(n_channels,1);
    mic_vh=[];
    
    figure(3), clf; set(gcf,'pos',[10 300 912 700]), hold on;
    scatter3(tips1(:,3),tips1(:,1),tips1(:,2),40,'.k')
    hold on
    
    %fill the recording channels into the 32 matrix 
    %we dont use all 32 chs, so we have match the 25 chs we use to the 32
    %matrix (2_0,2_2,.....5_7)
    %use the info in mic_info.mat file to do this
    disp("Mic IDs extracted from: " + mic_info_file)
    load(fullfile(cd, mic_info_file),'mic_ids','channels','mic_labels')

    %remove letters and underscores
    mic_labels(1:end-2) = regexprep(regexprep(mic_labels(1:end-2),'[a-zA-Z\s]',''),'_','');

    labelled_mics = round(cellfun(@str2num, mic_labels(1:3:n_mics*n_markers))/10);
    mic_idx = interp1(mic_ids(:,1),mic_ids(:,2),labelled_mics);
    
    %fill in nans where appropriate
    mic_loc = mic_nan_convert(mic_loc,mic_idx,n_channels);

    mic_vec = mic_nan_convert(mic_vec,mic_idx,n_channels);

    text(mic_loc(:,3),mic_loc(:,1),mic_loc(:,2),num2str((1:length(mic_loc(:,1)))'),'color','magenta')
    text(mic_loc(:,3),mic_loc(:,1),mic_loc(:,2),"    ("+channels+")",'color','black')
    view(75,-20)
    axis equal, grid on

    savepath = [save_folder 'mic_info_' date_data '_' trl_num '.mat'];
    
    save(savepath,...
      'mic_loc','mic_vec','mic_gain','mic_vh');

end


