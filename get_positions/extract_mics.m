function savepath = extract_mics(date_data,trl_num,no_calib, landmark_folder, save_folder, mic_info_file)
%modified by S.C.and GhostRider (NF) with help from JS
    
    %SETTING DEFAULT PARAMETERS
    if ~exist('no_calib','var')
         % parameter does not exist, so default it to something
          no_calib = 1;
    end
% 
%     if ~exist('landmark_folder','var')
%          % parameter does not exist, so default it to something
%           landmark_folder = fullfile(cd,'../../Data/Data/Vicon_Data/landmark/');
%        end
    if ~exist('landmark_folder','var')
        raw_path = fullfile(cd,'../../Data/Data/Vicon_Data/landmark/');
        landmark_folder = char(java.io.File(raw_path).getCanonicalPath());
    end



   

    if ~exist('save_folder','var')
         % parameter does not exist, so default it to something
          save_folder = fullfile(cd,'../../Data/Preprocessing/Vicon/Mic_Position/');
    end
%%
    if ~exist('mic_info_file','var')
         % parameter does not exist, so default it to something
        %get all the files in the mic_info folder

        files = dir(fullfile(cd,'mic_info_files', '*.mat'));
        mic_info_dates = zeros(size(files));

        for k = 1:numel(files)
            mic_info_dates(k) = str2double((files(k).name(end-11:end-4)));
        end

        diff = mic_info_dates - str2double(date_data); diff(diff>0) = Inf;
        [~, idx] = min(abs(diff));

%         mic_info_file = strjoin(['mic_info_files/mic_info_',string(mic_info_dates(idx)),'.mat'],'');
mic_info_file = ['mic_info_files/mic_info_' num2str(mic_info_dates(idx)) '.mat'];
%%
    end
    
    %work in working data
%     fn=[landmark_folder date_data '/' trl_num(1:end-2) '/' trl_num '.c3d'];
%    fn = fullfile(landmark_folder, date_data, [trl_num, '.c3d']);
    fn = fullfile(landmark_folder, date_data, 'landmark', [trl_num, '.c3d']);

    [point_array, ~, ~, ~, ~, ~] = lc3d(fn);
    point_names=cellfun(@(c) c.name,point_array,'uniformoutput',0);
    disp("Loading C3D from: " + fn)

    
    disp("Vicon label names extracted from: " + mic_info_file)
    load(fullfile(cd,mic_info_file),'mic_labels');
    
    mic_pos=nan(length(mic_labels),3);
    for mm=1:length(mic_labels)
      imic=find(strcmp(point_names,mic_labels{mm}));
      if ~isempty(imic)
        pmic=point_array{imic}.traj./1e3;
        indx=find(pmic(:,1),1);
        mic_pos(mm,:)=pmic(indx,:);
      end
    end
    
    figure(2); clf; set(gcf,'pos',[10 300 912 700]), hold on;
    cols='rgb';
    if no_calib
      for mm=1:3
        scatter3(mic_pos(mm:3:end,3),mic_pos(mm:3:end,1),mic_pos(mm:3:end,2),...
          cols(mm))
      end
    else
      for mm=1:3
        scatter3(mic_pos(mm:3:end-4,3),mic_pos(mm:3:end-4,1),mic_pos(mm:3:end-4,2),...
          cols(mm))
      end
      scatter3(mic_pos(end-3:end,3),mic_pos(end-3:end,1),mic_pos(end-3:end,2),'k')
      tip_indx=[(1:3:32*3)+2 32*3+2 32*3+4];
      tips=mic_pos(tip_indx,:);
      text(tips(:,3),tips(:,1),tips(:,2),num2str((1:length(tips(:,1)))'))
    end
    view(3)
    axis equal, grid on  
    view(75,-20)
    
    savepath = [save_folder 'mic_pos_' date_data '_' trl_num '.mat'];
    
    save(savepath,'mic_pos')

end

