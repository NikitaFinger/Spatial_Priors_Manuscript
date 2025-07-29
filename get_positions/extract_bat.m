function newfname = extract_bat(date,batID,mic_path,train_or_test,cond,vicon_trial_dir,bat_pos_dir,plotter)

    %Modified by Shivam for prism use case on Mar 13
    if ~exist('plotter','var')
        plotter = true;
    end
    
    if ~exist('vicon_trial_dir','var')
        % % directory does not exist, so default it to something
        vicon_trial_dir = fullfile(cd,strjoin({'../../Data/Data/Vicon_Data/',train_or_test,'/',batID,'/',date,'/trial_'},''));
    end
    
    [fname,pname]=uigetfile('*.c3d',[],vicon_trial_dir);
    
    if ~exist('landmark_folder','var')
        % % saving path does not exist, so default it to something
        bat_pos_dir = fullfile(cd,'../../Data/Preprocessing/Vicon/Bat_Position/');
    end
    
    D=strsplit(pname,pname(end));
    dateindx = contains(D,date(1:4));
    datename=D{dateindx};
    
    trialnum=regexprep(fname,'trial_','');
    trialnum=regexprep(trialnum,'.c3d','');
    
    newfname=[bat_pos_dir strjoin({batID,trialnum,datename,cond,'bat_pos'},'_') '.mat'];
    newfname=regexprep(newfname,'__','_');
    
    [point_array, frame_rate, ~, ~, ~, ~] = lc3d([pname,fname]);
    
    if isempty(point_array)
      return;
    end
    point_names=cellfun(@(c) c.name,point_array,'uniformoutput',0);
    
    %extracting
    markers={'Tip','Left','Right'};
    
    bat_pos={};
    for mm=1:length(markers)
      lab=~cellfun(@isempty,strfind(point_names,markers{mm}));
      if ~isempty(find(lab, 1))
        bat_pos{mm}=point_array{lab}.traj./1e3;
        bat_pos{mm}(bat_pos{mm}==0)=nan;
      end
    end
%     
%% 
% %% Added NF -02240205 interpolate position if xyz coordinates missing use if missing coordinates
% 
%    if isempty(bat_pos) || ~isempty(find(cellfun('isempty', bat_pos), 1))
%         return
%     end
% 
%     % Extrapolate missing XYZ coordinates using linear interpolation
%     for dim = 1:3 % Loop over X, Y, and Z dimensions
%         nan_indices = find(isnan(bat_pos{dim})); % Find NaN indices in the current dimension
%         if ~isempty(nan_indices)
%             for j = 1:length(nan_indices)
%                 nan_idx = nan_indices(j);
%                 before_nan_idx = find(isfinite(bat_pos{dim}(1:nan_idx)), 1, 'last');
%                 after_nan_idx = find(isfinite(bat_pos{dim}(nan_idx:end)), 1);
%                 if ~isempty(before_nan_idx) && ~isempty(after_nan_idx)
%                     % Linear interpolation
%                     interpolated_value = interp1([before_nan_idx, after_nan_idx + nan_idx - 1], ...
%                         [bat_pos{dim}(before_nan_idx), bat_pos{dim}(after_nan_idx + nan_idx - 1)], nan_idx);
%                     bat_pos{dim}(nan_idx) = interpolated_value;
%                 end
%             end
%         end
%     end
    %% 
    
    %creating a plane
    frames_w_alldata=find(isfinite(bat_pos{1}(:,1)) & isfinite(bat_pos{2}(:,1))...
      & isfinite(bat_pos{3}(:,1)));
    
    flA=[0 0 0];
    flB=[1 0 1];
    flC=[-1 0 1];
    nnfl=cross(flB-flA,flC-flA);
    
    nn=nan(size(bat_pos{1}));
    th_roll=nan(size(bat_pos{1},1),1);
    head_vec=nan(size(bat_pos{1}));
    for fr=frames_w_alldata'
      %three points on the triangle
      A=bat_pos{1}(fr,:);
      B=bat_pos{2}(fr,:);
      C=bat_pos{3}(fr,:);
      
      %get roll vector from plane
      nn(fr,:)=cross(C-A,B-A);
      th_roll(fr)=atan2(norm(cross(nn(fr,:),nnfl)),dot(nn(fr,:),nnfl));
      
      %direction vector
      midptBC=(B+C)./2;
      head_vec(fr,:)=A-midptBC;
    end
    
    %plotting
    if plotter
      cols='rgb';
      
      figure(1), clf; set(gcf,'pos',[10 40 520 480]); hold on;
      for mm=1:length(markers)
        scatter3(bat_pos{mm}(:,3),bat_pos{mm}(:,1),bat_pos{mm}(:,2),15,[cols(mm) '.'],'MarkerFaceAlpha',0.5)
      end
      
      text(bat_pos{mm}(frames_w_alldata(1),3),...
        bat_pos{mm}(frames_w_alldata(1),1),...
        bat_pos{mm}(frames_w_alldata(1),2),'Start');
      text(bat_pos{mm}(frames_w_alldata(end),3),...
        bat_pos{mm}(frames_w_alldata(end),1),...
        bat_pos{mm}(frames_w_alldata(end),2),'End');
        
      load(mic_path,'mic_loc')
      scatter3(mic_loc(:,3),mic_loc(:,1),mic_loc(:,2),40,'.k');

      xlabel('Z');
      ylabel('X');
      zlabel('Y');

      xlim([0 inf]);

      view(45,-50)
      axis equal, grid on
    end
    
    if ~isempty(bat_pos)
      if exist(newfname,'file')
        b=questdlg(['File ' newfname ' already exists. Overwrite?'], 'Overwrite?', ...
          'Yes','No','No');
        switch b
          case 'No'
            return
        end
      end
      save(newfname,'bat_pos','markers',...
        'nn','th_roll','head_vec','frame_rate')
      disp(['Saved ' newfname])
    end
end