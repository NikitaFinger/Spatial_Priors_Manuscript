function newfname = extract_bat(date,batID,mic_path,train_or_test,cond,vicon_trial_dir,bat_pos_dir,plotter)
    %Modified by Nikita (Feb 2025) 
    % Modified by Shivam for prism use case on Mar 13
    % Modified by ChatGPT to only fill in missing tip positions (bat_pos{1})
    % for frames where all three markers (Tip, Left, and Right) are missing.
    %
    % For gaps (of up to 100 frames) where all markers are missing, the missing
    % tip positions are filled using cubic spline interpolation. Only the
    % extrapolated tip positions are overlaid in black while any originally
    % present tip marker positions remain plotted in their original color.
    
    if ~exist('plotter','var')
        plotter = true;
    end
    
    if ~exist('vicon_trial_dir','var')
        % Directory does not exist, so default it to something.
        vicon_trial_dir = fullfile(cd,strjoin({'../../Data/Data/Vicon_Data/',train_or_test,'/',batID,'/',date,'/trial_'},''));
    end
    
    [fname, pname] = uigetfile('*.c3d', [], vicon_trial_dir);
    
    if ~exist('bat_pos_dir','var')
        % Saving path does not exist, so default it to something.
        bat_pos_dir = fullfile(cd,'../../Data/Preprocessing/Vicon/Bat_Position/');
    end
    
    D = strsplit(pname, pname(end));
    dateindx = contains(D, date(1:4));
    datename = D{dateindx};
    
    trialnum = regexprep(fname, 'trial_', '');
    trialnum = regexprep(trialnum, '.c3d', '');
    
    newfname = [bat_pos_dir strjoin({batID, trialnum, datename, cond, 'bat_pos'}, '_') '.mat'];
    newfname = regexprep(newfname, '__', '_');
    
    [point_array, frame_rate, ~, ~, ~, ~] = lc3d([pname, fname]);
    
    if isempty(point_array)
        return;
    end
    point_names = cellfun(@(c) c.name, point_array, 'uniformoutput', 0);
    
    % ----- Extract markers (Tip, Left, Right) -----
    markers = {'Tip','Left','Right'};
    bat_pos = {};
    for mm = 1:length(markers)
        lab = ~cellfun(@isempty, strfind(point_names, markers{mm}));
        if ~isempty(find(lab, 1))
            bat_pos{mm} = point_array{lab}.traj ./ 1e3;
            bat_pos{mm}(bat_pos{mm} == 0) = nan;
        end
    end
    
    if isempty(bat_pos) || ~isempty(find(cellfun('isempty', bat_pos), 1))
        return;
    end
    
    % ----- Compute roll and head vector using frames where all markers are present -----
    frames_w_alldata = find( isfinite(bat_pos{1}(:,1)) & isfinite(bat_pos{2}(:,1)) & isfinite(bat_pos{3}(:,1)) );
    
    flA = [0 0 0];
    flB = [1 0 1];
    flC = [-1 0 1];
    nnfl = cross(flB-flA, flC-flA);
    
    nn      = nan(size(bat_pos{1}));
    th_roll = nan(size(bat_pos{1},1),1);
    head_vec= nan(size(bat_pos{1}));
    for fr = frames_w_alldata'
        % Three points on the triangle:
        A = bat_pos{1}(fr,:);
        B = bat_pos{2}(fr,:);
        C = bat_pos{3}(fr,:);
        
        nn(fr,:) = cross(C-A, B-A);
        th_roll(fr) = atan2(norm(cross(nn(fr,:), nnfl)), dot(nn(fr,:), nnfl));
        
        midptBC = (B + C) ./ 2;
        head_vec(fr,:) = A - midptBC;
    end
    
    % ----- Plot the markers and microphone locations -----
    if plotter
        cols = 'rgb';
        figure(1), clf; set(gcf, 'pos', [10 40 520 480]); hold on;
        % Plot each marker in its original color.
        for mm = 1:length(markers)
            scatter3(bat_pos{mm}(:,3), bat_pos{mm}(:,1), bat_pos{mm}(:,2), 15, [cols(mm) '.'], 'MarkerFaceAlpha', 0.5);
        end
        
        % Mark the start and end using the Tip marker.
        text(bat_pos{1}(frames_w_alldata(1),3), bat_pos{1}(frames_w_alldata(1),1), bat_pos{1}(frames_w_alldata(1),2), 'Start');
        text(bat_pos{1}(frames_w_alldata(end),3), bat_pos{1}(frames_w_alldata(end),1), bat_pos{1}(frames_w_alldata(end),2), 'End');
        
        load(mic_path, 'mic_loc');
        scatter3(mic_loc(:,3), mic_loc(:,1), mic_loc(:,2), 40, '.k');
        
        xlabel('Z'); ylabel('X'); zlabel('Y');
        xlim([0 inf]);
        view(45,-50);
        axis equal, grid on;
    end
    
    % ----- Extrapolate missing tip marker (bat_pos{1}) positions -----
    %
    % Identify frames where all three markers are missing.
    % Only fill in the tip position for these frames (if gap <= max_gap frames).
    % We also record which frames are filled in so we can overlay them in black.
    
    all_missing = isnan(bat_pos{1}(:,1)) & isnan(bat_pos{2}(:,1)) & isnan(bat_pos{3}(:,1));
    missing_frames = find(all_missing);
    max_gap = 150;  % maximum consecutive missing frames to fill
    filled_idx = false(size(bat_pos{1},1),1);  % logical index for extrapolated frames
    
    % Use tip marker (bat_pos{1}) for interpolation from valid frames.
    valid_idx = find(~isnan(bat_pos{1}(:,1)));
    
    if ~isempty(missing_frames)
        gap_start = missing_frames(1);
        gap_end   = gap_start;
        for i = 2:length(missing_frames)
            if missing_frames(i) == gap_end + 1
                gap_end = missing_frames(i);
            else
                if (gap_end - gap_start + 1) <= max_gap
                    % Find the nearest valid tip frame indices before and after the gap.
                    idx_before = find(valid_idx < gap_start, 1, 'last');
                    idx_after  = find(valid_idx > gap_end, 1, 'first');
                    if isempty(idx_before) && isempty(idx_after)
                        % No valid neighbors: cannot interpolate.
                    elseif isempty(idx_before)
                        % Only neighbor after the gap available.
                        bat_pos{1}(gap_start:gap_end, :) = repmat(bat_pos{1}(valid_idx(idx_after), :), gap_end-gap_start+1, 1);
                        filled_idx(gap_start:gap_end) = true;
                    elseif isempty(idx_after)
                        % Only neighbor before the gap available.
                        bat_pos{1}(gap_start:gap_end, :) = repmat(bat_pos{1}(valid_idx(idx_before), :), gap_end-gap_start+1, 1);
                        filled_idx(gap_start:gap_end) = true;
                    else
                        % Both neighbors exist: use cubic spline interpolation.
                        t_valid = [valid_idx(idx_before), valid_idx(idx_after)];
                        t_gap   = gap_start:gap_end;
                        for dim = 1:3
                            bat_pos{1}(gap_start:gap_end, dim) = interp1(t_valid, bat_pos{1}(t_valid, dim), t_gap, 'spline');
                        end
                        filled_idx(gap_start:gap_end) = true;
                    end
                end
                % Start a new gap.
                gap_start = missing_frames(i);
                gap_end   = gap_start;
            end
        end
        % Process the final gap.
        if (gap_end - gap_start + 1) <= max_gap
            idx_before = find(valid_idx < gap_start, 1, 'last');
            idx_after  = find(valid_idx > gap_end, 1, 'first');
            if isempty(idx_before) && isempty(idx_after)
                % No valid neighbors.
            elseif isempty(idx_before)
                bat_pos{1}(gap_start:gap_end, :) = repmat(bat_pos{1}(valid_idx(idx_after), :), gap_end-gap_start+1, 1);
                filled_idx(gap_start:gap_end) = true;
            elseif isempty(idx_after)
                bat_pos{1}(gap_start:gap_end, :) = repmat(bat_pos{1}(valid_idx(idx_before), :), gap_end-gap_start+1, 1);
                filled_idx(gap_start:gap_end) = true;
            else
                t_valid = [valid_idx(idx_before), valid_idx(idx_after)];
                t_gap   = gap_start:gap_end;
                for dim = 1:3
                    bat_pos{1}(gap_start:gap_end, dim) = interp1(t_valid, bat_pos{1}(t_valid, dim), t_gap, 'spline');
                end
                filled_idx(gap_start:gap_end) = true;
            end
        end
    end
    
    % ----- Overlay only the extrapolated tip positions in black -----
    if plotter && any(filled_idx)
        scatter3(bat_pos{1}(filled_idx,3), bat_pos{1}(filled_idx,1), bat_pos{1}(filled_idx,2), ...
            20, 'k', 'filled');
    end
    
    % ----- Save the results -----
    if ~isempty(bat_pos)
        if exist(newfname, 'file')
            b = questdlg(['File ' newfname ' already exists. Overwrite?'], 'Overwrite?', 'Yes', 'No', 'No');
            if strcmp(b, 'No')
                return;
            end
        end
        save(newfname, 'bat_pos', 'markers', 'nn', 'th_roll', 'head_vec', 'frame_rate');
        disp(['Saved ' newfname]);
    end
end
