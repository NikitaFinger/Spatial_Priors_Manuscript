function createBeamPatternVideo(freq_desired, side_view, save_movie)
    addpath 'Z:\Nikita\Prism_Experiment\Analysis\beampattern_analysis\beampattern_preprocessing'
    
    a = [];
    vars = who;
    clear_vars = setdiff(vars, {'freq_desired', 'side_view', 'save_movie', 'bat_type'});
    clear(clear_vars{:});

    if ~exist('save_movie', 'var')
        save_movie = 1;
    end
    vid_frate = 12; %was at 12

    if ~exist('side_view', 'var')
        side_view = 0; %if not side view, then top view
    end

    if ~exist('freq_desired', 'var')
        freq_desired = 35; %khz
    end

    side_view = 0;
    save_movie = 1;

    mic_proc_dir = 'Z:\Nikita\Prism_Experiment\Data\Preprocessing\Mic\Beam_Output\';

    % Prompt user to select a trial file
    [filename, path] = uigetfile([mic_proc_dir '*.mat'], 'Select Trial File');
    
    if filename == 0
        % User canceled file selection, exit function
        return;
    end
    
    trials = dir(fullfile(path, filename));

    % Specify the output directory for the video
    outputDir = 'Z:\Nikita\Prism_Experiment\Data\Preprocessing\Mic\Beam_AV\';

    mic_data_dir = 'Z:\Nikita\Prism_Experiment\Data\Preprocessing\Mic\Mic_Data_Detect\';

    use_interp_beamshape = 1;
    interp_method = 'rb_natural';

    d2 = 1; % 2d data only
    if d2
        gaussfitfcn = 1;
    end

    frames_limit_hard = 1;

    plot_gaussian_beampattern = 1;
    plot_mics = 1;
    plot_mic_nums = 1;
    flatten_deg = 30;
    diag = 0;

    for tt=1  %:length(trials)
        bpp = load([mic_proc_dir trials(tt).name]);
        if exist('checked', 'var') && checked 
            if exist([mic_proc_dir trials(tt).name(1:end-4) '_checked.mat'], 'file')
                bpp_checked = load([mic_proc_dir trials(tt).name(1:end-4) '_checked.mat']);
                bpp.proc.chk_good_call = bpp_checked.proc.chk_good_call;
                bpp.proc.ch_ex = bpp_checked.proc.ch_ex;
            else
                continue
            end
        end
        bat = bpp.track.track_smooth;

        bpp.proc.call_psd_dB_comp_re20uPa_withbp(cellfun(@isempty, bpp.proc.call_psd_dB_comp_re20uPa_withbp)) = {NaN(1,63)};

        if save_movie
            out_fname = trials(tt).name(1:end-21);
            if side_view
                out_fname = [out_fname '_side'];
            else
                out_fname = [out_fname '_top'];
            end
            if use_interp_beamshape
                out_fname = [out_fname '_interp'];
            else
                out_fname = [out_fname '_raw'];
            end
            out_fname = [out_fname '_' num2str(freq_desired)];
            v = VideoWriter([outputDir out_fname '.mp4'], 'MPEG-4');
            
            v.FrameRate = vid_frate;
            open(v);
        end

        calls_with_track = bpp.mic_data.call_idx_w_track;
        call_track_locs = round([bpp.mic_data.call.call_start_idx] ...
            / bpp.mic_data.fs * bpp.track.fs);

        goodch = 1:length(bpp.mic_loc(:,1)); %alternatively load in the good ch from the ch_ex var
        mic_num = 1:bpp.mic_data.num_ch_in_file;
        
        close all;
        figure(1); set(gcf, 'pos', [10 40 900 900], 'color', 'w')
        set(gca, 'position', [.07 .06 .88 .89], 'units', 'normalized')
        hold on;

        if frames_limit_hard
            frames = call_track_locs(calls_with_track(find(bpp.proc.chk_good_call, 1))): ...
                call_track_locs(calls_with_track(find(bpp.proc.chk_good_call, 1, 'last')));
            
            % Limit to speed > 1.5 m/s
            speed = calc_speed(bat, bpp.track.fs, 0);
            speed = [nan; speed];
            good_speed = find(speed > 0.1);
            good_speed(good_speed < frames(1)) = [];
            good_speed(good_speed > frames(end)) = [];

            frames = max(good_speed(1), frames(1)):min(good_speed(end), frames(end));
        else
            frames = find(isfinite(bat(:,1)));
            % Limiting the frames plotted in animation to around the time of the first
            % and last good calls
            frames = max(frames(1), ...
                call_track_locs(calls_with_track(find(bpp.proc.chk_good_call, 1)))-10): ...
                min(call_track_locs(calls_with_track(find(bpp.proc.chk_good_call, 1, 'last')))+10, ...
                frames(end));
        end

        % Initial plot setup
        pb = plot3(bat(frames(1), 1), bat(frames(1), 2), bat(frames(1), 3), '-k', 'linewidth', 2);

        for fr = frames
            % Update the plot for the current frame
            set(pb, 'XData', bat(frames(1):fr, 1), 'YData', bat(frames(1):fr, 2), 'ZData', bat(frames(1):fr, 3));

            % Update axis limits based on the data in the current frame
            currentMax = max(bat(frames(1):fr, :), [], 1) + 0.1;  % Adding a small buffer
            currentMin = min(bat(frames(1):fr, :), [], 1) - 0.1;
            axis([currentMin(1) currentMax(1) currentMin(2) currentMax(2) currentMin(3) currentMax(3)]);
            
            % Maintain consistent aspect ratio
            axis equal;

            if save_movie
                writeVideo(v, getframe(gcf));
            end
            drawnow;
        end

        if save_movie
            close(v);  % Close the video writer object
        end
    end
end
