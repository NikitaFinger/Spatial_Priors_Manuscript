%This is a function modified from code by Ben Falk and Wu Jung Lee to create beam aim
%animations by Nikita Finger for use in the Beam Pattern Pipeline
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

    frames_limit_hard = 0;

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
            v = VideoWriter(['Z:\Nikita\Prism_Experiment\Data\Preprocessing\Mic\Beam_AV\' out_fname '.mp4'], 'MPEG-4');
            
            v.FrameRate = vid_frate;
            open(v);
        end

        calls_with_track = bpp.mic_data.call_idx_w_track;
        call_track_locs = round([bpp.mic_data.call.call_start_idx] ...
            / bpp.mic_data.fs * bpp.track.fs);

        goodch = 1:length(bpp.mic_loc(:,1)); %alternatively load in the good ch from the ch_ex var
        mic_num = 1:bpp.mic_data.num_ch_in_file;
        
        % Add locations of the  , 15 L
%        bpp.mic_loc(2,:) = [-0.05505141448974609 	2.323274658203125	1.2506273193359374]; % Original
%         bpp.mic_loc(3,:) = [-0.18473699951171876	2.323274658203125	1.2506273193359374]; % Moved perch location

%         Locations for 30 L
%        bpp.mic_loc(2,:) = [-0.05505141448974609 2.35066943359375 1.2496038818359374]; % Original
%        bpp.mic_loc(3,:) = [-0.44,	2.35066943359375,	1.24859375000000]; % Moved perch location

%       original
        
%         bpp.mic_loc(2,:) = [-0.05505141448974609 2.35066943359375 1.2496038818359374];

        close all;
        figure(1); set(gcf, 'pos', [10 40 900 900], 'color', 'w')
        set(gca, 'position', [.07 .06 .88 .89], 'units', 'normalized')
        hold on;
        % Commenting out the microphone plotting
        % pb_mics = plot3(bpp.mic_loc(goodch, 1), bpp.mic_loc(goodch, 2), bpp.mic_loc(goodch, 3), ...
        %     'ok', 'markerfacecolor', 'k');

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

        pb = plot3(bat(frames(1):frames(end), 1), ...
            bat(frames(1):frames(end), 2), ...
            bat(frames(1):frames(end), 3), '-k', 'linewidth', 2);

        bpp.mic_loc(3,:) = [-0.18473699951171876	2.323274658203125	1.2506273193359374];
        % Plot the new and moved perch locations
%         plot3(bpp.mic_loc(2, 1), bpp.mic_loc(2, 2), bpp.mic_loc(2, 3), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none', 'LineStyle', '--', 'MarkerSize', 12); % New perch (black dashed outline) % Original (red, dashed)
        plot3(bpp.mic_loc(3, 1), bpp.mic_loc(3, 2), bpp.mic_loc(3, 3), 'og', 'MarkerSize', 12, 'MarkerFaceColor', 'g'); % Moved perch
         % Add locations of the  , 15 L
    %    bpp.mic_loc(2,:) = [-0.05505141448974609 	2.323274658203125	1.2506273193359374]; % Original
%         bpp.mic_loc(3,:) = [-0.18473699951171876	2.323274658203125	1.2506273193359374]; % Moved perch location

%         Locations for 30 L
%        bpp.mic_loc(2,:) = [-0.05505141448974609 2.35066943359375 1.2496038818359374]; % Original
%        bpp.mic_loc(3,:) = [-0.44,	2.35066943359375,	1.24859375000000]; % Moved perch location

%       original
        
        if side_view
            view(0, 0)
        else
            view(2)
        end
        
% %         Grid 
%         axis equal, grid on;
%         aa = axis;
%         set(gca, 'fontsize', 22, ...
%             'XTick', aa(1):.4:aa(2), 'xticklabel', (aa(1):.4:aa(2))-aa(1), ...
%             'YTick', aa(3):.4:aa(4), 'yticklabel', (aa(3):.4:aa(4))-aa(3))
%         if side_view
%             title_name = 'Side View';
%         else
%             title_name = 'Top View';
%         end
%         title(title_name, 'fontsize', 24)
%         box on;
%         aa = axis;

% Grid and title setting - no grid or title 
axis equal;
aa = axis;
% set(gca, 'fontsize', 22, ...
%     'XTick', aa(1):.4:aa(2), 'xticklabel', (aa(1):.4:aa(2))-aa(1), ...
%     'YTick', aa(3):.4:aa(4), 'yticklabel', (aa(3):.4:aa(4))-aa(3))
box on;
aa = axis;


        delete(pb);  
        % Commenting out the microphone deletion
        % if ~plot_mics
        %     delete(pb_mics);
        % end

        call_fr = 0; pm = []; pmic = []; BX = []; BY = [];
        for fr = frames(1):frames(end)
            pb = plot3(bat(frames(1):fr, 1), bat(frames(1):fr, 2), bat(frames(1):fr, 3), ...
                '-', 'linewidth', 2, 'color', [.2 .2 .2]);

            % If a call is present on the frame and selected as 'good call' we plot the beam pattern
            good_call = find(bpp.proc.chk_good_call);
            call_present = ismember(call_track_locs(calls_with_track(good_call)), fr);

            call_indx_1 = find(call_present);
            call_indx = good_call(call_indx_1);
            if ~isempty(call_indx) && ismember(fr, frames)
                delete(pm); pm = [];

                if ~use_interp_beamshape
                    delete(pmic); pmic = [];
                end

                call_fr = fr;
                call_dB = nan(1, bpp.mic_data.num_ch_in_file);
                for iM = mic_num
                    freq = bpp.proc.call_freq_vec{call_indx, iM};
                    [~, fidx] = min(abs(freq - freq_desired * 1e3));
                    if isempty(bpp.proc.call_psd_dB_comp_re20uPa_withbp{call_indx, iM})
                        call_dB(iM) = 0;
                    else
                        call_dB(iM) = bpp.proc.call_psd_dB_comp_re20uPa_withbp{call_indx, iM}(fidx);
                    end
                end

                ch_ex_sig = find(isnan(call_dB)); % Low quality channel from call extraction function
                ch_good_loc = ~isnan(bpp.mic_loc(:,1))';  % Check for mics without location data

                angle_notnanidx = ~ismember(mic_num, ch_ex_sig) & ch_good_loc;
                goodch1 = goodch(angle_notnanidx);  % Xiaoyan

                mic_vec = bpp.mic_loc(goodch1, :) - ...
                    repmat(bat(fr, :), size(bpp.mic_loc(goodch1, :), 1), 1); % Vectors from bat to all mics pointing to bat
                [az, el] = ...
                    cart2sph(mic_vec(:, 1), mic_vec(:, 2), mic_vec(:, 3)); % Transform Cartesian to spherical coordinates

                % Doing the interpolation
                if use_interp_beamshape
                    [azq, elq] = meshgrid(min(az):pi/180:max(az), min(el):pi/180:max(el));
                    if strcmp(interp_method, 'rb_natural')  % Natural neighbor interpolation
                        vq = griddata(az, el, call_dB(angle_notnanidx), azq, elq, 'natural');
                    elseif strcmp(interp_method, 'rb_rbf')  % Radial basis function interpolation
                        vq = rbfinterp([azq(:)'; elq(:)'], rbfcreate([az(:)'; el(:)'], call_dB(angle_notnanidx), 'RBFFunction', 'multiquadrics'));
                        vq = reshape(vq, size(azq));
                    end
                    maxref = max(call_dB(angle_notnanidx));
                    vq_norm = vq - maxref;

                    % Find indices within measured polygon
                    k = boundary(az, el, 0);  % Outer boundary of all measured points
                    [in, on] = inpolygon(azq, elq, az(k), el(k));
                    in_smpl_poly = in | on;
                    clear in on
                    vq(~in_smpl_poly) = NaN;
                    vq_norm(~in_smpl_poly) = NaN;

                    [~, II] = max(vq_norm(:));
                    [I, J] = ind2sub(size(vq_norm), II);

                    elqm = elq;
                    elqm(isnan(vq_norm)) = NaN;
                    azqm = azq;
                    azqm(isnan(vq_norm)) = NaN;
                end

                % Plotting the beam direction Red =[228 26 28] ./ 255;,
                % blue=[0, 109, 219] ./ 255; green=[34, 139, 34] ./ 255;
                if bpp.proc.chk_good_call(call_indx)
                    beam_col = [0, 109, 219] ./ 255;
                else
                    beam_col = [0, 109, 219] ./ 255;
                end

                if d2
                    if use_interp_beamshape
                        if side_view
                            indx_flat = find(azqm <= 1/180*pi & azqm >= -1/180*pi); % Pull out indexes within 1 degree el in the interpolated data
                            angle_midline = elqm(indx_flat);
                        else
                            indx_flat = find(elqm <= 30/180*pi & elqm >= -30/180*pi); % Pull out indexes within 1 degree el in the interpolated data
                            angle_midline = azqm(indx_flat);
                        end
                        vq2d = vq_norm(indx_flat);
                        vq2d = vq2d - min(vq2d);
                        norm_vq2d = vq2d ./ max(vq2d);

                        [~, imax2d] = max(norm_vq2d);

                        [sigma, ~] = gaussfit(angle_midline, vq2d');
                        thdirfit = angle_midline(imax2d);

                        [bx, by] = pol2cart(thdirfit, .3);  % Transform polar to Cartesian coordinates.
                    else
                        if side_view
                            indx_flat = find(az <= flatten_deg/180*pi & az >= -flatten_deg/180*pi); % Pull out indexes within 1 degree el in the interpolated data
                            angle_midline = el(indx_flat);
                        else
                            indx_flat = find(el <= flatten_deg/180*pi & el >= -flatten_deg/180*pi); % Pull out indexes within 1 degree el in the interpolated data
                            angle_midline = az(indx_flat);
                        end
                        vq_norm = call_dB - max(call_dB);
                        vq2d = vq_norm(indx_flat);
                        vq2d = vq2d - min(vq2d);
                        norm_vq2d = vq2d' ./ max(abs(vq2d'));

                        if numel(angle_midline) >= 5
                            if gaussfitfcn
                                [sigma, thdirfit] = gaussfit(angle_midline, vq2d');
                            else
                                [thdirfit, sigma] = normfit(angle_midline, vq2d);
                            end
                        else
                            thdirfit = nan;
                        end

                        [bx, by] = pol2cart(thdirfit, .3);
                    end

                    if side_view
                        plot3([bat(fr, 1) bat(fr, 1) + bx], ...
                            [bat(fr, 2) bat(fr, 2)], ...
                            [bat(fr, 3) bat(fr, 3) + by], ...
                            '-', 'linewidth', 4, 'color', beam_col);
                    else
                        plot([bat(fr, 1) bat(fr, 1) + bx], [bat(fr, 2) bat(fr, 2) + by], ...
                            '-', 'linewidth', 4, 'color', beam_col);
                    end
                    BX = [BX; bx]; BY = [BY; by];
                else
                    [bx, by, bz] = sph2cart(azqm(I, J), elqm(I, J), .2);
                    plot3([bat(fr, 1) bat(fr, 1) + bx], [bat(fr, 2) bat(fr, 2) + by], ...
                        [bat(fr, 3) bat(fr, 3) + bz], ...
                        '-', 'linewidth', 2, 'color', beam_col);
                end

                % Plotting the beam pattern
                pm = [];
                if d2
                    if plot_gaussian_beampattern
                        % Create gaussian
                        thvals = -pi: 4/180*pi :pi;
                        yvals = 1/(sqrt(2*pi)*sigma) * exp( -(thvals - thdirfit).^2 / (2*sigma^2));

                        [Bx, By] = pol2cart(thvals, yvals ./ max(yvals) * .35); % Transform polar to Cartesian coordinates.
                        symbol = '';
                        % Plot
                        if side_view
                            pm = plot3(repmat(bat(fr, 1), size(Bx)) + Bx, ...
                                repmat(bat(fr, 2), size(By)) + zeros(size(By)), ...
                                repmat(bat(fr, 3), size(By)) + By, ...
                                ['-' symbol], 'linewidth', 2, 'color', [.4 .4 .4]);
                        else
                            pm = plot(repmat(bat(fr, 1), size(Bx)) + Bx, ...
                                repmat(bat(fr, 2), size(By)) + By, ...
                                ['-' symbol], 'linewidth', 2, 'color', [.4 .4 .4]);
                        end
                    else
                        [~, isort] = sort(angle_midline);
                        [Bx, By] = pol2cart(angle_midline(isort), norm_vq2d(isort) * .3);
                        if use_interp_beamshape
                            symbol = '';
                        else
                            symbol = '+';
                        end
                        if side_view
                            pm = plot3(repmat(bat(fr, 1), size(Bx)) + Bx, ...
                                repmat(bat(fr, 2), size(By)) + zeros(size(By)), ...
                                repmat(bat(fr, 3), size(By)) + By, ...
                                ['-' symbol], 'linewidth', 1, 'color', [.4 .4 .4]);
                        else
                            pm = plot(repmat(bat(fr, 1), size(Bx)) + Bx, repmat(bat(fr, 2), size(By)) + By, ...
                                ['-' symbol], 'linewidth', 1, 'color', [.4 .4 .4]);
                        end
                    end
                else
                    mic_vec = bpp.mic_loc(goodch, :) - ...
                        repmat(bat(fr, :), size(bpp.mic_loc(goodch, :), 1), 1);
                    norm_mic_vec = mic_vec ./ ...
                        repmat(sqrt((mic_vec(:, 1).^2 + mic_vec(:, 2).^2 + mic_vec(:, 3).^2)), 1, 3);
                    I_mic_vec = norm_mic_vec .* repmat(call_dB' - min(call_dB), 1, 3);
                    norm_I_mic_vec = I_mic_vec ./ max(call_dB - min(call_dB)) * .2; % Normalizing dB's and scaling to fit

                    if use_interp_beamshape % Not sure this is working correctly...
                        [bx_int, by_int, bz_int] = sph2cart(azqm, elqm, .2);
                        nan_indx = isnan(bx_int) | isnan(by_int) | isnan(bz_int);
                        bx_int(nan_indx) = [];
                        by_int(nan_indx) = [];
                        bz_int(nan_indx) = [];

                        XX = [bat(fr, 1); bat(fr, 1) + bx_int'];
                        YY = [bat(fr, 2); bat(fr, 2) + by_int'];
                        ZZ = [bat(fr, 3); bat(fr, 3) + bz_int'];
                    else
                        XX = [bat(fr, 1) + norm_I_mic_vec(:, 1)];
                        YY = [bat(fr, 2) + norm_I_mic_vec(:, 2)];
                        if ~d2
                            ZZ = [bat(fr, 3) + norm_I_mic_vec(:, 3)];
                        end
                    end

                    XYZ = convexHull(delaunayTriangulation([XX YY ZZ]));

                    hullFacets = convexHull(delaunayTriangulation([XX YY ZZ])); 
                    pm = trisurf(hullFacets, XX, YY, ZZ);
                    set(pm, 'FaceColor', [.3 .3 .3], 'FaceAlpha', 0.2, 'edgecolor', [.5 .5 .5]);

                end

                if ~use_interp_beamshape && plot_mics
                    pmic = plot3(bpp.mic_loc(indx_flat, 1), bpp.mic_loc(indx_flat, 2), ...
                        bpp.mic_loc(indx_flat, 3), 'og', 'markerfacecolor', 'g');
                end

                if diag
                    if d2
                        figure(3); set(gcf, 'pos', [2238 322 560 420]); clf
                        [~, isort] = sort(angle_midline);
                        plot(angle_midline(isort), vq2d(isort))

                        hold on;
                        plot(thdirfit, max(vq2d), '+r')

                        if plot_gaussian_beampattern
                            plot(thvals, yvals * 130 * 2)
                        end
                    end

                    figure(1);
                end

                if ~save_movie
                    pause(.25)
                end
            end

%             axis(aa);
            if save_movie
                writeVideo(v, getframe(gcf));
            else
                drawnow
            end

            if fr - call_fr > 15
                delete(pm)
                pm = [];
            end
        end

        if save_movie
            close(v);  % Close the video writer object

            % Construct the output file name
            S = strsplit(trials(tt).name, '_');
            trl_indic = strjoin(S(1:4), '_');
            aud_fn = [mic_data_dir regexprep(trials(tt).name, '_bp_proc', '')];

            % Check if the audio file exists
            if exist(aud_fn, 'file')
                disp(['Loading audio file: ' aud_fn]);
                wav_data = load(aud_fn);

                % Crop the audio
                samp1 = round(frames(1) / (bpp.track.fs) * wav_data.fs);
                samp2 = min(round(frames(end) / (bpp.track.fs) * wav_data.fs), length(wav_data.sig));

                % Ratio of signal to noise
                norm_sig = wav_data.sig ./ repmat(max(abs(wav_data.sig)), size(wav_data.sig, 1), 1);
                noise = median(abs(norm_sig));
                [~, ch] = min(noise);

                y = wav_data.sig(samp1:samp2, ch);
                
                % Write the cropped audio to a WAV file
                wav_file = ['Z:\Nikita\Prism_Experiment\Data\Preprocessing\Mic\Beam_Output\' trials(tt).name '.wav'];
                disp(['Writing WAV file: ' wav_file]);
                audiowrite(wav_file, y ./ max(abs(y)), round(wav_data.fs / (bpp.track.fs / vid_frate)));

                
                
                % Convert audio to AAC with improved quality
                aac_file = ['Z:\Nikita\Prism_Experiment\Data\Preprocessing\Mic\Beam_Output\' trials(tt).name '.aac'];
                command = ['C:\ffmpeg\bin\ffmpeg -y -i ' wav_file ' -c:a aac -b:a 48k ' aac_file];
                disp(['Converting WAV to AAC: ' command]);
                [status, cmdout] = system(command);
                
                if status == 0
                    disp('WAV to AAC conversion successful.');
                else
                    disp('Error converting WAV to AAC:');
                    disp(cmdout);
                end
                
                % Delete the temporary WAV file
                if exist(wav_file, 'file')
                    delete(wav_file);
                else
                    disp(['WAV file not found for deletion: ' wav_file]);
                end
                
                % Merge audio and video with improved quality and noise reduction
                mp4_file_input = ['Z:\Nikita\Prism_Experiment\Data\Preprocessing\Mic\Beam_AV\' out_fname '.mp4'];
                mp4_file_output = ['Z:\Nikita\animate_beam_dirs\' out_fname '.mp4'];
                command = ['C:\ffmpeg\bin\ffmpeg -y -i ' mp4_file_input ' -i ' aac_file ...
                    ' -c:v libx264 -crf 18 -b:v 2000k -vf "hqdn3d" -pix_fmt yuv420p -c:a copy -bsf:a aac_adtstoasc ' mp4_file_output];
                disp(['Merging audio and video: ' command]);
                [status, cmdout] = system(command);
                
                % Check if the merge was successful
                if status == 0
                    disp('Audio and video merged successfully.');
                else
                    disp('Error merging audio and video:');
                    disp(cmdout);
                end
                
                % Delete the temporary AAC file
                if exist(aac_file, 'file')
                    delete(aac_file);
                else
                    disp(['AAC file not found for deletion: ' aac_file]);
                end
                
                
          
            % Only if audio isn't present - still create animation
            if ~exist(['Z:\Nikita\animate_beam_dirs\' out_fname '.mp4'], 'file')
                mp4_file_input = ['Z:\Nikita\Prism_Experiment\Data\Preprocessing\Mic\Beam_AV\' out_fname '.mp4'];
                command = ['C:\ffmpeg\bin\ffmpeg -y -i ' mp4_file_input ' -c:v libx264 -crf 20 -pix_fmt yuv420p ' ...
                    'Z:\Nikita\animate_beam_dirs\' out_fname '_no_audio.mp4'];
                disp(['Creating video without audio: ' command]);
                [status, cmdout] = system(command);

                % Check if the video creation was successful
                if status == 0
                    disp('Video created successfully without audio.');
                else
                    disp('Error creating video without audio:');
                    disp(cmdout);
                end
            end

            % Delete the original MP4 file if it exists
            if exist(['Z:\Nikita\Prism_Experiment\Data\Preprocessing\Mic\Beam_AV\' out_fname '.mp4'], 'file')
                delete(['Z:\Nikita\Prism_Experiment\Data\Preprocessing\Mic\Beam_AV\' out_fname '.mp4']);
            else
                disp(['MP4 file not found: ' mp4_file_input]);
            end
        end
    end
end
