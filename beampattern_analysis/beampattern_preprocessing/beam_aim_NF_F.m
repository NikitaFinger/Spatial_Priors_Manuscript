function createBeamPatternVideo(freq_desired, side_view, save_movie)
    addpath 'Z:\Nikita\Prism_Experiment\Analysis\beampattern_analysis\beampattern_preprocessing\'
    
    a = [];
    vars = who;
    clear_vars = setdiff(vars, {'freq_desired', 'side_view', 'save_movie', 'bat_type'});
    clear(clear_vars{:});

    if ~exist('save_movie', 'var')
        save_movie = 1;
    end
    vid_frate = 12;

    if ~exist('side_view', 'var')
        side_view = 0; %if not side view, then top view
    end

    if ~exist('freq_desired', 'var')
        freq_desired = 35; %khz
    end

    side_view = 0;
    save_movie = 1;

    mic_proc_dir = 'Z:\Nikita\Prism_Experiment\Data\Preprocessing\Mic\Beam_Output\';

    trials = dir([mic_proc_dir 'batS_6_20230630_095042_mic_data_bp_proc*']);

    close all;

    mic_data_dir = 'Z:\Nikita\Prism_Experiment\Data\Preprocessing\Mic\Mic_Data_Detect';

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

    for tt = 1 % :length(trials)
       bpp=load([mic_proc_dir trials(tt).name]);
  if exist('checked','var') && checked 
    if exist([mic_proc_dir trials(tt).name(1:end-4) '_checked.mat'],'file')
      bpp_checked=load([mic_proc_dir trials(tt).name(1:end-4) '_checked.mat']);
      bpp.proc.chk_good_call = bpp_checked.proc.chk_good_call;
      bpp.proc.ch_ex=bpp_checked.proc.ch_ex;
    else
      continue
    end
  end
  bat=bpp.track.track_smooth;
  
  bpp.proc.call_psd_dB_comp_re20uPa_withbp(cellfun(@isempty,bpp.proc.call_psd_dB_comp_re20uPa_withbp)) = {NaN(1,63)};

    end

    if save_movie
        % Save the video and audio if needed
    end
end
