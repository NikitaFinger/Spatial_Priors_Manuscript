%allows freq_desired to be passed into script as a variable
addpath 'Y:\Spatial Memory Project (Xiaoyan and Nikita)\beampattern_analysis\Working Code\'
%  addpath 'Y:\Spatial Memory Project (Xiaoyan and Nikita)\beampattern_analysis\Working Code\'
% addpath '/Volumes/data1/Spatial Memory Project (Xiaoyan and Nikita)/beampattern_analysis/beampattern_processing-master'

a=[];
vars=who;
clear_vars=setdiff(vars,{'freq_desired','side_view','save_movie','bat_type'});
clear(clear_vars{:});

if ~exist('save_movie','var')
  save_movie=1;
end
vid_frate=12;  

if ~exist('side_view','var')
  side_view=0; %if not side view, then top view
end

if ~exist('freq_desired','var')
  freq_desired=35; %khz
end

 
  side_view=0;
  save_movie=1;
% if exist('bat_type','var') && strcmp(bat_type,'rousettus')
  mic_proc_dir='Z:\Nikita\Prism_Experiment\Data\Preprocessing\Mic\Beam_Output\';
%   mic_proc_dir='Y:\Spatial Memory Project (Xiaoyan and Nikita)\beampattern_analysis\Working Code\proc_output_final\';
%    mic_proc_dir='/Volumes/data1/Spatial Memory Project (Xiaoyan and Nikita)/beampattern_analysis/beampattern_processing-master/bp_proc_main_folder/proc_output/';
  
trials=dir([mic_proc_dir 'batS_6_20230630_095042_mic_data_bp_proc*']);
% else
%   mic_proc_dir='..\proc_output_eptesicus_new\';
%   trials=dir([mic_proc_dir 'eptesicus_20150824_*_mic_data_bp_proc.mat']);
%   checked=1;
% end

close all;

mic_data_dir='Z:\Nikita\Prism_Experiment\Data\Preprocessing\Mic\Mic_Data_Detect';
% mic_data_dir='Y:\Spatial Memory Project (Xiaoyan and Nikita)\beampattern_analysis\Working Code\mic_data\';

% mic_proc_dir='..\proc_output_eptesicus_new\';
% trials=dir([mic_proc_dir 'eptesicus_20150824_*_mic_data_bp_proc.mat']);
% checked=1;

use_interp_beamshape=1;
% interp_method='rb_natural';  % natural neighbor interpolation
interp_method='rb_natural';  % radial basis function interpolation'

d2=1; %2d data only
if d2
  gaussfitfcn=1; %1 to use fit function of matlab to determine beam direction
end

frames_limit_hard=1;

plot_gaussian_beampattern=1;
plot_mics=1;
plot_mic_nums=1;
flatten_deg=30;
diag=0;

for tt=1  %:length(trials)
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
  
  if save_movie
    out_fname=trials(tt).name(1:end-21);
    if side_view
      out_fname=[out_fname '_side'];
    else
      out_fname=[out_fname '_top'];
    end
    if use_interp_beamshape
      out_fname=[out_fname '_interp'];
    else
      out_fname=[out_fname '_raw'];
    end
    out_fname=[out_fname '_' num2str(freq_desired)];
    v=VideoWriter(['Z:\Nikita\Prism_Experiment\Data\Preprocessing\Mic\Beam_AV' out_fname '.avi'],'uncompressed avi');
%     v=VideoWriter(['Y:\Spatial Memory Project (Xiaoyan and Nikita)\beampattern_analysis\Working Code\beam_aim_plots\' out_fname '.avi'],'uncompressed avi');
   
    v.FrameRate=vid_frate;
    open(v);
  end
  
  calls_with_track = bpp.mic_data.call_idx_w_track;
  call_track_locs = round([bpp.mic_data.call.call_start_idx]...
    /bpp.mic_data.fs*bpp.track.fs);
%   For trials with only 4 seconds of recorded audio.
%     call_track_locs = round([bpp.mic_data.call.call_start_idx]...
%     /bpp.mic_data.fs*bpp.track.fs) + 400;

   goodch=1:length(bpp.mic_loc(:,1)); %alternatively load in the good ch from the ch_ex var
  mic_num = 1:bpp.mic_data.num_ch_in_file;
  
  %%% add the original (2) and testing (3) perch location of the mic
  %%% location array from mic pos file:
%   6/28/2022
    %bpp.mic_loc(2,:)=[-0.1563638153076172, 2.271655517578125, 1.2348447265625]; %original perch location 6/28/2022
    %-0.18473699951171876	2.323274658203125	1.2506273193359374
    %-0.18473699951171876	2.323274658203125	1.2506273193359374
    %bpp.mic_loc(2,:)=[-0.8388746337890625	2.994169677734375	1.2176396484375];
    %bpp.mic_loc(2,:)=[-0.023652957916259765 2.3611552734375 1.24591162109375]



    %bpp.mic_loc(2,:)=[-0.18473699951171876	2.323274658203125	1.2506273193359374]; % 15L

% 7/01/2022
    %bpp.mic_loc(2,:)=[-0.023652957916259765 2.3611552734375 1.24591162109375]; %original 7/01/2022





% bpp.mic_loc(3,:)=[-0.267540313720703,2.27989038085938,1.24859375000000];   %testing perch location   15L dark
%   bpp.mic_loc(3,:)=[-0.44,	2.57671435546875,	1.24859375000000];   %testing perch location   30L
 
%   bpp.mic_loc(3,:)=[-0.184737,	2.32327,	1.25063];   %testing perch location   15L light
%  bpp.mic_loc(3,:)=[0.000488172352313995	2.57671435546875	1.2458001708984374]; 

  close all;
  figure(1); set(gcf,'pos',[10 40 900 900],'color','w')
  set(gca,'position',[.07 .06 .88 .89],'units','normalized')
  hold on;
  pb_mics=plot3(bpp.mic_loc(goodch,1),bpp.mic_loc(goodch,2),bpp.mic_loc(goodch,3),...
    'ok','markerfacecolor','k');
%% change original loc to blue and testing loc to red
%   hold on;
%   plot3(bpp.mic_loc(2,1),bpp.mic_loc(2,2),bpp.mic_loc(2,3),...
%     'ok','markerfacecolor','b','Marker','pentagram','MarkerSize',18);
%   hold on;
%   plot3(bpp.mic_loc(3,1),bpp.mic_loc(3,2),bpp.mic_loc(3,3),...
%     'ok','markerfacecolor','g','Marker','pentagram','MarkerSize',18);

%   if plot_mic_nums
%     text(bpp.mic_loc(goodch,1),bpp.mic_loc(goodch,2),bpp.mic_loc(goodch,3),...
%       num2str(mic_num'))
%   end

  if frames_limit_hard
    frames = call_track_locs( calls_with_track( find(bpp.proc.chk_good_call,1) ) ) :...
      call_track_locs( calls_with_track( find(bpp.proc.chk_good_call,1,'last') ) );

    %limit to speed > 1.5 m/s
    speed = calc_speed(bat,bpp.track.fs,0);
    speed = [nan; speed];
    good_speed = find(speed > 0.1);
    good_speed (good_speed < frames(1))=[];
    good_speed (good_speed > frames(end))=[];

    frames = max(good_speed(1),frames(1)) : min(good_speed(end),frames(end));
  else
    frames=find(isfinite(bat(:,1)));
    %limiting the frames plotted in animation to around the time of the first
    %and last good calls
    frames=max(frames(1),...
      call_track_locs(calls_with_track( find(bpp.proc.chk_good_call,1) ) )-10):...
      min(call_track_locs(calls_with_track(find(bpp.proc.chk_good_call,1,'last')))+10,...
      frames(end));
  end
  
  
  pb=plot3(bat(frames(1):frames(end),1),...
    bat(frames(1):frames(end),2),...
    bat(frames(1):frames(end),3),'-k','linewidth',2);
  if side_view
    view(0,0)
  else
    view(2)
  end
  axis equal, grid on;
  aa=axis;
  set(gca,'fontsize',22,...
    'XTick',aa(1):.4:aa(2),'xticklabel',(aa(1):.4:aa(2))-aa(1),...
    'YTick',aa(3):.4:aa(4),'yticklabel',(aa(3):.4:aa(4))-aa(3))
  if side_view
    title_name='Side View';
  else
    title_name='Top View';
  end
  title(title_name,'fontsize',24)
  box on;
  aa=axis;
  
  delete(pb);  
  if ~plot_mics
    delete(pb_mics);
  end
  
  call_fr=0; pm=[]; pmic=[]; BX=[]; BY=[];
  for fr=frames(1): frames(end)
   
    pb=plot3(bat(frames(1):fr,1),bat(frames(1):fr,2),bat(frames(1):fr,3),...
      '-','linewidth',2,'color',[.2 .2 .2]);
    %     call_indx= find(proc.call_loc_on_track_interp<=fr,1,'last');
    
    %if a call is present on the frame and selected as 'good call' we plot the beam pattern
    good_call = find(bpp.proc.chk_good_call);
    call_present=ismember(call_track_locs(calls_with_track(good_call)),fr);
    
    call_indx_1=find(call_present);
    call_indx = good_call(call_indx_1);
    if ~isempty(call_indx) && ismember(fr,frames)
      delete(pm); pm=[];
      
      if ~use_interp_beamshape
        delete(pmic); pmic=[];
      end
      
      call_fr=fr;
      call_dB = nan(1,bpp.mic_data.num_ch_in_file);
      for iM=mic_num
        freq = bpp.proc.call_freq_vec{call_indx,iM};
        [~,fidx] = min(abs(freq-freq_desired*1e3));
        if isempty (bpp.proc.call_psd_dB_comp_re20uPa_withbp{call_indx,iM})
            call_dB(iM) = 0;
        else
        call_dB(iM) = bpp.proc.call_psd_dB_comp_re20uPa_withbp{call_indx,iM}(fidx);
        end
      end
      
      ch_ex_sig = find(isnan(call_dB)); % low quality channel from call extraction function
      ch_good_loc = ~isnan(bpp.mic_loc(:,1))';  % Check for mics without location data
      
      angle_notnanidx = ~ismember(mic_num,ch_ex_sig) & ch_good_loc;
      goodch1 = goodch(angle_notnanidx);  %xiaoyan
      
%       mic_to_bat_angle = squeeze(bpp.proc.mic_to_bat_angle(call_indx,:,:));
%       az = mic_to_bat_angle(angle_notnanidx,1);
%       el = mic_to_bat_angle(angle_notnanidx,2);
      
%       %create my own angles for bat to mic because we are plotting on
%       cartesian coordinates, not on head centered plane
      
      mic_vec=bpp.mic_loc(goodch1,:)-...
        repmat(bat(fr,:),size(bpp.mic_loc(goodch1,:),1),1); %vectors from bat to all mics pointing to bat
      [az,el] =...
        cart2sph(mic_vec(:,1),mic_vec(:,2),mic_vec(:,3)); %Transform Cartesian to spherical coordinates
      
      % doing the interpolation:
      if use_interp_beamshape
%         [azq,elq] = meshgrid(linspace(min(az),max(az),150),linspace(min(el),max(el),150));
        [azq,elq] = meshgrid(min(az):pi/180:max(az),min(el):pi/180:max(el));
        if strcmp(interp_method,'rb_natural')  % natural neighbor interpolation
          vq = griddata(az,el,call_dB(angle_notnanidx),azq,elq,'natural');
        elseif strcmp(interp_method,'rb_rbf')  % radial basis function interpolation
          vq = rbfinterp([azq(:)';elq(:)'],rbfcreate([az(:)';el(:)'],call_dB(angle_notnanidx),'RBFFunction','multiquadrics'));
          vq = reshape(vq,size(azq));
        end
        maxref = max(call_dB(angle_notnanidx));
        vq_norm = vq-maxref;

        % Find indices within measured polygon
        k = boundary(az,el,0);  % outer boundary of all measured points
        [in,on] = inpolygon(azq,elq,az(k),el(k));
        in_smpl_poly = in|on;
        clear in on
        vq(~in_smpl_poly) = NaN;
        vq_norm(~in_smpl_poly) = NaN;

        [~,II] = max(vq_norm(:));
        [I,J]=ind2sub(size(vq_norm),II);

        elqm = elq;
        elqm(isnan(vq_norm)) = NaN;
        azqm = azq;
        azqm(isnan(vq_norm)) = NaN;
      end
      
%       head_dir_xyz=bpp.head_aim.head_aim_int(...
%         bpp.track.call_loc_idx_on_track_interp(call_indx),:);
%       [h_az,h_el]=cart2sph(head_dir_xyz(1),head_dir_xyz(2),head_dir_xyz(3));
      
      %plotting the beam direction
      if bpp.proc.chk_good_call(call_indx)
        beam_col=[0, 109, 219]./255;
      else
        beam_col=[228 26 28]./255;
      end
      
      if d2
        if use_interp_beamshape
          if side_view
            indx_flat=find( azqm <= 1/180*pi & azqm >= -1/180*pi ); %pull out indexes within 1 degree el in the interpolated data
            angle_midline = elqm(indx_flat);
          else
            indx_flat=find( elqm <= 30/180*pi & elqm >= -30/180*pi ); %pull out indexes within 1 degree el in the interpolated data
            angle_midline = azqm(indx_flat);
          end
          vq2d = vq_norm(indx_flat);
          vq2d = vq2d - min(vq2d);
          norm_vq2d = vq2d ./max(vq2d);

          [~,imax2d]=max(norm_vq2d);
          
          [sigma,~]=gaussfit(angle_midline,vq2d');
          thdirfit = angle_midline(imax2d);

          [bx,by]=pol2cart(thdirfit,.3);  %Transform polar to Cartesian coordinates.
        else
          if side_view
            indx_flat=find( az <= flatten_deg/180*pi & az >= -flatten_deg/180*pi ); %pull out indexes within 1 degree el in the interpolated data
            angle_midline = el(indx_flat);
          else
            indx_flat=find( el <= flatten_deg/180*pi & el >= -flatten_deg/180*pi ); %pull out indexes within 1 degree el in the interpolated data
            angle_midline = az(indx_flat);
          end
          vq_norm=call_dB-max(call_dB);
          vq2d = vq_norm(indx_flat);
          vq2d = vq2d - min(vq2d);
          norm_vq2d = vq2d'./max(abs(vq2d'));
          
          if numel(angle_midline)>=5
            if gaussfitfcn
              [sigma, thdirfit]=gaussfit(angle_midline,vq2d');
            else
  %           f = fit(th,call_dB','gauss1');
  %           thdirfit=f.b1;
              [thdirfit,sigma]=normfit(angle_midline,vq2d);
            end
          else
            thdirfit=nan;
          end

          [bx,by]=pol2cart(thdirfit,.3);
        end
        
        if side_view
          plot3([bat(fr,1) bat(fr,1)+bx],...
            [bat(fr,2) bat(fr,2)],...
            [bat(fr,3) bat(fr,3)+ by],...
            '-','linewidth',4,'color',beam_col);
        else
          plot([bat(fr,1) bat(fr,1)+bx],[bat(fr,2) bat(fr,2)+ by],...
            '-','linewidth',4,'color',beam_col);
        BX=[BX; bx]; BY=[BY; by];
        end
      else
        [bx,by,bz]=sph2cart(azqm(I,J),elqm(I,J),.2);
        plot3([bat(fr,1) bat(fr,1)+bx],[bat(fr,2) bat(fr,2)+ by],...
          [bat(fr,3) bat(fr,3)+ bz],...
          '-','linewidth',2,'color',beam_col);
      end
      
      %plotting the beampattern
      pm=[];
      if d2
        if plot_gaussian_beampattern
          %create gaussian
          thvals=-pi: 4/180*pi :pi;
          yvals = 1/(sqrt(2*pi)*sigma)*exp( -(thvals - thdirfit).^2 / (2*sigma^2));
          
          [Bx,By]=pol2cart(thvals,yvals./max(yvals)*.35); %Transform polar to Cartesian coordinates.
          symbol='';
          %plot
          if side_view
            pm=plot3(repmat(bat(fr,1),size(Bx))+Bx,...
              repmat(bat(fr,2),size(By))+zeros(size(By)),...
              repmat(bat(fr,3),size(By))+By,...
              ['-' symbol],'linewidth',2,'color',[.4 .4 .4]);
          else
            pm=plot(repmat(bat(fr,1),size(Bx))+Bx,...
              repmat(bat(fr,2),size(By))+By,...
              ['-' symbol],'linewidth',2,'color',[.4 .4 .4]);
          end
          
        else
          [~,isort]=sort(angle_midline);
          [Bx,By]=pol2cart(angle_midline(isort),norm_vq2d(isort).*.3);
          if use_interp_beamshape
            symbol='';
          else
            symbol='+';
          end
          if side_view
            pm=plot3(repmat(bat(fr,1),size(Bx))+Bx,...
              repmat(bat(fr,2),size(By))+zeros(size(By)),...
              repmat(bat(fr,3),size(By))+By,...
              ['-' symbol],'linewidth',1,'color',[.4 .4 .4]);
          else
            pm=plot(repmat(bat(fr,1),size(Bx))+Bx,repmat(bat(fr,2),size(By))+By,...
              ['-' symbol],'linewidth',1,'color',[.4 .4 .4]);
          end
        end
      else
        mic_vec=bpp.mic_loc(goodch,:)-...
          repmat(bat(fr,:),size(bpp.mic_loc(goodch,:),1),1);
        norm_mic_vec=mic_vec./...
          repmat(sqrt((mic_vec(:,1).^2+mic_vec(:,2).^2+mic_vec(:,3).^2)),1,3);
        I_mic_vec=norm_mic_vec.*repmat(call_dB'-min(call_dB),1,3);
        norm_I_mic_vec=I_mic_vec./max(call_dB-min(call_dB)).*.2; %normalizing dB's and scaling to fit
%         [az,el,Rho_I]=cart2sph(norm_I_mic_vec(:,1),norm_I_mic_vec(:,2),norm_I_mic_vec(:,3));

        if use_interp_beamshape %not sure this is working correctly...
          [bx_int,by_int,bz_int]=sph2cart(azqm,elqm,.2);
          nan_indx=isnan(bx_int) | isnan(by_int) | isnan(bz_int);
          bx_int(nan_indx)=[];
          by_int(nan_indx)=[];
          bz_int(nan_indx)=[];

          XX=[bat(fr,1); bat(fr,1)+bx_int'];
          YY=[bat(fr,2); bat(fr,2)+by_int'];
          ZZ=[bat(fr,3); bat(fr,3)+bz_int'];
        else
          XX=[bat(fr,1)+norm_I_mic_vec(:,1)];
          YY=[bat(fr,2)+norm_I_mic_vec(:,2)];
          if ~d2
            ZZ=[bat(fr,3)+norm_I_mic_vec(:,3)];
          end
        end
        
        XYZ = convexHull(delaunayTriangulation([XX YY ZZ]));
        
        hullFacets = convexHull(delaunayTriangulation([XX YY ZZ])); 
        pm=trisurf(hullFacets,XX,YY,ZZ);
        set(pm,'FaceColor',[.3 .3 .3],'FaceAlpha',0.2,'edgecolor',[.5 .5 .5]);
%         colormap gray
%         light('Position',[50 -15 29])
%         lighting phong
%         shading interp
        
%         pm=scatter3(XX,YY,ZZ,[],call_dB,'o','filled');
%         colorbar
      end
      
      if ~use_interp_beamshape && plot_mics
        pmic=plot3(bpp.mic_loc(indx_flat,1),bpp.mic_loc(indx_flat,2),...
          bpp.mic_loc(indx_flat,3),'og','markerfacecolor','g');
      end
      
      
      if diag
        if d2
          figure(3); set(gcf,'pos',[2238 322 560 420]); clf
          [~,isort]=sort(angle_midline);
          plot(angle_midline(isort),vq2d(isort))
          
          hold on;
          plot(thdirfit,max(vq2d),'+r')
          
          if plot_gaussian_beampattern
            plot(thvals,yvals*130*2)
          end
        end
        
%         xp = linspace(min(th),max(th),50);
%         if gaussfitfcn
%           yp = 1/(sqrt(2*pi)* sigma ) * exp( - (xp-thdirfit).^2 / (2*sigma^2));
%         else
%           %         yp=feval(f,xp)';
%           yp = 1/(sqrt(2*pi)* sigma ) * exp( - (xp-thdirfit).^2 / (2*sigma^2));
%         end
%         
%         figure(2); clf;
%         polar(th,Rho_I,'+');
%         hold on;
%         plot(thdir,max(call_dB),'+k')
%         plot(xp,yp.*max(call_dB).*2,'r')
%         
%         plot(thdirfit,max(call_dB),'+g')
        
        figure(1);
      end
      
      %       [XD,YD]=pol2cart(xp,yp/max(yp)/2);
      %       pf=plot(XD+bat(fr,1), YD+bat(fr,2),'-','color',[.4 .4 .4]);
      if ~save_movie
        pause(.25)
      end
    end
    
    axis(aa);
    if save_movie
%       set(gca, 'gridcolor','k','linewidth',1.5,'GridAlpha',.25)
      writeVideo(v,getframe(gcf));
    else
      drawnow
    end
    
%     delete(pb);   %keep the bat fly path, xiaoyan
    if fr-call_fr > 15
      delete(pm)
      pm=[];
    end
   
  end
  
  if save_movie
    close(v);
    
    S=strsplit(trials(tt).name,'_');
    trl_indic = strjoin(S(1:4),'_');
  
    aud_fn=[mic_data_dir regexprep(trials(tt).name,'_bp_proc','')];
    if exist(aud_fn,'file')
      wav_data=load(aud_fn);
      %crop the audio
      samp1=round(frames(1)/(bpp.track.fs)*wav_data.fs);
      samp2=min(round(frames(end)/(bpp.track.fs)*wav_data.fs),length(wav_data.sig));
      
      %ratio of signal to noise
      norm_sig = wav_data.sig./ repmat(max(abs(wav_data.sig)),size(wav_data.sig,1),1);
      noise=median(abs(norm_sig));
      [~,ch]=min(noise);
      
      y=wav_data.sig(samp1:samp2,ch);
%       audiowrite(['Y:\Spatial Memory Project (Xiaoyan and Nikita)\beampattern_analysis\Working Code\proc_output\' trials(tt).name '.wav'],y./max(abs(y)),...
      audiowrite(['Z:\Nikita\Prism_Experiment\Data\Preprocessing\Mic\Beam_Output\' trials(tt).name '.wav'],y./max(abs(y)),...
        round(wav_data.fs/(bpp.track.fs/vid_frate)) )
      
      %convert audio
      [~,~]=system(['C:\video_tools\ffmpeg_git\bin\ffmpeg -y ' ...
        '-i Z:\Spatial Memory Project (Xiaoyan and Nikita)\beampattern_analysis\Working Code\proc_output\' trials(tt).name '.wav ' ...
        '-c:a libfdk_aac -b:a 48k -cutoff 18000 ' ...
        'Z:\Spatial Memory Project (Xiaoyan and Nikita)\beampattern_analysis\Working Code\proc_output\' trials(tt).name '.aac ']);
      delete(['Z:\Spatial Memory Project (Xiaoyan and Nikita)\beampattern_analysis\Working Code\proc_output\' trials(tt).name '.wav'])
      
      %merge audio and video
      [~,~]=system(['C:\video_tools\ffmpeg_64\bin\ffmpeg -y '...
        '-i Z:\Spatial Memory Project (Xiaoyan and Nikita)\beampattern_analysis\Working Code\proc_output\' out_fname '.avi '...
        '-i Z:\Spatial Memory Project (Xiaoyan and Nikita)\beampattern_analysis\Working Code\proc_output\' trials(tt).name '.aac ' ...
        '-c:v libx264 -crf 20 -pix_fmt yuv420p -c:a copy -bsf:a aac_adtstoasc ' ...
        '..\animate_beam_dirs\' out_fname '.mp4']);
      delete(['Z:\Spatial Memory Project (Xiaoyan and Nikita)\beampattern_analysis\Working Code\proc_output\' trials(tt).name '.aac'])
    end
    
    %only if audio isn't present - still create animation
    if ~exist(['..\animate_beam_dirs\' out_fname '.mp4'],'file')
      [~,~]=system(['C:\video_tools\ffmpeg_64\bin\ffmpeg -y ' ...
        '-i Z:\Spatial Memory Project (Xiaoyan and Nikita)\beampattern_analysis\Working Code\proc_output\' out_fname '.avi '...
        '-c:v libx264 -crf 20 -pix_fmt yuv420p ' ...
        '..\animate_beam_dirs\' out_fname '_no_audio.mp4']);
    end
    delete(['Z:\Nikita\Prism_Experiment\Data\Preprocessing\Mic\Beam_Output\' out_fname '.avi'])
  end
end
% fprintf('Size of call_psd_dB_comp_re20uPa_withbp: %s\n', mat2str(size(bpp.proc.call_psd_dB_comp_re20uPa_withbp)));
% fprintf('Value of fidx: %d\n', fidx);
% fprintf('Value of call_indx: %d\n', call_indx);
% fprintf('Value of iM: %d\n', iM);
