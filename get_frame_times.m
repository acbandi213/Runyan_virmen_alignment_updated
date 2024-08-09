function [alignment_info] = get_frame_times(imaging_base_path, sync_base_path, microscope_id, channel_number,plot_on)


%%%**warnings: do not have spurious folders/files called 'TSeries' in
%%%imaging_base_path. We are assuming 2 imaging channels.

%%%%%inputs
%%%imaging_base_path is the folder where the subfolders of tifs are stored
%%%e.g. imaging_base_path = '/Volumes/Runyan3/Noelle/2P/2P LC/gc700_20220809/'

%%%sync_base_path is the folder where the synchronization data is stored.
%%%e.g. sync_base_path = '/Volumes/Runyan/Code/Imaging and Behavior Alignment/Sample Data/nf_gc700_20220809/'

%%%microscope_id: '2pplus' or 'investigator'

%%%channel_number is the synch channel where slow galvo info is recorded

%%%plot_on: 0 if no plotting, 1 if plotting to check frame detection

%%%%%%terminology
%%%'t_series': a single imaging acquisition
%%%'acquisition_number': a pair of imaging and synchronizations taken
%%%during a single t-series


%%%%outputs
%%alignment_info(acquisition_number).frames_times is 1xnframes vector, the time of each imaging frame in "sync" time (i.e.
%%wavesurfer or pclamp time).



%%special cases to consider in getting frame times:
%%1) photostimulation pulseslonger than 10 or so ms, can get black tifs
%%with no galvo movement, so no frame will be detected in the sync file.

%%2) galvo continues running for some period of time after acquiring the
%%total number of tifs




%%%%run through each acquisition, load the sync file that matches each
%%%%t-series
if length(dir([sync_base_path '*.abf']))>0
    sync_dir = dir([sync_base_path '*.abf']);
    num_syncs = length(sync_dir);
    is_pclamp = 1;
else
    sync_dir = dir([sync_base_path '*.h5']);
    num_syncs = length(sync_dir);
    is_pclamp = 0;
end

%%%%%make sure that we have the same number of synch and t-series
num_tseries = length(dir([imaging_base_path '*TSeries*']));
if num_syncs==num_tseries
    num_acqs = num_syncs;
else
    %%%display an error message
    error('Number of Tseries does not match Number of Syncs')
end

%%get the number of tifs for each TSeries
imaging_dir=dir([imaging_base_path '*TSeries*']);
for i=1:length(imaging_dir)
    %     num_tifs(i)=length(dir([imaging_base_path imaging_dir(i).name '/*.tif']))/2;
    num_tifs(i)=length(dir([imaging_base_path imaging_dir(i).name '/*Ch2*'])); %%%to do: this is slow
end

%%%%TO DO: decide on the best way to force a verification that the sync and
%%%%imaging files are properly matched up. For now printing a list of paired TSeries and Sync files so that they can
%%%%be visually inspected to be sure they match.
for acq_number = 1:num_acqs
    imaging_dir(acq_number).name
    sync_dir(acq_number).name
    
end

%%%%TO DO: get the imaging frame rate from the meta data. Find a way that
%%%%works for all microscopes (ideally), and all matlab versions.


for acq_number = 1:num_acqs
    alignment_info(acq_number).imaging_id = imaging_dir(acq_number).name;
    alignment_info(acq_number).sync_id = sync_dir(acq_number).name;
    %%%load the sync files in order, get slow galvo channel
    if is_pclamp==1;
        [sync_data,sync_sampling_interval,header] = abfload([sync_base_path sync_dir(acq_number).name]);
        sync_sampling_rate = 1/sync_sampling_interval*1e6; %%%convert sampling interval (in us) to sampling rate in hz
        galvo_signal = sync_data(:,channel_number);
        alignment_info(acq_number).sync_sampling_rate = sync_sampling_rate;
    else %%%%%TO DO: add same functionality for h5/wavesurfer files
        
        data = ws.loadDataFile([sync_base_path sync_dir(acq_number).name]);
        fields = fieldnames(data);
        sweep_id = fields{2};
        sync_data = eval(['data.' sweep_id '.analogScans']);
        galvo_signal = sync_data(:,channel_number);
        alignment_info(acq_number).sync_sampling_rate = data.header.AcquisitionSampleRate;
        sync_sampling_rate = data.header.AcquisitionSampleRate;
    end
    
    %%%Normalize the galvo signal so that the same thresholds can be used
    %%%across zooms, etc?
    
    galvo_signal_norm = galvo_signal./max(galvo_signal);
    %%%%identify the frame times using the slow galvo signal
    [~,frame_times]=findpeaks(galvo_signal_norm,'MinPeakHeight',.3,'MinPeakProminence',0.1);
    %     [~,frame_times]=findpeaks(abs(diff(galvo_signal_norm)),'MinPeakHeight',.1);
    %     [~,frame_times]=findpeaks(galvo_signal_norm,'MinPeakHeight',.3,'MinPeakDistance',sync_sampling_rate/35);  %%frame times are in "sync time", which depends on the sync_sample_rate
    
    scanning_amplitude = mean(galvo_signal(frame_times));  %%%in future, could get this by reading in meta data from imaging to get zoom
    good_frame_times = find(galvo_signal(frame_times)<scanning_amplitude*1.1 & galvo_signal(frame_times)>.95*scanning_amplitude); %%%getting ride of artifacts before and after photostims
    frame_times = frame_times(good_frame_times);
    imaging_frame_rate = 1/(mode(diff(frame_times))/sync_sampling_rate); %%overall frame rate, in hz, using sync frame times. TO DO: compare to imaging metadata.
    
    frame_intervals_sec = diff(frame_times).*(1/sync_sampling_rate);
    frame_rate = round(1./frame_intervals_sec);  %%ongoing frame rate in hz
    
    
    if plot_on
        figure(87)
        clf
        plot(galvo_signal)
        hold on
        plot(frame_times,galvo_signal(frame_times),'*r')
        title('Galvo Signal and Defined Frame Times')
        set(gca,'fontsize',12)
        
        
        figure(88)
        clf
        hist(diff(frame_times))
        title('Unique frame intervals - if seeing range of values then check alignment carefully')
        set(gca,'fontsize',12)
        %pause
    end
    %%%frames_to_delete=find(diff(locs)<300);
    %%%locs(frames_to_delete)=nan;
    %%%frame_times=locs(~isnan(locs));
    
    
    %%%%possible way to deal with gaps in galvo scanning due to long
    %%%%photostims
    long_intervals = find(diff(frame_times)>(mode(diff(frame_times)*1.5)));
    long_intervals_for_bad_frames_calc = long_intervals;
    if length(long_intervals)>0
        figure(87)
        hold on
        plot(frame_times(long_intervals),galvo_signal(frame_times(long_intervals)),'*c')
        
        temp_frame_times = [];
        temp_frame_times(1:long_intervals(1)) = frame_times(1:long_intervals(1));
        temp_bad_frames = [];
        for li = 1:length(long_intervals);
            
            
            
            frame_intervals = [];
            last_good_frame_time = frame_times(long_intervals(li));  %%currently, sometimes this is still a partial frame....
            next_good_frame_time = frame_times(long_intervals(li)+1);
            num_frames_to_fill = round((next_good_frame_time-last_good_frame_time)/sync_sampling_rate*imaging_frame_rate)-1; %%was using "floor" rather than "round", thinking that would be the right way, but round seems to be working better.
            
            %%make those frame times
            frame_intervals(2:num_frames_to_fill+1) = 1/imaging_frame_rate*sync_sampling_rate; %%frame rate in sync time units   frame_intervals_sec(2:num_tifs(acq_number));
            frame_intervals(1) = last_good_frame_time;
            fill_frame_times = cumsum(frame_intervals);
            fill_frame_times = fill_frame_times(2:end);
            temp_frame_times = cat(2,temp_frame_times,fill_frame_times);
            if li==length(long_intervals)
                temp_frame_times = cat(2,temp_frame_times,frame_times(long_intervals(li)+1:end)');
            else
                temp_frame_times = cat(2,temp_frame_times,frame_times(long_intervals(li)+1:long_intervals(li+1))');
            end
            %         temp_frame_times(long_intervals(li)+1:long_intervals(li)+num_frames_to_fill) = fill_frame_times;
            %         temp_frame_times(long_intervals(li)+num_frames_to_fill+1:long_intervals(li+1)) = frame_times(long_intervals(li)+1:long_intervals(li+1));  %%need to get the right frame times
            temp_bad_frames = cat(2,temp_bad_frames,long_intervals_for_bad_frames_calc(li)+1:long_intervals_for_bad_frames_calc(li)+num_frames_to_fill);
            long_intervals_for_bad_frames_calc(li:end) = long_intervals_for_bad_frames_calc(li:end)+num_frames_to_fill;  %%%correct bad frames definition for the number of tifs that were missing from previous photostim period. Originally, I thought that we also needed to update frame times, but the galvo-defined frame times become farther and farther apart from the tif-defined frame number.
            figure(87)
            plot(fill_frame_times,galvo_signal(fill_frame_times),'ok')
            title('Galvo Signal-Defined Plus Filled in Frames During Photostims')
            %pause
        end
        
        frame_times = temp_frame_times;
        alignment_info(acq_number).bad_frames = temp_bad_frames;
        
        
        %%%%%%end possible way to deal with gaps due to photostim
    end
    %%%%want to make sure that we have detected the correct number of frames.
    %%%%%%options: use meta data or use the number of tifs in the imaging
    %%%%%%folder. For now using number of tifs, detected above.
    if length(frame_times)>num_tifs(acq_number)  %%this is only the case for standard acquisitions without the long photostims that are missing scans
        frame_times = frame_times(1:num_tifs(acq_number));
    end
    
    alignment_info(acq_number).frame_times = frame_times;
    
    %%%%%%%%%Not using, too much jitter
    %     %%%%alternate option: use syncdata to define the time of the first frame.
    %     %%%%Use the metadata to get precise frame rate. Assign frame times in
    %     %%%%"synch time". Then compare these frame times to the slow galvo signal,
    %     %%%%this allows us to define some types of bad frames (during longer photostim pulses) like above
    %     %%%%because slow galvo is static. Also compare these as a sanity check.
    %     start_time = frame_times(1);  %%%for the option of defining all frame times based on the start and end time
    %     frame_intervals = zeros(1,num_tifs(acq_number));
    %     frame_intervals(2:end) = 1/imaging_frame_rate*sync_sampling_rate; %%frame rate in sync time units   frame_intervals_sec(2:num_tifs(acq_number));
    %     frame_intervals(1) = start_time;
    %     alt_frame_times = cumsum(frame_intervals);
    %
    %     if plot_on && length(frame_times) == length(alt_frame_times)
    %         figure(888)
    %         clf
    %         plot(frame_times,alt_frame_times,'o')
    %         xlabel('Frame times defined by galvo signal')
    %         ylabel('Frame times defined by first galvo frame and then by frame rate')
    %         set(gca,'fontsize',12)
    %         axis square
    %     elseif plot_on && length(frame_times) == length(alt_frame_times)
    %         figure(888)
    %         clf
    %         plot(galvo_signal)
    %         hold on
    %         plot(alt_frame_times,galvo_signal(alt_frame_times),'*r')
    %
    %     end
    %     if length(unique(frame_rate))>1
    %         %%%need to use the alternate alignment method
    %         alignment_info(acq_number).frame_times = alt_frame_times;
    %     else
    %         alignment_info(acq_number).frame_times = frame_times;
    %     end
     %%%Pausing here so user can check the frame times. Zoom in, especially on areas with photostimulations to be sure that frame definitions are correct and do not include artifacts
    
end






