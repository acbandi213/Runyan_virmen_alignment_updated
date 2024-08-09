function [sounds_per_file,all_sounds_together] = binarize_sounds_simple(virmen_it,sound_condition_array, trial_its,sound_info,sound_st,file_trial_ids,alignment_info)
all_sounds_together = []; %used for passive data

%add up all frames in the prior alignment structures (so things line up
%with suite2p where all frames are concatenated together)
frame_lengths = cellfun(@length,{alignment_info.frame_times});

for file = 1:length(virmen_it)

    within_trials_all = file_trial_ids(file,1):file_trial_ids(file,2);
    
    sound_trigger_time = virmen_it(file).sound_trigger_time; %need to figure out indexing here
    
    
    [a,b] =find(sound_trigger_time(1) - [sound_condition_array(file).VR_sounds{:,3}]<0,1,'first');
    first_sound_trial = b-1;
    end_trial_time = virmen_it(file).it_times(trial_its.end_trial_its(within_trials_all)-virmen_it(file).actual_it_values(1)+1);
    first_onset = [];

sounds_per_file(file).onsets = [sound_st(file).file.onset];
sounds_per_file(file).offsets = [sound_st(file).file.offset];
binary_sounds = zeros(1,(max(sounds_per_file(file).offsets)+sound_info.sync_sampling_rate));
binary_sounds_frames = [];
min_distance = mode(diff(alignment_info(file).frame_times))+(2*alignment_info(file).sync_sampling_rate/1000); %minimum distance between frames in digidata time
binary_sounds_imaging_time = zeros(1,length(alignment_info(file).frame_times));

for p = 1:length(sounds_per_file(file).offsets)
    %binarize signal based on digidata time
    binary_sounds(sounds_per_file(file).onsets(p):sounds_per_file(file).offsets(p)) = 1; %in digidata time
    %find onset and offset frames and binarize based on their positions
    binary_sounds_frames(p,:) =  [find_closest_frames(alignment_info(file).frame_times,sounds_per_file(file).onsets(p),min_distance),find_closest_frames(alignment_info(file).frame_times,sounds_per_file(file).offsets(p),min_distance)];
    new_sound_st(file).frames(p,:) =  binary_sounds_frames(p,:);
%     new_sound_st(file).corr_frames(p,:) = new_sound_st(file).frames(p,:)+previous_sum;

    if all(~isnan(binary_sounds_frames(p,:)))
        binary_sounds_imaging_time(binary_sounds_frames(p,1):binary_sounds_frames(p,2)) = 1;
    end
end



figure(557);clf;
hold on
title(strcat('Frame times sounds in file # ', num2str(file)));
plot(binary_sounds_imaging_time,'-k')
hold off
pause(2)

sounds_per_file(file).binary_digidata_times = binary_sounds;
sounds_per_file(file).binary_frame_times = binary_sounds_imaging_time;


end
