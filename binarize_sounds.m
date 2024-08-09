function [sounds_per_file,all_sounds_together] = binarize_sounds(virmen_it,sound_condition_array, trial_its,sound_info,sound_st,file_trial_ids,alignment_info)
overall_diff = [];

%add up all frames in the prior alignment structures (so things line up
%with suite2p where all frames are concatenated together)
frame_lengths = cellfun(@length,{alignment_info.frame_times});

for file = 1:length(virmen_it)

%     within_trials = virmen_it(file).file_trial_id_start:virmen_it(file).file_trial_id_end;
    within_trials_all = file_trial_ids(file,1):file_trial_ids(file,2);
    
    %sound_trigger_its_diff = virmen_it(file).it_times(trial_its.sound_trigger_its(within_trials_all)-virmen_it(file).actual_it_values(1)+1)-virmen_it(file).it_times(trial_its.sound_trigger_its(within_trials_all)-virmen_it(file).actual_it_values(1)+1-1);
%     sound_trigger = trial_its.sound_trigger_its(within_trials_all)-virmen_it(file).actual_it_values(1)+1;
    sound_trigger_time = virmen_it(file).sound_trigger_time; %need to figure out indexing here
    
    sound_onsets = cell(1,1);
    sound_offsets= cell(1,1);
    [a,b] =find(sound_trigger_time(1) - [sound_condition_array(file).VR_sounds{:,3}]<0,1,'first');
    first_sound_trial = b-1;
    end_trial_time = virmen_it(file).it_times(trial_its.end_trial_its(within_trials_all)-virmen_it(file).actual_it_values(1)+1);
    first_onset = [];

    % adding to make nice structure similar to passive
    previous_sum = sum(frame_lengths(1:file-1)); %keep track of all previous frames

    relevant_trials = 1:length([sound_st(file).file.true_condition]);%find([sound_st(file).file.trial_num]>=within_trials_all(1) & [sound_st(file).file.trial_num]<=within_trials_all(end));
    cond = [sound_st(file).file.true_condition]';
    new_sound_st(file).condition = cond(relevant_trials);
    og_trial = [sound_st(file).file.trial_num]';
    og_trial_file(file,:) = og_trial(end);
    new_sound_st(file).og_trial_num = og_trial(relevant_trials);
    if file == 1
        new_sound_st(file).trial_num = og_trial;
    else
        new_sound_st(file).trial_num = og_trial+sum(og_trial_file(file-1:-1:1,:));
    end

    for s = 1:length(sound_trigger_time) %within_trials_all
            onset = [];
            
            %if nan determine onset and estimate where sounds might be
            %it before gap is about ~180ms from sound onset
            %on for 1 sec off for 200ms
            
            if ~isnan(virmen_it(file).mean_sound_distance)
                onset = sound_trigger_time(s)+virmen_it(file).mean_sound_distance; %conver to ms
            else
                onset = sound_trigger_time(s);
            end
            if onset < end_trial_time(s)
                sound_onsets{s,:} = round(onset:sound_info.distance_within_sounds+sound_info.sound_duration:end_trial_time(s));
                sound_offsets{s,:} = round(onset+sound_info.sound_duration:sound_info.distance_within_sounds+sound_info.sound_duration:end_trial_time(s));
                if size(sound_onsets{s,:},2)>size(sound_offsets{s,:},2)% if last sound is unifinsehd
                    sound_offsets{s,:} = [sound_offsets{s,:},round(end_trial_time(s))];
                end
                first_onset = [first_onset,sound_onsets{s,1}(1,1)];
            end

    end
    file
%     overall_diff = [overall_diff, [[sound_condition_array(file).VR_sounds{1:21,2}]-difference_test(1:21)]];
    
sounds_per_file(file).onsets = [sound_onsets{:,:}];%[[sound_st(file).file.onset],[sound_onsets{:,:}]];
sounds_per_file(file).offsets = [sound_offsets{:,:}];%[[sound_st(file).file.offset],[sound_offsets{:,:}]];

%to deal with trials where sound started before maze (fixed in may 2023)
    og_weird = [];
    if ismember('weird_trial',fields(sound_condition_array(file)))
        % need to figure out indexing 
        og_weird = sound_condition_array(file).weird_trial;
        sounds_per_file(file).weird_trial = og_weird+within_trials_all(1)-first_sound_trial-1;

        dist1 = [];
        for s = 1:length(first_onset)
            temp = min(abs([sound_condition_array(file).VR_sounds{~isnan([sound_condition_array(file).VR_sounds{:,2}]),2}] - first_onset(s)));
            dist1 = [dist1,temp];
        end
        
        %find weird trials that might be missed (happens if it's the very
        %first trial of the file)--- assumes there will be more weird
        %trials within the entire session
        temp_weird = [];
        for d = 1:length(dist1)
            if dist1(d) > 0.03*sound_info.sync_sampling_rate %more than a frame apart
                temp_weird = [temp_weird,d+within_trials_all(1)-first_sound_trial];
                og_weird = [og_weird,d+first_sound_trial];
            end
        end
        og_weird = unique(og_weird);
        sounds_per_file(file).weird_trial = unique([sounds_per_file(file).weird_trial,temp_weird]);
    end
    
    

binary_sounds = zeros(1,(max(sounds_per_file(file).offsets)+sound_info.sync_sampling_rate));
for p = 1:length(sounds_per_file(file).offsets)
    binary_sounds(sounds_per_file(file).onsets(p):sounds_per_file(file).offsets(p)) = 1;
end

%visualize the binary sounds
ex_data = abfload(strcat(virmen_it(file).directory));
figure(555);clf;
title(strcat('Predicted binary sounds in file # ', num2str(file)));
hold on;cc = plot(rescale(ex_data(:,sound_info.spkr_channel_number(1)),0,1),'-b');
if length(sound_info.spkr_channel_number) > 1
    dd = plot(rescale(ex_data(:,sound_info.spkr_channel_number(2)),0,1),'-m');
    %legend([ cc dd  a(1) ],'Speaker 1','Speaker 2', 'binary sounds')

end
a = plot(binary_sounds,'-k');
legend([ cc  a(1) ],'Speaker 1', 'binary sounds')

%plot(rescale(ex_data(:,sound_info.spkr_channel_number(3)),-1,0),'-r');
if length(sound_info.spkr_channel_number)>2
    plot(rescale(ex_data(:,sound_info.spkr_channel_number(3)),0,1),'-r')
end
if ~isempty(og_weird)
    plot([sound_condition_array(file).VR_sounds{og_weird,2}],1,'*c')
end
hold off;
pause

%organize with actually detected sounds
sound_onsets = cell(1,1);
sound_offsets= cell(1,1);
for s = 1:length(sound_trigger_time) %
    if s+first_sound_trial < length(sound_condition_array(file).VR_sounds) && isnan(sound_condition_array(file).VR_sounds{s+first_sound_trial,3}) 
        onset = [];
        %if nan determine onset and estimate where sounds might be
        %it before gap is about ~180ms from sound onset
        %on for 1 sec off for 200ms
        if ~isnan(virmen_it(file).mean_sound_distance)
                onset = sound_trigger_time(s)+virmen_it(file).mean_sound_distance; %conver to ms
        else
            onset = sound_trigger_time(s);
        end
        if onset < end_trial_time(s)

            sound_onsets{s,:} = round(onset:sound_info.distance_within_sounds+sound_info.sound_duration:end_trial_time(s));
            sound_offsets{s,:} = round(onset+sound_info.sound_duration:sound_info.distance_within_sounds+sound_info.sound_duration:end_trial_time(s));
            if size(sound_onsets{s,:},2)>size(sound_offsets{s,:},2)% if last sound is unifinsehd
                sound_offsets{s,:} = [sound_offsets{s,:},round(end_trial_time(s))];
            end
        end
    end
end
sounds_per_file(file).onsets = [[sound_st(file).file.onset],[sound_onsets{:,:}]];
sounds_per_file(file).offsets = [[sound_st(file).file.offset],[sound_offsets{:,:}]];
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
    new_sound_st(file).corr_frames(p,:) = new_sound_st(file).frames(p,:)+previous_sum;

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

dist = [];
for s = 1:length(first_onset)
    temp = min(abs([sound_condition_array(file).VR_sounds{~isnan([sound_condition_array(file).VR_sounds{:,2}]),2}] - first_onset(s)));
    temp = temp(find(temp<.1*sound_info.sync_sampling_rate));
    dist = [dist,temp];
end
sounds_per_file(file).distance_to_true = dist*1000/sound_info.sync_sampling_rate; %convert to ms

end

%make structure similar to passive data!
all_sounds_together = [];
f = fieldnames(new_sound_st);
for i = 1:length(f)
    all_sounds_together.(f{i}) = cat(1,new_sound_st(1:length(sound_st)).(f{i}));%[new_sound_st(1:length(sound_st)).(f{i})];%{[new_sound_st(1:2).(f{i})]};
end

 mean_al = [];var_al =[];for i = 1:length(virmen_it); mean_al = [mean_al,mean(sounds_per_file(i).distance_to_true)];var_al = [var_al,var(sounds_per_file(i).distance_to_true)];end
 figure(556);clf; scatter(mean_al,var_al,[],'markeredgecolor','k','linewidth',1.5);xlabel('mean distance from true sound in ms'); ylabel('variance in  ms'); title('distance between true and predicted sound onsets')


