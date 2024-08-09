function [sound_outputs,trialConditions, condition_onset_array_all]=find_spkr_output_task(server,mouse,date,alignment_info,spkr_channel_number,string,distance_between_sounds,threshold_spk,distance_within_sounds) 
%pc=1 if windows, any other number if mac
cd(strcat(server,'\Connie\RawData\',num2str(mouse),'\wavesurfer\',num2str(date)));
sync_dir = dir(strcat('*',string,'*.abf'));
num_files = length(sync_dir);
file_ind = 0;

for file = 1:num_files
    file_ind = file_ind +1
    [sync_data,sync_sampling_interval,~]  = abfload(sync_dir(file).name);
    sync_sampling_rate = 1/sync_sampling_interval*1e6;
    sync_data=sync_data';
    sync_data=double(sync_data);
    rawSounds=sync_data(spkr_channel_number,:);
    min_threshold = 0.0007; %finds peaks above this threshold
    freq_range = [0.2 2];
    [pure_tones] = bandpass_sound(normr(rawSounds),freq_range,distance_between_sounds,min_threshold,sync_sampling_rate);
%     reward_tones = bandpass(rawSounds(1,:),[0.2 2],sync_sampling_rate); %trying to obtain pure tones within sounds
%     reward_sounds = find(abs(reward_tones)>40);
%     for i  = 1:length(rawSounds);rawSoundsCombined(i) =max(abs(rawSounds(:,i)));end
  
    frames_times{file} = alignment_info(file).frame_times; 
    for i=1:length(spkr_channel_number)
        rawSounds(i,:)=(rawSounds(i,:)-mean(rawSounds(i,:)))/range(rawSounds(i,:));
        rawSounds(i,:)=abs(rawSounds(i,:));
        rawSounds(i,:)=smooth(rawSounds(i,:),0.025*sync_sampling_rate,'moving');
    end
    baseline = min(rawSounds,[],'all');
    rawSounds = rawSounds - (baseline*10);

    norm_sound = normr(rawSounds);
    
    %figure(1);for i=1:length(spkr_channel_number); subplot(length(spkr_channel_number),1,i); plot(norm_sound(i,:)); end; movegui(gcf,'center')
    s=0;
    sounds1=[];
    %threshold_spk = 1.2e-4;%0.0001; %was doing 0.0003 works with nat sounds but not reward/incorrect sounds
    for i=2:length(rawSounds)-2
        if max(norm_sound(:,i))> threshold_spk && max(norm_sound(:,i-1))< threshold_spk%find where sounds occur across all speakers
            s=s+1;
%         elseif max(norm_sound(:,i))< threshold_spk && max(norm_sound(:,i-1))> threshold_spk && max(norm_sound(:,i+1))> threshold_spk 
%             % to deal with really low values during reward/incorrect (it's
%             % usually one value between two good ones)
%             s=s+1;
        end
            sounds1(i)=s;
    end  
    x=find(max(norm_sound)<threshold_spk);
    sounds1(x)=0;
    sounds=zeros(size(sounds1));
    
    norm_sounds=norm_sound;

% getting rid of trial sounds that are not fully done!
% assumes the distance between sounds is the same (~200ms) and the
% sounds are ~1s in duration!
diff_sounds = diff(sounds1);
[pks,locs] = findpeaks(abs(diff_sounds));
sound_onset = []; sound_offset = [];
sound_onset = locs(find(diff_sounds(locs)>0));
sound_offset = locs(find(diff_sounds(locs)<0));
    % Print information for debugging
fprintf('Total detected onsets: %d\n', length(sound_onset));
fprintf('Total detected offsets: %d\n', length(sound_offset));
 % Find unpaired onsets and offsets
unpaired_onset = setdiff(sounds1(sound_onset + 1), sounds1(sound_offset)); % onset not in offset
unpaired_offset = setdiff(sounds1(sound_offset), sounds1(sound_onset + 1)); % offset not in onset

onset_index = setdiff(sounds1(sound_onset+1),unpaired_onset);
offset_index = setdiff(sounds1(sound_offset),unpaired_offset);
ss = 0; sound_pairs = [];true_sound_pairs = [];
for s = 1:min([length(onset_index),length(offset_index)])%min([length(sound_onset),length(sound_offset)])

    sound_pairs(s,:) = [sound_onset(sounds1(sound_onset+1) == onset_index(s));sound_offset(sounds1(sound_offset) == offset_index(s))];
        if sound_pairs(s,1) < size(norm_sounds,2) && sound_pairs(s,2) < size(norm_sounds,2) %unsure why this would be bigger than array but does happen
            ss = ss+1;
            true_sound_pairs(ss,:) = sound_pairs(s,:);
        else
            sound_pairs(s,:) = nan;
        end
end
sound_duration = [0.99*sync_sampling_rate,1.1*sync_sampling_rate];%[0.98*sync_sampling_rate,1.02*sync_sampling_rate];

% Print the valid pairs for debugging
fprintf('Number of valid pairs: %d\n', length(true_sound_pairs));
disp(true_sound_pairs);

%distance_between = 0.22; %distance between sounds within a trial
    
difference = [true_sound_pairs(:,2) - true_sound_pairs(:,1)]; 
all_trial_sounds = [];
all_trial_sounds = true_sound_pairs(find(difference >sound_duration(1) & difference < sound_duration(2)),:); %sounds that are outside limits of sound duration
% adding code to also include sounds that are cut off early
count = 0; unfinished_sounds = [];
unfinished_sounds = setdiff(1:length(difference),find(difference >sound_duration(1) & difference < sound_duration(2)));
if ~isempty(unfinished_sounds)
    for es = 1:length(unfinished_sounds)
        extra_sound = unfinished_sounds(es);
        if extra_sound > 1 && extra_sound<unfinished_sounds(end) && [sound_pairs(extra_sound,1) - sound_pairs(extra_sound-1,2)] < distance_within_sounds && [sound_pairs(extra_sound,1) - sound_pairs(extra_sound-1,2)] > (distance_within_sounds-(distance_within_sounds*.15))
            count = count+1;
             %use the diff to see if you can find offset of unfinished sound
            sound_derivative = diff(max(norm_sound(:,sound_pairs(extra_sound,1):sound_pairs(extra_sound,2))));
            % was using first zero but it is not always accurate
            below_change =  find(diff(sound_derivative)==0);%find(abs(sound_derivative)<1e-6);
            %find first one that is close to other ones
            
            new_offset = find(diff(sound_derivative)==0,1,'first');%find(sound_derivative==0,1,'first'); 
            if abs(sound_derivative(new_offset+1))> 1e-6
                new_offset = below_change(2);
                fprintf(num2str(sound_pairs(extra_sound,1)));
                unfinished_sounds(es)
            end
            unfinished_sounds_toadd(count,:) = [sound_pairs(extra_sound,1),(new_offset+sound_pairs(extra_sound,1)-1)];
        end
    end
    all_trial_sounds = sort([all_trial_sounds; unfinished_sounds_toadd]);
end

condition_array = []; onset_array = []; offset_array= []; tt = 0;sound_struc = {};
  for t = 1:length(all_trial_sounds)
        current_sound = all_trial_sounds(t,1):all_trial_sounds(t,2);
        
            if size(norm_sounds,1) == 2 %eventually need to adapt to more speakers
                if any(norm_sounds(1,current_sound)>threshold_spk) && max(norm_sounds(2,current_sound))<threshold_spk %speaker 1
                    sounds(current_sound)=1;
    
                elseif any(norm_sounds(2,current_sound)>threshold_spk) && max(norm_sounds(1,current_sound))<threshold_spk %spekear 2
                    sounds(current_sound)=2;
                end 
            end
         if sounds(all_trial_sounds(t,1)+2)>0 && sounds(all_trial_sounds(t,2))>0
             tt = tt+1;
            sound_struc(tt).onset = all_trial_sounds(t,1);
            sound_struc(tt).offset = all_trial_sounds(t,2);
            sound_struc(tt).condition = sounds(all_trial_sounds(t,1)+2);
            condition_array = [condition_array,sound_struc(tt).condition];
            onset_array = [onset_array,all_trial_sounds(t,1)];
            offset_array = [offset_array,all_trial_sounds(t,2)];
        end
  end


j=0;
for ii=1:length(onset_array)
    if sound_struc(ii).condition > 0
        j = j+1;
        %FIND FRAMES THAT MATCH THIS
        [x,y] = min(abs(sound_struc(j).onset - frames_times{1,file}));
        sound_struc(j).onset_frames = y;
        [x2,y2] = min(abs(sound_struc(j).offset - frames_times{1,file}));
        sound_struc(j).offset_frames = y2;
%         sound_struc(j).onset_frames_ad = y+frame_added;
%         sound_struc(j).offset_frames_ad = y2+frame_added; %adding number of frames from previous files to concatenate
    end
end

numSounds = length(onset_array);
numConditions = length(unique(condition_array));
timeThreshold = distance_between_sounds;
expectedDistance = distance_within_sounds; %distance_between*sync_sampling_rate; 

% Initialize cell arrays to store trial information for each condition
trialsPerCondition = cell(numConditions, 1);
groupNum = 0; % Initialize group number
onsettimes= [];
condition_group_array = [];
% Iterate over each sound
for i = 1:numSounds
    onsetTime = onset_array(i);
    offsetTime = offset_array(i);
    condition = condition_array(i);
    % Initialize or get the current condition's trials
    if isempty(trialsPerCondition{condition})
        trialsPerCondition{condition} = [onsetTime, offsetTime];
    else
        lastTrial = trialsPerCondition{condition}(end, :);
        lastOffsetTime = lastTrial(2);
        
        % Check if the current sound can be added to the last trial
        if abs(lastOffsetTime - onsetTime) <= expectedDistance
            % Extend the last trial
            trialsPerCondition{condition}(end, 2) = offsetTime;
        else
            % Create a new trial
            trialsPerCondition{condition} = [trialsPerCondition{condition}; [onsetTime, offsetTime]];
        end
    end
    % add to sound struct
    if i > 1 && onsetTime - offset_array(i - 1) <= distance_between_sounds
        % Assign the same group number as the previous sound
        sound_struc(i).trial_num = sound_struc(i - 1).trial_num;
        
    else
        % Increment the group number and assign to the current sound
        groupNum = groupNum + 1;
        sound_struc(i).trial_num = groupNum;
        onsettimes = [onsettimes,onsetTime];
        condition_group_array(groupNum, 2) = onsettimes(1);
    end
    
    % Store condition and group number information
    condition_group_array(groupNum, 1) = condition;
    
    condition_group_array(groupNum, 3) = offsetTime;
    condition_group_array(groupNum, 4) = groupNum;
    onsettimes=[];
end



%% Display the resulting trials
figure(111);clf
hold on
plot(sounds, 'b');
plot(rescale(norm_sound(1,:),0,2));
plot(rescale(norm_sound(2,:),0,2));
hold off

figure(110);clf
subplot(2,1,1)
% Plot the sound_array as a binary plot
plot(sounds, 'b');
ylim([-0.5, 2.5]);
yticks([0, 1]);
yticklabels({'No Sound', 'Sound'});
xlabel('Time');
title('Sound Array');

hold on;

% Plot the trials for each condition
for c = 1:numConditions
    currentConditionTrials = trialsPerCondition{c};
    if ~isempty(currentConditionTrials)
        for t = 1:size(currentConditionTrials, 1)
            trial = currentConditionTrials(t, :);
            rectangle('Position', [trial(1), c-0.4, trial(2)-trial(1), 0.8], 'FaceColor', [1, 0, 0, 0.3], 'EdgeColor', 'none');
        end
    end
end

% Adjust y-axis and labels for the conditions
yticks(1:numConditions);
yticklabels({'Condition 1', 'Condition 2'}); % Add labels as needed

hold off;
subplot(2,1,2)
hold on
plot(norm_sound(1,:));
plot(norm_sound(2,:));
hold off
xlabel('Normalized sounds');
title('Sound Array');
legend('Condition 1', 'Condition 2')
%%

sound_outputs(file).file = sound_struc;
trialConditions(file).file = trialsPerCondition;
condition_onset_array_all(file).file = condition_group_array;
pause



end