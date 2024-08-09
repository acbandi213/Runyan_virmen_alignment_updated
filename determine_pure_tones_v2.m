function [condition_onset_array_all_final] = determine_pure_tones_v2(pure_tone_signal,sync_sampling_rate,correct_length,incorrect_length,condition_onset_array_all)
min_distance_between_pure_tones = 9;%in seconds
pos_peaks = []; pos_peaks2 = [];

%new method to find peaks using zeros
zero_values = find(abs(diff(pure_tone_signal(1,:))) == 0);
binarized_zeros = heaviside(rescale(diff(zero_values),0,1)*-1); %
[~,b] = findpeaks(abs(diff(binarized_zeros)));
[distance,location] = findpeaks(diff(b));

correct_length = correct_length*sync_sampling_rate;
incorrect_length = incorrect_length*sync_sampling_rate;
time_after = [0.1*sync_sampling_rate,1*sync_sampling_rate]; %examples were 0.188 0.192 0.193 (incorrect), .523 .664 0.687 correct
target_differences = [correct_length-correct_length*.04:correct_length*.04+correct_length];
target_differences_incorrect = [incorrect_length-incorrect_length*.02:incorrect_length+incorrect_length*.02];

% corr_zero = []; incorr_zero = [];
% for val = 1:length(distance); 
%     if ismember(distance(val),target_differences)
%         corr_zero = [corr_zero,zero_values(b(location(val)))];
%     elseif ismember(distance(val),target_differences_incorrect)
%         incorr_zero = [incorr_zero,zero_values(b(location(val)))];
%     end
% end
count =0; pure_tones_trial2=[];
for trial = 1:length(distance)
    if ismember(distance(trial),target_differences)
        count = count+1;
        pure_tones_trial2(count,:) = [1,zero_values(b(location(trial))), zero_values(b(location(trial)))+correct_length,count];
    elseif ismember(distance(trial),target_differences_incorrect)
        count = count+1;
        pure_tones_trial2(count,:) = [0,zero_values(b(location(trial))), zero_values(b(location(trial)))+incorrect_length,count];
    end
    count
end

%organize data
order1 =[] ; order2 =[];
for t = 1:length(pure_tones_trial2)
    difference = [[condition_onset_array_all.VR_sounds{:,3}] - pure_tones_trial2(t,2)];
    index_diff = [find([[condition_onset_array_all.VR_sounds{:,3}] - pure_tones_trial2(t,2)] < 0)];
    if ~isempty(index_diff) && abs(difference(index_diff(end))) < (sync_sampling_rate* 1.5) %has to be within a second after the sound is played
        order1 = [order1, t]; %good index - sound that can be paired with ITI sound
    else
        order2 =[order2, t]; %nan index
    end
        
end

condition_onset_array_all_final.ITI_sounds = pure_tones_trial2;
%adding this so arrays are the same size!
if length(order1) < length(condition_onset_array_all.VR_sounds) && condition_onset_array_all.VR_sounds{end,3} > pure_tones_trial2(end,2)
    combined_list(order1,:) = [condition_onset_array_all.VR_sounds(1:end-1,:)]; %if there is a sound not paired with pure tune at the end!
    if ~isempty(order2)
        combined_list(order2,:) = {nan};
    end
elseif length(order1) < length(condition_onset_array_all.VR_sounds)  
    combined_list(order1,:) = [condition_onset_array_all.VR_sounds(order1,:)]; %if there is a sound not paired with pure tune at the end!
    if ~isempty(order2)
        combined_list(order2,:) = {nan};
    end
else
    combined_list(order1,:) = condition_onset_array_all.VR_sounds;
    if ~isempty(order2)
        combined_list(order2,:) = {nan};
    end
end

if ~isempty(pure_tones_trial2)
    condition_onset_array_all_final.VR_sounds = combined_list;
end


%%
combined_pure_tones = pure_tones_trial2;
figure(112);clf; hold on; title('ITI Pure tones with NaN')
hold on
plot(diff(pure_tone_signal(1,:)));
plot((pure_tone_signal(1,:)));
plot(combined_pure_tones(:,2),0,'*c')
if length(find(isnan(combined_pure_tones(:,1))))>0  
    
    a = plot([condition_onset_array_all.VR_sounds{(find(isnan(combined_pure_tones_ordered(:,1)))),3}],0,'*m');
end
hold off
legend('Pure Tones Diff','Pure Tones regular','Start Pure Tone','Last sound near NaN')

%make sure they match what you are saying
final_diff = combined_pure_tones(:,3) - combined_pure_tones(:,2);
if all(combined_pure_tones(find(combined_pure_tones(:,3) - combined_pure_tones(:,2)>target_differences_incorrect(1)-1),1) == 0) && all(combined_pure_tones(find(combined_pure_tones(:,3) - combined_pure_tones(:,2)<target_differences_incorrect(1)-1),1) == 1)
    fprintf('Values verified to be correct\n')
else
    fprintf("----------------Values don't match the expected correct/incorrect lengths----------------\n")
    keyboard
end
