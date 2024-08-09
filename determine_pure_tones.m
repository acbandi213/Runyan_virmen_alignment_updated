function [condition_onset_array_all,pure_tones_trial] = determine_pure_tones(pure_tone_signal,sync_sampling_rate,correct_length,incorrect_length,condition_onset_array_all)
% pos_peaks = []; pos_peaks = [];
% for s = 1:size(pure_tone_signal,1)
%     [~,pos_peaks2] = findpeaks(abs(diff(pure_tone_signal(s,:))),'MinPeakHeight', 0.005);
%     pos_peaks = [pos_peaks,pos_peaks2];
% end
% correct_length = correct_length*sync_sampling_rate;
% incorrect_length = incorrect_length*sync_sampling_rate;
% time_after = [0.1*sync_sampling_rate,1*sync_sampling_rate]; %examples were 0.188 0.192 0.193 (incorrect), .523 .664 0.687 correct
% target_differences = [correct_length-correct_length*.025:correct_length*.025+correct_length];
% target_differences_incorrect = [incorrect_length-incorrect_length*.025:incorrect_length+incorrect_length*.025];
% 
% combos = nchoosek(1:length(pos_peaks),2);
% correct_pure_tones = [];incorrect_pure_tones = [];trial_options=[];trial_options2=[]; possible_values = [];pure_tones_trial = [];
% all_differences = abs(pos_peaks(combos(:,2)) - pos_peaks(combos(:,1)));
% in_diff = find(all_differences>target_differences_incorrect(1) & all_differences < target_differences_incorrect(end));
% cor_diff = find(all_differences>target_differences(1) & all_differences < target_differences(end));
% %  pure_tones_trial = zeros(size(condition_onset_array_all.VR_sounds,1),4);
% good_diff= [ in_diff,cor_diff];
% for t = 1:size(condition_onset_array_all.VR_sounds,1)
%     
%     for ii = 1:length(good_diff)
%     
%     possible_values = [];
% 
%             %find differences between each peak value
%             i = good_diff(ii);
%             difference = all_differences(i);
%             % check if difference is equal to target differences
%             if ismember(i,in_diff) && min([pos_peaks(combos(i,1)), pos_peaks(combos(i,2))]) - condition_onset_array_all.VR_sounds(t,3)<time_after(2) && min([pos_peaks(combos(i,1)), pos_peaks(combos(i,2))]) - condition_onset_array_all.VR_sounds(t,3)>time_after(1)
%                 %incorrect_pure_tones = [incorrect_pure_tones;sort([pos_peaks(combos(i,1)), pos_peaks(combos(i,2))])];
% 
%                 trial_options = [trial_options;sort([pos_peaks(combos(i,1)), pos_peaks(combos(i,2))])];
%                 pure_tones_trial(t,:) = [0,trial_options(1,:),t];
%                  
%                 incorrect_pure_tones(t,:) = [trial_options(1,:)];
%             
%             elseif ismember(i,cor_diff) && min([pos_peaks(combos(i,1)), pos_peaks(combos(i,2))]) - condition_onset_array_all.VR_sounds(t,3)<time_after(2) && min([pos_peaks(combos(i,1)), pos_peaks(combos(i,2))]) - condition_onset_array_all.VR_sounds(t,3)>time_after(1) && size(pure_tones_trial,1)==(t-1)
%                 %correct_pure_tones = [correct_pure_tones;sort([pos_peaks(combos(i,1)), pos_peaks(combos(i,2))])];
%                 trial_options2 = [trial_options2;sort([pos_peaks(combos(i,1)), pos_peaks(combos(i,2))])];
%                 pure_tones_trial(t,:) = [1,trial_options2(1,:),t];
%                 correct_pure_tones(t,:) = [trial_options2(1,:)];
%             %end
%             end
%             
%     end
%     trial_options=[];
%     trial_options2=[];
% 
%     %different method of finding pure tones without peaks
%     %relies on the fact that between pure tones the signal =0
%     if size(pure_tones_trial,1) < t 
%         if t == size(condition_onset_array_all.VR_sounds,1)
%             pure_tones_trial(t,:) = [nan,nan,nan,t];
%         else
%             dist = find(abs(diff(pure_tone_signal(1,condition_onset_array_all.VR_sounds(t,3)+1:condition_onset_array_all.VR_sounds(t+1,3))))>0);
%             distance = diff(dist);
%             a = find(distance>target_differences(1) & distance < target_differences(end));
%             a2 = find(distance>target_differences_incorrect(1) & distance < target_differences_incorrect(end));
%             if ~isempty(a)
%                 pure_tones_trial(t,:) = [1,dist(a(1))+condition_onset_array_all.VR_sounds(t,3)+1, dist(a(1)+1)+condition_onset_array_all.VR_sounds(t,3)+1,t];
%             elseif ~isempty(a2)
%                 pure_tones_trial(t,:) = [0,dist(a2(1))+condition_onset_array_all.VR_sounds(t,3)+1, dist(a2(1)+1)+condition_onset_array_all.VR_sounds(t,3)+1,t];
%             else
%                 pure_tones_trial(t,:) = [nan,nan,nan,t];
%             end
%         end
%         
%     end
%     condition_onset_array_all.ITI_sounds(t,:) = pure_tones_trial(t,:);
% end
% 
% figure(112);clf; hold on; title('ITI Pure tones with NaN')
% hold on
% plot(diff(pure_tone_signal(1,:)));
% plot((pure_tone_signal(1,:)));
% if length(find(isnan(pure_tones_trial(:,1))))>0  
%     
%     plot(condition_onset_array_all.VR_sounds([find(isnan(pure_tones_trial(:,1)))],3),0,'*m');
% end
% hold off
% legend('Pure Tones Diff','Pure Tones regular','Last sound near NaN')

%% 3rd try
%% change code to adjust when there are speakers missing but still want the ITI
min_distance_between_pure_tones = 9;%in seconds
pos_peaks = []; pos_peaks2 = [];

%reprocessing signal
%binary_pure_signal = heaviside(abs(pure_tone_signal)*-1);
for s = 1:size(pure_tone_signal,1)
    [~,pos_peaks2] = findpeaks(abs(diff(pure_tone_signal(s,:))),'MinPeakHeight', 0.006);%findpeaks(abs(diff(binary_pure_signal(s,:))),'MinPeakHeight', 0.005);%
    pos_peaks = [pos_peaks,pos_peaks2];
end
pos_peaks = unique(pos_peaks);
correct_length = correct_length*sync_sampling_rate;
incorrect_length = incorrect_length*sync_sampling_rate;
time_after = [0.1*sync_sampling_rate,1*sync_sampling_rate]; %examples were 0.188 0.192 0.193 (incorrect), .523 .664 0.687 correct
target_differences = [correct_length-correct_length*.04:correct_length*.04+correct_length];
target_differences_incorrect = [incorrect_length-incorrect_length*.02:incorrect_length+incorrect_length*.02];

combos = nchoosek(1:length(pos_peaks),2);
pure_tones = [];trial_options=[];trial_options2=[]; possible_values = [];pure_tones_trial = [];
pure_tones2 = [];pure_tones2 = [];pure_tones_trial2 = [];

all_differences = abs(pos_peaks(combos(:,2)) - pos_peaks(combos(:,1))); 
in_diff = find(all_differences>target_differences_incorrect(1) & all_differences < target_differences_incorrect(end));
cor_diff = find(all_differences>target_differences(1) & all_differences < target_differences(end));
%  pure_tones_trial = zeros(size(condition_onset_array_all.VR_sounds,1),4);

% write code to get rid of correct peaks if the incorrect peaks are really
% nearby
corr_peakss = pos_peaks(combos(cor_diff,2));
incorr_peakss = pos_peaks(combos(in_diff,2));
corr_to_eliminate = [];
for ii = 1:length(incorr_peakss)
    corr_to_eliminate=[corr_to_eliminate,find(abs(corr_peakss - incorr_peakss(ii))<.2*sync_sampling_rate)]; %if peaks are within 200ms of each other probably inside the same ITI window
end
corr_to_eliminate = unique(corr_to_eliminate);
cor_diff(corr_to_eliminate) = [];

%new method to find peaks using zeros
zero_values = find(abs(diff(pure_tone_signal(s,:))) == 0);
binarized_zeros = heaviside(rescale(diff(zero_values),0,1)*-1); %
[~,b] = findpeaks(abs(diff(binarized_zeros)));
[distance,location] = findpeaks(diff(b));

corr_zero = []; incorr_zero = [];
for val = 1:length(distance); 
    if ismember(distance(val),target_differences)
        corr_zero = [corr_zero,zero_values(b(location(val)))];
    elseif ismember(distance(val),target_differences_incorrect)
        incorr_zero = [incorr_zero,zero_values(b(location(val)))];
    end
end

good_diff=sort([ in_diff,cor_diff]); 
for t = 1:size(condition_onset_array_all.VR_sounds,1)
    
    if t == size(condition_onset_array_all.VR_sounds,1)
            for ii = 1:length(good_diff)
                possible_values = [];
                            %find differences between each peak value
                    i = good_diff(ii);
                    difference = all_differences(i);
                    % check if difference is equal to target differences
                    if min([pos_peaks(combos(i,1)), pos_peaks(combos(i,2))]) < length(pure_tone_signal) && min([pos_peaks(combos(i,1)), pos_peaks(combos(i,2))]) - condition_onset_array_all.VR_sounds{t,3}>time_after(1)
                        trial_options = [trial_options;sort([pos_peaks(combos(i,1)), pos_peaks(combos(i,2))])];
                        if ismember(i,cor_diff) && trial_options(1,1) == pos_peaks(combos(i,1)) && trial_options(1,2) == pos_peaks(combos(i,2))
                            pure_tones_trial(t,:) = [1,trial_options(1,:),t];
                        elseif ismember(i,in_diff) && trial_options(1,1) == pos_peaks(combos(i,1)) && trial_options(1,2) == pos_peaks(combos(i,2))
                            pure_tones_trial(t,:) = [0,trial_options(1,:),t];
                        end  

                        pure_tones = [pure_tones;trial_options(1,:)];
                        for tt = 2:size(trial_options,1)
                            if size(trial_options,1)>1 && trial_options(tt,1) - max(pure_tones(:,2))> min_distance_between_pure_tones*sync_sampling_rate %if there are multiple peaks and they are at least 1 sec apart
                                if ismember(i,cor_diff)
                                    pure_tones_trial2 = [pure_tones_trial2;1,trial_options(tt,:),t];
                                else
                                    pure_tones_trial2 = [pure_tones_trial2;0,trial_options(tt,:),t];
                                end
                                pure_tones = [pure_tones;trial_options(tt,:)];
                            end
                        end
                    end
                    
            end
            trial_options=[];
            trial_options2=[];

    else
        
          for ii = 1:length(good_diff)
            
%          if t == 11 
%             t
%         end
            possible_values = [];
            

                    %find differences between each peak value
                    i = good_diff(ii);
                    %[pos_peaks(combos(i,1)), pos_peaks(combos(i,2))]
                    difference = all_differences(i);
                    % check if difference is equal to target differences
                    if min([pos_peaks(combos(i,1)), pos_peaks(combos(i,2))]) < condition_onset_array_all.VR_sounds{t+1,3} && min([pos_peaks(combos(i,1)), pos_peaks(combos(i,2))]) - condition_onset_array_all.VR_sounds{t,3}>time_after(1)
                        trial_options = [trial_options;sort([pos_peaks(combos(i,1)), pos_peaks(combos(i,2))])];
                        
                        if ismember(i,cor_diff) && trial_options(1,1) == pos_peaks(combos(i,1)) && trial_options(1,2) == pos_peaks(combos(i,2))
                            pure_tones_trial(t,:) = [1,trial_options(1,:),t];
                        elseif ismember(i,in_diff) && trial_options(1,1) == pos_peaks(combos(i,1)) && trial_options(1,2) == pos_peaks(combos(i,2))
                            pure_tones_trial(t,:) = [0,trial_options(1,:),t];
                        end    
                        pure_tones = [pure_tones;trial_options(1,:)];
                        %this loop to deal with multiple peaks inside pure
                        %tone - it is likely peaks could find a correct
                        %even though it is incorrect just bc there are so
                        %many peaks
                        for tt = 2:size(trial_options,1)
                            if size(trial_options,1)>1 && trial_options(tt,1) - max(pure_tones(:,2))> min_distance_between_pure_tones*sync_sampling_rate %if there are multiple peaks and they are at least 1 sec apart
                                if ismember(i,cor_diff)
                                    pure_tones_trial2 = [pure_tones_trial2;1,trial_options(tt,:),t];
                                else
                                    pure_tones_trial2 = [pure_tones_trial2;0,trial_options(tt,:),t];
                                end
                                pure_tones = [pure_tones;trial_options(tt,:)];
                            end
                            
                        end
                        
                    end
                    
            end
            trial_options=[];
            trial_options2=[];

    end
    
    %different method of finding pure tones without peaks
    %relies on the fact that between pure tones the signal =0
    if size(pure_tones_trial,1) < t 
        if t == size(condition_onset_array_all.VR_sounds,1)
            pure_tones_trial(t,:) = [nan,nan,nan,t];
        else
            dist = find(abs(diff(pure_tone_signal(1,condition_onset_array_all.VR_sounds{t,3}+1:condition_onset_array_all.VR_sounds{t+1,3})))>0);
            distance = diff(dist);
            a = find(distance>target_differences(1) & distance < target_differences(end));
            if ~isempty(a)
            a = distance(a);
            end
            a2 = find(distance>target_differences_incorrect(1) & distance < target_differences_incorrect(end));
            if ~isempty(a2)
            a2 = distance(a2);
            end
            all_diff = [a;a2];
            if ~isempty(all_diff)
                for tt = 1:length(all_diff)
                    if ismember(all_diff (tt),a)
                        pure_tones_trial(t,:) = [1,dist(find(distance == all_diff(1)))+condition_onset_array_all.VR_sounds{t,3}+1, dist(find(distance == all_diff(1))+1)+condition_onset_array_all.VR_sounds{t,3}+1,t];
                    else
                        pure_tones_trial(t,:) = [0,dist(find(distance == all_diff(1)))+condition_onset_array_all.VR_sounds{t,3}+1, dist(find(distance == all_diff(1))+1)+condition_onset_array_all.VR_sounds{t,3}+1,t];
                    end
                    if length(all_diff)>1 && dist(find(distance == all_diff(tt)))+condition_onset_array_all.VR_sounds{t,3}+1 -  dist(find(distance == all_diff(tt-1))+1)+condition_onset_array_all.VR_sounds{t,3}+1> min_distance_between_pure_tones*sync_sampling_rate %if there are multiple peaks and they are at least 1 sec apart
                       if ismember(all_diff (tt),a)
                            pure_tones_trial2(t,:) = [1,dist(find(distance == all_diff(tt)))+condition_onset_array_all.VR_sounds{t,3}+1, dist(find(distance == all_diff(tt+1)))+condition_onset_array_all.VR_sounds{t,3}+1,t];
                        else
                            pure_tones_trial2(t,:) = [0,dist(find(distance == all_diff(tt)))+condition_onset_array_all.VR_sounds{t,3}+1, dist(find(distance == all_diff(tt+1)))+condition_onset_array_all.VR_sounds{t,3}+1,t];
                        end
                        
                    end
                end
            else
                pure_tones_trial(t,:) = [nan,nan,nan,t];
            end
        end
        
    end
%         condition_onset_array_all.ITI_sounds(t,:) = pure_tones_trial(t,:);
end
% fine tune pure tones
% [~,pure_tone_peaks] = findpeaks(abs(actual_sounds(1,:)),'minpeakheight',0.1);
pure_tones_together = [pure_tones_trial;pure_tones_trial2];
% for p = 1:length(pure_tones_together)
%     temp = find(pure_tone_peaks>pure_tones_together(p,2) & pure_tone_peaks<pure_tones_together(p,3));
%     if ~isempty(temp)
%         onset_offset(p,:) = [pure_tone_peaks(temp(1)), pure_tone_peaks(temp(end))];
%     else
%         onset_offset(p,:) = [nan, nan];
%     end
% end
% pure_tones_together(:,2:3) = onset_offset;
%organize data
[combined_pure_tones,order] = sortrows(pure_tones_together,2);
[combined_pure_tones_ordered,order2] = sortrows(pure_tones_together,4);
condition_onset_array_all.ITI_sounds = combined_pure_tones_ordered;
%adding this so arrays are the same size!
nan_array = cell(size(pure_tones_trial2,1),size(pure_tones_trial2,2)+1);
nan_array(1:size(pure_tones_trial2,1),1:size(pure_tones_trial2,2)+2) = {nan};
combined_list = [condition_onset_array_all.VR_sounds;nan_array];%nan(size(pure_tones_trial2))

if ~isempty(pure_tones_trial2)
    condition_onset_array_all.VR_sounds = combined_list(order2,:);
end


% % %% adding code to deal with extra pure tones
% % condition_onset_array_all_corrected = condition_onset_array_all;
% %  to_add = [];
% % if ~isempty(pure_tones_trial2)
% %     to_add = find(pure_tones_trial2(:,2));
% % end
% % tt=0; t2 = 0;
% % for t = 1:(size(condition_onset_array_all.VR_sounds,1))
% %     t2 = t2+1;
% %     tt = tt+1;
% %     if ~any(t == to_add) % has to be t2 = 9
% %         condition_onset_array_all_corrected.VR_sounds(tt,:) = condition_onset_array_all.VR_sounds(t,:);
% %         condition_onset_array_all_corrected.ITI_sounds(tt,:) = pure_tones_trial(t,:);
% %     else 
% %         t2=t2+1;
% %         tt = tt+1;%find((t  == to_add));
% %         condition_onset_array_all_corrected.ITI_sounds(tt,:) = pure_tones_trial2(t,:); %tt has to be 11
% %         condition_onset_array_all_corrected.ITI_sounds(tt,4) = tt;
% %         condition_onset_array_all_corrected.VR_sounds(tt,:) = [nan,nan,nan,tt];
% %     end
% %     
% % end
%%
% condition_onset_array_all_corrected = condition_onset_array_all;
% 
% if ~isempty(pure_tones_trial2)
%     to_add = find(pure_tones_trial2(:,2));
% end
% %define insertion positions
% insertion_positions = to_add;
% %initialize shift to keep track of the shift due to the insertions
% shift = 0; count=0;
% %perform insertions at specified positions
% for i =1:length(insertion_positions)
%     count = count+1;
%     insertion_positions = to_add(i)+shift;
%     %split data into two parts
%     first_part = 1:insertion_positions;
%     second_part = insertion_positions+1:size(condition_onset_array_all_corrected.VR_sounds,1);
%     condition_onset_array_all_corrected.VR_sounds = [condition_onset_array_all_corrected.VR_sounds(first_part,:);[nan,nan,nan,insertion_positions];condition_onset_array_all_corrected.VR_sounds(second_part,:)];
%     condition_onset_array_all_corrected.ITI_sounds =[condition_onset_array_all_corrected.ITI_sounds(first_part,:);pure_tones_trial2(to_add(i),:);condition_onset_array_all_corrected.ITI_sounds(second_part,:)];
% 
%     shift = shift+count;
% end
%%
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





