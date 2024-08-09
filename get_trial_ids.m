function [file_trial_ids,file_estimated_trial_info_updated] = get_trial_ids(file_matching_trials,file_estimated_trial_info,alignment_info,sync_base_path,task_info)
%%%Output: file_trial_inds(starting full trial across all files,ending full
%%%trial across all files,starting full trial in this file,ending full
%%%trial in this file,yes/no preivous ITI within imaging,yes/no next trial
%%%within imaging)
file_trial_ids = zeros(size(file_matching_trials,1),6);
for file = 1:size(file_matching_trials,1)
%% 1) find the first trial that would have imaging data and start with this trial if it is a full trial

trial_id = file_estimated_trial_info(file).trial_id; %within each file
trial_id(2,:) = [file_matching_trials(file,1):file_matching_trials(file,1)+trial_id(end)-1]; %relative to all other files

%use frame times to determine which trials should be included
frame_start_trial = find(file_estimated_trial_info(file).start_trials_digidata_time - alignment_info(file).frame_times(1)>0);
frame_end_ITI = find(file_estimated_trial_info(file).end_iti_digidata_time - alignment_info(file).frame_times(end)<0);

all_trials = min(frame_start_trial):1:max(frame_end_ITI);
excluded_trials = setdiff(1:length(trial_id),all_trials);% sort(union(find(isnan(frame_start_trial)),find(isnan(frame_end_ITI)))); %no frames from start trial to end ITI within a trial
included_trials = setdiff(all_trials,excluded_trials);

start_trial(1,1) = trial_id(2,find(trial_id(1,:) == included_trials(1)))+1; %+1 bc it is based on the previous ITI indicating there was at least one trial before
start_trial(1,2) = trial_id(1,find(trial_id(1,:) == included_trials(1))); %keep original number for plotting

if ~isempty(trial_id(2,find(trial_id(1,:) == included_trials(end)))) %empty when the number is bigger than max
    end_trial(1,1) = trial_id(2,find(trial_id(1,:) == included_trials(end)));
    end_trial(1,2) = trial_id(1,find(trial_id(1,:) == included_trials(end)));
else
    end_trial(1,1) = trial_id(2,end);
    end_trial(1,2) = trial_id(1,end);
end
file_trial_ids(file,1:4) = [start_trial(1,1),end_trial(1,1),start_trial(1,2),end_trial(1,2)];
file_estimated_trial_info_updated = file_estimated_trial_info;

%include iti before start trial and full maze after end iti if they are
%within imaging limits
frame_end_maze= [];
frame_start_trial_iti = find(file_estimated_trial_info(file).start_iti_digidata_time(start_trial(1,2)) - alignment_info(file).frame_times(1)>0);
if length(file_estimated_trial_info(file).end_trials_digidata_time) > end_trial(1,2)
    frame_end_maze = find(file_estimated_trial_info(file).end_trials_digidata_time(end_trial(1,2)+1) - alignment_info(file).frame_times(end)<0);
end
if ~isempty(frame_start_trial_iti)
    file_trial_ids(file,5) = [1];
end
if ~isempty(frame_end_maze)
    file_trial_ids(file,6) = [1];
end

%load example file to plot it
ex_data = abfload([sync_base_path alignment_info(file).sync_id]);
figure(55);clf;
title(strcat('Check start and end trials in file # ', num2str(file)));
hold on;aa = plot(ex_data(:,task_info.channel_number(1)));bb = plot(ex_data(:,task_info.channel_number(2)),'color',[0.7 0.7 0.7]);  cc = plot(rescale(ex_data(:,task_info.channel_number(3)),-1,0),'-b');

a = plot(file_estimated_trial_info(file).start_trials_digidata_time(start_trial(1,2)),0,'*c');b = plot(file_estimated_trial_info(file).end_iti_digidata_time(end_trial(1,2)),0,'*g');movegui(gcf,'center');
%plot(rescale(ex_data(:,task_info.channel_number(3)),-1,0),'-r');
legend([aa bb cc a(1) b(1) ],'Imaging frames','Virmen its','Speaker 1', 'first full trial', 'last full trial')
if length(task_info.channel_number)>3
    dd = plot(rescale(ex_data(:,task_info.channel_number(4)),-1,0),'-m');
elseif length(task_info.channel_number)>4
    plot(rescale(ex_data(:,task_info.channel_number(5)),-1,0),'-r')
end
hold off;

pause


end
%figure(); hold on;plot(ex_data(:,6)); plot(rescale(ex_data(:,4),-1,0));plot(possible_it_times,0,'*c');hold off; movegui(gcf,'center');

