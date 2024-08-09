function [trial_its,trial_its_time] = virmen_it_rough_estimation_adjusted(data,sampling_rate,it_times,it_values)
% use virmen iterations to get rough estimate of which iterations are
% specific trial events based on stereotyped gaps between them!
% also find iteration where mouse crosses 50 in the y position!- assumes a
% threshold of 350- takes more than these # iterations to do a trial

iterations_in_time = it_times;


%look at iti to determine isITI and also start of trial
start_iti_its = find(diff(data.data(9,:))>0);
start_iti_its = start_iti_its + 2; %+2 bc it seems one earlier than it actually should be
end_iti_its = find(diff(data.data(9,:))<0);
start_trial_its = end_iti_its + 3;
start_trial_its = sort([start_trial_its,1]); %making sure I have an it for the very first trial
end_trial_its = start_iti_its - 1;
%get complete trials
min_trials = min([length(end_trial_its),length(start_trial_its)]);
if length(start_trial_its) > min_trials
    start_trial_its(end) = [];
end


%check with virmen iterations to see if the gaps are in the right place
gaps = find(diff(iterations_in_time) > .7*sampling_rate); %might need to adjust this threshold
gap_iterations = it_values(gaps);
iteration_threshold = 8; %3 was before - used to move ITI iterations to the next gap

for g = 1:length(gaps)
    
    [val,current_gap] = min(abs(gap_iterations(g) - end_trial_its));

    if val < iteration_threshold %&& gap_iterations(g) - end_trial_its(current_gap) >= 0
        end_trial_its(current_gap) = gap_iterations(g);
    end
end

gaps = find(diff(iterations_in_time) > .07*sampling_rate & diff(iterations_in_time) < .7*sampling_rate); %greater than 80ms apart
gap_iterations = it_values(gaps);
for g = 1:length(gaps)
   
    [val,current_gap] = min(abs(gap_iterations(g) - end_iti_its));

    if val < iteration_threshold && gap_iterations(g) - end_iti_its(current_gap) >= 0
        end_iti_its(current_gap) = gap_iterations(g);
    end
end

start_iti_its = zeros(1,length(start_iti_its));
start_iti_its = end_trial_its(1:length(start_iti_its)) +1;

start_trial_its = zeros(1,length(start_trial_its)-1); %one less to add 1 to the front
start_trial_its = sort([0,end_iti_its(1:(length(start_trial_its)-1))])+1;

%check to make sure start_trial_its has a one in front
if isempty(find(start_trial_its == 1))
    start_trial_its = sort([1,start_trial_its]);
elseif length(find(start_trial_its == 1)) > 1
    start_trial_its = unique(start_trial_its );
end

%get complete trials
min_trials = min([length(end_trial_its),length(start_trial_its)]);
if length(start_trial_its) > min_trials
    start_trial_its(end) = [];
end

%add to strcuture
trial_its.end_trial_its = end_trial_its;
trial_its.start_iti_its = start_iti_its;
trial_its.end_iti_its = end_iti_its;
trial_its.start_trial_its = start_trial_its;

%% looking at iteration timing info! try to tie it back to the big sterotyped gaps that I see!
it_times = data.data(1,:).* 86400; %convert to seconds
it_time_gaps = diff(it_times);

%look at iti to determine isITI and also start of trial in iteration times
%did not use this really?
end_trial_its_time = it_times(find(it_time_gaps>.8));
start_iti_its_time = it_times(find(it_time_gaps>.8)+1);
incorrect_its_time = it_times(find(it_time_gaps >.8 & it_time_gaps < .95));
correct_its_its_time = it_times(find(it_time_gaps > .95));
temp_its = find(it_time_gaps > .25 & it_time_gaps < .51)+1;
start_trial_its_time = it_times([1,temp_its]); %adding one for very first trial

trial_its_time.start_trial_its_time = start_trial_its_time;
trial_its_time.end_trial_its_time = end_trial_its_time;
trial_its_time.start_iti_its_time = start_iti_its_time;
trial_its_time.correct_its_its_time = correct_its_its_time;
trial_its_time.incorrect_its_time = incorrect_its_time;

%look at sound onset first
temp_onset = sort([find(round(floor(data.data(3,:))) == 50),find(round(floor(data.data(3,:))) == 51),find(round(floor(data.data(3,:))) == 52)]);
sound_trigger_its = [];
% keep only the first one and delete the rest
for t = 1:min([length(trial_its.start_trial_its),length(trial_its.end_trial_its)])
    temp = temp_onset(find(temp_onset >trial_its.start_trial_its(t) & temp_onset <trial_its.end_trial_its(t),1,'first'));
    sound_trigger_its = [sound_trigger_its,temp];
end

% %find ones that are too close together and delete them
% temp = find(diff(temp_onset)< threshold);
% temp_onset(temp+1) = [];
% sound_trigger_its = temp_onset;
% %make sure there is one trigger per trial! finding repeats getting rid of
% %anything after first one
% temp_sound_trigger = [];
% for t = 1:min([length(trial_its.start_trial_its),length(trial_its.end_trial_its)])
%     find(abs(trial_its.start_trial_its(t) - sound_trigger_its) == min(abs(trial_its.start_trial_its(t) - sound_trigger_its))); 
%     temp_sound_trigger = [temp_sound_trigger, t];
% end
% [U, I] = unique(temp_sound_trigger, 'first');
% x = 1:length(temp_sound_trigger);
% x(I) = [];
% sound_trigger_its(x) = [];% getting rid of repeats after first

trial_its.sound_trigger_its = sound_trigger_its;