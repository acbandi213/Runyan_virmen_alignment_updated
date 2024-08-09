%% provide all inputs
info.mousename = 'HA11-1R';%;
info.mouse = info.mousename;
info.date = '2023-05-05'; %;
info.server = 'V:'; %/Volumes/Runyan5
runyan5 = "V:";
runyan4 = 'W:';
data_base = 'CBHA11-1R_230505';%;
info.experimenter_name = 'Connie';
info.sync_base_path = [ info.server '/' info.experimenter_name '/RawData/' info.mousename '/wavesurfer/' info.date '/'];
info.virmen_base = [info.server '/' info.experimenter_name '/RawData/' info.mousename '/virmen/' data_base ];
info.imaging_base_path=[info.server '/' info.experimenter_name '/RawData/' info.mousename '/' info.date '/'];
info.save_path = [info.server '/' info.experimenter_name '/ProcessedData/' info.mousename '/' info.date '/VR/'];
info.processed_path = [info.server '/' info.experimenter_name '/ProcessedData/' info.mousename '/' info.date '/'];

% info.sync_base_path = [ info.server '/Akhil/' info.mousename '/' info.date '/layer4/wavesurfer/'];
% info.virmen_base = [ info.server '/Akhil/' info.mousename '/' info.date '/layer4/virmen/' data_base];
% info.imaging_base_path=[info.server '/Akhil/' info.mousename '/' info.date '/layer4/'];
% info.save_path = [info.server '/Akhil/ProcessedData/' info.mousename '/' info.date '/layer4/'];
% info.processed_path = [info.server '/Akhil/' info.mousename '/' info.date '/layer4/spikes/'];

info.is_stim_dataset = 1; 
% give data inputs!
info.galvo_channel = 7;
info.virmen_channel = 6;
info.vr_sync_string = 'VR'; %string within VR files to look for so it doesn't read other unnecessary files

sound_info = {};
sound_info.spkr_channel_number = [4,8];%[4,5,8];
sound_info.speaker_ids = [1,2];%[1,2,4]; 
sound_info.mult_spkr = 0; %if multiple speakers are used in a single trial (8 locs)
%load conditions per speaker in runyan 5 info.server
load(strcat(runyan5,'/Connie/condition_per_speaker'));
%load('/Volumes/Runyan5/Connie/condition_per_speaker.mat');
sound_info.condition_per_speaker = conditions_per_speaker;

task_info.correct = 3; %correct ITI time in sec (distance between end trial and end ITI is ITI +1 sec due to virmen bug)
task_info.incorrect = 5; %incorrect ITI time in sec
task_info.min = 3.5;% minimum time in sec to complete a trial (without ITI)
task_info.channel_number = [info.galvo_channel,info.virmen_channel,sound_info.spkr_channel_number];

%reward info from runyan 4 server
load(strcat(runyan4,'/Connie/extra_tests/2023-12-11/mdl_pure_sol.mat'));
load(strcat(runyan4,'/Connie/extra_tests/2023-12-11/mdl_end_trial_sol.mat'));

code_folder = uigetdir; 
addpath(genpath(code_folder));
% cd('C:\Code\Github\Alternative_virmen_alignment');
% addpath(genpath('C:\Code\Align_signals_imaging'));
%% load imaging data
load(strcat(info.processed_path,'dff.mat'));
load(strcat(info.processed_path,'deconv/deconv.mat'));
% dff = zeros(10,15e4);
% deconv = dff;
%% load virmen data
if isfile(strcat(info.virmen_base, '_Cell_2.mat'))
    data = load(strcat(info.virmen_base, '_2.mat'));
    dataCell = load(strcat(info.virmen_base, '_Cell_2.mat'));
elseif isfile(strcat(info.virmen_base, '_Cell_1.mat'))
    data = load(strcat(info.virmen_base, '_1.mat'));
    dataCell = load(strcat(info.virmen_base, '_Cell_1.mat'));
else
    data = load(strcat(info.virmen_base, '.mat'));
    dataCell = load(strcat(info.virmen_base, '_Cell.mat'));
end
fprintf(['number virmen trials: ' num2str(length(dataCell.dataCell)) '\n'])
%% get frame times of all files in this folder
mkdir(strcat(info.server,'/', info.experimenter_name, '/ProcessedData/',num2str(info.mouse),'/',num2str(info.date)))
cd(strcat(info.server,'/', info.experimenter_name, '/ProcessedData/',num2str(info.mouse),'/',num2str(info.date)));
if isfile("alignment_info.mat")
    load("alignment_info.mat");
else
    [alignment_info] = get_frame_times(info.imaging_base_path, info.sync_base_path, [], info.galvo_channel,1,[],[]); %7 is res galvo channel in investigator
end
save ('alignment_info','alignment_info');
%% find sound info for each file
sound_info.sync_sampling_rate = alignment_info(1).sync_sampling_rate;
sound_info.distance_between_sounds = 2*sound_info.sync_sampling_rate ;%min distance between sounds in digidata units in task is ~4 seconds between reward sound and start of trial- (passive at 10k was about 45000)
sound_info.distance_within_sounds = 0.2*sound_info.sync_sampling_rate; %for task should be 200
sound_info.sound_duration = 1*sound_info.sync_sampling_rate;%[0.99*sound_info.sync_sampling_rate,1.1*sound_info.sync_sampling_rate];
sound_info.correct = .250; %correct_trial_ITI_length in seconds
sound_info.incorrect = .40; %incorrect_trial_ITI_length in seconds
sound_info.smoothing_factor = 15; %almost always 15 sometimes 20

sound_info.unique_detection_threshold = [];%list specific file and threshold wanted [file#1,threshold1; file#2,threshold2]
sound_info.detection_threshold = 0.450;%for 1k (0.45)between 0.4 and 0.5 (0.5 gets rid of more noise) - for some 10k 0.8 (one file #8 in HA10-1L\2023-03-24)
sound_info.corrected_iti = [];%only needed when no thresholds work [file#,trial#,correctorno,start,end] 

sound_info.iti_tone_version = 1;% use version 2 if you have nice data

%if you don't want ITI info (need at least 2 speakers) make sound_info.correct = [];

% include string that is in the name of the abf files you want to read
[sound_st, sound_trials, sound_condition_array] = find_spkr_output_task_new(info,alignment_info,info.vr_sync_string,sound_info,sound_info.iti_tone_version);

%% get trial info using the virmen files!
for tr = 1:length(dataCell.dataCell)
    trial_info(tr).correct = dataCell.dataCell{1,tr}.result.correct;
    trial_info(tr).condition = dataCell.dataCell{1,tr}.maze.condition;
    if info.is_stim_dataset == 1
        trial_info(tr).is_stim = dataCell.dataCell{1,tr}.maze.is_stim_trial;
    end
end

%% get digidata iteration locations and difference between them

digidata_its = get_digidata_iterations(info.sync_base_path,info.vr_sync_string, info.virmen_channel);

%% find its in the data that best match the its for each trial dividing files into trials that match them (IF ITERATIONS ARE WEIRD USE THIS)

[file_estimated_trial_info,file_matching_trials,sound_condition_array] = match_trialsperfile(digidata_its, trial_info,sound_condition_array,task_info); %uses ITI sounds

%try this one if above does not work
%[file_estimated_trial_info,file_matching_trials,sound_condition_array] = match_trialsperfile_v2(digidata_its, trial_info,sound_condition_array,task_info,data); %uses ITI sounds if you have positive peaks

%use this one if you have pulse sounds
% [file_estimated_trial_info,file_matching_trials,sound_condition_array] = match_trialsperfile_pulse_maze(digidata_its, trial_info,sound_condition_array,task_info,data); %changes to ITI gap distance for pulse maze


% find the start and end trials that are within the imaging frames!// also
% puts trials into context of all other trials (file_digidata_trial_info)
[file_trial_ids,file_digidata_trial_info] = get_trial_ids(file_matching_trials,file_estimated_trial_info,alignment_info,info.sync_base_path,task_info);

%% shift iterations in time until they match positive peaks or first trial iteration (IF ITERATIONS ARE WEIRD USE THIS!!)
[virmen_it,trial_its,sound_condition_array] = shift_sync_data(data,file_trial_ids,digidata_its,file_estimated_trial_info,sound_condition_array,task_info);

% IF ITERATIONS ARE GOOD CAN USE THIS AND SKIP CODE ABOVE
% [virmen_it,trial_its,sound_condition_array] = get_virmen_iterations_and_times_digidata_positive_peaks(info.sync_base_path,info.virmen_channel,info.vr_sync_string,sound_condition_array,data);

%% binarize trial sounds and determine sounds for trials without speaker
 % (assumes all sounds are the same distance apart)
[sounds_per_file,vr_sound_frames] = binarize_sounds(virmen_it,sound_condition_array, trial_its,sound_info,sound_st,file_trial_ids,alignment_info);

% IF ITERATIONS ARE GOOD CAN USE THIS AND SKIP CODE ABOVE - use this if you have all the sounds and are happy with onsets/offsets
% [sounds_per_file,vr_sound_frames] = binarize_sounds_simple(virmen_it,sound_condition_array, trial_its,sound_info,sound_st,file_trial_ids,alignment_info); %can ignore vr_sound_frames
%% determine reward location 
[rewards_per_file,reward_loc_pure_frames,reward_loc_end_trial,reward_loc_pure] = find_reward(virmen_it,sound_condition_array,mdl_end_trial_sol,mdl_pure_sol,file_trial_ids,trial_its,digidata_its,trial_info,alignment_info);

%%  align virmen data!
%(dff,deconv,virmen_aq,alignment_info,data,dataCell,trial_its,stimulus_info,reward_info)
imaging = align_virmen_data(dff,deconv,virmen_it,alignment_info,data,dataCell,trial_its,sounds_per_file,reward_loc_pure_frames);
selected_fields = [1,6,7]; %1 y position, 6 reward, 7 ITI %inside movememnt_in_imaging_time
plot_random_trials_alignment (imaging,selected_fields);
%% save data!
mkdir(info.save_path)
cd(info.save_path)
save('alignment_variables','task_info','sound_info','sounds_per_file','virmen_it','sound_st', 'sound_trials', 'sound_condition_array','reward_loc_pure_frames','trial_its','file_trial_ids','file_matching_trials','digidata_its','info');
save('imaging','imaging');
save('vr_sound_frames',"vr_sound_frames"); 
%% redo imaging (loads necesarry info and reruns aling_virmen_data)

redo_imaging(info,0);