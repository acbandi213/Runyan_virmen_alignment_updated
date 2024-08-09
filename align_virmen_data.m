function imaging = align_virmen_data(dff,deconv,virmen_aq,alignment_info,data,dataCell,trial_its,stimulus_info,reward_info)
%z-score all dff (assumes neurons x time)
z_dff=zscore(dff,0,2);

imaging ={};
file_ind = 1; % equivalent to imaging_trials in match_virmen_imaging - wavesurfer files 
previous_frames_sum = 0; previous_frames_temp = 0; previous_frames = 0;
turning_threshold = 0.1; %this is used to decide when the mouse has started to turn looking at x_position for that trial
weird_trials = [];

for vr_trial = 1:length(dataCell.dataCell)-1%1:length(dataCell.dataCell)-1 % virmen trials within each acquisition (ex: VR_stim01)
    %new method to get iterations
    start_it = trial_its.start_trial_its(vr_trial); 
    end_it = trial_its.end_trial_its(vr_trial); 
    iti_start_it = trial_its.start_iti_its(vr_trial);  
    iti_end_it = trial_its.end_iti_its(vr_trial); 

    % for each trial in this folder/acquisition find the start and end iteration #
    imaging(vr_trial).start_it = start_it;
    imaging(vr_trial).end_it = end_it;
    imaging(vr_trial).iti_start_it = iti_start_it;
    imaging(vr_trial).iti_end_it = iti_end_it;

    %virmen trial info! for all trials even without imaging
    virmen_trial_info(vr_trial).correct = dataCell.dataCell{1,vr_trial}.result.correct; %correct or not for that trial
    virmen_trial_info(vr_trial).left_turn = dataCell.dataCell{1,vr_trial}.result.leftTurn; %left turn or not for that trial
    virmen_trial_info(vr_trial).condition = dataCell.dataCell{1,vr_trial}.maze.condition; %stimulus condition for that trial
    if isfield(dataCell.dataCell{1, vr_trial}.maze,'is_stim_trial')%stim_dataset == 1
       virmen_trial_info(vr_trial).is_stim_trial = [dataCell.dataCell{1, vr_trial}.maze.is_stim_trial];
    end 
    imaging(vr_trial).virmen_trial_info = virmen_trial_info(vr_trial);
   % for each trial making sure you are within bounds of this tseries folder's iteration numbers else go to the next folder!
    if start_it>= virmen_aq(file_ind).actual_it_values(1)&& end_it <  virmen_aq(file_ind).actual_it_values(end) 
        
        output_data = {}; % file_frame_data equivalent to output_data from CR code
        output_data.frame_times = alignment_info(file_ind).frame_times; % frame times from res galvo signal in digidata time
        output_data.iteration_times = virmen_aq(file_ind).it_times; % it times from virmen iterations in digidata time
        output_data.iteration_ids = virmen_aq(file_ind).actual_it_values; % it_ids from virmen iterations
        
%             %digidata info
%             imaging(vr_trial).start_digidata = output_data.iteration_times(find(output_data.iteration_ids==start_it));
%             imaging(vr_trial).iti_end_digidata = output_data.iteration_times(find(output_data.iteration_ids==iti_end_it));


          if  iti_end_it <= output_data.iteration_ids(end) && any(output_data.iteration_times(find(output_data.iteration_ids==start_it)) >= output_data.frame_times(1)) ...
                  && any(output_data.iteration_times(find(output_data.iteration_ids==end_it))<=alignment_info(file_ind).frame_times(end))
                    %this if statement gets rid of trials if not within the frame limit at the very start since finding min will say multiple trials start and end in frame 1

                    %now that we know this is a good trial with imaging data we can get more info
                  %movement in virmen time using virmen iterations for that
                  %trial! starting with start it ending with the last iteration of the iti (keep in mind that reward period is not
                  %accounted for in virmen time)
                    movement_in_virmen_time(vr_trial).y_position = data.data(3,start_it:iti_end_it);
                    movement_in_virmen_time(vr_trial).x_position = data.data(2,start_it:iti_end_it);
                    movement_in_virmen_time(vr_trial).view_angle = data.data(4,start_it:iti_end_it);
                    movement_in_virmen_time(vr_trial).x_velocity = data.data(5,start_it:iti_end_it);
                    movement_in_virmen_time(vr_trial).y_velocity = data.data(6,start_it:iti_end_it);
                    movement_in_virmen_time(vr_trial).cw = data.data(7,start_it:iti_end_it); %current world
                    movement_in_virmen_time(vr_trial).is_reward = data.data(8,start_it:iti_end_it);
                    movement_in_virmen_time(vr_trial).in_ITI = data.data(9,start_it:iti_end_it);
            
                    imaging(vr_trial).movement_in_virmen_time = movement_in_virmen_time(vr_trial);
                    %virmen_trial_info(vr_trial).maze_length = round(max(movement_in_virmen_time(1).y_position),-2);%y position goes to like 310.9 so rounding down

            %finding the closest frames that match the start_it and iti_end_it for this trial
            frame_ids = [];
            for its = start_it:iti_end_it
                %[~,frame] = min(abs(output_data.frame_times - output_data.iteration_times(find(output_data.iteration_ids==its))));
                frame = find(output_data.frame_times - output_data.iteration_times(find(output_data.iteration_ids==its))>0,1);
                frame_ids = [frame_ids,frame];
            end
            imaging(vr_trial).frame_id = frame_ids;% virmen iterations that fit within this trial including ITI
            maze_start_frame = imaging(vr_trial).frame_id(1);
            iti_end_frame = imaging(vr_trial).frame_id(end);
   
            %finding the closest frames per trial epoch
            frame_ids = [];
            for its = start_it:end_it
                %[~,frame] = min(abs(output_data.frame_times - output_data.iteration_times(find(output_data.iteration_ids==its))));
                frame = find(output_data.frame_times - output_data.iteration_times(find(output_data.iteration_ids==its))>0,1);
                frame_ids = [frame_ids,frame];
            end
            imaging(vr_trial).frame_id_events.maze = frame_ids;

            frame_ids = [];
            for its = iti_start_it:iti_end_it
                %[~,frame] = min(abs(output_data.frame_times - output_data.iteration_times(find(output_data.iteration_ids==its))));
                frame = find(output_data.frame_times - output_data.iteration_times(find(output_data.iteration_ids==its))>0,1);
                frame_ids = [frame_ids,frame];
            end
            imaging(vr_trial).frame_id_events.iti = frame_ids;

            % include frames for reward period between end of trial and
            % start of iti
            imaging(vr_trial).frame_id_events.reward = imaging(vr_trial).frame_id_events.maze(end)+1:imaging(vr_trial).frame_id_events.iti(1)-1;
            %imaging(vr_trial).frame_id = sort([imaging(vr_trial).frame_id,imaging(vr_trial).frame_id_events.reward]);   

            %use frames to get the neural activity - start from first frame
            %of trial and end with last frame of iti (includes reward
            %period) frames included maze_start_frame:iti_end_frame

            %have to add previous frames from previous folders to make sure things match up
            %e.g. frame #1 in t-series folder 2 might actually be frames 10001
            imaging(vr_trial).dff = dff(:,maze_start_frame+previous_frames_sum:iti_end_frame+previous_frames_sum);
            imaging(vr_trial).z_dff = z_dff(:,maze_start_frame+previous_frames_sum:iti_end_frame+previous_frames_sum);
            imaging(vr_trial).deconv = deconv(:,maze_start_frame+previous_frames_sum:iti_end_frame+previous_frames_sum);
            
            imaging(vr_trial).relative_frames = maze_start_frame+previous_frames_sum:iti_end_frame+previous_frames_sum;
            
            imaging(vr_trial).file_num = file_ind;
            
            this_stimulus = [];
            %add stimulus information if it exits!
            if ~isempty(stimulus_info)
                %use frames included to find which sounds are in this trial
                this_stimulus = stimulus_info(file_ind).binary_frame_times(maze_start_frame:iti_end_frame);
            end
            
            this_reward_onset = [];
            this_pure_tone = [];
            %add reward information if it exists!
            if ~isempty(reward_info) && virmen_trial_info(vr_trial).correct == 1
                %use frames included to find which sounds are in this trial
                this_reward_onset = reward_info(vr_trial,1);
                this_pure_tone = reward_info(vr_trial,2):reward_info(vr_trial,3);
            else
                this_pure_tone = reward_info(vr_trial,2):reward_info(vr_trial,3); %onset to offset frames of pure tone
            end

            
            %initialize variables to be same size as imaging frames for
            %this trial
            this_y_position = zeros(1,size(imaging(vr_trial).dff,2));
            this_x_position = zeros(1,size(imaging(vr_trial).dff,2));
            this_x_velocity = zeros(1,size(imaging(vr_trial).dff,2));
            this_y_velocity = zeros(1,size(imaging(vr_trial).dff,2));
            this_view_angle = zeros(1,size(imaging(vr_trial).dff,2));
            this_reward = zeros(1,size(imaging(vr_trial).dff,2));
            this_ITI = zeros(1,size(imaging(vr_trial).dff,2));
            this_pure_tones = zeros(1,size(imaging(vr_trial).dff,2));
                
            frame_indices = [];
            maze_frames =[];
            iti_frames = [];
            reward_frames = [];
            for frame = 1:size(imaging(vr_trial).dff,2)
                %%%here finding the index of the virmen iteration that matches this frame id, then takes the value only at that FIRST iteration that matches the imaging frame id. 
                % Could in the future decide to bin across all virmen iterations that occur during an imaging frame
                ind = find(imaging(vr_trial).frame_id-min(imaging(vr_trial).frame_id)+1==frame,1,'first');  
                if ~isempty(ind) 
                    % align data based on current trial's frames! getting the first one that matches
                    this_y_position(frame) = movement_in_virmen_time(vr_trial).y_position(ind); %this is carolines variable's names
                    this_y_velocity(frame) = movement_in_virmen_time(vr_trial).y_velocity(ind);
                    this_x_velocity(frame) = movement_in_virmen_time(vr_trial).x_velocity(ind);
                    this_x_position(frame) = movement_in_virmen_time(vr_trial).x_position(ind);
                    this_view_angle(frame) = movement_in_virmen_time(vr_trial).view_angle(ind);
                    this_ITI(frame) = movement_in_virmen_time(vr_trial).in_ITI(ind);

                    % ---------ADJUSTING REWARD--------------
                    %adjusting reward to the frame of onset (which is
                    %incorrect in the virmen time reward)
                    this_reward(frame) = 0;%movement_in_virmen_time(vr_trial).is_reward(ind);
                    
                    if ~isempty(this_stimulus)
                        this_stimulus(frame) = this_stimulus(frame);
                    end  
                    frame_indices(frame) = ind;

                    % define what event the frame belongs to
                    if ismember(imaging(vr_trial).frame_id(ind),imaging(vr_trial).frame_id_events.maze)
                        maze_frames = [maze_frames,frame];
                    else
                        iti_frames = [iti_frames,frame];
                    end
                    
                else 
                    % if there are no frames for that iteration
                    %happens between trials, between end of trial and start
                    %iti and in the middle of sounds and before first sound
                    ind2 = find(imaging(vr_trial).frame_id_events.reward-min(imaging(vr_trial).frame_id_events.maze)+1==frame,1,'first');
                    % INCLUDE REWARD ONSET FRAME
                    if ~isempty(this_reward_onset) && frame == this_reward_onset-min(imaging(vr_trial).frame_id_events.maze)+1
                        this_reward(frame) = 1;
                    end
                    % INCLUDE PURE TONES!
                    if ~isempty(this_pure_tone) && ismember(frame,this_pure_tone-min(imaging(vr_trial).frame_id_events.maze)+1)
                        this_pure_tones(frame) = 1; 
                    end
                    if ~isempty(ind2) 
                        reward_frames = [reward_frames,frame];
                    elseif imaging(vr_trial).frame_id_events.maze(1)-min(imaging(vr_trial).frame_id)+1 <= frame && imaging(vr_trial).frame_id_events.maze(end)-min(imaging(vr_trial).frame_id)+1 >= frame
                        maze_frames = [maze_frames,frame];
                        this_pure_tones(frame) = nan;
                        this_reward(frame) = nan;
                    else
                        iti_frames = [iti_frames,frame];
                    end
                     
                    this_y_position(frame) = NaN;
                    this_y_velocity(frame) = NaN;
                    this_x_velocity(frame) = NaN;
                    this_x_position(frame) = NaN;
                    this_view_angle(frame) = NaN;
                    %this_reward(frame) = NaN;
                    this_ITI(frame) = NaN;
                    if ~isempty(this_stimulus)
                        this_stimulus(frame) = NaN;
                    end
                    frame_indices(frame) = NaN;
                    
                end
                
            end

            movement_in_imaging_time.y_position = this_y_position;
            movement_in_imaging_time.x_position = this_x_position;
            movement_in_imaging_time.y_velocity = this_y_velocity;
            movement_in_imaging_time.x_velocity = this_x_velocity;
            movement_in_imaging_time.view_angle = this_view_angle;
            movement_in_imaging_time.is_reward = this_reward;
            movement_in_imaging_time.in_ITI = this_ITI;
            movement_in_imaging_time.pure_tones = this_pure_tones;
            movement_in_imaging_time.maze_frames = maze_frames;
            movement_in_imaging_time.iti_frames = iti_frames;
            movement_in_imaging_time.reward_frames = reward_frames;

            % check if ITI frames are after gap!
            %ones before NANs of reward period!
            if isnan(this_ITI(find(this_ITI == 1,1,'first')+1))
                figure(55);clf;plot(this_ITI);
                error('Error. n\Check ITI iteration location');
            end
            
            %deal with weird trials where sound plays in the wrong place
            last_val = movement_in_virmen_time(vr_trial).y_position(end);
            if last_val<round(max(movement_in_virmen_time(vr_trial).y_position))-9 && last_val>0 %could need to be readjusted
                weird_trials = [weird_trials,vr_trial+1];
            end

            if ~isempty(this_stimulus)
                movement_in_imaging_time.stimulus = this_stimulus;
            end

            %keep track of which indices in virmen time go with imaging time!
            movement_in_imaging_time.frame_indices = frame_indices; 
            
            %finding when the mouse starts to turn at the end of maze (value is
            %greater/equal to 0.1)
            turn_frame = find(abs (movement_in_imaging_time.x_position(2:end)) >= turning_threshold,1,'first');
            movement_in_imaging_time.turn_frame = turn_frame; %maybe put in a different place??   
            
            %save movement in imaging frames
            imaging(vr_trial).movement_in_imaging_time = movement_in_imaging_time; 
            
            %good trial: the entire trial is usable
            imaging(vr_trial).good_trial = 1;
          else %if not within the frame limit (of start_it and end_it) go to the next trial
            vr_trial = vr_trial +1 
          end
    elseif start_it < virmen_aq(file_ind).actual_it_values(1) % if any wavesurfer iterations happen before virmen trial, go to next trial
        vr_trial = vr_trial+1
    elseif file_ind<length(virmen_aq) %%%% if the wavesurfer iterations are not within virmen iterations (start_it and end_it )move to the next wavesurfer file
        file_ind = file_ind +1
        % this if statement is adding frames at the start of each folder based on how many folders were before it
        if file_ind == 1
            previous_frames_temp = 0;
        else
            previous_frames_temp = length(alignment_info(file_ind-1).frame_times);
        end
        previous_frames = [previous_frames,previous_frames_temp];
        previous_frames_sum = sum(previous_frames);
 
    end
end

%empty trials if they are considered weird
if ~isempty(weird_trials)
    weird_trials = unique([weird_trials,[stimulus_info(:).weird_trial]]);
    fn = fieldnames(imaging);
    for w = 1:length(weird_trials)
%         for k = 6:length(fn)
%             imaging(weird_trials(w)).(fn{k})= [];
%         end
        imaging(weird_trials(w)).good_trial = []; %keep info in case want to use ITI or for something else
    end
end
end

