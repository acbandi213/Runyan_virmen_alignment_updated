function [possible_alignment] = determine_shift(file_trial_ids,sound_condition_array, trial_its, file_digidata_trial_info,digidata_its,data)

for file = 1:length(digidata_its)
%% 2) aling data. Start by looking at first trial then shift the real iterations around until they make sense    
    start_trial_number = file_trial_ids(file,1);
    end_trial_number = file_trial_ids(file,2);

    %start_trial_its is the iti #
    %possible iterations IDs within file
    possible_iterations = trial_its.start_trial_its(start_trial_number): trial_its.end_iti_its(end_trial_number);
    
    % first say that the first iteration time in the first full trial is in the
    % position of the first iteration in the trua dataset (virmen data)
    file_trial_index = find([file_digidata_trial_info(file).estimated_digidata.trial_id] == start_trial_number);
    possible_digidata_time_start = file_digidata_trial_info(file).estimated_digidata(file_trial_index).start_trial_digidata_time; 
    %it_time_gaps
    %use time difference and assume first time difference to be exactly zero
    it_times_this_file = [data.data(1,possible_iterations)-data.data(1,possible_iterations(1))].* 86400;% difference with first iteration for each iteration in seconds!
    it_times_this_file = it_times_this_file*digidata_its(file).sync_sampling_rate;
    
    %start by assuming that the estimated digidata time alings perfectly with
    %true data- here is where I figure out how to shift over time
    possible_it_times = it_times_this_file + possible_digidata_time_start;

    % CHECKPOINTS!
    %1) see if there are any positive peaks- use this info for alignment if
    %it's there use this for alignment and compare with previous value!
    pos_pk_id=[];
    if ~isempty(digidata_its(file).pos_loc)
        pos_locs = digidata_its(file).pos_loc;
        pos_peaks = digidata_its(file).pos_pks;
        % determine the ID number based on the peak magnitude
        fprintf(strcat('#',num2str(file),' file has positive peaks!\n'));
        for pk = 1:length(pos_peaks)
            if pos_peaks(pk) == 10
                pos_pk_id(pk,:) = 1;
            else
                pos_pk_id(pk,:) = pos_peaks(pk)*1e5;
            end
        end   
        % shift things according to the peaks! use the first value if there
        % are multiple peaks
        if pos_pk_id(1) ==1
            possible_digidata_time_start_wpeaks = pos_locs(1);
            possible_it_times = it_times_this_file + possible_digidata_time_start_wpeaks;
        else
            %determine by how much to shift things by
            current_loc = it_times_this_file(pos_pk_id(1)-possible_iterations(1));
            %check how far away it is with previous shift
            if abs(current_loc - pos_locs(1)) > 0.0060 * digidata_its(file).sync_sampling_rate
                fprintf(strcat('#',num2str(file),' positive peak location does not make sense. Will use positive peak for alignment!\n'));
                estimated_shift = pos_locs(1) - current_loc;
                possible_it_times = it_times_this_file + possible_digidata_time_start+estimated_shift;

            end
        end
    end
    
    %2) see if sound is around iteration number 7 from small gap
    sound_onsets_speakers = [sound_condition_array(file).VR_sounds{:,2}]; %iterations_in_time(trial_its.start_trial_its(start_trial_number));
%     sound_onsets_iterations = trial_its.sound_trigger_its(find(trial_its.sound_trigger_its > possible_iterations(1) & trial_its.sound_trigger_its < possible_iterations(end)))+7; 
%     possible_sound_onsets = possible_it_times(sound_onsets_iterations-possible_iterations(1));

    sound_onsets_iterations = trial_its.sound_trigger_its(find(trial_its.sound_trigger_its > trial_its.start_trial_its(start_trial_number) & trial_its.sound_trigger_its <trial_its.end_iti_its(end_trial_number)))+6; %sound happens within 7 iterations
    possible_sound_onsets = possible_it_times(sound_onsets_iterations);

    for s = 1:length(possible_sound_onsets)
        difference_sounds = min(abs(possible_sound_onsets(s) - sound_onsets_speakers));
        if difference_sounds < 0.012 * digidata_its(file).sync_sampling_rate
            fprintf('Sound distances make sense!\n');
        else
            fprintf('Sound distances do not make sense!\n');
        end
    end

    possible_alignment(file).it_times = possible_it_times;
    possible_alignment(file).sound_onsets = possible_sound_onsets;
    ex_data = abfload(strcat(digidata_its(file).directory));
    figure(file);clf; hold on;plot(ex_data(:,7));plot(ex_data(:,6)); plot(rescale(ex_data(:,4),-1,0));plot(rescale(ex_data(:,8),-1,0));plot(possible_it_times,0,'*c');plot(possible_sound_onsets,0,'*r');hold off; movegui(gcf,'center');
end
%figure(); hold on;plot(ex_data(:,6)); plot(rescale(ex_data(:,4),-1,0));plot(possible_it_times,0,'*c');hold off; movegui(gcf,'center');

