function [acquisitions,updated_trial_its,sound_condition_array]  = get_virmen_iterations_and_times_digidata_positive_peaks(base, virmen_channel_number,string,sound_condition_array,data,file_trial_ids,file_estimated_trial_info)
file_ind = 0;

%initialize variables
possible_it_times = [];
possible_iterations = [];

cd(base);
z = dir(strcat('*',string,'*.abf'));
num_files = length(z);

[updated_trial_its] = virmen_it_rough_estimation(data);

for n = 1:num_files  %%%%still need to deal with iteration #1
  %keep n file_ind z virmen_channel_number acquisitions string
  %if contains(z(n).name,string)==1
  file_ind = n;
  [sync_data,sync_sampling_interval,~]= abfload(z(n).name);
    [sync_data,sync_sampling_interval,~]= abfload(z(n).name);
    sync_sampling_rate = 1/sync_sampling_interval*1e6;
    temp = sync_data(:,virmen_channel_number);
    time = 1:length(temp).*1/sync_sampling_rate;
    temp = double(temp);
    temp = temp-median(temp);
    
    [pks,locs] = findpeaks(abs(temp),'MINPEAKHEIGHT',0.09,'MinPeakDistance',5);
    temp = temp./abs(mean(temp(locs)));
    %temp = temp./abs(mean(temp(locs)))*1e3; %do I need to do this???
    %%%if it's the first file, ignore Virmen initializing? Very negative
    %%%initial pulse
    it_values = temp(locs); %%%wavesurfer values at pulse times
    it_times = locs;  %%%in wavesurfer units
    
    marked_its = find(it_values>0);
    if isempty(marked_its)
        disp('Error! Current file does not have positive values. File is probably too short')
        n=n+1;
%         file_ind = file_ind +1;
    else
    
    %findpeaks could find peaks right next to each other this code will
    %find the min of a group and use that as the marked_its
    %it will get rid of extra neighbouring detected peaks and account for
    %this by subtracting it from the marked_its
    group = {}; count = 1; 
    temp2=[];
    temp3=[];
    if length(marked_its)>1
        for m = 1:length(marked_its)-1
            if marked_its(m)+2>marked_its(m+1) %if current value +2 is greater than current value+1
                temp2 = [temp2,marked_its(m)];
                temp3 = min(temp2);
                group(count).group = temp3;
                group(count).temp2 = temp2;
            else
                temp2 = [temp2,marked_its(m)];
                group(count).temp2 = temp2;
                group(count).group = min(group(count).temp2);
                temp2=[];
                count= count+1;
            end
        end
    %dealing with if you have a last value with no group
    if marked_its(length(marked_its)) - marked_its(length(marked_its)-1) >10 %10 value could change
        group(length(group)+1).group = marked_its(length(marked_its));
        group(length(group)).temp2 = marked_its(length(marked_its));
    end
    marked_its_old = marked_its; 
    extra_num = {};
    marked_its=[];
    marked_its_shifted = [];
        for i =1:length(group)
            extra_num(i).group = length(group(i).temp2)-1;
            if i == 1 %no value in front of first one so this is the true value
               marked_its(i) = group(i).group;
               marked_its_shifted(i)=group(i).group;
            else
            temp5= [];
            temp3 = group(i).group;
                for j = i-1:-1:1
                temp4 = sum(extra_num(j).group);
                temp5 = [temp5, temp4];
                end
            marked_its(i) = temp3 - sum(temp5);
            marked_its_shifted(i)= temp3;
            end
        end
    extra_its = ~ismember(marked_its_old,marked_its_shifted);
    it_times(marked_its_old(extra_its)) = [];
    it_values = temp(it_times);
    end
    actual_it_values = zeros(1,length(it_values));
    
    %if first 10000?
    if round(it_values(marked_its(1)),0)==10 %round(it_values(marked_its(1)),1)==10 was set to 1 instead of 0
        actual_it_values(marked_its(1)) = 1; %might need to make zero to deal with negative iteration which is not part of the true iterations
         
    else
        actual_it_values(marked_its(1)) = round(((it_values(marked_its(1))*1e5)/1e4))*1e4;
       
        next_value = actual_it_values((marked_its(1))); %%%changed cb 2/18/22
        
        next_its = marked_its(1)-1:-1:1;
        %%%count back from first known landmark,
        for it = 1:length(next_its)
            actual_it_values(next_its(it)) = next_value-1;
            next_value =  actual_it_values(next_its(it));
        end
    end
      %%%and forward from first....
    previous_value = actual_it_values((marked_its(1)))-1
    if length(marked_its)>1
        for m = 2:length(marked_its) 
           if m ==2 
                for it = marked_its(m-1):marked_its(m)-1 %% changed from it = marked_its(m-1):marked_its(m)-1
                    actual_it_values(it) = previous_value+1;
                    previous_value = actual_it_values(it);
                end
           end

           if round(previous_value+1,-4)==round(((it_values(marked_its(m))*1e5)/1e4))*1e4
               if m <length(marked_its) 
                   for it = marked_its(m):marked_its(m+1)-1%length(it_values)%marked_its(m+1)-1%marked_its(m):marked_its(m+1) %have to figure out numbers
                        actual_it_values(it) = previous_value+1;
                        previous_value = actual_it_values(it);
                   end
               else
                    for it = marked_its(m):length(it_values)%%have to figure out numbers
                        actual_it_values(it) = previous_value+1;
                        previous_value = actual_it_values(it);
                    end
               end
            else
                if m <length(marked_its) 
                   for it = marked_its(m):marked_its(m+1)-1%length(it_values)%marked_its(m+1)-1%marked_its(m):marked_its(m+1) %have to figure out numbers
                        actual_it_values(it) = previous_value+1;%round(((it_values(marked_its(m))*1e5)/1e4))*1e4+1;
                        previous_value = actual_it_values(it);
                   end
               else
                    for it = marked_its(m):length(it_values)%%have to figure out numbers
                        actual_it_values(it) = previous_value+1;%round(((it_values(marked_its(m))*1e5)/1e4))*1e4+1;
                        previous_value = actual_it_values(it);
                    end
               end
           end
        end
    else
        for it = marked_its(1):length(actual_it_values) % marked_its(1)+1:length(actual_it_values)
            actual_it_values(it) = previous_value+1;
            previous_value = actual_it_values(it);
        end
        actual_it_values(length(actual_it_values))
    end

    %CHECK TO MAKE YOU ARE NOT SKIPPING ANY ITERATIONS!!
    if ~isempty(find(diff(actual_it_values)>1)) || ~isempty(find(diff(actual_it_values)==0))
        error('Error. n\Some iterations are missing');
    end
    %
    % get trial ids! used in other code so we have to do it
    % frame_start_trial = find(file_estimated_trial_info(file).start_trials_digidata_time - alignment_info(n).frame_times(1)>0);
    % frame_end_ITI = find(file_estimated_trial_info(file).end_iti_digidata_time - alignment_info(n).frame_times(end)<0);

    if ~isempty(file_trial_ids)
        start_trial_number = file_trial_ids(file_ind,1);%+file_trial_ids(file,3)-2; 
        end_trial_number = file_trial_ids(file_ind,2);
    
        %adjust iterations from virmen using actual recorded iterations! to
        %figure out true start iti and end trial iterations and keep track of
        %them for each file
        [trial_its,~] = virmen_it_rough_estimation_adjusted(data,sync_sampling_rate,it_times,actual_it_values); 
        if start_trial_number == 1
            updated_trial_its.end_trial_its(start_trial_number:end_trial_number+1) = trial_its.end_trial_its(start_trial_number:end_trial_number+1);
            updated_trial_its.end_iti_its(start_trial_number:end_trial_number+1) = trial_its.end_iti_its(start_trial_number:end_trial_number+1);
            updated_trial_its.start_trial_its(start_trial_number:end_trial_number+1) = trial_its.start_trial_its(start_trial_number:end_trial_number+1);
            updated_trial_its.start_iti_its(start_trial_number:end_trial_number+1) = trial_its.start_iti_its(start_trial_number:end_trial_number+1);
        elseif end_trial_number == length(trial_its)
            updated_trial_its.end_trial_its(start_trial_number-1:end_trial_number) = trial_its.end_trial_its(start_trial_number-1:end_trial_number);
            updated_trial_its.end_iti_its(start_trial_number-1:end_trial_number) = trial_its.end_iti_its(start_trial_number-1:end_trial_number);
            updated_trial_its.start_trial_its(start_trial_number-1:end_trial_number) = trial_its.start_trial_its(start_trial_number-1:end_trial_number);
            updated_trial_its.start_iti_its(start_trial_number-1:end_trial_number) = trial_its.start_iti_its(start_trial_number-1:end_trial_number);
        else
            updated_trial_its.end_trial_its(start_trial_number-1:end_trial_number+1) = trial_its.end_trial_its(start_trial_number-1:end_trial_number+1);
            updated_trial_its.end_iti_its(start_trial_number-1:end_trial_number+1) = trial_its.end_iti_its(start_trial_number-1:end_trial_number+1);
            updated_trial_its.start_trial_its(start_trial_number-1:end_trial_number+1) = trial_its.start_trial_its(start_trial_number-1:end_trial_number+1);
            updated_trial_its.start_iti_its(start_trial_number-1:end_trial_number+1) = trial_its.start_iti_its(start_trial_number-1:end_trial_number+1);
        end
        
        
        %assign iterations within file limits!
        possible_iterations = updated_trial_its.start_trial_its(start_trial_number):updated_trial_its.end_iti_its(end_trial_number); %limit iterations to ones within imaging frames// this is iterations ids
        possible_it_locs = find(ismember(actual_it_values,possible_iterations)); %locations of iterations within limits
    else
                % possible_iterations = updated_trial_its.start_trial_its(start_trial_number):updated_trial_its.end_iti_its(end_trial_number); %limit iterations to ones within imaging frames// this is iterations ids
        possible_it_locs = actual_it_values; %locations of iterations within limits

    end


    %include iti before first trial and full maze trial at the end if possible?
%     if file_trial_ids(file_ind,5) == 1 && file_trial_ids(file_ind,6) == 1
%         possible_iterations = trial_its.start_iti_its(start_trial_number-1):trial_its.end_trial_its(end_trial_number+1); %limit iterations to ones within imaging frames// this is iterations ids
%         possible_it_locs = find(ismember(actual_it_values,possible_iterations)); %locations of iterations within limits
%     elseif file_trial_ids(file_ind,5) == 1
%         possible_iterations = trial_its.start_iti_its(start_trial_number-1):trial_its.end_iti_its(end_trial_number); %limit iterations to ones within imaging frames// this is iterations ids
%         possible_it_locs = find(ismember(actual_it_values,possible_iterations)); %locations of iterations within limits
%     elseif file_trial_ids(file_ind,6) == 1
%         possible_iterations = trial_its.start_trial_its(start_trial_number):trial_its.end_trial_its(end_trial_number+1); %limit iterations to ones within imaging frames// this is iterations ids
%         possible_it_locs = find(ismember(actual_it_values,possible_iterations)); %locations of iterations within limits
%     else
%         possible_iterations = trial_its.start_trial_its(start_trial_number):trial_its.end_trial_its(end_trial_number); %limit iterations to ones within imaging frames// this is iterations ids
%         possible_it_locs = find(ismember(actual_it_values,possible_iterations)); %locations of iterations within limits
%     end

    %for dealing with weird zeros before acqusition counter for first file

    if ~isempty(sound_condition_array)
        if file_ind ==1
            if length(find(actual_it_values==0)) >0
                no_zeros = find(actual_it_values>0);
                no_zeros = [no_zeros(1)-1,no_zeros];
                acquisitions(file_ind).actual_it_values = actual_it_values(possible_it_locs);%(no_zeros);
        %         acquisitions(file_ind).it_times = it_times(no_zeros); %unsure if I need this
                acquisitions(file_ind).it_times = it_times(possible_it_locs);
                acquisitions(file_ind).directory = z(file_ind).name;
                %save sound onsets
                sound_onsets_speakers = [sound_condition_array(file_ind).VR_sounds{:,2}]; 
                nan_ind = find(isnan(sound_onsets_speakers));
                sound_onsets_speakers = sound_onsets_speakers(~isnan(sound_onsets_speakers));
                possible_sound_onsets = sound_onsets_speakers;
                acquisitions(file_ind).sound_trigger_time = possible_sound_onsets;
            else
                acquisitions(file_ind).actual_it_values = actual_it_values(possible_it_locs);
                acquisitions(file_ind).it_times = it_times(possible_it_locs);
                acquisitions(file_ind).directory = z(file_ind).name;
                %save sound onsets
                sound_onsets_speakers = [sound_condition_array(file_ind).VR_sounds{:,2}]; 
                nan_ind = find(isnan(sound_onsets_speakers));
                sound_onsets_speakers = sound_onsets_speakers(~isnan(sound_onsets_speakers));
                possible_sound_onsets = sound_onsets_speakers;
                acquisitions(file_ind).sound_trigger_time = possible_sound_onsets;
            
            end
        else
            acquisitions(file_ind).actual_it_values = actual_it_values(possible_it_locs);
            acquisitions(file_ind).it_times = it_times(possible_it_locs);
            acquisitions(file_ind).directory = z(file_ind).name;
        
            %save sound onsets
            sound_onsets_speakers = [sound_condition_array(file_ind).VR_sounds{:,2}]; 
            nan_ind = find(isnan(sound_onsets_speakers));
            sound_onsets_speakers = sound_onsets_speakers(~isnan(sound_onsets_speakers));
            possible_sound_onsets = sound_onsets_speakers;
            acquisitions(file_ind).sound_trigger_time = possible_sound_onsets;
    
        end
    else
                if file_ind ==1
            if length(find(actual_it_values==0)) >0
                no_zeros = find(actual_it_values>0);
                no_zeros = [no_zeros(1)-1,no_zeros];
                acquisitions(file_ind).actual_it_values = actual_it_values;%(no_zeros);
        %         acquisitions(file_ind).it_times = it_times(no_zeros); %unsure if I need this
                acquisitions(file_ind).it_times = it_times;
                acquisitions(file_ind).directory = z(file_ind).name;
                % %save sound onsets
                % sound_onsets_speakers = [sound_condition_array(file_ind).VR_sounds{:,2}]; 
                % nan_ind = find(isnan(sound_onsets_speakers));
                % sound_onsets_speakers = sound_onsets_speakers(~isnan(sound_onsets_speakers));
                % possible_sound_onsets = sound_onsets_speakers;
                % acquisitions(file_ind).sound_trigger_time = possible_sound_onsets;
            else
                acquisitions(file_ind).actual_it_values = actual_it_values;
                acquisitions(file_ind).it_times = it_times;
                acquisitions(file_ind).directory = z(file_ind).name;
                %save sound onsets
                % sound_onsets_speakers = [sound_condition_array(file_ind).VR_sounds{:,2}]; 
                % nan_ind = find(isnan(sound_onsets_speakers));
                % sound_onsets_speakers = sound_onsets_speakers(~isnan(sound_onsets_speakers));
                % possible_sound_onsets = sound_onsets_speakers;
                % acquisitions(file_ind).sound_trigger_time = possible_sound_onsets;
                % 
            end
        else
            acquisitions(file_ind).actual_it_values = actual_it_values;
            acquisitions(file_ind).it_times = it_times;
            acquisitions(file_ind).directory = z(file_ind).name;
        
            %save sound onsets
            % sound_onsets_speakers = [sound_condition_array(file_ind).VR_sounds{:,2}]; 
            % nan_ind = find(isnan(sound_onsets_speakers));
            % sound_onsets_speakers = sound_onsets_speakers(~isnan(sound_onsets_speakers));
            % possible_sound_onsets = sound_onsets_speakers;
            % acquisitions(file_ind).sound_trigger_time = possible_sound_onsets;
    
        end
    end
    figure(2);clf; hold on;plot(temp); plot(it_times(marked_its),0,'*c'); hold off
    
   % figure(4);clf; hold on; plot(sync_data(:,virmen_channel_number)); plot(acquisitions(n).it_times,0,'*c');plot(acquisitions(n).it_times(find(ismember(acquisitions(n).actual_it_values,updated_trial_its.start_iti_its(file_trial_ids(n,1):file_trial_ids(n,2))))),0,'*r');plot(acquisitions(n).it_times(find(ismember(acquisitions(n).actual_it_values,updated_trial_its.end_trial_its(file_trial_ids(n,1):file_trial_ids(n,2))))),0,'*b');plot(acquisitions(n).it_times(find(ismember(acquisitions(n).actual_it_values,updated_trial_its.start_trial_its(file_trial_ids(n,1):file_trial_ids(n,2))))),0,'*g');hold off
    pause()
    end
    
end

figure(3);clf;hold on;title('virmen iteration values across files'); for i= 1:length(acquisitions) ; plot(acquisitions(i).it_times,acquisitions(i).actual_it_values);end;hold off
%figure(3);clf;hold on; for i= 1:length(acquisitions) ; plot(acquisitions(i).it_times,acquisitions(i).actual_it_values);end;hold off
%figure (3); clf;for i = 1:8; subplot(8,1,i); hold on; plot(ex_file_long(:,i)'); hold off; end