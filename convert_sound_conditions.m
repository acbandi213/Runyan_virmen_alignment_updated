function updated_condition_array = convert_sound_conditions(condition_array,conditions_per_speaker,speaker_ids)
%speaker_ids = [1,2,4];
total_speakers = 1:4;
updated_condition_array={};
if length(speaker_ids ) == 3
    for t = 1:length(condition_array)
            %CONDITION NAN
            if isnan(condition_array(1,t))
                id = setdiff(total_speakers,speaker_ids);
                updated_condition_array{1,t}= setdiff(find(conditions_per_speaker(:,id)),[find(conditions_per_speaker(:,speaker_ids(1)));find(conditions_per_speaker(:,speaker_ids(2)));find(conditions_per_speaker(:,speaker_ids(3)))]);
            
            %CONDITION1
            elseif condition_array(1,t) == 1 && ismember(2,speaker_ids) %if speaker 1 and 2 are present then the multisounds are detected separately
                id = speaker_ids(1);
                updated_condition_array{1,t}= setdiff(find(conditions_per_speaker(:,id)),[find(conditions_per_speaker(:,speaker_ids(2)));find(conditions_per_speaker(:,speaker_ids(3)))]);
            elseif condition_array(1,t) == 1 && ~ismember(2,speaker_ids) %if speaker 1 and 2 is not present then the multisounds are not detected separately
                id =  speaker_ids(1);
                updated_condition_array{1,t}= find(conditions_per_speaker(:,id));
            
             %CONDITION2
            elseif condition_array(1,t) == 3 && ismember(1,speaker_ids) %if speaker 1 and 2 are present then the multisounds are detected separately
                id =  speaker_ids(3);
                updated_condition_array{1,t}= setdiff(find(conditions_per_speaker(:,id)),[find(conditions_per_speaker(:,speaker_ids(1)));find(conditions_per_speaker(:,speaker_ids(2)))]);
            elseif condition_array(1,t) == 3 && ~ismember(1,speaker_ids) %if speaker 1 and 2 is not present then the multisounds are not detected separately
                id =  speaker_ids(3);
                updated_condition_array{1,t}= find(conditions_per_speaker(:,id));
            
            %CONDITION3 minus whatever sound they share with speaker_ids(1)
            elseif speaker_ids(2) == 2 && condition_array(1,t) == 2 && ismember(1,speaker_ids) %if speaker 3 and 4 are present then the multisounds are detected separately
                id =  speaker_ids(2);
                updated_condition_array{1,t}= setdiff(find(conditions_per_speaker(:,id)),[find(conditions_per_speaker(:,speaker_ids(1)));find(conditions_per_speaker(:,speaker_ids(3)))]);
           
                %if other speaker plugged in then this is what I need
           elseif speaker_ids(2) == 3 && condition_array(1,t) == 2 && ismember(4,speaker_ids) %if speaker 3 and 4 are present then the multisounds are detected separately
                id =  speaker_ids(2);
                updated_condition_array{1,t}= setdiff(find(conditions_per_speaker(:,id)),[find(conditions_per_speaker(:,speaker_ids(1)));find(conditions_per_speaker(:,speaker_ids(3)))]);
            
    
            %CONDITION4 - speakers next to each other refers to a single
            %condition
            elseif condition_array(1,t) == 4 && all(ismember([1,2],speaker_ids)) %if speaker  1 and 2 are present then the multisounds are detected separately
                id = [1,2];
                updated_condition_array{1,t}= find(conditions_per_speaker(:,id(1))&conditions_per_speaker(:,id(2)));
            elseif condition_array(1,t) == 4 && all(ismember([3,4],speaker_ids)) %if speaker  1 and 2 are present then the multisounds are detected separately
                id = [3,4];
                updated_condition_array{1,t}= find(conditions_per_speaker(:,id(1))&conditions_per_speaker(:,id(2)));
            end
    
    end
% else
%     %FIXING LOGIC FOR PASSIVE 
%     for t = 1:length(condition_array)
%             %CONDITION NAN
%             if condition_array(1,t) == 1 && ismember(2,speaker_ids) %if speaker 1 and 2 are present then the multisounds are detected separately
%                 id = speaker_ids(1);
%                 updated_condition_array{1,t}= setdiff(find(conditions_per_speaker(:,id)),[find(conditions_per_speaker(:,speaker_ids(2)));find(conditions_per_speaker(:,speaker_ids(3)))]);    
%              %CONDITION2
%             elseif condition_array(1,t) == 3 && ismember(1,speaker_ids) %if speaker 1 and 2 are present then the multisounds are detected separately
%                 id =  speaker_ids(3);
%                 updated_condition_array{1,t}= setdiff(find(conditions_per_speaker(:,id)),[find(conditions_per_speaker(:,speaker_ids(1)));find(conditions_per_speaker(:,speaker_ids(2)))]);
%             elseif condition_array(1,t) == 3 && ~ismember(1,speaker_ids) %if speaker 1 and 2 is not present then the multisounds are not detected separately
%                 id =  speaker_ids(3);
%                 updated_condition_array{1,t}= find(conditions_per_speaker(:,id));
%             
%             %CONDITION3 minus whatever sound they share with speaker_ids(1)
%             elseif speaker_ids(2) == 2 && condition_array(1,t) == 2 && ismember(1,speaker_ids) %if speaker 3 and 4 are present then the multisounds are detected separately
%                 id =  speaker_ids(2);
%                 updated_condition_array{1,t}= setdiff(find(conditions_per_speaker(:,id)),[find(conditions_per_speaker(:,speaker_ids(1)));find(conditions_per_speaker(:,speaker_ids(3)))]);
%            
%                 %if other speaker plugged in then this is what I need
%            elseif speaker_ids(2) == 3 && condition_array(1,t) == 2 && ismember(4,speaker_ids) %if speaker 3 and 4 are present then the multisounds are detected separately
%                 id =  speaker_ids(2);
%                 updated_condition_array{1,t}= setdiff(find(conditions_per_speaker(:,id)),[find(conditions_per_speaker(:,speaker_ids(1)));find(conditions_per_speaker(:,speaker_ids(3)))]);
%             
%     
%             %CONDITION4 - speakers next to each other refers to a single
%             %condition
%             elseif condition_array(1,t) == 4 && all(ismember([1,2],speaker_ids)) %if speaker  1 and 2 are present then the multisounds are detected separately
%                 id = [1,2];
%                 updated_condition_array{1,t}= find(conditions_per_speaker(:,id(1))&conditions_per_speaker(:,id(2)));
%             elseif condition_array(1,t) == 4 && all(ismember([3,4],speaker_ids)) %if speaker  1 and 2 are present then the multisounds are detected separately
%                 id = [3,4];
%                 updated_condition_array{1,t}= find(conditions_per_speaker(:,id(1))&conditions_per_speaker(:,id(2)));
%             end
%     
%     end
end