function frames = find_closest_frames(frame_times,data,min_distance)
frames = zeros(1,length(data));
count = 0;
for datapoint = data
    count = count+1;
    temp = find(frame_times - datapoint>0,1);
    if frame_times(temp) - datapoint < min_distance
        frames(count) = temp;
    else
        frames(count) =nan;
    end
end