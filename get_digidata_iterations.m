%% %loop through each file
function digidata_its = get_digidata_iterations(sync_base_path,string, virmen_channel)

cd(sync_base_path);
sync_dir = dir(strcat('*',string,'*.abf'));
num_files = length(sync_dir);
file_ind = 0;
for file = 1:num_files
    file_ind = file_ind +1
    [sync_data,sync_sampling_interval,~]  = abfload(sync_dir(file).name);
    sync_sampling_rate = 1/sync_sampling_interval*1e6;
    sync_data=sync_data';
    sync_data=double(sync_data);
    [pks,locs] = findpeaks(abs(sync_data(virmen_channel,:)),'MINPEAKHEIGHT',0.09,'MinPeakDistance',5); %might need to change distance depending on rate
    it_gaps = diff(locs);
    %find any positive peaks in the file!
    [pos_pks,pos_locs] = findpeaks((sync_data(virmen_channel,:)),'MINPEAKHEIGHT',0.09,'MinPeakDistance',5); %might need to change distance depending on rate
    % organize per file
    digidata_its(file_ind).locs = locs;
    digidata_its(file_ind).it_gaps = it_gaps;
    digidata_its(file_ind).sync_sampling_rate = sync_sampling_rate;
    digidata_its(file_ind).directory = strcat(sync_base_path,sync_dir(file).name);
    digidata_its(file_ind).pos_loc = pos_locs;
    digidata_its(file_ind).pos_pks = round(pos_pks,1);

   
end