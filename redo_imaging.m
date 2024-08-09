function redo_imaging(info,redoornot)
% Load necessary variables
basepath = strcat(info.server,'\Connie\ProcessedData\',info.mouse,'\',info.date,'\');
cd(basepath);
load('dff.mat');
load('alignment_info.mat');
cd([basepath '\deconv']);
load('deconv.mat');


cd([basepath '\VR']);
load('alignment_variables.mat');

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

if redoornot == 1
    %align data
    imaging = align_virmen_data(dff,deconv,virmen_it,alignment_info,data,dataCell,trial_its,sounds_per_file,reward_loc_pure_frames);
    
    % save new imaging structure
    mkdir(info.save_path)
    cd(info.save_path)
    save('imaging','imaging');
end