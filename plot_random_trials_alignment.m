function plot_random_trials_alignment (imaging,selected_fields)
empty_trials = find(cellfun(@isempty,{imaging.good_trial}));
good_trials =  setdiff(1:length(imaging),empty_trials);
fprintf(['number of full trials: ' num2str(length(good_trials)) '\n'])
sel_trials = randperm(length(good_trials),9);
selected_trials = good_trials(sel_trials);
fieldname = fieldnames(imaging(good_trials(1)).movement_in_imaging_time);
colors = jet(length(selected_fields)+5);
figure(1010);clf;
tiledlayout(3,3,"TileSpacing","compact");
for t = 1:9
    imaging_cell = struct2cell(imaging(selected_trials(t)).movement_in_imaging_time);
    nexttile
    hold on
    for f = 1:length(selected_fields)
        if max(abs(imaging_cell{selected_fields(f),1})) > 1
            plot(rescale(imaging_cell{selected_fields(f),1},0,1),'color',colors(f,:));
        else
            plot(imaging_cell{selected_fields(f),1},'color',colors(f,:));
        end
    end
    legend(fieldname{selected_fields});
    title(num2str(selected_trials(t)))
    hold off
end
set(gcf,'position', [1,1,1500,900])
