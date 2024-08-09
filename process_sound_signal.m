function [binary_sound_signal,pure_tone_signal] = process_sound_signal(sound_signal,threshold,sync_sampling_rate, smoothing_val)
% process sound signal by taking the max and min of the signal
% scaling ampltude nonlinearly multiplied by scale_value

smooth_val = 0.006*sync_sampling_rate; % usually 0.005 works well
xmin = movmin(sound_signal,smooth_val);%ordfilt2(x, 1, true(value_n));
xmax = movmax(sound_signal,smooth_val);%ordfilt2(x, value_n*value_n, true(value_n));
%     figure(1);clf; 
%     hold on;plot(x,LineWidth=0.05);plot(xmin);plot(xmax); hold off
scale_value = 5;
%x1 = tanh(value*(upperenv-lowerenv)); %scale amplitude nonlinearly
scaled_sound_signal = tanh(scale_value*(xmax-xmin)); %scale amplitude nonlinearly

value2 =  smoothing_val; %0.015
smoothed_sound_signal = conv2(scaled_sound_signal,value2);%gaussian filter
% threshold =0.45;%

%sounds at 0, anything else 1
binary_sound_signal = heaviside(threshold-smoothed_sound_signal);

%nexttile
hold on; plot(smoothed_sound_signal);plot(binary_sound_signal); hold off

pure_tone_signal = binary_sound_signal.*sound_signal;
% figure(3);clf; plot(pure_tone_signal)


