function dif = trial_its_difference (trial_its2,data)
%find the difference to adjust its so that first iteration after ITI is
%after gap!
[trial_its,~] = virmen_it_rough_estimation(data);

dif.end_trial_its = trial_its.end_trial_its - trial_its2.end_trial_its;
dif.end_iti_its = trial_its.end_iti_its - trial_its2.end_iti_its;
dif.start_trial_its = trial_its.start_trial_its - trial_its2.start_trial_its;
dif.start_iti_its = trial_its.start_iti_its - trial_its2.start_iti_its;