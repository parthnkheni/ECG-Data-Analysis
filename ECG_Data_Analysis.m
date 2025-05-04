%Final_Project 
%1
load('Resting.mat');
resttime = b1(:, 1);
restecg = b1(:, 2);

load('Exercise.mat');
exertime = b1(:, 1);
exerecg = b1(:, 2);

load('BoxBreathing.mat');
boxtime = b1(:, 1);
boxecg = b1(:, 2);
disp('Data loaded successfully.');
%%
%2
figure;
plot(resttime(1:end), restecg(1:end), 'LineWidth', 2.5); 
title('Resting ECG', 'FontSize', 16);
xlabel('Time (s)', 'FontSize', 16);
ylabel('ECG (mV)', 'FontSize', 16);
xlim([10, 20]);
ylim([-0.2, 0.8]); 
set(gca, 'FontSize', 16);


figure;
plot(exertime(1:end), exerecg(1:end), 'LineWidth', 2.5);
title('Exercise ECG', 'FontSize', 16);
xlabel('Time (s)', 'FontSize', 16);
ylabel('ECG (mV)', 'FontSize', 16);
xlim([20, 30]);
ylim([-0.4, 0.8]); 
set(gca, 'FontSize', 16);



figure;
plot(boxtime(1:end), boxecg(1:end), 'LineWidth', 2.5);
title('Box Breathing ECG', 'FontSize', 16);
xlabel('Time (s)', 'FontSize', 16);
ylabel('ECG (mV)', 'FontSize', 16);
xlim([10, 20]);
ylim([-0.2, 0.8]); 
set(gca, 'FontSize', 16);

%subplot code

figure;

subplot(3, 1, 1);
plot(resttime(1:end), restecg(1:end), 'LineWidth', 1.5);
title('Resting ECG', 'FontSize', 16);
xlabel('Time (s)', 'FontSize', 16);
ylabel('ECG (mV)', 'FontSize', 16);
xlim([120, 130]);
ylim([-0.2, 0.8]);
set(gca, 'FontSize', 16);

subplot(3, 1, 2);
plot(exertime(1:end), exerecg(1:end), 'LineWidth', 1.5);
title('Exercise ECG', 'FontSize', 16);
xlabel('Time (s)', 'FontSize', 16);
ylabel('ECG (mV)', 'FontSize', 16);
xlim([120, 130]);
ylim([-0.4, 0.8]);
set(gca, 'FontSize', 16);

subplot(3, 1, 3);
plot(boxtime(1:end), boxecg(1:end), 'LineWidth', 1.5);
title('Box Breathing ECG', 'FontSize', 16);
xlabel('Time (s)', 'FontSize', 16);
ylabel('ECG (mV)', 'FontSize', 16);
xlim([120, 130]);
ylim([-0.2, 0.8]);
set(gca, 'FontSize', 16);

sgtitle('ECG Plots Under Different Conditions', 'FontSize', 18);



%%
%3
help findpeaks;
%%
%4
[rest_peaks, rest_locs] = findpeaks(restecg);

[exercise_peaks, exercise_locs] = findpeaks(exerecg);

[box_peaks, box_locs] = findpeaks(boxecg);


figure;
subplot(2, 1, 1);
plot(resttime, restecg, 'LineWidth', 1.5);
title('Resting', 'FontSize', 16);
xlabel('Time (s)', 'FontSize', 14);
ylabel('mV', 'FontSize', 14);
xlim([120, 130]);

subplot(2, 1, 2);
plot(resttime, restecg, 'LineWidth', 1.5);
hold on;
plot(resttime(rest_locs), rest_peaks, 'v', 'MarkerFaceColor', 'yellow', 'MarkerEdgeColor', 'red');
hold off;
title('Resting with Peaks', 'FontSize', 16);
xlabel('Time (s)', 'FontSize', 14);
ylabel('mV', 'FontSize', 14);
xlim([120, 130]);


figure;
subplot(2, 1, 1);
plot(exertime, exerecg, 'LineWidth', 1.5);
title('Exercise ', 'FontSize', 16);
xlabel('Time (s)', 'FontSize', 14);
ylabel('mV', 'FontSize', 14);
xlim([120, 130]);

subplot(2, 1, 2);
plot(exertime, exerecg, 'LineWidth', 1.5);
hold on;
plot(exertime(exercise_locs), exercise_peaks, 'v', 'MarkerFaceColor', 'yellow', 'MarkerEdgeColor', 'red');
hold off;
title('Exercise with Peaks', 'FontSize', 16);
xlabel('Time (s)', 'FontSize', 14);
ylabel('mV', 'FontSize', 14);
xlim([120, 130]);


figure;
subplot(2, 1, 1);
plot(boxtime, boxecg, 'LineWidth', 1.5);
title('Box Breathing ', 'FontSize', 16);
xlabel('Time (s)', 'FontSize', 14);
ylabel('mV', 'FontSize', 14);
xlim([120, 130]);

subplot(2, 1, 2);
plot(boxtime, boxecg, 'LineWidth', 1.5);
hold on;
plot(boxtime(box_locs), box_peaks, 'v', 'MarkerFaceColor', 'yellow', 'MarkerEdgeColor', 'red');
hold off;
title('Box Breathing with Peaks', 'FontSize', 16);
xlabel('Time (s)', 'FontSize', 14);
ylabel('mV', 'FontSize', 14);
xlim([120, 130]);
%%
%5
MinPeakHeight_rest = 0.4;

[rest_peaks, rest_locs] = findpeaks(restecg, 'MinPeakHeight', MinPeakHeight_rest);

MinPeakHeight_exercise = 0.4;

[exercise_peaks, exercise_locs] = findpeaks(exerecg, 'MinPeakHeight', MinPeakHeight_exercise);

MinPeakHeight_BoxBreathing = 0.4;

[Box_peaks, Box_locs] = findpeaks(boxecg, 'MinPeakHeight', MinPeakHeight_BoxBreathing);


figure;
subplot(2, 1, 1);
plot(resttime, restecg, 'LineWidth', 1.5);
title('Resting ECG', 'FontSize', 16);
xlabel('Time (s)', 'FontSize', 14);
ylabel('mV', 'FontSize', 14);
xlim([120, 130]);

subplot(2, 1, 2);
plot(resttime, restecg, 'LineWidth', 1.5);
hold on;
plot(resttime(rest_locs), rest_peaks, 'v', 'MarkerFaceColor', 'yellow', 'MarkerEdgeColor', 'red');
hold off;
title('Resting ECG with R-Peaks', 'FontSize', 16);
xlabel('Time (s)', 'FontSize', 14);
ylabel('mV', 'FontSize', 14);
xlim([120, 130]);



figure;
subplot(2, 1, 1);
plot(exertime, exerecg, 'LineWidth', 1.5);
title('Exercise ECG', 'FontSize', 16);
xlabel('Time (s)', 'FontSize', 14);
ylabel('mV', 'FontSize', 14);
xlim([120, 130]);

subplot(2, 1, 2);
plot(exertime, exerecg, 'LineWidth', 1.5);
hold on;
plot(exertime(exercise_locs), exercise_peaks, 'v', 'MarkerFaceColor', 'yellow', 'MarkerEdgeColor', 'red');
hold off;
title('Exercise ECG with R-Peaks', 'FontSize', 16);
xlabel('Time (s)', 'FontSize', 14);
ylabel('mV', 'FontSize', 14);
xlim([120, 130]);



figure;
subplot(2, 1, 1);
plot(boxtime, boxecg, 'LineWidth', 1.5);
title('Box Breathing ECG', 'FontSize', 16);
xlabel('Time (s)', 'FontSize', 14);
ylabel('mV', 'FontSize', 14);
xlim([120, 130]);

subplot(2, 1, 2);
plot(boxtime, boxecg, 'LineWidth', 1.5);
hold on;
plot(boxtime(Box_locs), Box_peaks, 'v', 'MarkerFaceColor', 'yellow', 'MarkerEdgeColor', 'red');
hold off;
title('Box breathing ECG with R-Peaks', 'FontSize', 16);
xlabel('Time (s)', 'FontSize', 14);
ylabel('mV', 'FontSize', 14);
xlim([120, 130]);


%%
%6

intervalsbtwR_rest = diff(rest_locs);
avepeakdistance_rest = mean(intervalsbtwR_rest);

figure;
subplot(3, 1, 1);
hold on;
for i = 2:length(rest_locs)-1
    startidx = max(1, round(rest_locs(i) - avepeakdistance_rest * 1/4));
    endidx = min(length(restecg), round(rest_locs(i) + avepeakdistance_rest * 3/4));
    percycle = restecg(startidx:endidx); % this is done this way even though it is not asked so it is eassier to compare the graphs it also provids a better visuale for the exercise ECG.
    percycle = percycle - mean(percycle);
    percycle = percycle / max(abs(percycle));
    plot(0:length(percycle)-1, percycle, 'Color', [1, 0, 0, 0.02]);
end
title('Resting ECG');
xlabel('Count');
ylabel('Normalized ECG (mV)');
ylim([-0.3,1.3])
grid on;
hold off;

intervalsbtwR_exercise = diff(exercise_locs);
avepeakdistance_exercise = mean(intervalsbtwR_exercise);

subplot(3, 1, 2);
hold on;
for i = 2:length(exercise_locs)-1
    startidx = max(1, round(exercise_locs(i) - avepeakdistance_exercise * 1/4));
    endidx = min(length(exerecg), round(exercise_locs(i) + avepeakdistance_exercise * 3/4));
    percycle = exerecg(startidx:endidx);
    percycle = percycle - mean(percycle);
    percycle = percycle / max(abs(percycle));
    plot(0:length(percycle)-1, percycle, 'Color', [0, 0, 1, 0.02]);
end
title('Exercise ECG');
xlabel('Count');
ylabel('Normalized ECG (mV)');
grid on;
hold off;

intervalsbtwR_box = diff(Box_locs);
avepeakdistance_box = mean(intervalsbtwR_box);

subplot(3, 1, 3);
hold on;
for i = 2:length(Box_locs)-1
    startidx = max(1, round(Box_locs(i) - avepeakdistance_box * 1/4));
    endidx = min(length(boxecg), round(Box_locs(i) + avepeakdistance_box * 3/4));
    percycle = boxecg(startidx:endidx);
    percycle = percycle - mean(percycle);
    percycle = percycle / max(abs(percycle));
    plot(0:length(percycle)-1, percycle, 'Color', [0, 1, 0, 0.02]);
end
title('Box Breathing ECG');
xlabel('Count');
ylabel('Normalized ECG (mV)');
grid on;
hold off;
%%
inPeakHeight_rest = 0.4;
MinPeakHeight_exercise = 0.4;

resttime_filtered = resttime(resttime >= 100 & resttime <= 180);
restecg_filtered = restecg(resttime >= 100 & resttime <= 180);
exertime_filtered = exertime(exertime >= 100 & exertime <= 180);
exerecg_filtered = exerecg(exertime >= 100 & exertime <= 180);

[~, rest_locs_filtered] = findpeaks(restecg_filtered, 'MinPeakHeight', inPeakHeight_rest);
[~, exercise_locs_filtered] = findpeaks(exerecg_filtered, 'MinPeakHeight', MinPeakHeight_exercise);

rest_intervals_by_hand = resttime_filtered(rest_locs_filtered(2:end)) - resttime_filtered(rest_locs_filtered(1:end-1));
exercise_intervals_by_hand = exertime_filtered(exercise_locs_filtered(2:end)) - exertime_filtered(exercise_locs_filtered(1:end-1));

rest_intervals_diff = diff(resttime_filtered(rest_locs_filtered));
exercise_intervals_diff = diff(exertime_filtered(exercise_locs_filtered));

rest_heart_rate_by_hand = 60 ./ rest_intervals_by_hand;
exercise_heart_rate_by_hand = 60 ./ exercise_intervals_by_hand;

rest_heart_rate_diff = 60 ./ rest_intervals_diff;
exercise_heart_rate_diff = 60 ./ exercise_intervals_diff;

figure;


subplot(2, 1, 1);
hold on;
histogram(rest_heart_rate_by_hand, 'Normalization', 'probability', 'BinEdges', 70:0.38:105, 'FaceColor', 'r', 'EdgeColor', 'k');
histogram(exercise_heart_rate_by_hand, 'Normalization', 'probability', 'BinEdges', 70:0.38:105, 'FaceColor', 'b', 'EdgeColor', 'k');
hold off;
title('Heart Rate (By Hand): Resting vs Exercise', 'FontSize', 16);
xlabel('Heart Rate (bpm)', 'FontSize', 14);
ylabel('Probability', 'FontSize', 14);
legend('Resting', 'Exercise');
xlim([70 105]);


subplot(2, 1, 2);
hold on;
histogram(rest_heart_rate_diff, 'Normalization', 'probability', 'BinEdges', 70:0.38:105, 'FaceColor', 'r', 'EdgeColor', 'k');
histogram(exercise_heart_rate_diff, 'Normalization', 'probability', 'BinEdges', 70:0.38:105, 'FaceColor', 'b', 'EdgeColor', 'k');
hold off;
title('Heart Rate (Using diff()): Resting vs Exercise', 'FontSize', 16);
xlabel('Heart Rate (bpm)', 'FontSize', 14);
ylabel('Probability', 'FontSize', 14);
legend('Resting', 'Exercise');
xlim([70 105]);

resting = [rest_heart_rate_by_hand; rest_heart_rate_diff];
exercises = [exercise_heart_rate_by_hand; exercise_heart_rate_diff];

group = [ones(size(resting)); 2 * ones(size(exercises))];

figure;
boxchart(group, [resting; exercises])
set(gca, 'XTickLabel', {'Resting', 'Exercise'})
xticks([1 2]);
ylabel('Heart Rate (bpm)', 'FontSize', 14);
title('Heart Rate Comparison: Resting vs Exercise', 'FontSize', 16);
%%
% 8
rpkmin = 0.4;  

[rntpk, locs_rnt] = findpeaks(restecg, 'MinPeakHeight', rpkmin);  

indices_rpk = locs_rnt;  

[tpk, locs_tpk] = findpeaks(restecg, 'MinPeakHeight', 0.1);  

valid_id = tpk < rpkmin;  
tpk = tpk(valid_id);  
locs_tpk = locs_tpk(valid_id);  

deltat = 1/200;
derivative_rpk = zeros(1, length(indices_rpk));

for i = 1:length(indices_rpk)
   start_idx = max(1, indices_rpk(i) - 10);  
   end_idx = min(length(restecg), indices_rpk(i) + 10);  
   sample_window = restecg(start_idx:end_idx);
   delta_amp = abs(diff(sample_window));  
   derivative_set = delta_amp / deltat;  
   derivative_rpk(i) = max(derivative_set);  
end

derivative_tpk = zeros(1, length(locs_tpk));

for i = 1:length(locs_tpk)
   start_idx = max(1, locs_tpk(i) - 10);  
   end_idx = min(length(restecg), locs_tpk(i) + 10);  
   sample_window = restecg(start_idx:end_idx);
   delta_amp = abs(diff(sample_window));  
   derivative_set = delta_amp / deltat;  
   derivative_tpk(i) = max(derivative_set);  
end

figure;
histogram(derivative_tpk, 'Normalization', 'probability', 'EdgeColor', 'yellow', 'FaceAlpha', 0.5, 'EdgeAlpha', 0.7, 'LineWidth', 1.5);
hold on;
histogram(derivative_rpk, 'Normalization', 'probability', 'EdgeColor', 'blue', 'FaceAlpha', 0.4, 'EdgeAlpha', 0.7, 'LineWidth', 1.5);
xlabel('Voltage Change (mV/s)', 'FontSize', 14);
title('Rate of Change of Voltage for R and T Phase', 'FontSize', 16);
legend('T', 'R');
hold off;

%%
%9

pds_re = diff(resttime(rest_locs));
pds_ex = diff(exertime(exercise_locs));

frequency_re = 1 ./ pds_re;
frequency_ex = 1 ./ pds_ex;

hrate_re = frequency_re * 60;
hrate_ex = frequency_ex * 60;

fprintf('For the two-sample t-test between resting and exercise heart rate\n');
[h1, p1] = ttest2(hrate_re, hrate_ex);
fprintf('The p-value is %.4f\n', p1);
if h1 == 1
   fprintf('Therefore, we can reject the null hypothesis\n');
else
   fprintf('Therefore, we cannot reject the null hypothesis\n');
end
fprintf('\n');

rest_period1_idx = resttime(rest_locs) >= 140 & resttime(rest_locs) <= 150;
rest_period2_idx = resttime(rest_locs) >= 150 & resttime(rest_locs) <= 160;

rest_period1 = hrate_re(rest_period1_idx);
rest_period2 = hrate_re(rest_period2_idx);

fprintf('For the two-sample t-test between two different periods of resting heart rate\n');
[h2, p2] = ttest2(rest_period1, rest_period2);
fprintf('The p-value is %.4f\n', p2);
if h2 == 1
   fprintf('Therefore, we can reject the null hypothesis\n');
else
   fprintf('Therefore, we cannot reject the null hypothesis\n');
end