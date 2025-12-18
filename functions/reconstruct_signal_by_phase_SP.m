function [reconstructed_signal, ideal_total_curve, quality, fig_handles, fit_results] = ...
    reconstruct_signal_by_phase_SP(t_ms, y_data, cardiac_phase, model_lv, model_rv, fit_results_input, target_phase)
%RECONSTRUCT_SIGNAL_BY_PHASE_SP Retrospective cardiac phase-alignment for HP-MRS.
%
%   Description:
%     Demodulates the measured pyruvate timecourse using a time-dependent 
%     modulation factor M(t) derived from ventricular volume V(phi). [cite: 1151]
%     Aligns the dynamic signal to a specific user-defined cardiac phase. [cite: 1153, 1214]
%
%   Inputs:
%     t_ms          - Acquisition time points in milliseconds (ms).
%     y_data        - Measured signal intensities (AIF).
%     cardiac_phase - Vector of recorded cardiac phases (phi from 0 to 1) per excitation. [cite: 1133]
%     model_lv/rv   - MATLAB cfit objects representing ventricular volume curves from 1H cine.
%     fit_results_input - Kinetic parameters estimated from the hemodynamic fit.
%     target_phase  - Target cardiac phase for reconstruction (e.g., 0.31 for end-systole). 
%
%   Outputs:
%     reconstructed_signal - Phase-corrected pyruvate timecourse aligned to target_phase.
%     ideal_total_curve    - Theoretical model curve without phase-induced fluctuations.

fit_results = fit_results_input;
quality = []; 

%% --- 1. Process and Visualize Volume Models ---
phase_points = linspace(0, 1, 501)';
rv_vol_raw = model_rv(phase_points);
lv_vol_raw = model_lv(phase_points);

% Normalize to Max=1 (ED volume)
[rv_max_vol, idx_rv_max] = max(rv_vol_raw);
[rv_min_vol, idx_rv_min] = min(rv_vol_raw);
rv_volume_profile_normalized = rv_vol_raw / rv_max_vol;

[lv_max_vol, idx_lv_max] = max(lv_vol_raw);
[lv_min_vol, idx_lv_min] = min(lv_vol_raw);
lv_volume_profile_normalized = lv_vol_raw / lv_max_vol;

% --- Plotting Volume Profiles with Annotations ---
fig_handles.volume_profile = figure('Name', 'Volume Profiles');
hold on;

% 1. Plot Curves
plot(phase_points, rv_volume_profile_normalized, 'b-', 'LineWidth', 2, 'DisplayName', 'Normalized RV Volume');
plot(phase_points, lv_volume_profile_normalized, 'g-', 'LineWidth', 2, 'DisplayName', 'Normalized LV Volume');

% 2. Identify ED (End-Diastole / Max) and ES (End-Systole / Min)
% RV Points
rv_ed_phase = phase_points(idx_rv_max); rv_ed_val = 1;
rv_es_phase = phase_points(idx_rv_min); rv_es_val = rv_volume_profile_normalized(idx_rv_min);

% LV Points
lv_ed_phase = phase_points(idx_lv_max); lv_ed_val = 1;
lv_es_phase = phase_points(idx_lv_min); lv_es_val = lv_volume_profile_normalized(idx_lv_min);

% 3. Add Markers and Text Labels
% RV Markers (Blue)
plot(rv_ed_phase, rv_ed_val, 'bo', 'MarkerFaceColor', 'b', 'HandleVisibility', 'off');
text(rv_ed_phase, rv_ed_val, ' RV ED', 'Color', 'b', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold', 'FontSize', 9);

plot(rv_es_phase, rv_es_val, 'bo', 'MarkerFaceColor', 'b', 'HandleVisibility', 'off');
text(rv_es_phase, rv_es_val, ' RV ES', 'Color', 'b', 'VerticalAlignment', 'top', 'FontWeight', 'bold', 'FontSize', 9);

% LV Markers (Dark Green for visibility)
dark_green = [0 0.5 0]; 
plot(lv_ed_phase, lv_ed_val, 'go', 'MarkerFaceColor', 'g', 'HandleVisibility', 'off');
text(lv_ed_phase, lv_ed_val, ' LV ED', 'Color', dark_green, 'VerticalAlignment', 'bottom', 'FontWeight', 'bold', 'FontSize', 9);

plot(lv_es_phase, lv_es_val, 'go', 'MarkerFaceColor', 'g', 'HandleVisibility', 'off');
text(lv_es_phase, lv_es_val, ' LV ES', 'Color', dark_green, 'VerticalAlignment', 'top', 'FontWeight', 'bold', 'FontSize', 9);

% Styling
title('Normalized Ventricular Volume vs. Cardiac Phase');
xlabel('Cardiac Phase (0-1)'); ylabel('Relative Volume (Normalized to ED=1)');
legend('show', 'Location', 'southwest'); grid on;
ylim([0.25 1.1]); % Add some headroom for the "ED" labels

%% --- 2. Calculate "Modulation Factor" based on the FULL signal model ---
t_sec = t_ms(:) / 1000;
predict_time = fit_results.Predict_curve_time;

% Get all components of the ideal signal (RV, LV, and Recirculation)
ideal_rv_curve_dynamic = fit_results.Predict_RV_FP - fit_results.C_fixed;
ideal_lv_curve_dynamic = fit_results.Predict_LV_FP - fit_results.C_fixed;
ideal_recirc_curve_dynamic = fit_results.Predict_Recirc_RV + fit_results.Predict_Recirc_LV;

% Interpolate to find the ideal signal values at the exact measurement time points
ideal_rv_at_t_sec = interp1(predict_time, ideal_rv_curve_dynamic, t_sec, 'linear', 'extrap');
ideal_lv_at_t_sec = interp1(predict_time, ideal_lv_curve_dynamic, t_sec, 'linear', 'extrap');
ideal_recirc_at_t_sec = interp1(predict_time, ideal_recirc_curve_dynamic, t_sec, 'linear', 'extrap');

% This is the ideal, phase-aligned total dynamic signal
ideal_total_dynamic_at_t_sec = ideal_rv_at_t_sec + ideal_lv_at_t_sec + ideal_recirc_at_t_sec;

% Get the volume scaling factors at each measurement time point
rv_vol_factor = interp1(phase_points, rv_volume_profile_normalized, cardiac_phase(:), 'linear', 'extrap');
lv_vol_factor = interp1(phase_points, lv_volume_profile_normalized, cardiac_phase(:), 'linear', 'extrap');

% Calculate the weighted average volume factor for the recirculation component
% Weight based on the amplitude ratio found in the fit
total_amp = fit_results.AmpRV + fit_results.AmpLV;
if total_amp == 0, total_amp = 1; end % Prevent div by zero
rv_weight = fit_results.AmpRV / total_amp;
lv_weight = fit_results.AmpLV / total_amp;
recirc_vol_factor = (rv_vol_factor * rv_weight) + (lv_vol_factor * lv_weight);

% This is the PREDICTED signal WITH cardiac phase modulation
predicted_modulated_signal = (ideal_rv_at_t_sec .* rv_vol_factor) + ...
                             (ideal_lv_at_t_sec .* lv_vol_factor) + ...
                             (ideal_recirc_at_t_sec .* recirc_vol_factor);

% The modulation factor is the ratio of the phase-modulated signal to the ideal signal
modulation_factor = ones(size(t_sec));
valid_indices = abs(ideal_total_dynamic_at_t_sec) > 1e-7; % Avoid division by zero
modulation_factor(valid_indices) = predicted_modulated_signal(valid_indices) ./ ideal_total_dynamic_at_t_sec(valid_indices);

%% --- 3. Perform Signal Reconstruction ---
% Determine the single modulation factor for the target phase
target_vol_factor_rv = interp1(phase_points, rv_volume_profile_normalized, target_phase);
target_vol_factor_lv = interp1(phase_points, lv_volume_profile_normalized, target_phase);
target_modulation_factor = rv_weight * target_vol_factor_rv + lv_weight * target_vol_factor_lv;

% Reconstruct: remove original modulation, apply target modulation
reconstructed_signal_dynamic = (y_data(:) - fit_results.C_fixed) ./ modulation_factor * target_modulation_factor;
reconstructed_signal = reconstructed_signal_dynamic + fit_results.C_fixed;

% The "ideal total curve" for comparison
ideal_total_curve.time_sec = predict_time;
ideal_total_curve.signal = fit_results.Predict_total; 

%% --- 4. Visualize the Reconstruction ---
fig_handles.reconstruction = figure('Name', 'Reconstruction Plot');
plot(t_sec, y_data, 'ko-', 'MarkerFaceColor', 'k', 'MarkerSize', 4, 'DisplayName', 'Measured pyruvate');
hold on;
plot(ideal_total_curve.time_sec, ideal_total_curve.signal, '--', 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 2, 'DisplayName', 'Ideal total fit');
plot(t_sec, reconstructed_signal, 'bo-', 'MarkerFaceColor', 'b', 'MarkerSize', 4, 'DisplayName', sprintf('Reconstructed (aligned to Phase %.2f)', target_phase));
xlabel('Time (seconds)');
ylabel('Signal Intensity');
legend('show', 'Location', 'northeast');
grid on;
end