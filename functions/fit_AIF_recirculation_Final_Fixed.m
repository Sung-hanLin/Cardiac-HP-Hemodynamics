function [results, quality, fig_handles] = fit_AIF_recirculation_Final_Fixed(t_ms, y_data, TarrivalRV_ini, early_data_endpoint_sec, full_fit_endpoint_sec)
%FIT_AIF_RECIRCULATION_FINAL_FIXED Multi-stage hemodynamic fitting for cardiac HP-MRS.
%
%   Description:
%     Implements a three-component model (RV first-pass, LV first-pass, and 
%     global recirculation) to deconvolve complex cardiac hemodynamics. [cite: 1189, 1222]
%     The model accounts for T1 decay and Gamma-variate AIF kinetics. [cite: 1101]
%
%   Inputs:
%     t_ms        - Vector of acquisition time points in milliseconds (ms).
%     y_data      - Measured [1-13C]pyruvate signal intensities (arbitrary units).
%     TarrivalRV_ini - Initial guess for RV bolus arrival time (seconds).
%     early_data_endpoint_sec - Time limit for the first-pass initial guess fitting (seconds).
%     full_fit_endpoint_sec  - Total time window for the final joint optimization (seconds).
%
%   Outputs:
%     results     - Struct containing estimated kinetic parameters (PTT, k, theta, scale). [cite: 1186]
%     quality     - Goodness-of-fit metrics (R-squared). [cite: 1190]
%     fig_handles - Handles to diagnostic plots (Fit composition, Residuals).
%
%   Reference: Lin et al., MRM (2025).

%% --- 0. Pre-processing ---
t_sec = t_ms(:) / 1000;
y_data = y_data(:);
C_fixed = mean(y_data(t_sec < TarrivalRV_ini));
y_data_subtracted = y_data - C_fixed;

% Optimization Options
options = optimoptions('lsqcurvefit', 'Algorithm', 'trust-region-reflective', 'Display', 'off', 'MaxFunctionEvaluations', 5000);

%% --- 1. Initial Guess: First Pass (Early Data) ---
early_indices = t_sec <= early_data_endpoint_sec;
t_early = t_sec(early_indices);
y_early_subtracted = y_data_subtracted(early_indices);

[p_fp_initial, ~] = get_first_pass_initial_guess(t_early, y_early_subtracted, TarrivalRV_ini);
y_fp_initial_full = dual_gamma_model(p_fp_initial, t_sec);
y_residual = y_data_subtracted - y_fp_initial_full;

%% --- 2. Initial Guess: Recirculation (Residual) ---
recirc_model_func = @(p_recirc, t) p_recirc(1) * interp1(t, y_fp_initial_full, t - p_recirc(2), 'linear', 0);
p0_recirc = [0.2, 22];
lb_recirc = [0, 15];
ub_recirc = [1.0, 45];
p_recirc_initial = lsqcurvefit(recirc_model_func, p0_recirc, t_sec, y_residual, lb_recirc, ub_recirc, options);

%% --- 3. Final Joint Optimization ---
p0_final = [p_fp_initial, p_recirc_initial];
% Bounds
lb_final = [0, 3.0, 2.5, 0, 6, TarrivalRV_ini-1, 5.0, 15, lb_recirc(1), lb_recirc(2)];
ub_final = [inf, 8.0, 3.9, inf, 12,  TarrivalRV_ini+1, 8.0, 67, ub_recirc(1), ub_recirc(2)];

final_indices = t_sec <= full_fit_endpoint_sec;
t_final = t_sec(final_indices);
y_final_subtracted = y_data_subtracted(final_indices);

full_model_func = @(p, t) full_recirculation_model(p, t);
final_fitted_params = lsqcurvefit(full_model_func, p0_final, t_final, y_final_subtracted, lb_final, ub_final, options);

%% --- 4. Process Results ---
results.AmpRV = final_fitted_params(1);
results.kRV = final_fitted_params(2);
results.thetaRV = final_fitted_params(3);
results.AmpLV = final_fitted_params(4);
results.kLV = final_fitted_params(5);
results.Time_arrival_RV = final_fitted_params(6);
results.delta_t = final_fitted_params(7);
results.T1_shared = final_fitted_params(8);
results.alpha = final_fitted_params(9);
results.t_offset = final_fitted_params(10);
results.thetaLV_fixed = 2.0;
results.C_fixed = C_fixed;
results.Time_arrival_LV = results.Time_arrival_RV + results.delta_t;

% Generate High-Res Curves
t_fine = linspace(min(t_sec), max(t_sec), 2000)';
p_fp_final = final_fitted_params(1:8);
alpha = final_fitted_params(9);
t_offset = final_fitted_params(10);

[~, y_rv_fine, y_lv_fine] = dual_gamma_model(p_fp_final, t_fine);
y_recirc_rv_fine = alpha * interp1(t_fine, y_rv_fine, t_fine - t_offset, 'linear', 0);
y_recirc_lv_fine = alpha * interp1(t_fine, y_lv_fine, t_fine - t_offset, 'linear', 0);
y_total_fine = y_rv_fine + y_lv_fine + y_recirc_rv_fine + y_recirc_lv_fine;

results.Predict_curve_time = t_fine;
results.Predict_total = y_total_fine + C_fixed;
results.Predict_RV_FP = y_rv_fine + C_fixed;
results.Predict_LV_FP = y_lv_fine + C_fixed;
results.Predict_Recirc_RV = y_recirc_rv_fine;
results.Predict_Recirc_LV = y_recirc_lv_fine;

% Quality Metrics
[y_fit_dynamic, ~, ~, ~, ~] = full_recirculation_model(final_fitted_params, t_sec);
residuals = y_data_subtracted - y_fit_dynamic;
SST = sum((y_data - mean(y_data)).^2);
SSE = sum(residuals.^2);
quality.R_squared = 1 - (SSE / SST);

% Calculate Peak Times
[~, Idx_RV] = max(y_rv_fine); [~, Idx_LV] = max(y_lv_fine);
results.T_peak_RV = t_fine(Idx_RV); results.T_peak_LV = t_fine(Idx_LV);
results.PTT_actual = results.T_peak_LV - results.T_peak_RV;
results.RatioLV = results.AmpLV / results.AmpRV;

%% --- 5. REPORTING ---
fprintf('\n============================================================\n');
fprintf('             AIF ANALYSIS REPORT (FINAL MERGED)\n');
fprintf('============================================================\n');
fprintf('--- [1] FIRST PASS METRICS ---\n');
fprintf('   RV Params     : k = %.2f, theta = %.2f\n', results.kRV, results.thetaRV);
fprintf('   LV Params     : k = %.2f, theta = 2.00 (Fixed)\n', results.kLV);
fprintf('   Peak Time     : RV @ %.2f s  |  LV @ %.2f s\n', results.T_peak_RV, results.T_peak_LV);
fprintf('   PTT           : %.2f s\n', results.PTT_actual);
fprintf('\n');
fprintf('--- [2] RECIRCULATION METRICS ---\n');
fprintf('   RV Params     : k = %.2f, theta = %.2f (Inherited)\n', results.kRV, results.thetaRV);
fprintf('   LV Params     : k = %.2f, theta = 2.00 (Inherited)\n', results.kLV);
fprintf('   Peak Time     : RV @ %.2f s  |  LV @ %.2f s\n', results.T_peak_RV + results.t_offset, results.T_peak_LV + results.t_offset);
fprintf('   PTT           : %.2f s\n', results.PTT_actual);
fprintf('   T-Offset      : %.2f s\n', results.t_offset);
fprintf('   Amplitude     : %.2f%% of First Pass\n', results.alpha * 100);
fprintf('\n');
fprintf('--- [3] GLOBAL FIT STATS ---\n');
fprintf('   R-Squared     : %.4f\n', quality.R_squared);
fprintf('   Amp Ratio     : LV is %.1f%% of RV\n', results.RatioLV * 100);
fprintf('============================================================\n');

%% --- 6. VISUALIZATION ---

% Figure 1: Model Fit & Composition (Master Plot)
fig_handles.fig1_fit_composition = figure('Name', 'Model Fit & Composition');
hold on;

% Transparent Overlays
t_poly = [t_fine; flipud(t_fine)];
zero_base = zeros(size(t_fine));
c_rv = [0 0.4470 0.7410];       % Blue
c_lv = [0.9290 0.6940 0.1250];  % Yellow/Orange
c_re_rv = [0.3010 0.7450 0.9330]; % Light Blue
c_re_lv = [0.4660 0.6740 0.1880]; % Green
c_total = [0.85 0.325 0.098];   % Red

h1 = fill(t_poly, [y_rv_fine; zero_base], c_rv, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
h2 = fill(t_poly, [y_lv_fine; zero_base], c_lv, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
h3 = fill(t_poly, [y_recirc_rv_fine; zero_base], c_re_rv, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
h4 = fill(t_poly, [y_recirc_lv_fine; zero_base], c_re_lv, 'FaceAlpha', 0.4, 'EdgeColor', 'none');

% Lines & Points
h_total = plot(t_fine, results.Predict_total - C_fixed, 'Color', c_total, 'LineWidth', 2.5);
h_raw = plot(t_sec, y_data_subtracted, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);

title(sprintf('Model Fit & Signal Composition (R^2=%.4f)', quality.R_squared));
xlabel('Time (s)'); ylabel('Signal Intensity (Baseline Corrected)');
grid on; xlim([0 max(t_sec)]);

legend([h_raw, h_total, h1, h2, h3, h4], ...
    'Raw Data', 'Total Fit', ...
    'RV First Pass', 'LV First Pass', 'RV Recirc', 'LV Recirc', ...
    'Location', 'northeast');

% Figure 2: Residual Plot (Restored for QC)
fig_handles.fig2_residuals = figure('Name', 'Residual Analysis');
subplot(2,1,1);
plot(t_sec, residuals, 'k.-', 'LineWidth', 1);
yline(0, 'r--');
title('Residuals (Data - Fit)'); xlabel('Time (s)'); ylabel('Difference'); grid on; xlim([0 max(t_sec)]);
subplot(2,1,2);
histogram(residuals, 20, 'FaceColor', [0.4 0.4 0.4]);
title('Residual Distribution'); xlabel('Residual Value'); ylabel('Count'); grid on;

end

%% --- NESTED FUNCTIONS ---
function [p_fp, h] = get_first_pass_initial_guess(t, y, TarrivalRV_ini)
h = [];
options = optimoptions('lsqcurvefit', 'Algorithm', 'levenberg-marquardt', 'Display', 'off');
p0_rv = [max(y), 4, 3, TarrivalRV_ini];
lb_rv = [0,   3.0, 2.5, TarrivalRV_ini-1];
ub_rv = [inf, 8.0, 6.0, TarrivalRV_ini+4];
fitted_params_rv = lsqcurvefit(@(p,t) single_curve_model(p,t,29), p0_rv, t, y, lb_rv, ub_rv, options);
y_fit_rv = single_curve_model(fitted_params_rv, t, 29);
y_residual_lv = y - y_fit_rv;
p0_lv = [max(y_residual_lv), 9.1, 6.48, 40];
lb_lv = [0,   6, 5.0, 20];
ub_lv = [inf, 12, 8.0, 100];
Time_arrival_RV_fixed = fitted_params_rv(4);
model_for_lv = @(p,t) single_curve_model(p, t, Time_arrival_RV_fixed, 2.0, true);
fitted_params_lv = lsqcurvefit(model_for_lv, p0_lv, t, y_residual_lv, lb_lv, ub_lv, options);
p_fp = [fitted_params_rv(1), fitted_params_rv(2), fitted_params_rv(3), fitted_params_lv(1), fitted_params_lv(2), fitted_params_rv(4), fitted_params_lv(3), fitted_params_lv(4)];
end

function [y_total, y_fp, y_rv, y_lv, y_recirc] = full_recirculation_model(p, t)
p_fp = p(1:8);
alpha = p(9);
t_offset = p(10);
[y_fp, y_rv, y_lv] = dual_gamma_model(p_fp, t);
y_recirc = alpha * interp1(t, y_fp, t - t_offset, 'linear', 0);
y_total = y_fp + y_recirc;
end

function [y_total, y_rv, y_lv] = dual_gamma_model(p, t)
p_rv_full = [p(1:3), p(6)];
p_lv_full = [p(4:5), p(7), p(8)];
T1_shared = p(8);
y_rv = single_curve_model(p_rv_full, t, T1_shared);
y_lv = single_curve_model(p_lv_full, t, p(6), 2.0, true);
y_total = y_rv + y_lv;
end

function y = single_curve_model(p, t, varargin)
T1 = varargin{1};
if nargin > 3, is_lv = varargin{3};
    if is_lv, Amp = p(1); k = p(2); delta_t = p(3); T_arr_RV = varargin{1}; T1 = p(4); theta = varargin{2}; T_arr = T_arr_RV + delta_t; end
else, Amp = p(1); k = p(2); theta = p(3); T_arr = p(4);
end
y = zeros(size(t)); idx = t >= T_arr;
if ~any(idx), return;  end; t_shifted = t(idx) - T_arr; t_shifted(t_shifted <= 0) = 1e-9;
gamma_part = (t_shifted.^(k-1)) .* exp(-t_shifted / theta); norm_factor = gamma(k) * theta^k; decay = exp(-t_shifted / T1);
y(idx) = Amp * (gamma_part / norm_factor) .* decay;
end