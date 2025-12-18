%% MAIN SCRIPT: Cardiac HP-13C MRS Reconstruction Pipeline
%  Project: Hemodynamic modeling and phase-resolved reconstruction of 
%           hyperpolarized cardiac 13C MRS using 1H cine and ECG timing.
%  Author: Sung-Han Lin, et al. (UT Southwestern Medical Center) [cite: 1039]
%  Journal: Magnetic Resonance in Medicine [cite: 1051]
%
%  Description:
%    This script performs the full analysis pipeline for in-vivo human datasets:
%    1. Loads HP-MRS pyruvate timecourses and subject-specific volume models.
%    2. Executes multi-stage hemodynamic fitting (RV, LV, and Recirculation).
%    3. Performs retrospective cardiac phase-alignment using 1H cine templates. [cite: 1214]
%    4. Exports analysis results and high-resolution diagnostic figures.
%
%  Required Data Files (should be in the same folder or path):
%    - Signal_restoration_inj2_1121.mat (Raw HP-MRS data)
%    - Fitting_models_02252025.mat (Ventricular volume cfit models)
%
%  MATLAB Toolbox Requirements:
%    - Optimization Toolbox & Curve Fitting Toolbox

clc; clear; close all;

% --- 1. Load Data ---
load('Signal_restoration_inj2_1121.mat');
load('Fitting_models_02252025.mat');

% Define Inputs
y_data = p_Pyr;
t_ms = Time_RF * 1000;
cardiac_phase = CardiacPhase_RF;

% User Parameters
TarrivalRV_ini = 7.0;           % (seconds)
early_data_endpoint_sec = 18.0; % (seconds) - for initial guess
full_fit_endpoint_sec = 70.0;   % (seconds) - range for final fit
target_phase = 0.838;           % (0-1) Target Phase

%% --- 2. Execute Fitting ---
fprintf('Starting AIF Analysis (Robust Multi-Stage Method)...\n');

% Call the modified function
[fit_results, quality, handles] = ...
    fit_AIF_recirculation_Final_Fixed(t_ms, y_data, TarrivalRV_ini, early_data_endpoint_sec, full_fit_endpoint_sec);

%% --- 3. Signal Reconstruction (Requirement #4) ---
fprintf('Starting Signal Reconstruction...\n');

% Call reconstruction (logic unchanged)
[reconstructed_signal, ideal_total_curve, ~, handles_recon, ~] = ...
    reconstruct_signal_by_phase_SP(t_ms, y_data, cardiac_phase, ...
    LV_fitted_02252025, RV_fitted_02252025, ...
    fit_results, target_phase);

% Combine handles
handles.fig3_reconstruction = handles_recon.reconstruction;
if isfield(handles_recon, 'volume_profile')
    handles.fig0_volume_profile = handles_recon.volume_profile;
end

%% --- 4. Save Outputs (Requirement #1) ---
output_basename = ['analysis_robust_', datestr(now, 'yyyymmdd_HHMMSS')];

% Save Workspace
save([output_basename, '.mat'], ...
    'fit_results', 'quality', 'reconstructed_signal', 't_ms', 'y_data', ...
    'cardiac_phase', 'target_phase');
fprintf('Workspace results saved to: %s.mat\n', output_basename);

% Save Figures
fprintf('Saving Figures...\n');
fn = fieldnames(handles);
for i = 1:numel(fn)
    if isgraphics(handles.(fn{i}))
        fig_name = [output_basename, '_', fn{i}];
        savefig(handles.(fn{i}), [fig_name, '.fig']);
        exportgraphics(handles.(fn{i}), [fig_name, '.png'], 'Resolution', 600);
        fprintf('   Saved: %s (.fig/.png)\n', fig_name);
    end
end
fprintf('Processing Complete.\n');