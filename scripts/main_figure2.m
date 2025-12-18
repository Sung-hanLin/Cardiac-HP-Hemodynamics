%% main_figure2.m
% Generates Figure 2: Hemodynamic simulation of HP pyruvate with cardiac phase misalignment.

clear; clc; close all;
addpath('functions'); 

%% 1. Parameters
dataFile = 'data/LVRV_curvefit_model_0613.mat';
if isfile(dataFile)
    load(dataFile, 'fittedmodel');
else
    error('Data file not found in data/ folder.');
end

TR = 2000;         [cite_start]% ms [cite: 61, 74]
NumTimeP = 80;     
startPhase = 0.48; [cite_start]% ES Target [cite: 76]
T1 = 30;           [cite_start]% s [cite: 64]

HRs = [45, 59, 73, 87, 101]; [cite_start]% Heart Rates for Figure 2 [cite: 70]
PTT_ref = 6.7; 
HR_ref = 59;

%% 2. Initialization
fig_handle = figure('Name', 'Figure 2: Simulation Results', 'Color', 'w', 'Position', [50, 50, 2400, 800]);
tiledlayout(2, 5, 'TileSpacing', 'normal', 'Padding', 'compact');

%% 3. Loop
for i = 1:length(HRs)
    current_HR = HRs(i);
    
    [cite_start]% PTT calculation based on HR [cite: 73]
    PTT = PTT_ref * (HR_ref / current_HR);
    RR_ms = (60 / current_HR) * 1000;
    phase_drift = (TR / RR_ms) - floor(TR / RR_ms); [cite_start]% [cite: 126]
    
    phases = zeros(NumTimeP, 1);
    phases(1) = startPhase;
    for k = 1:NumTimeP-1
        next_ph = phases(k) + phase_drift;
        phases(k+1) = next_ph - floor(next_ph);
    end
    
    results_sim = simulate_hemodynamics(TR, NumTimeP, PTT, T1);
    signal_ideal = results_sim.input_function_Total; 
    vols = fittedmodel(phases);
    norm_vols = vols ./ fittedmodel(startPhase);
    signal_unaligned = signal_ideal .* norm_vols; [cite_start]% Apply volume modulation [cite: 76, 117]
    
    % --- Top Row (Signal) ---
    nexttile(i); hold on;
    p1 = plot(1:NumTimeP, signal_ideal, 'k-', 'LineWidth', 2.0);
    p2 = plot(1:NumTimeP, signal_unaligned, 'b--', 'LineWidth', 1.5);
    title(sprintf('HR=%d bpm (PTT=%.1fs)', current_HR, PTT));
    ylim([0, 1.8]); grid off; box off;
    if i == 1, ylabel('Signal intensity (a.u.)'); legend([p1, p2], {'Aligned','Unaligned'}, 'Box','off'); end
    set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'XTickLabel', []);

    % --- Bottom Row (Phase) ---
    nexttile(i + 5); hold on;
    yline(startPhase, 'r:', 'LineWidth', 2.0);
    plot(1:NumTimeP, phases, '-o', 'Color', 'b', 'MarkerSize', 4);    
    title(sprintf('\\Delta\\phi = %.2f', phase_drift));
    ylim([0, 1]); grid off; box off;
    xlabel('Time Points (TR)');
    if i == 1, ylabel('Cardiac Phase (\phi)'); end
    set(gca, 'FontSize', 14, 'FontWeight', 'bold');
end

% Final aesthetic adjustment
all_text = findall(gcf, '-property', 'FontName');
set(all_text, 'FontName', 'Arial');