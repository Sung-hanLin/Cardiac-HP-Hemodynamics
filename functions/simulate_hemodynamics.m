function [results] = simulate_hemodynamics(TR_ms, NumTimeP, PTT, T1)
% SIMULATE_HEMODYNAMICS: Models HP [1-13C]pyruvate kinetics in cardiac ventricles.
%
%   Description:
%     Generates synthetic cardiac pyruvate timecourses based on Gamma-variate AIF
%     and pulmonary transit time (PTT) for digital phantom simulations. [cite: 1100, 1101]
%
%   Inputs:
%     TR_ms      - Repetition time in milliseconds (ms).
%     NumTimeP   - Total number of time points (TRs) to simulate.
%     PTT        - Pulmonary Transit Time in seconds (s). [cite: 1108]
%     T1         - T1 relaxation time of [1-13C]pyruvate in vivo (s). [cite: 1101]
%
%   Outputs:
%     results    - Struct containing 'input_function_Total' (RV + LV signals).
    
    experiment.kPL = 0.02; 
    experiment.std_noise = 0.005;
    acq.TR = TR_ms/1000;
    acq.N = NumTimeP;
    t = (1:acq.N)*acq.TR;

    % --- Right Ventricle (RV) Modeling ---
    Time_arrival_RV = 12; 
    Tp_arrival_RV = Time_arrival_RV/acq.TR;
    k_RV = 2; theta_RV = 1.5;
    input_function_RV = zeros(length(t),1);
    for tt = 1:acq.N
        if tt >= Tp_arrival_RV
            [cite_start]% Gamma-variate function accounting for T1 decay [cite: 64, 67]
            input_function_RV(tt,1) = (1/(gamma(k_RV)*theta_RV^k_RV)) * ...
                ((tt-Tp_arrival_RV)^(k_RV-1)) * exp(-(tt-Tp_arrival_RV)/theta_RV) * ...
                exp(-(tt-Tp_arrival_RV)/T1); 
        end
    end
    input_function_RV = input_function_RV/experiment.std_noise;
    input_function_RV = input_function_RV/max(input_function_RV,[],'all'); 

    % --- Left Ventricle (LV) Modeling ---
    Time_arrival_LV = Time_arrival_RV + PTT;
    Tp_arrival_LV = Time_arrival_LV/acq.TR;
    k_LV = 4; theta_LV = 1.5;
    input_function_LV = zeros(length(t),1);
    for tt = 1:acq.N
        if tt >= Tp_arrival_LV
            input_function_LV(tt,1) = (1/(gamma(k_LV)*theta_LV^k_LV)) * ...
                ((tt-Tp_arrival_LV)^(k_LV-1)) * exp(-(tt-Tp_arrival_LV)/theta_LV) * ...
                exp(-(tt-Tp_arrival_LV)/T1); 
        end
    end
    input_function_LV = 1.6 * input_function_LV/experiment.std_noise;
    input_function_LV = input_function_LV/max(input_function_LV,[],'all'); 
    
    results.input_function_Total = input_function_RV + input_function_LV;
    results.input_function_Total = results.input_function_Total / max(results.input_function_Total,[],'all');
end