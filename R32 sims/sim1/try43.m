% --- Import base workspace ---
modelWorkspace = get_param('sim1', 'ModelWorkspace');

% --- Inputs ---
T_set_C = 27;               % Target temperature (°C)
T_cond_C = 55;              % Condensation temperature (°C)
T_evap_C = 10;              % Evaporation temperature (°C)
T_env_C = 30;               % Enviroenment air temperature (°C)
subcooling_K = 8.3;         % Subcooling (K)
superheat_K = 11.1;         % Superheat (K)
mass_flow_rate = 57.4/3600; % Mass flow rate (kg/s)
refrigerant = 'R32';        % Refrigerant
discharge_tube_dia = 0.3188976378; % discharge tube diameter (in inches)
discharge_tube_thickness = 0.01; % discharge tube thickness (in inches)
suction_tube_dia = 0.3858267717;   % suction tube diameter (in inches)
suction_tube_thickness = 0.01; % discharge tube thickness (in inches)
tonnage = 1;            % tonnage of cooling
compressor_RPM = 5400;      % Compressor RPM

fin_efficiency = 0.75;

% condensor coil inputs
cond_coil_height = 0.4;    % Condensor coil height (in m)
cond_coil_width = 0.6;    % Condensor coil width (in m)
cond_coil_transverse_pitch = 25.4; % Transverse tube pitch (perpendicular to flow direction) (in mm)
cond_coil_longitudinal_pitch = 22; % Longitudinal tube pitch (along flow direction) (in mm)
cond_circuits_n = 4 ; 
flow_air_cond = 1000; % CPM, example for air side (adjust as needed)
pitch_cond = 1.5/1000 ; % pitch of fins (in m)

% evaporator coil inputs
evap_coil_height = 0.15;    % evaporator coil height (in m)
evap_coil_width = 0.7;    % evaporator coil width (in m)
evap_coil_transverse_pitch = 25.4; % Transverse tube pitch (perpendicular to flow direction) (in mm)
evap_coil_longitudinal_pitch = 22; % Longitudinal tube pitch (along flow direction) (in mm)
evap_circuits_n = 4 ;
flow_air_evap = 245; % CFM, example for air side (adjust as needed)
pitch_evap = 1.5 / 1000; % pitch of fins (in m)


% --- CoolProp Setup ---
import py.CoolProp.CoolProp;
CP = @py.CoolProp.CoolProp.PropsSI;


% Convert temperatures to Kelvin
T_cond_K = T_cond_C + 273.15;
T_evap_K = T_evap_C + 273.15;

% Get condenser and evaporator pressures
P_cond = CP('P', 'T', T_cond_K, 'Q', 0, refrigerant);   % Pa
P_evap = CP('P', 'T', T_evap_K, 'Q', 1, refrigerant);   % Pa

% evaporator
T_evap_subcooled_K = T_cond_K - subcooling_K;  % Subcooled temperature in Kelvin
h_evap = CP('H', 'T', T_evap_subcooled_K, 'P', P_cond, refrigerant); % J/kg

% Compressor inlet state (after superheat)
T_inlet_K = T_evap_K + superheat_K;
h1 = CP('H', 'T', T_inlet_K, 'P', P_evap, refrigerant); % J/kg
s1 = CP('S', 'T', T_inlet_K, 'P', P_evap, refrigerant); % J/kg/K

% Isentropic outlet state (same entropy, condenser pressure)
h2s = CP('H', 'P', P_cond, 'S', s1, refrigerant);       % J/kg


T_outlet_K = CP('T', 'P', P_cond, 'S', s1, refrigerant);
T_outlet_K = T_outlet_K + 21;

% Actual compressor outlet enthalpy
h2 = CP('H', 'T', T_outlet_K, 'P', P_cond, refrigerant); % J/kg

% Calculate isentropic efficiency
eta_isentropic = (h2s - h1) / (h2 - h1);  % should be less than 1

% Calculate logs
logP = log(P_cond / P_evap);
logT = log(T_outlet_K / T_inlet_K);

% Formula for polytropic exponent n
poly_n = logP / (logP - logT);

T_cond_inlet_K = T_outlet_K;
T_cond_outlet_K = T_cond_K + subcooling_K;
T_evap_inlet_K = T_evap_K;
T_evap_outlet_K = T_inlet_K;

% --- Calculate Cycle Performance Parameters ---

% Unit conversion for enthalpies (J/kg to kJ/kg)
h1_kJ_kg = h1 / 1000;
h2_kJ_kg = h2 / 1000;
h_evap_kJ_kg = h_evap / 1000; % Enthalpy at evaporator inlet

% Evaporator Cooling Capacity (kW) - Heat absorbed by refrigerant
Q_evap = mass_flow_rate * (h1_kJ_kg - h_evap_kJ_kg);

% Compressor Work Input (kW) - Energy consumed by compressor
W_comp = mass_flow_rate * (h2_kJ_kg - h1_kJ_kg);

% Condenser Heat Rejection Rate (kW) - Heat rejected by refrigerant
Q_cond = mass_flow_rate * (h2_kJ_kg - h_evap_kJ_kg); % Assuming condenser outlet is at h_evap

% Coefficient of Performance (dimensionless) - COP = Cooling Output / Work Input
COP = Q_evap / W_comp;

% Expected cooling capacity from the tonnage (kW)
tonnage_kW = tonnage * 3.517; % Convert tonnage to kW

% Condensor
% --- Inputs ---
tube_dia_cond_m = discharge_tube_dia * 0.0254; % tube diameter (inches to meters)
discharge_tube_thickness_m = discharge_tube_thickness * 0.0254; % tube diameter (inches to meters)

% For refrigerant (already defined)
mass_flow_refrigerant = mass_flow_rate; % kg/s

% Fluid properties at average temperatures (adjust if needed)
T_air_cond_avg_K = (T_env_C + T_cond_C) / 2 + 273.15;
T_refrigerant_cond_avg_K = (T_cond_inlet_K + T_cond_outlet_K) / 2;

% Get air properties
rho_air_cond = CP('D', 'T', T_air_cond_avg_K, 'P', 101325, 'Air');
mass_flow_air_cond = (flow_air_cond * rho_air_cond) * ((0.3048)^3) / 60; % convert to kg/s
mu_air_cond = CP('V', 'T', T_air_cond_avg_K, 'P', 101325, 'Air'); % dynamic viscosity (Pa.s)
k_air_cond = CP('L', 'T', T_air_cond_avg_K, 'P', 101325, 'Air'); % thermal conductivity (W/mK)
Cp_air_cond = CP('C', 'T', T_air_cond_avg_K, 'P', 101325, 'Air'); % J/kgK

% Get refrigerant properties
P_avg_refrigerant_cond = (P_cond + P_evap) / 2;
rho_refrigerant_cond = CP('D', 'T', T_refrigerant_cond_avg_K, 'P', P_avg_refrigerant_cond, refrigerant);
mu_refrigerant_cond = CP('V', 'T', T_refrigerant_cond_avg_K, 'P', P_avg_refrigerant_cond, refrigerant); % Pa.s
k_refrigerant_cond = CP('L', 'T', T_refrigerant_cond_avg_K, 'P', P_avg_refrigerant_cond, refrigerant); % W/mK
Cp_refrigerant_cond = CP('C', 'T', T_refrigerant_cond_avg_K, 'P', P_avg_refrigerant_cond, refrigerant); % J/kgK

% --- Air side calculations ---
A_air_cond = cond_coil_width * cond_coil_height;
v_air_cond = ( mass_flow_air_cond / double(rho_air_cond) ) / A_air_cond; % m/s

Re_air_cond = (double(rho_air_cond) * v_air_cond * tube_dia_cond_m) / double(mu_air_cond);
Pr_air_cond = (double(Cp_air_cond) * double(mu_air_cond)) / double(k_air_cond);
Nu_air_cond = 0.037 * (Re_air_cond^0.8) * (Pr_air_cond^0.4);

h_air_cond = (Nu_air_cond * double(k_air_cond)) / pitch_cond;

% --- Refrigerant side calculations ---
A_ref_cond = pi * (tube_dia_cond_m^2) / 4;
v_refrigerant_cond = ( mass_flow_refrigerant / double(rho_refrigerant_cond) ) / A_ref_cond; % m/s

Re_ref_cond = (double(rho_refrigerant_cond) * v_refrigerant_cond * tube_dia_cond_m) / double(mu_refrigerant_cond);
Pr_ref_cond = (double(Cp_refrigerant_cond) * double(mu_refrigerant_cond)) / double(k_refrigerant_cond);
Nu_ref_cond = 0.023 * (Re_ref_cond^0.8) * (Pr_ref_cond^0.4);

h_refrigerant_cond = (Nu_ref_cond * double(k_refrigerant_cond)) / tube_dia_cond_m;

% --- Overall heat transfer coefficient ---
h_coeff_cond = 1 / ( (1 / h_air_cond) + (1 / h_refrigerant_cond) );

% Evaporator
% --- Inputs ---
tube_dia_evap_m = suction_tube_dia * 0.0254; % tube diameter (inches to meters)
suction_tube_thickness_m = suction_tube_thickness * 0.0254; % tube diameter (inches to meters)

% For refrigerant (already defined)
mass_flow_refrigerant = mass_flow_rate; % kg/s

% Fluid properties at average temperatures (adjust as needed)
T_air_evap_avg_K = (T_env_C + T_evap_C) / 2 + 273.15;
T_refrigerant_evap_avg_K = (T_evap_inlet_K + T_evap_outlet_K) / 2;

% Get air properties
rho_air_evap = CP('D', 'T', T_air_evap_avg_K, 'P', 101325, 'Air');
mass_flow_air_evap = (flow_air_evap * rho_air_evap) * ((0.3048)^3) / 60; % convert to kg/s
mu_air_evap = CP('V', 'T', T_air_evap_avg_K, 'P', 101325, 'Air'); % dynamic viscosity (Pa.s)
k_air_evap = CP('L', 'T', T_air_evap_avg_K, 'P', 101325, 'Air'); % thermal conductivity (W/mK)
Cp_air_evap = CP('C', 'T', T_air_evap_avg_K, 'P', 101325, 'Air'); % J/kgK

% Get refrigerant properties
P_avg_refrigerant_evap = (P_cond + P_evap) / 2;
rho_refrigerant_evap = CP('D', 'T', T_refrigerant_evap_avg_K, 'P', P_avg_refrigerant_evap, refrigerant);
mu_refrigerant_evap = CP('V', 'T', T_refrigerant_evap_avg_K, 'P', P_avg_refrigerant_evap, refrigerant); % Pa.s
k_refrigerant_evap = CP('L', 'T', T_refrigerant_evap_avg_K, 'P', P_avg_refrigerant_evap, refrigerant); % W/mK
Cp_refrigerant_evap = CP('C', 'T', T_refrigerant_evap_avg_K, 'P', P_avg_refrigerant_evap, refrigerant); % J/kgK

% --- Air side calculations ---
A_air_evap = evap_coil_width * evap_coil_height;
v_air_evap = mass_flow_air_evap / double(rho_air_evap) / A_air_evap; % m/s

Re_air_evap = (double(rho_air_evap) * v_air_evap * tube_dia_evap_m) / double(mu_air_evap);
Pr_air_evap = (double(Cp_air_evap) * double(mu_air_evap)) / double(k_air_evap);
Nu_air_evap = 0.037 * (Re_air_evap^0.8) * (Pr_air_evap^0.4);

h_air_evap = (Nu_air_evap * double(k_air_evap)) / pitch_evap;

% --- Refrigerant side calculations ---
A_ref_evap = pi * (tube_dia_evap_m^2) / 4;
v_refrigerant_evap = mass_flow_refrigerant / double(rho_refrigerant_evap) / A_ref_evap; % m/s

Re_ref_evap = (double(rho_refrigerant_evap) * v_refrigerant_evap * tube_dia_evap_m) / double(mu_refrigerant_evap);
Pr_ref_evap = (double(Cp_refrigerant_evap) * double(mu_refrigerant_evap)) / double(k_refrigerant_evap);
Nu_ref_evap = 0.023 * (Re_ref_evap^0.8) * (Pr_ref_evap^0.4);

h_refrigerant_evap = (Nu_ref_evap * double(k_refrigerant_evap)) / tube_dia_evap_m;

% --- Overall heat transfer coefficient ---
h_coeff_evap = 1 / ( (1 / h_air_evap) + (1 / h_refrigerant_evap) );

% --- Required surface area (m²) for condensor ---
delta_T_cond = abs(T_env_C - T_cond_C);
A_required_cond = (Q_cond * 1000) / (h_coeff_cond * delta_T_cond);

total_length_coil_cond = A_required_cond/(2 * pi * (tube_dia_cond_m)/2);
total_length_coil_cond = round(total_length_coil_cond);

n_sr_cond = (cond_coil_height * 1000) / cond_coil_transverse_pitch;
n_sr_cond = round(n_sr_cond);

n_r_cond = total_length_coil_cond / (n_sr_cond * cond_coil_width);
n_r_cond = round(n_r_cond);

cond_coil_length = n_r_cond * cond_coil_longitudinal_pitch / 1000;

cond_fin_area = cond_coil_length * cond_coil_height;
total_cond_fin_area = cond_fin_area * (cond_coil_width / pitch_cond) / fin_efficiency;

% --- Required surface area (m²) for evaporator ---
delta_T_evap = T_env_C - T_evap_C;
A_required_evap = (Q_evap * 1000) / (h_coeff_evap * delta_T_evap);

total_length_coil_evap = A_required_evap/(2 * pi * (tube_dia_evap_m)/2);
total_length_coil_evap = round(total_length_coil_evap);

% n_r_cond =
n_sr_evap = (evap_coil_height * 1000) / evap_coil_transverse_pitch;
n_sr_evap = round(n_sr_evap);

n_r_evap = total_length_coil_evap / (n_sr_evap * evap_coil_width);
n_r_evap = round(n_r_evap);

evap_coil_length = n_r_evap * evap_coil_longitudinal_pitch / 1000;

evap_fin_area = evap_coil_length * evap_coil_height;
total_evap_fin_area = evap_fin_area * (evap_coil_width / pitch_evap) / fin_efficiency;

%=========add code here==========
T_refrigerant_cond_outlet = T_cond_K - subcooling_K;
rho_refrigerant_cond_outlet = CP('D', 'T', T_refrigerant_cond_outlet, 'P', P_cond, refrigerant);
T_refrigerant_cond_inlet = T_outlet_K;
rho_refrigerant_cond_inlet = CP('D', 'T', T_refrigerant_cond_outlet, 'P', P_cond, refrigerant);
T_refrigerant_evap_outlet = T_inlet_K;
rho_refrigerant_evap_outlet = CP('D', 'T', T_refrigerant_evap_outlet, 'P', P_evap, refrigerant);
T_refrigerant_evap_inlet = T_cond_K - subcooling_K;
rho_refrigerant_evap_inlet = CP('D', 'T', T_refrigerant_evap_inlet, 'P', P_evap, refrigerant);


% --- assign inputs ---
assignin(modelWorkspace, 'T_set', T_set_C);
assignin(modelWorkspace, 'T_cond', T_cond_C);
assignin(modelWorkspace, 'T_evap', T_evap_C);
assignin(modelWorkspace, 'T_env', T_env_C);
assignin(modelWorkspace, 'P_cond', P_cond);
assignin(modelWorkspace, 'H_cond', h2_kJ_kg);
assignin(modelWorkspace, 'P_evap', P_evap);
assignin(modelWorkspace, 'H_evap', h_evap_kJ_kg);
assignin(modelWorkspace, 'Subcooling', subcooling_K);
assignin(modelWorkspace, 'Superheat', superheat_K);
assignin(modelWorkspace, 'ton_of_cooling', tonnage);
assignin(modelWorkspace, 'mass_flow_rate', mass_flow_rate);
assignin(modelWorkspace, 'fin_efficiency', fin_efficiency);

% condensor coil
assignin(modelWorkspace, 'n_r_cond', n_r_cond);
assignin(modelWorkspace, 'n_sr_cond', n_sr_cond);
assignin(modelWorkspace, 'cond_circuits_n', cond_circuits_n);
assignin(modelWorkspace, 'total_length_coil_cond', total_length_coil_cond);
assignin(modelWorkspace, 'tube_dia_cond_m', tube_dia_cond_m);
assignin(modelWorkspace, 'tube_thickness_cond_m', discharge_tube_thickness_m);
assignin(modelWorkspace, 'cond_coil_longitudinal_pitch', cond_coil_longitudinal_pitch);
assignin(modelWorkspace, 'cond_coil_transverse_pitch', cond_coil_transverse_pitch);
assignin(modelWorkspace, 'cond_coil_length', cond_coil_length);
assignin(modelWorkspace, 'cond_coil_width', cond_coil_width);
assignin(modelWorkspace, 'cond_coil_height', cond_coil_height);
assignin(modelWorkspace, 'total_cond_fin_area', total_cond_fin_area);
assignin(modelWorkspace, 'A_air_cond', A_air_cond);
assignin(modelWorkspace, 'rho_cond_outlet', rho_refrigerant_cond_outlet);
assignin(modelWorkspace, 'rho_cond_inlet', rho_refrigerant_cond_inlet);
assignin(modelWorkspace, 'flow_air_cond', flow_air_cond);

% evaporator coil
assignin(modelWorkspace, 'n_r_evap', n_r_evap);
assignin(modelWorkspace, 'n_sr_evap', n_sr_evap);
assignin(modelWorkspace, 'evap_circuits_n', evap_circuits_n);
assignin(modelWorkspace, 'total_length_coil_evap', total_length_coil_evap);
assignin(modelWorkspace, 'tube_dia_evap_m', tube_dia_evap_m);
assignin(modelWorkspace, 'tube_thickness_evap_m', suction_tube_thickness_m);
assignin(modelWorkspace, 'evap_coil_longitudinal_pitch', evap_coil_longitudinal_pitch);
assignin(modelWorkspace, 'evap_coil_transverse_pitch', evap_coil_transverse_pitch);
assignin(modelWorkspace, 'evap_coil_length', evap_coil_length);
assignin(modelWorkspace, 'evap_coil_width', evap_coil_width);
assignin(modelWorkspace, 'evap_coil_height', evap_coil_height);
assignin(modelWorkspace, 'total_evap_fin_area', total_evap_fin_area);
assignin(modelWorkspace, 'A_air_evap', A_air_evap);
assignin(modelWorkspace, 'rho_evap_outlet', rho_refrigerant_evap_outlet);
assignin(modelWorkspace, 'rho_evap_inlet', rho_refrigerant_evap_inlet);
assignin(modelWorkspace, 'flow_air_evap', flow_air_evap);


assignin(modelWorkspace, 'eta_isentropic', eta_isentropic);
assignin(modelWorkspace, 'polytropic_const', poly_n);

% % --- running model for ISEER calculations ---
% 
% % --- BEFORE the loop, pre‑allocate storage ---
% T_outdoor = 24:1:43;
% N = numel(T_outdoor);
% 
% % Cell to hold each sim’s output
% simOutResults = cell(N,1);
% 
% % OPTIONAL: if you’ve previously run sims in this MATLAB session,
% % clear any leftover logsout in the base workspace
% if evalin('base','exist(''logsout'',''var'')')
%     evalin('base','clear logsout')
% end
% 
% % --- RUN the sims, one by one ---
% for i = 1:N
%     % assign your bin hours & env temp as before
%     assignin(modelWorkspace,'bin_hrs', bin_hours(i));
%     assignin(modelWorkspace,'T_env',  T_outdoor(i));
% 
%     % run the model; return workspace outputs
%     simOutResults{i} = sim('sim1', ...
%                            'ReturnWorkspaceOutputs','on', ...
%                            'SaveOutput','on');
% 
%     % now immediately clear the logsout handle so the next sim is fresh
%     % (if you rely on logsout rather than simOut.Results)
%     if evalin('base','exist(''logsout'',''var'')')
%         evalin('base','clear logsout')
%     end
% end
% 
% % --- AFTER the loop: extract & compute ---
% cooling_energy_array  = zeros(N,1);
% comp_work_array       = zeros(N,1);
% 
% for i = 1:N
%     out = simOutResults{i};
% 
%     % pull your time series objects
%     q_evap_ts  = out.logsout.get('torque').Values;
%     w_comp_ts  = out.logsout.get('Compressor.mechancial_power').Values;
% 
%     % integrate
%     cooling_energy_array(i) = trapz(q_evap_ts.Time,  q_evap_ts.Data);
%     comp_work_array(i)       = trapz(w_comp_ts.Time,  w_comp_ts.Data);
% end
% 
% % Grand totals & ISEER
% total_cooling_energy  = sum(cooling_energy_array);
% total_compressor_work = sum(comp_work_array);
% ISEER = total_cooling_energy / total_compressor_work;
% 
% % (Optional) clear simOutResults if you want a completely “clean” workspace
% clear simOutResults
