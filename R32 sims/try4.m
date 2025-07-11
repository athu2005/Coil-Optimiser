%% 

% --- Import base workspace ---
modelWorkspace = get_param('sim1', 'ModelWorkspace');

%% 

% --- Inputs ---
T_cond_C = 55;              % Condensation temperature (°C)
T_evap_C = 10;              % Evaporation temperature (°C)
T_env_C = 30;               % Enviroenment air temperature (°C)
subcooling_K = 8.3;         % Subcooling (K)
superheat_K = 11.1;         % Superheat (K)
mass_flow_rate = 57.4/3600; % Mass flow rate (kg/s)
refrigerant = 'R32';        % Refrigerant
discharge_tube_dia = 0.375; % discharge tube diameter (in inches)
suction_tube_dia = 0.375;   % suction tube diameter (in inches)
tonnage = 1.25 ;            % tonnage of cooling
compressor_RPM = 5400;      % Compressor RPM


%% 

% --- CoolProp Setup ---
import py.CoolProp.CoolProp;
CP = @py.CoolProp.CoolProp.PropsSI;


% Convert temperatures to Kelvin
T_cond_K = T_cond_C + 273.15;
T_evap_K = T_evap_C + 273.15;

% Get condenser and evaporator pressures
P_cond = CP('P', 'T', T_cond_K, 'Q', 0, refrigerant);   % Pa
P_evap = CP('P', 'T', T_evap_K, 'Q', 1, refrigerant);   % Pa

% Compressor inlet state (after superheat)
T_inlet_K = T_evap_K + superheat_K;
h1 = CP('H', 'T', T_inlet_K, 'P', P_evap, refrigerant); % J/kg
s1 = CP('S', 'T', T_inlet_K, 'P', P_evap, refrigerant); % J/kg/K

% Isentropic outlet state (same entropy, condenser pressure)
h2s = CP('H', 'P', P_cond, 'S', s1, refrigerant);       % J/kg


T_outlet_K = CP('T', 'P', P_cond, 'S', s1, refrigerant);
T_outlet_K = T_outlet_K + 30;

% Actual compressor outlet enthalpy
h2 = CP('H', 'T', T_outlet_K, 'P', P_cond, refrigerant); % J/kg

% Calculate isentropic efficiency
eta_isentropic = (h2s - h1) / (h2 - h1);  % should be less than 1


%% 

%=========add code here==========



% --- assign inputs ---
assignin(modelWorkspace, 'T_cond', T_cond_C);
assignin(modelWorkspace, 'T_evap', T_evap_C);
assignin(modelWorkspace, 'T_env', T_env_C);
assignin(modelWorkspace, 'P_cond', P_cond);
assignin(modelWorkspace, 'P_evap', P_evap);
assignin(modelWorkspace, 'Subcooling', subcooling_K);
assignin(modelWorkspace, 'Superheat', superheat_K);
assignin(modelWorkspace, 'mass_flow_rate', mass_flow_rate);

assignin(modelWorkspace, 'eta_isentropic', eta_isentropic);


%% 
