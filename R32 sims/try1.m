%--- Inputs ---
T_cond_C = 55;       % Condensation temperature (Celsius)
T_evap_C = 10;       % Evaporation temperature (Celsius)
subcooling_K = 8.3;     % Subcooling (Kelvin)
superheat_K = 11.1;      % Superheat (Kelvin)
mass_flow_rate = 56.4/3600; % Mass flow rate (kg/s)
refrigerant = 'R32';  % Refrigerant (string)

%--- CoolProp Setup ---
import py.CoolProp.CoolProp;
CP = @py.CoolProp.CoolProp.PropsSI;

%--- Temperature Conversions ---
T_cond_K = T_cond_C + 273.15;
T_evap_K = T_evap_C + 273.15;

%--- 1. Thermodynamic Properties at Key Points ---

%--- Evaporator Outlet (Compressor Inlet) ---
T_comp_in_K = T_evap_K + superheat_K;
P_comp_in = CP('P', 'T', T_evap_K, 'Q', 1, refrigerant); % Pressure at saturation
h_comp_in = CP('H', 'P', P_comp_in, 'T', T_comp_in_K, refrigerant);
s_comp_in = CP('S', 'P', P_comp_in, 'T', T_comp_in_K, refrigerant);

%--- Condenser Outlet (Expansion Valve Inlet) ---
T_cond_out_K = T_cond_K - subcooling_K;
P_cond_out = CP('P', 'T', T_cond_K, 'Q', 1, refrigerant); % Pressure at saturation
h_cond_out = CP('H', 'P', P_cond_out, 'T', T_cond_out_K, refrigerant);

%--- Compressor Outlet ---
P_comp_out = P_cond_out; % Assuming condenser and compressor outlet pressures are the same
h_comp_out_isentropic = CP('H', 'P', P_comp_out, 'S', s_comp_in, refrigerant);

%Estimate actual compressor outlet enthalpy (using an assumed isentropic efficiency)
%This is a simplification, in a real system you would get this from the compressor model
eta_isentropic = 0.8; % Assumed isentropic efficiency
h_comp_out_actual = h_comp_in + (h_comp_out_isentropic - h_comp_in) / eta_isentropic;

%--- 2. Coil Calculations ---

%--- Evaporator ---

%Evaporator inlet enthalpy (assuming saturated liquid at the expansion valve outlet)
h_evap_in = CP('H', 'P', P_comp_in, 'Q', 0, refrigerant);

%Refrigerant Heat Absorption
Q_evap_refrig = mass_flow_rate * (h_comp_in - h_evap_in);

fprintf('Evaporator Refrigerant Heat Absorption: %.4f W\n', Q_evap_refrig);

%----debug-----
fprintf('\nSanity Check:\n');
fprintf('h_evap_in: %.2f\n', h_evap_in);
fprintf('h_comp_in: %.2f\n', h_comp_in);
fprintf('h_comp_out_isentropic: %.2f\n', h_comp_out_isentropic);
fprintf('h_comp_out_actual: %.2f\n', h_comp_out_actual);
fprintf('h_cond_out: %.2f\n\n', h_cond_out);

fprintf('Q_evap: %.2f\n', Q_evap_refrig);
fprintf('W_compressor_actual: %.2f\n', W_compressor_actual);
fprintf('Q_cond: %.2f\n', Q_cond_refrig);
fprintf('Q_evap + W_comp = %.2f\n', Q_evap_refrig + W_compressor_actual);





%--- Condenser ---

%Refrigerant Heat Rejection
Q_cond_refrig = mass_flow_rate * (h_comp_out_actual - h_cond_out);

fprintf('Condenser Refrigerant Heat Rejection: %.4f W\n', Q_cond_refrig);

%--- 3. Compressor Work ---

%Compressor Work (Actual)
W_compressor_actual = mass_flow_rate * (h_comp_out_actual - h_comp_in);

%Compressor Work (Isentropic)
W_compressor_isentropic = mass_flow_rate * (h_comp_out_isentropic - h_comp_in);

fprintf('Compressor Work (Actual): %.4f W\n', W_compressor_actual);
fprintf('Compressor Work (Isentropic): %.4f W\n', W_compressor_isentropic);

%--- 4. Isentropic Efficiency ---

%Calculate Isentropic Efficiency
eta_compressor = W_compressor_isentropic / W_compressor_actual;

fprintf('Compressor Isentropic Efficiency: %.4f\n', eta_compressor);

%--- 5. Coefficient of Performance (COP) ---

%Coefficient of Performance (COP)
COP = Q_evap_refrig / W_compressor_actual;

fprintf('COP: %.4f\n', COP);