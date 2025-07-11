classdef ThermoCycleModel < matlab.System
    % ThermoCycleModel: A System block to compute cooling cycle parameters
    % and assign all internal variables to the Model Workspace

    % Public, tunable properties
    properties
        % --- Inputs ---
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
        
        fin_efficiency = 0.75; % fin efficiency
        

        cond_coil_height = 0.4;    % Condensor coil height (in m)
        cond_coil_width = 0.6;    % Condensor coil width (in m)
        cond_coil_transverse_pitch = 25.4 % Transverse tube pitch (perpendicular to flow direction) (in mm)
        cond_coil_longitudinal_pitch = 22 % Longitudinal tube pitch (along flow direction) (in mm)
        cond_circuits_n = 4 ; 
        flow_air_cond = 1000; % CPM, example for air side (adjust as needed)
        pitch_cond = 1.5/1000 ; % pitch of fins (in m)
        

        evap_coil_height = 0.15;    % evaporator coil height (in m)
        evap_coil_width = 0.7;    % evaporator coil width (in m)
        evap_coil_transverse_pitch = 25.4 % Transverse tube pitch (perpendicular to flow direction) (in mm)
        evap_coil_longitudinal_pitch = 22 % Longitudinal tube pitch (along flow direction) (in mm)
        evap_circuits_n = 4 ;
        flow_air_evap = 245; % CFM, example for air side (adjust as needed)
        pitch_evap = 1.5 / 1000; % pitch of fins (in m)
    end

    methods (Access = protected)

        function setupImpl(obj)
            % One-time setup
            import py.CoolProp.CoolProp;
            obj.CP = @py.CoolProp.CoolProp.PropsSI;
        end

        function [COP, Q_evap, Q_cond, eta_isentropic] = stepImpl(obj)

            % Constants
            obj.refrigerant = obj.refrigerant;
            obj.CP = obj.CP;

            % Get model workspace
            modelName = bdroot(gcs); % Get top-level model name
            modelWorkspace = get_param(modelName, 'ModelWorkspace');

            % Convert temperatures to Kelvin
            T_cond_K = obj.T_cond_C + 273.15;
            T_evap_K = obj.T_evap_C + 273.15;
            
            % Get condenser and evaporator pressures
            P_cond = obj.CP('P', 'T', T_cond_K, 'Q', 0, obj.refrigerant);   % Pa
            P_evap = obj.CP('P', 'T', T_evap_K, 'Q', 1, obj.refrigerant);   % Pa
            
            % evaporator
            T_evap_subcooled_K = T_cond_K - obj.subcooling_K;  % Subcooled temperature in Kelvin
            h_evap = obj.CP('H', 'T', T_evap_subcooled_K, 'P', P_cond, obj.refrigerant); % J/kg
            
            % Compressor inlet state (after superheat)
            T_inlet_K = T_evap_K + obj.superheat_K;
            h1 = obj.CP('H', 'T', T_inlet_K, 'P', P_evap, obj.refrigerant); % J/kg
            s1 = obj.CP('S', 'T', T_inlet_K, 'P', P_evap, obj.refrigerant); % J/kg/K
            
            % Isentropic outlet state (same entropy, condenser pressure)
            h2s = obj.CP('H', 'P', P_cond, 'S', s1, obj.refrigerant);       % J/kg
            
            
            T_outlet_K = obj.CP('T', 'P', P_cond, 'S', s1, obj.refrigerant);
            T_outlet_K = T_outlet_K + 21;
            
            % Actual compressor outlet enthalpy
            h2 = obj.CP('H', 'T', T_outlet_K, 'P', P_cond, obj.refrigerant); % J/kg
            
            % Calculate isentropic efficiency
            eta_isentropic = (h2s - h1) / (h2 - h1);  % should be less than 1
            
            % Calculate logs
            logP = log(P_cond / P_evap);
            logT = log(T_outlet_K / T_inlet_K);
            
            % Formula for polytropic exponent n
            poly_n = logP / (logP - logT);
            
            T_cond_inlet_K = T_outlet_K;
            T_cond_outlet_K = T_cond_K + obj.subcooling_K;
            T_evap_inlet_K = T_evap_K;
            T_evap_outlet_K = T_inlet_K;
            
            % --- Calculate Cycle Performance Parameters ---
            
            % Unit conversion for enthalpies (J/kg to kJ/kg)
            h1_kJ_kg = h1 / 1000;
            h2_kJ_kg = h2 / 1000;
            h_evap_kJ_kg = h_evap / 1000; % Enthalpy at evaporator inlet
            
            % Evaporator Cooling Capacity (kW) - Heat absorbed by refrigerant
            Q_evap = obj.mass_flow_rate * (h1_kJ_kg - h_evap_kJ_kg);
            
            % Compressor Work Input (kW) - Energy consumed by compressor
            W_comp = obj.mass_flow_rate * (h2_kJ_kg - h1_kJ_kg);
            
            % Condenser Heat Rejection Rate (kW) - Heat rejected by refrigerant
            Q_cond = obj.mass_flow_rate * (h2_kJ_kg - h_evap_kJ_kg); % Assuming condenser outlet is at h_evap
            
            % Coefficient of Performance (dimensionless) - COP = Cooling Output / Work Input
            COP = Q_evap / W_comp;
            
            % Expected cooling capacity from the tonnage (kW)
            tonnage_kW = obj.tonnage * 3.517; % Convert tonnage to kW
            
            % Condensor
            % --- Inputs ---
            tube_dia_cond_m = obj.discharge_tube_dia * 0.0254; % tube diameter (inches to meters)
            discharge_tube_thickness_m = obj.discharge_tube_thickness * 0.0254; % tube diameter (inches to meters)
            
            % For refrigerant (already defined)
            mass_flow_refrigerant = obj.mass_flow_rate; % kg/s
            
            % Fluid properties at average temperatures (adjust if needed)
            T_air_cond_avg_K = (obj.T_env_C + obj.T_cond_C) / 2 + 273.15;
            T_refrigerant_cond_avg_K = (T_cond_inlet_K + T_cond_outlet_K) / 2;
            
            % Get air properties
            rho_air_cond = obj.CP('D', 'T', T_air_cond_avg_K, 'P', 101325, 'Air');
            mass_flow_air_cond = (obj.flow_air_cond * rho_air_cond) * ((0.3048)^3) / 60; % convert to kg/s
            mu_air_cond = obj.CP('V', 'T', T_air_cond_avg_K, 'P', 101325, 'Air'); % dynamic viscosity (Pa.s)
            k_air_cond = obj.CP('L', 'T', T_air_cond_avg_K, 'P', 101325, 'Air'); % thermal conductivity (W/mK)
            Cp_air_cond = obj.CP('C', 'T', T_air_cond_avg_K, 'P', 101325, 'Air'); % J/kgK
            
            % Get refrigerant properties
            P_avg_refrigerant_cond = (P_cond + P_evap) / 2;
            rho_refrigerant_cond = obj.CP('D', 'T', T_refrigerant_cond_avg_K, 'P', P_avg_refrigerant_cond, obj.refrigerant);
            mu_refrigerant_cond = obj.CP('V', 'T', T_refrigerant_cond_avg_K, 'P', P_avg_refrigerant_cond, obj.refrigerant); % Pa.s
            k_refrigerant_cond = obj.CP('L', 'T', T_refrigerant_cond_avg_K, 'P', P_avg_refrigerant_cond, obj.refrigerant); % W/mK
            Cp_refrigerant_cond = obj.CP('C', 'T', T_refrigerant_cond_avg_K, 'P', P_avg_refrigerant_cond, obj.refrigerant); % J/kgK
            
            % --- Air side calculations ---
            A_air_cond = obj.cond_coil_width * obj.cond_coil_height;
            v_air_cond = ( mass_flow_air_cond / double(rho_air_cond) ) / A_air_cond; % m/s
            
            Re_air_cond = (double(rho_air_cond) * v_air_cond * tube_dia_cond_m) / double(mu_air_cond);
            Pr_air_cond = (double(Cp_air_cond) * double(mu_air_cond)) / double(k_air_cond);
            Nu_air_cond = 0.037 * (Re_air_cond^0.8) * (Pr_air_cond^0.4);
            
            h_air_cond = (Nu_air_cond * double(k_air_cond)) / obj.pitch_cond;
            
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
            tube_dia_evap_m = obj.suction_tube_dia * 0.0254; % tube diameter (inches to meters)
            suction_tube_thickness_m = obj.suction_tube_thickness * 0.0254; % tube diameter (inches to meters)
            
            % For refrigerant (already defined)
            mass_flow_refrigerant = obj.mass_flow_rate; % kg/s
            
            % Fluid properties at average temperatures (adjust as needed)
            T_air_evap_avg_K = (obj.T_env_C + obj.T_evap_C) / 2 + 273.15;
            T_refrigerant_evap_avg_K = (T_evap_inlet_K + T_evap_outlet_K) / 2;
            
            % Get air properties
            rho_air_evap = obj.CP('D', 'T', T_air_evap_avg_K, 'P', 101325, 'Air');
            mass_flow_air_evap = (obj.flow_air_evap * rho_air_evap) * ((0.3048)^3) / 60; % convert to kg/s
            mu_air_evap = obj.CP('V', 'T', T_air_evap_avg_K, 'P', 101325, 'Air'); % dynamic viscosity (Pa.s)
            k_air_evap = obj.CP('L', 'T', T_air_evap_avg_K, 'P', 101325, 'Air'); % thermal conductivity (W/mK)
            Cp_air_evap = obj.CP('C', 'T', T_air_evap_avg_K, 'P', 101325, 'Air'); % J/kgK
            
            % Get refrigerant properties
            P_avg_refrigerant_evap = (P_cond + P_evap) / 2;
            rho_refrigerant_evap = obj.CP('D', 'T', T_refrigerant_evap_avg_K, 'P', P_avg_refrigerant_evap, obj.refrigerant);
            mu_refrigerant_evap = obj.CP('V', 'T', T_refrigerant_evap_avg_K, 'P', P_avg_refrigerant_evap, obj.refrigerant); % Pa.s
            k_refrigerant_evap = obj.CP('L', 'T', T_refrigerant_evap_avg_K, 'P', P_avg_refrigerant_evap, obj.refrigerant); % W/mK
            Cp_refrigerant_evap = obj.CP('C', 'T', T_refrigerant_evap_avg_K, 'P', P_avg_refrigerant_evap, obj.refrigerant); % J/kgK
            
            % --- Air side calculations ---
            A_air_evap = obj.evap_coil_width * obj.evap_coil_height;
            v_air_evap = mass_flow_air_evap / double(rho_air_evap) / A_air_evap; % m/s
            
            Re_air_evap = (double(rho_air_evap) * v_air_evap * tube_dia_evap_m) / double(mu_air_evap);
            Pr_air_evap = (double(Cp_air_evap) * double(mu_air_evap)) / double(k_air_evap);
            Nu_air_evap = 0.037 * (Re_air_evap^0.8) * (Pr_air_evap^0.4);
            
            h_air_evap = (Nu_air_evap * double(k_air_evap)) / obj.pitch_evap;
            
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
            delta_T_cond = abs(obj.T_env_C - obj.T_cond_C);
            A_required_cond = (Q_cond * 1000) / (h_coeff_cond * delta_T_cond);
            
            total_length_coil_cond = A_required_cond/(2 * pi * (tube_dia_cond_m)/2);
            total_length_coil_cond = round(total_length_coil_cond);
            
            n_sr_cond = (obj.cond_coil_height * 1000) / obj.cond_coil_transverse_pitch;
            n_sr_cond = round(n_sr_cond);
            
            n_r_cond = total_length_coil_cond / (n_sr_cond * obj.cond_coil_width);
            n_r_cond = round(n_r_cond);
            
            cond_coil_length = n_r_cond * obj.cond_coil_longitudinal_pitch / 1000;
            
            cond_fin_area = cond_coil_length * obj.cond_coil_height;
            total_cond_fin_area = cond_fin_area * (obj.cond_coil_width / obj.pitch_cond) / obj.fin_efficiency;
            
            % --- Required surface area (m²) for evaporator ---
            delta_T_evap = obj.T_env_C - obj.T_evap_C;
            A_required_evap = (Q_evap * 1000) / (h_coeff_evap * delta_T_evap);
            
            total_length_coil_evap = A_required_evap/(2 * pi * (tube_dia_evap_m)/2);
            total_length_coil_evap = round(total_length_coil_evap);
            
            % n_r_cond =
            n_sr_evap = (obj.evap_coil_height * 1000) / obj.evap_coil_transverse_pitch;
            n_sr_evap = round(n_sr_evap);
            
            n_r_evap = total_length_coil_evap / (n_sr_evap * obj.evap_coil_width);
            n_r_evap = round(n_r_evap);
            
            evap_coil_length = n_r_evap * obj.evap_coil_longitudinal_pitch / 1000;
            
            evap_fin_area = evap_coil_length * obj.evap_coil_height;
            total_evap_fin_area = evap_fin_area * (obj.evap_coil_width / obj.pitch_evap) / obj.fin_efficiency;
            
            %=========add code here==========
            T_refrigerant_cond_outlet = T_cond_K - obj.subcooling_K;
            rho_refrigerant_cond_outlet = obj.CP('D', 'T', T_refrigerant_cond_outlet, 'P', P_cond, obj.refrigerant);
            T_refrigerant_cond_inlet = T_outlet_K;
            rho_refrigerant_cond_inlet = obj.CP('D', 'T', T_refrigerant_cond_outlet, 'P', P_cond, obj.refrigerant);
            T_refrigerant_evap_outlet = T_inlet_K;
            rho_refrigerant_evap_outlet = obj.CP('D', 'T', T_refrigerant_evap_outlet, 'P', P_evap, obj.refrigerant);
            T_refrigerant_evap_inlet = T_cond_K - obj.subcooling_K;
            rho_refrigerant_evap_inlet = obj.CP('D', 'T', T_refrigerant_evap_inlet, 'P', P_evap, obj.refrigerant);


            %===============================
            % Assign to Model Workspace
            %===============================
            % --- assign inputs ---
            assignin(modelWorkspace, 'T_cond', obj.T_cond_C);
            assignin(modelWorkspace, 'T_evap', obj.T_evap_C);
            assignin(modelWorkspace, 'T_env', obj.T_env_C);
            assignin(modelWorkspace, 'P_cond', P_cond);
            assignin(modelWorkspace, 'H_cond', h2_kJ_kg);
            assignin(modelWorkspace, 'P_evap', P_evap);
            assignin(modelWorkspace, 'H_evap', h_evap_kJ_kg);
            assignin(modelWorkspace, 'Subcooling', obj.subcooling_K);
            assignin(modelWorkspace, 'Superheat', obj.superheat_K);
            assignin(modelWorkspace, 'ton_of_cooling', obj.tonnage);
            assignin(modelWorkspace, 'mass_flow_rate', obj.mass_flow_rate);
            assignin(modelWorkspace, 'fin_efficiency', obj.fin_efficiency);
            
            % condensor coil
            assignin(modelWorkspace, 'n_r_cond', n_r_cond);
            assignin(modelWorkspace, 'n_sr_cond', n_sr_cond);
            assignin(modelWorkspace, 'cond_circuits_n', obj.cond_circuits_n);
            assignin(modelWorkspace, 'total_length_coil_cond', total_length_coil_cond);
            assignin(modelWorkspace, 'tube_dia_cond_m', tube_dia_cond_m);
            assignin(modelWorkspace, 'tube_thickness_cond_m', discharge_tube_thickness_m);
            assignin(modelWorkspace, 'cond_coil_longitudinal_pitch', obj.cond_coil_longitudinal_pitch);
            assignin(modelWorkspace, 'cond_coil_transverse_pitch', obj.cond_coil_transverse_pitch);
            assignin(modelWorkspace, 'cond_coil_length', cond_coil_length);
            assignin(modelWorkspace, 'cond_coil_width', obj.cond_coil_width);
            assignin(modelWorkspace, 'cond_coil_height', obj.cond_coil_height);
            assignin(modelWorkspace, 'total_cond_fin_area', total_cond_fin_area);
            assignin(modelWorkspace, 'A_air_cond', A_air_cond);
            assignin(modelWorkspace, 'rho_cond_outlet', rho_refrigerant_cond_outlet);
            assignin(modelWorkspace, 'rho_cond_inlet', rho_refrigerant_cond_inlet);
            assignin(modelWorkspace, 'flow_air_cond', obj.flow_air_cond);
            
            % evaporator coil
            assignin(modelWorkspace, 'n_r_evap', n_r_evap);
            assignin(modelWorkspace, 'n_sr_evap', n_sr_evap);
            assignin(modelWorkspace, 'evap_circuits_n', obj.evap_circuits_n);
            assignin(modelWorkspace, 'total_length_coil_evap', total_length_coil_evap);
            assignin(modelWorkspace, 'tube_dia_evap_m', tube_dia_evap_m);
            assignin(modelWorkspace, 'tube_thickness_evap_m', suction_tube_thickness_m);
            assignin(modelWorkspace, 'evap_coil_longitudinal_pitch', obj.evap_coil_longitudinal_pitch);
            assignin(modelWorkspace, 'evap_coil_transverse_pitch', obj.evap_coil_transverse_pitch);
            assignin(modelWorkspace, 'evap_coil_length', evap_coil_length);
            assignin(modelWorkspace, 'evap_coil_width', obj.evap_coil_width);
            assignin(modelWorkspace, 'evap_coil_height', obj.evap_coil_height);
            assignin(modelWorkspace, 'total_evap_fin_area', total_evap_fin_area);
            assignin(modelWorkspace, 'A_air_evap', A_air_evap);
            assignin(modelWorkspace, 'rho_evap_outlet', rho_refrigerant_evap_outlet);
            assignin(modelWorkspace, 'rho_evap_inlet', rho_refrigerant_evap_inlet);
            assignin(modelWorkspace, 'flow_air_evap', obj.flow_air_evap);
            
            
            assignin(modelWorkspace, 'eta_isentropic', eta_isentropic);
            assignin(modelWorkspace, 'polytropic_const', poly_n);

            % Optionally, you can log a message:
            % disp('[ThermoCycleModel] Variables updated in Model Workspace.');
        end

        function resetImpl(~)
            % No states to reset
        end

        function num = getNumInputsImpl(~)
            num = 0; % []
        end

        function num = getNumOutputsImpl(~)
            num = 4; % [COP, Q_evap, Q_cond, eta_isentropic]
        end

        function [o1, o2, o3, o4] = getOutputSizeImpl(~)
            o1 = [1 1]; o2 = [1 1]; o3 = [1 1]; o4 = [1 1];
        end

        function [o1, o2, o3, o4] = getOutputDataTypeImpl(~)
            o1 = 'double'; o2 = 'double'; o3 = 'double'; o4 = 'double';
        end

        function [o1, o2, o3, o4] = isOutputComplexImpl(~)
            o1 = false; o2 = false; o3 = false; o4 = false;
        end

        function [o1, o2, o3, o4] = isOutputFixedSizeImpl(~)
            o1 = true; o2 = true; o3 = true; o4 = true;
        end

    end

    properties (Access = private)
        CP % CoolProp handle
    end
end
