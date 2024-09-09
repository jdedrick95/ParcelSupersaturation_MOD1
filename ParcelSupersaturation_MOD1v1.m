%% ParcelSupersaturation_MOD1v1.m
% 
% PURPOSE: Calculates parcel-based supersaturations from the Quasi-Steady
% state approximation, and a lognormal, single-mode, internally-mixed
% aerosol population (Abdul-Razzak et al., 1998). 
% 
% INPUTS: None.
% 
% OUTPUTS: None. 
% 
% AUTHOR: Jeramy Dedrick
%         Scripps Institution of Oceanography, La Jolla, CA
%         June 4, 2024
%
%
clear all; close all; clc


%% Set Up Parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Aerosol & Microphysics %%%%%%%%%%%%%%%%%%%%%%%%%%%

N_drop     = 100e6;   % droplet number (m^-3)
rbar_drop  = 10e-6;   % mean drop radius (meters)
w          = 0.5;     % updraft speed (m/s)
N_aero     = 200e6;   % aerosol mode number concentration (m^-3)
r_aero     = 0.06e-6; % aerosol mode radius (m)
sigma_aero = 1.5;     % aerosol mode width (unitless)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Aerosol Hygroscopicity %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% composition densities (Russotto et al., 2013)
rho_NH4SO4   = 1.769;
rho_sea_salt = 2.17;
rho_org      = 1.00;
rho_BC       = 1.7;

% hygroscopicity parameters (Ghan et al., 2001)
B_NH4SO4 = 0.51;
B_NaCl   = 1.16;
B_org    = 0.14;
B_BC     = 5e-7;

% mass fractions of compositions (values based on clean marine accumulation 
% mode, Dedrick et al., 2024; https://doi.org/10.1029/2023GL105798)
f_NH4SO4   = 0.30;
f_sea_salt = 0.60;
f_org      = 0.10;
f_BC       = 0.00;

% B (Solute Effect) Internal Mixture of
% NH4SO4, Sea salt, and Org, Gahn et al., 2001 (EQ. 1)
hygroscopicity_mix_a = sum([(B_NH4SO4 .* f_NH4SO4 ./ rho_NH4SO4);
                          (B_NaCl .* f_sea_salt ./ rho_sea_salt);
                          (B_org .* f_org ./ rho_org);
                          (B_BC .* f_BC ./ rho_BC)]);

hygroscopicity_mix_b = sum([(f_NH4SO4 ./ rho_NH4SO4);
                          (f_sea_salt  ./ rho_sea_salt); 
                          (f_org ./ rho_org); ...
                          (f_BC ./ rho_BC)]);

hygroscopicity = hygroscopicity_mix_a ./ hygroscopicity_mix_b;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Thermodynamic Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%

lv       = 4.5e4;      % latent heat of vaporization, J/mol (AMS Glossary)
cp       = 1005.7;     % specific heat of air at constant pressure, J/kg/K (AMS Glossary)
Ma       = 28.97e-3;   % molar mass of air, kg/mol (https://www.engineeringtoolbox.com/molecular-mass-air-d_679.html)
Mw       = 18.02e-3;   % molar mass of water, kg/mol (https://www.burgettstown.k12.pa.us/cms/lib/PA01000191/Centricity/Domain/155/Chapter%207%20section%203%20web%20notes.pdf)
R        = 8.3144621;  % universal gas constant, J/mol/K (AMS Glossary)
Dv       = 2.55e-5;    % (https://doi.org/10.5194/acp-19-639-2019) molecular diffusion constant (water in air), cm^2/s (https://www.thermopedia.com/content/696/)
g        = 9.81;       % gravity, m/s^2
kv       = 2.40e-2;    % coefficient of thermal conductivity of air, J/m/s/K (S&P 2016)
rho_l    = 1e3;        % density of water, kg/m^3 
cp       = cp .* (Ma); % converts cp to J/mol/K units
sigma_sa = 0.076;      % surface tension of solution wrt air at 0 deg C, kg/s^2 (Ghan et al., 2001) 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Cloud Thermodynamics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T  = 15;       % deg C
P  = 1000;     % hPa

T  = T + 273.15; % K 
P  = P.*(1e2);   % Pascal

% saturation vapor pressure based on Clasius-Clapeyron (August–Roche–Magnus formula wiki https://en.wikipedia.org/wiki/Clausius%E2%80%93Clapeyron_relation#August%E2%80%93Roche%E2%80%93Magnus_formula)
es = 6.1094 .* exp( (17.625 .* (T - 273.15)) ./ ...
                    ((T - 273.15) + 243.04) );
es = es .* 100; % convert to Pascals



%% Calculate Supersaturation

% EQ. (A2) Yang et al., 2019
Q1 = (((lv)./(cp.*T)) - 1) .* ...
           ((Ma.*g)./(R.*T));

% EQ. (A3) Yang et al., 2019)
Q2 = ((lv.^2)./(Mw.*cp.*P.*T)) + ...
           ((R.*T)./(Mw.*es));

% EQ. (A4) Yang et al., 2019
G  = ((rho_l.*R.*T)./(Mw.*Dv.*es)) + ...
     (((rho_l.*lv)./(Mw.*kv.*T)) .* ...
     ((lv./(R.*T)) - 1));

% EQ. (A1) Yang et al., 2019
A = (Q1 .* G) ./ ...
    (4.*pi.*rho_l.*Q2);

% quasi-steady-state supersaturation, EQ. (2) Yang et al., 2019)
s_qs = (A .* w) ./ (N_drop .* rbar_drop);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Kohler Terms  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A (Kelvin Effect) Ghan et al. (2001, EQ. 5)
A_kohler  = (2.*Mw.*sigma_sa)./(R.*T.*rho_l);

% B (Solute Effect)
B_kohler  = hygroscopicity;

% Critical supersaturation Abdul-Razzak et al., 1998 (EQ. 8)
Sm        = (2./sqrt(B_kohler)) .* ...
            ((A_kohler./(3.*r_aero)).^(3/2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% alpha term, EQ. (11) Abdul-Razzak et al., 1998
alpha = ((g.*Mw.*lv)./(cp.*R.*(T.^2))) - ...
              (g.*Ma)./(R.*T);

% gamma term, EQ. (12) Abdul-Razzak et al., 1998
gamma = ((R.*T)./(es.*Mw)) + ...
        ((Mw.*lv.^2)./(cp.*P.*Ma.*T));

% f1 parameter, EQ. (28) Abdul-Razzak et al., 1998
f1 = 1.5 .* exp(2.25 .* (log(sigma_aero).^2));

% f2 parameter, EQ. (29) Abdul-Razzak et al., 1998
f2 = 1 + (0.25 .* log(sigma_aero));

% eta term, EQ. (22) Abdul-Razzak et al., 1998
eta = (((alpha.*w)./((1./G))).^(3/2)) ./ ...
      (2.*pi.*rho_l.*gamma.*N_aero);

% zeta term, EQ. (23) Abdul-Razzak et al., 1998
zeta = (2/3) .* ...
       sqrt((alpha.*w)./((1./G))) .* ...
       (A_kohler);

% Smax, EQ. (31) Abdul-Razzak et al., 1998
Smax_AR98 = Sm ./ ...
            ((f1.*((zeta./eta).^(3/2))) + ...
            (f2.*(((Sm.^2)./(eta)).^(3/4)))).^(1/2);

% smallest activated aerosol size, EQ. (32) Abdul-Razzak et al., 1998
a_c = r_aero .* ...
      ( (f1 .* (zeta/eta)^(3/2)) + ...
        (f2 .* (((Sm.^2)/eta)^(3/4))) ); % meters
a_c = a_c .* 1e6; % micrometers


% convert fraction to percent
s_qs      = s_qs.* 100;
Smax_AR98 = Smax_AR98 .* 100;



disp(strcat('Quasi-Steady State=', num2str(round(s_qs,3)), '%'))
disp(strcat('AR98=', num2str(round(Smax_AR98,3)), '%'))
disp(strcat('smallest activated aerosol radius= ', num2str(round(a_c*1000,2)), ' nm'))

