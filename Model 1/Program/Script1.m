addpath ./readonly
load('readonly/CellData.mat');
load('readonly/CellModel.mat');

% Estimating SOC based on OCV only
socEstimate1 = SOCfromOCVtemp(voltage,25,model);

% SOC Estimation Plot
t = 0:length(socEstimate1)-1; t = t/60;
plot(t,100*socEstimate1, t,100*soc);
title('True SOC and voltage-based estimate');
xlabel('Time (min)'); ylabel('SOC and estimate (%)');
axis([0 500 0 100]); legend('SOC estimate','True SOC'); grid on;

sqrt(mean((soc-socEstimate1).^2))
%% 

R0=getParamESC('R0param',25,model);

% Estimating Soc based on OCV with series resistance compensation
socEstimate2 = SOCfromOCVtemp(voltage+current*R0,25,model);

% SOC Estimation Plot
t = 0:length(socEstimate2)-1; t = t/60;
plot(t,100*socEstimate2, t,100*soc);
title('True SOC and voltage-based estimate');
xlabel('Time (min)'); ylabel('SOC and estimate (%)');
axis([0 500 0 100]); legend('SOC estimate','True SOC'); grid on;

sqrt(mean((soc-socEstimate2).^2))
