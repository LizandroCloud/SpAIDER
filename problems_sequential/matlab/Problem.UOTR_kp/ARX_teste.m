% Clear workspace and close figures
clear; clc; close all;

%% Parameters
Ts = 0.1; % Sample time
t = (0:Ts:10)'; % Time vector

%% Generate Input and Output Signals
% Generating input signal (sine wave with noise)
inputSignal = sin(t) + 0.1 * randn(size(t)); 

% Generating output signal (filtered input signal with noise)
outputSignal = filter([1], [1, -1.5, 0.7], inputSignal) + 0.1 * randn(size(t)); 

%% Create iddata Object
data = iddata(outputSignal, inputSignal, Ts);

%% Visualize Input and Output Signals
figure;
subplot(2,1,1);
plot(t, inputSignal);
title('Input Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,1,2);
plot(t, outputSignal);
title('Output Signal');
xlabel('Time (s)');
ylabel('Amplitude');

%% Identify ARX Model
% Define the orders for the ARX model
na = 2; % Order of the autoregressive part
nb = 2; % Order of the exogenous input part
nk = 0; % Input delay

% Estimate the ARX model
model = arx(data, [na nb nk]);

%% Prepare Future Inputs for Prediction
% Define future inputs for prediction
u_future = sin((t(end)+Ts):Ts:(t(end)+5))'; % Future input signal

% Create a combined dataset for prediction
% Combine past output data and future input data
% Extend the time vector for future predictions
t_future = (t(end)+Ts):Ts:(t(end)+5);
data_future = iddata([], u_future, Ts); 

% Use the last output to initialize predictions
% Create an empty output vector for predictions
y_future = zeros(size(u_future));

% Predict future outputs using the identified model
% Iterate to predict future outputs step by step
for i = 1:length(u_future)
    % Predict the next output based on the previous outputs and future input
    if i == 1
        % Use the last known output for the first prediction
        y_future(i) = predict(model, iddata(outputSignal(end), inputSignal(end), Ts));
    else
        % Use the previous predicted output for subsequent predictions
        y_future(i) = predict(model, iddata(y_future(i-1), u_future(i), Ts));
    end
end

%% Visualize Predictions
figure;
subplot(2,1,1);
plot(t_future, y_future, 'r');
title('Predicted Future Output');
xlabel('Time (s)');
ylabel('Predicted Amplitude');

subplot(2,1,2);
plot(t, outputSignal, 'b');
hold on;
plot(t_future, y_future, 'r--');
title('Output Signal and Predictions');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Actual Output', 'Predicted Future Output');
grid on;
