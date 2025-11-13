clc
%% Step 1
% Constants
DeltaX = ((2.9985-0.02)/(128-1));
Fs = 16000;     % sampling frequency

% Filter numbers
n = (1:128)';  

% Formulas
fp = 8000 .* 10.^(-0.667 .* (128 - n + 1) .* DeltaX);
Qp = 5 + ((5 / (max(n) - 1)) .* (n - 1));
BWp = fp ./ Qp;
thetap = 2 .* pi .* fp ./ Fs;
rp = 1 - (BWp ./ Fs) .* pi;
b1 = 2 .* rp .* cos(thetap);
b2 = rp .^ 2;

% Combine into table
FilterTable = table(n, fp, Qp, BWp, thetap, rp, b1, b2, ...
    'VariableNames', {'Filter_no','fp','Qp','BWp','thetap','rp','b1','b2'});

% Display
disp(FilterTable)

%% Step 2

numFilters = height(FilterTable);
impRes = cell(numFilters,1);
freqRes = cell(numFilters,1);
freqVectors = cell(numFilters,1);

freqSelection = cell(2,1);

for k = 1:numFilters
    b1k = FilterTable.b1(k);
    b2k = FilterTable.b2(k);
    
    % (1 - z^-2) / (1 - b1*z^-1 + b2*z^-2)
    b = [1 0 -1];
    a = [1 -b1k b2k];
    
    % truncation
    len = 160;
    
    [h, t] = impz(b, a);

    if length(h) < 160
        h = [h; zeros(160-length(h), 1)];
    elseif length(h) > 160
        h = h(1:160);
    end

    impRes{k} = h;

    [H,f] = freqz(b,a, len, Fs);
    freqRes{k} = H;
    freqVectors{k} = f;

    if k == 45
        freqSelection{1} = freqz(h,1, len, Fs);
    elseif k == 65
        freqSelection{2} = freqz(h, 1, len, Fs);
    end

 end

figure
subplot(2,3,1)
plot(impRes{45})
xlim([0 200])

subplot(2,3,4)
plot(impRes{65})
xlim([0 200])

subplot(2,3,2)
plot(freqVectors{100}, 20*log10(abs(freqRes{100})))
subplot(2,3,5)
plot(freqVectors{100}, unwrap(angle(freqRes{100})))

length(impRes{45})

subplot(2,3,[3,6])
plot(freqVectors{45}, 20*log10(abs(freqSelection{1})))
hold on
plot(freqVectors{45}, 20*log10(abs(freqRes{45})))
plot(freqVectors{65}, 20*log10(abs(freqSelection{2})))
plot(freqVectors{65}, 20*log10(abs(freqRes{65})))
hold off



