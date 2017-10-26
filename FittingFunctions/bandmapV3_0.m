function [fitresult] = bandmapV3_0(u, eachplot)
% peakVisibilityV6_2 without fitting satellite peaks

a=-log(u);
F = a;

% % rescaling of the figure, done automatically in Winview
% in = -log(double(a(:,:,1)));
% a = imrotate(in,-4);
% F = a(800:1025,385:535);
% % figure;
% % imagesc(F);

r = size(F(:,1),1); % the y-image size
index = (1:r)'; % a column vector enumerating the pixels of the image

% Taking a sum to use as the main fitting data
projection = sum(F,2);



%Low pass filtering for thermal fit
fourier = fft(projection); % fourier transform
filter = 80; % filter order, the higher the number the sharper the cutoff is, 50 seems acceptable I wouldn't play with this number
lowpass = real(ifft(ifftshift(fftshift(fourier).*(exp((-linspace(-r/2,r/2,r).^2)/((r/2)/filter)))')));

[lowpassvalue,lowpasscenter] = max(lowpass); % finding the amplitude and center for inital guesses for the fitter

% Now thermal peak fitting for the first guess
[xData0, yData0] = prepareCurveData(index,lowpass);

% Set up fittype and options for the thermal fraction
f0 = fittype('(a*exp(-(1/2)*((x-b)/c)^2))+d','independent','x','dependent','y');
opts = fitoptions(f0);
opts.Display = 'Off';

% fitting paramaters format, alphabetical: [a b c d]
opts.Lower = [0 0 0 -10];
opts.StartPoint = [lowpassvalue lowpasscenter r/2 0];
opts.Upper = [2*max(projection) r 2*r 10];

% Fit model to data.
[peakfit0,~] = fit(xData0,yData0,f0,opts);

fitcoeff0 = coeffvalues(peakfit0); % coefficients used as initial guesses for the final fitter
% curve that will serve as an initial thermal component guess
peaks0 = feval(peakfit0,index);

% guesses for the fittter based on what we've done so far
Cthermal = 1.0.*fitcoeff0(2);
Othermal = 1.0.*fitcoeff0(4);
[Asf1,Csf1] = max(projection);

if fitcoeff0(3) < 15
    Athermal = 0.05*fitcoeff0(1);
    Wthermal = 20;
else
    Athermal = 0.7.*fitcoeff0(1);
    Wthermal = 1.0.*fitcoeff0(3);
    Asf1 = Asf1 - 0.7.*peaks0(Csf1);
end

Wsf1 = Wthermal/3;

% Now fit central and thermal, but not sattelites
[xData1, yData1] = prepareCurveData(index,projection);

% Set up fittype and options for the sf peaks
f1 = fittype('(a*exp(-(1/2)*((x-b3)/c)^2))+d+(a5*exp(-(1/2)*((x-b3)/c5)^2))','independent','x','dependent','y');
opts = fitoptions(f1);
opts.Display = 'Off';

% format: [a a5 b3 c c5 d]
opts.Lower = [0  0 Csf1-r/4 Wthermal/3 Wsf1/3 -5*abs(Othermal)];
opts.StartPoint = [Athermal Asf1 Csf1 Wthermal Wsf1 Othermal];
opts.Upper = [2*max(projection) 2*max(projection) Csf1+r/4 3*Wthermal 1.5*Wsf1 5*abs(Othermal)];

% Fit model to data.
[peakfit1,~] = fit(xData1,yData1,f1,opts);
fitcoeff1 = coeffvalues(peakfit1);

%figure(12), plot(projection); hold on; plot(peakfit1); hold off

peaks1 = feval(peakfit1,index); % the fitted thermal plus sf fraction
resid1 = (projection-peaks1); % fitting residual

Athermal = fitcoeff1(1);
Asf1 = fitcoeff1(2);
Csf1 = fitcoeff1(3);
Wthermal = fitcoeff1(4);
Wsf1 = fitcoeff1(5);
Othermal = fitcoeff1(6);



thermal = Athermal*exp(-(1/2)*((index-Csf1)./Wthermal).^2);
CenterPeak = Asf1*exp(-(1/2)*((index-Csf1)./Wsf1).^2);

baseValue = min(thermal+Othermal);
%condensateFraction = (Ca*Cc + (S1a + S2a)*Sc)/((Ca*Cc + (S1a + S2a)*Sc)+Ta*Tc+offset*r);
condensateFraction = sum(CenterPeak)/sum(CenterPeak+thermal+Othermal-baseValue);

fitresult = [condensateFraction, Wsf1, Wthermal];

if eachplot == 1
    figure;
    subplot(1,2,1);
    plot(index,projection,index,thermal+Othermal,index,CenterPeak,index,resid1-5)
    subplot(1,2,2);
    imagesc(F);
    axis tight
    title(['Coherent fraction: ',num2str(condensateFraction)]);
end

end
