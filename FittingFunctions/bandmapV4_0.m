function [fitresult] = bandmapV4_0(u, eachplot)
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
f0 = fittype('B * exp(cos(pi*(x-center)/xt)/T) * (heaviside((x-center)+xt)-heaviside((x-center)-xt)) + d','independent','x','dependent','y');
opts = fitoptions(f0);
opts.Display = 'Off';

% fitting paramaters format: [B T center d xt]
opts.Lower = [0 0 0 -10 0];
opts.StartPoint = [lowpassvalue 1/10 lowpasscenter 0 r/2];
opts.Upper = [2*max(projection) inf r 10 inf];

% Fit model to data.
[peakfit0,~] = fit(xData0,yData0,f0,opts);

% plot(peakfit0, 1:r,projection)

fitcoeff0 = coeffvalues(peakfit0); % coefficients used as initial guesses for the final fitter

% guesses for the fittter based on what we've done so far
BGuess = fitcoeff0(1);
TGuess = fitcoeff0(2);
centerGuess = fitcoeff0(3);
offsetGuess = fitcoeff0(4);
xtGuess = fitcoeff0(5);


% if fitcoeff0(3) < 15
%     Athermal = 0.05*fitcoeff0(1);
%     Wthermal = 20;
% else
%     Athermal = 0.7.*fitcoeff0(1);
%     Wthermal = 1.0.*fitcoeff0(3);
%     Asf1 = Asf1 - 0.7.*peaks0(Csf1);
% end

% Wsf1 = Wthermal/3;

% Now fit central and thermal
[xData1, yData1] = prepareCurveData(index,projection);

% Set up fittype and options for the final fit
f1 = fittype('A * (1-((x-center)/x0)^2)^2 * (heaviside((x-center)+x0)-heaviside((x-center)-x0)) + B * exp(cos(pi*(x-center)/xt)/T) * (heaviside((x-center)+xt)-heaviside((x-center)-xt)) + d','independent','x','dependent','y');
opts = fitoptions(f1);
opts.Display = 'Off';

% format: [A B T center d x0 xt]
opts.Lower = [0 0 0 0 -5*abs(offsetGuess) xtGuess/10 xtGuess/3];
opts.StartPoint = [max(projection) BGuess TGuess centerGuess offsetGuess xtGuess/4 xtGuess];
opts.Upper = [2*max(projection) 2*max(projection) inf r 5*abs(offsetGuess) xtGuess 3*xtGuess];

% % format: [a a5 b3 c c5 d]
% opts.Lower = [0  0 Csf1-r/4 Wthermal/3 Wsf1/3 -5*abs(Othermal)];
% opts.StartPoint = [Athermal Asf1 Csf1 Wthermal Wsf1 Othermal];
% opts.Upper = [2*max(projection) 2*max(projection) Csf1+r/4 3*Wthermal 1.5*Wsf1 5*abs(Othermal)];

% Fit model to data.
[peakfit1,~] = fit(xData1,yData1,f1,opts);
fitcoeff1 = coeffvalues(peakfit1);

%figure(12), plot(projection); hold on; plot(peakfit1); hold off

peaks1 = feval(peakfit1,index); % the fitted thermal plus sf fraction
resid1 = (projection-peaks1); % fitting residual

AFit = fitcoeff1(1);
BFit = fitcoeff1(2);
TFit = fitcoeff1(3);
centerFit = fitcoeff1(4);
offsetFit = fitcoeff1(5);
x0Fit = fitcoeff1(6);
xtFit = fitcoeff1(7);

thermal = BFit * exp(cos(pi*(index-centerFit)/xtFit)/TFit).*(heaviside((index-centerFit)+xtFit)-heaviside((index-centerFit)-xtFit));
coherentPeak = AFit * (1-((index-centerFit)/x0Fit).^2).^2.* (heaviside((index-centerFit)+x0Fit)-heaviside((index-centerFit)-x0Fit));

baseValue = min(thermal+offsetFit);
%condensateFraction = (Ca*Cc + (S1a + S2a)*Sc)/((Ca*Cc + (S1a + S2a)*Sc)+Ta*Tc+offset*r);
condensateFraction = sum(coherentPeak)/sum(coherentPeak+thermal+offsetFit-baseValue);

fitresult = [condensateFraction, x0Fit, xtFit];

if eachplot == 1
    figure;
    subplot(1,2,1);
    plot(index,projection,index,thermal+offsetFit,index,coherentPeak,index,resid1-5)
    subplot(1,2,2);
    imagesc(F);
    axis tight
    title(['Coherent fraction: ',num2str(condensateFraction)]);
end

end
