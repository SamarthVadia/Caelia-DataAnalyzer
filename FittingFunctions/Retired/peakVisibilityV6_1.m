function condensateFraction = peakVisibilityV6_1(u, eachplot)
% peakVisibilityV6_1 outputs the condensate fraction of a time-of-flight
% image of a superfluid with some thermal component. It assumes that a, the 
% input, is a two-dimensional matrix with  values proportional to the 
% optical density, cropped to a region of  interest that is just outside 
% the superfluid peaks.
% This first projects the input onto the y-axis and then fits four
% Gaussians to the resulting 1D plot. The Gaussians are associated with the
% thermal fraction, main coherent peak, and two sattelite peaks. Since we
% have already integrated twice, the number of atoms in a given peak with
% height A and width W is proportional to A*W.
% The visibility is defined as the fitted number of atoms under the
% coherent peaks divided by the total fitted number of atoms.

% This version first fits to just a thermal background, then a thermal peak
% with the central coherent peak, and finally the whole set of four peaks.
% It uses the results from each fit as the starting values for the next
% one.
% Note that it also uses the size for some upper and lower limits for the
% fit.
% We have removed the cuts that force the condensate fraction to zero in known
% cases where the fitter tends to fail.

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
filter = 50; % filter order, the higher the number the sharper the cutoff is, 50 seems acceptable I wouldn't play with this number
lowpass = 2*real(ifft(fourier.*(1-exp(-((index-r/2)./(r/2)).^filter)))); % low pass filtering the lineout
[lowpassvalue,lowpasscenter] = max(lowpass); % finding the amplitude and center for inital guesses for the fitter

% Now thermal peak fitting for the first guess
[xData0, yData0] = prepareCurveData(index,lowpass);

% Set up fittype and options for the thermal fraction
f0 = fittype('(a*exp(-(1/2)*((x-b)/c)^2))+d','independent','x','dependent','y');
opts = fitoptions(f0);
opts.Display = 'Off';

% fitting paramaters format, alphabetical: [a b c d]
opts.Lower = [0 0 0 -100];
opts.StartPoint = [lowpassvalue lowpasscenter r/4 0];
opts.Upper = [2*max(projection) r 10*r 100];

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

peaks1 = feval(peakfit1,index); % the fitted thermal plus sf fraction
resid1 = (projection-peaks1); % fitting residual

Athermal = fitcoeff1(1);
Asf1 = fitcoeff1(2);
Csf1 = fitcoeff1(3);
Wthermal = fitcoeff1(4);
Wsf1 = fitcoeff1(5);
Othermal = fitcoeff1(6);

% Run a peak finder on the residuals for initial values for the sattelite
% peaks.

% using lowpass filtered data to take derivatives and find the peak positions
smoothed = 2*real(ifft(fft(resid1).*(1-exp(-((index-r/2)./(6*r/16)).^filter))));
smoothdiff = diff(smoothed);

peakguess = find(smoothdiff(1:r-2).*smoothdiff(2:r-1) < 0); % finding all zero crossings
[~, peakindex] = sort(resid1(peakguess),'descend'); % sorted list of values and indicies of peaks in the peakguesses list

if length(peakindex) >= 2
    Csfdiff = abs((peakguess(peakindex(1))-peakguess(peakindex(2)))/2);
    Asf2 = projection(min([round(Csf1+Csfdiff),r-1]));
    Asf3 = projection(max([round(Csf1-Csfdiff) 1]));
else
    Csfdiff = r/4;
    Asf2 = 0;
    Asf3 = 0;
end
Wsf2 = Wsf1/2;


% Now fit all peaks
[xData2, yData2] = prepareCurveData(index,projection);

% Set up fittype and options for the sf peaks
f2 = fittype('(a*exp(-(1/2)*((x-b3)/c)^2))+d+(a1*exp(-(1/2)*((x-(b3+b1))/c1)^2))+(a2*exp(-(1/2)*((x-(b3-b1))/c1)^2))+(a5*exp(-(1/2)*((x-b3)/c5)^2))','independent','x','dependent','y');
opts = fitoptions(f2);
opts.Display = 'Off';

% format: [a a1 a2 a5 b1 b3 c c1 c5 d]
opts.Lower = [0 0 0 0 r/6 Csf1-r/4 Wthermal/3 Wsf2/3 Wsf1/3 -5*abs(Othermal)];
opts.StartPoint = [Athermal Asf2 Asf3 Asf1 Csfdiff Csf1 Wthermal Wsf2 Wsf1 Othermal];
opts.Upper = [2*max(projection) 2*max(projection) 2*max(projection) 2*max(projection) r/2 Csf1+r/4 3*Wthermal 1.5*Wsf2 1.5*Wsf1 5*abs(Othermal)];

% Fit model to data.
[peakfit2,gof] = fit(xData2,yData2,f2,opts);
fitcoeff2 = coeffvalues(peakfit2);

peaks2 = feval(peakfit2,index); % the fitted thermal plus sf fraction
resid2 = (projection-peaks2); % fitting residual

Ta = fitcoeff2(1);
S1a = fitcoeff2(2);
S2a = fitcoeff2(3);
Ca = fitcoeff2(4);
Sb = fitcoeff2(5);
Cb = fitcoeff2(6);
Tc = fitcoeff2(7);
Sc = fitcoeff2(8);
Cc = fitcoeff2(9);
offset = fitcoeff2(10);

thermal = Ta*exp(-(1/2)*((index-Cb)./Tc).^2);
SattPeaks = S1a*exp(-(1/2)*((index-(Cb+Sb))./Sc).^2)+S2a*exp(-(1/2)*((index-(Cb-Sb))./Sc).^2);
CenterPeak = Ca*exp(-(1/2)*((index-Cb)./Cc).^2);

condensateFraction = (Ca*Cc + (S1a + S2a)*Sc)/((Ca*Cc + (S1a + S2a)*Sc)+Ta*Tc);

% % cuts for 0 condensate fraction
% if (max([S1a S2a])/Ca > 3/4)
%     condensateFraction = 0;
% 
% elseif (max(resid1) > Ca)
%     condensateFraction = 0;
% 
% elseif gof.rsquare < 0.7
%     condensateFraction = 0;
% 
% elseif (Sb < 2*(Cc + Sc))
%     condensateFraction = 0;
% 
% elseif (max(resid1) > 2*max([S1a S2a]))
%     condensateFraction = 0;
% 
% elseif (round(Cb - Sb) < 1 || round(Cb + Sb) > r)
%     condensateFraction = 0;
% 
% elseif (max([resid1(round(Cb - Sb)) resid1(round(Cb + Sb))]) > max([S1a S2a]))
%     condensateFraction = 0;
% 
% elseif (Sc > 2*Cc || Sc < Cc/3)
%     condensateFraction = 0;
% end

if eachplot == 1
    figure;
    subplot(1,2,1);
    plot(index,peaks2,index,projection,index,thermal+offset,index,SattPeaks,index,CenterPeak,index,resid2-5)
    subplot(1,2,2);
    imagesc(F);
    axis tight
    title(['Coherent fraction: ',num2str(condensateFraction)]);
end

end






