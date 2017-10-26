function [ave_vis] = peakVisV5_0(u, eachplot)

% rescaling of the figure, done automatically in Winview
% F = -log(double(a(670:770,540:640,1)));
% G = -log(double(a(670:770,540:640,1)))';
% figure;imagesc(F)

a=-log(u);
F = a';
G = a;


%finding the center and subtracting background
r = size(F(:,1),1); % the y-image size
s = size(G(:,1),1); % the x-image size

indexy = (1:r)'; % a column vector enumerating the pixels of the lineout
indexx = (1:s)';

% prefourier = fft(F)


COMx = zeros(1,s); % Preallocating space for ycom
COMy = zeros(1,r); % for xcom

for n=1:r;
    COMx(n) = round(sum(F(n,:).*indexx')/sum(F(n,:),2));
end
for n=1:s;
    COMy(n) = round(sum(F(:,n).*indexy)/sum(F(:,n),1));
end

centerx = round(sum(COMx(1,(round(s/2)-3):(round(s/2)+3)))/size(COMx(1,(round(s/2)-3):(round(s/2)+3)),2));

background = (sum(F(1,:))+sum(F(:,1))+sum(F(size(F,1),:))+sum(F(:,size(F,2))))/(2*size(F,2)+2*size(F,1)); % sum of all the edge pixels to estimate the image background offset
F = F - background; %Subtracting the average background

% Taking a lineout to use as the main fitting data
lineouty = sum(F(:,max(1,(centerx-7)):(centerx+7)),2)/15; % 1D lineout of the center 5 pixels in a column vector format for both x and y directions

% Moving to fourier space for instituting some pre-fitting and cuts
fouriery = fft(lineouty); % fourier transform

%Low pass filtering for thermal fit
filter = 50; % filter order, the higher the number the sharper the cutoff is, 12 seems acceptable I wouldn't play with this number
lowpassy = 2*real(ifft(fouriery.*(1-exp(-((indexy-r/2)./(r/2)).^filter)))); % low pass filtering the lineout

[lowpassvaluey,lowpasscentery] = max(lowpassy); % finding the amplitude and center for inital guesses for the fitter

ycount = sum(lowpassy,1); % counting the atom number in our lineout, will be used as a cut for visibility data where no atom number = 0 visibility

expectedcounty = background*r*10; % how many atoms we expect to see just from background noise, is the comparison point for the atom number cut



% Now thermal peak fitting for the first guess in the y direction
[xData0, yData0] = prepareCurveData(indexy,lowpassy);

% Set up fittype and options for the thermal fraction
f0t = fittype('(a*exp(-(1/2)*((x-b)/c)^2))+d','independent','x','dependent','y');
opts = fitoptions(f0t);
opts.Display = 'Off';
opts.Lower = [0 0 0 -100];
% fitting paramaters format, alphabetical: [a b c d]
opts.StartPoint = [lowpassvaluey lowpasscentery r/4 0];
opts.Upper = [2*max(max(F)) r 10*r 100];
% Fit model to data.
[peakfit0,~] = fit(xData0,yData0,f0t,opts);
fitcoeff0 = coeffvalues(peakfit0); % coefficients used as initial guesses for the final fitter
% curve that will serve as an initial thermal component guess
peaks0 = feval(peakfit0,indexy);

%now using a lowpass filtered lineout to take derivatives and find the peak positions
smoothedy = 2*real(ifft(fouriery.*(1-exp(-((indexy-r/2)./(6*r/16)).^filter))));
smoothdiffy = diff(smoothedy);

peakguessy = find(smoothdiffy(1:r-2).*smoothdiffy(2:r-1) < 0); % finding all zero crossings

valleyguessy = zeros(length(peakguessy),1); %preallocating the size of the valleyguess vector
% loop to filter all maxima from minima
for iter = 1:length(peakguessy);
    if smoothedy(max(1,peakguessy(iter)-1))>smoothedy(peakguessy(iter)+1)
        valleyguessy(iter) = peakguessy(iter);
        peakguessy(iter) = 1;
    elseif smoothedy(max(1,peakguessy(iter)-1))<smoothedy(peakguessy(iter)+1)
        valleyguessy(iter) = 1;
    end
end
% peakguess = peakguess(peakguess~=0); % vector of peak postions
% valleyguess = valleyguess(valleyguess~=0); % vector of valley positions
[~, peakindexy] = sort(lineouty(peakguessy),'descend'); % sorted list of values and indicies of peaks in the peakguesses list


% inital guesses for the fittter based on what we've done so far
% if fitcoeff0(3) < 15
if fitcoeff0(3) < r/7
    Athermaly = 0.05*fitcoeff0(1);
    Wthermaly = 20;
    Asf1y = lineouty(peakguessy(peakindexy(1)));
    if length(peakindexy) >= 3
        Asf2y = lineouty(peakguessy(peakindexy(2)));
        Asf3y = lineouty(peakguessy(peakindexy(3)));
    else
        Asf2y = 0;
        Asf3y = 0;
    end
else
    Athermaly = 0.7.*fitcoeff0(1);
    Wthermaly = 1.0.*fitcoeff0(3);
    Asf1y = lineouty(peakguessy(peakindexy(1))) - 0.7.*peaks0(peakguessy(peakindexy(1)));
    if length(peakindexy) >= 3
        Asf2y = lineouty(peakguessy(peakindexy(2))) - 0.7.*peaks0(peakguessy(peakindexy(2)));
        Asf3y = lineouty(peakguessy(peakindexy(3))) - 0.7.*peaks0(peakguessy(peakindexy(3)));
    else
        Asf2y = 0;
        Asf3y = 0;
    end
end
Cthermaly = 1.0.*fitcoeff0(2);
Othermaly = 1.0.*fitcoeff0(4);
% Csfdiffy = 32;
% Csfdiffy = sum(diff(sort(peakguessy(peakindexy(1:5)))))/4;
Csf1y = peakguessy(peakindexy(1));
if length(peakindexy) >= 3
    Csfdiffy = abs((peakguessy(peakindexy(2))-peakguessy(peakindexy(3)))/2);
else
    Csfdiffy = 30;
end
% Wsf2y = 6; %sum(diff(sort(peakguess(peakindex(1:5)))))/25;
% Wsf1y = 10; %sum(diff(sort(peakguess(peakindex(1:5)))))/25;
Wsf2y = r/10;
Wsf1y = r/10;

% Now full SF peak fitting
[xData1, yData1] = prepareCurveData(indexy,lineouty);

% Set up fittype and options for the sf peaks
f1t = fittype('(a*exp(-(1/2)*((x-b)/c)^2))+d+(a1*exp(-(1/2)*((x-(b3+b1))/c1)^2))+(a2*exp(-(1/2)*((x-(b3-b1))/c1)^2))+(a5*exp(-(1/2)*((x-b3)/c5)^2))','independent','x','dependent','y');
opts = fitoptions(f1t);
opts.Display = 'Off';
% format: [a a1 a2 a5 b b1 b3 c c1 c5 d]
opts.Lower = [0 0 0 0 3*Cthermaly/4 14 0.8*Csf1y Wthermaly/3 Wsf2y/3 Wsf1y/3 -1];
opts.StartPoint = [Athermaly Asf2y Asf3y Asf1y Cthermaly Csfdiffy Csf1y Wthermaly Wsf2y Wsf1y Othermaly];
opts.Upper = [200 200 200 200 5*Cthermaly/4 40 1.2*Csf1y 3*Wthermaly 1.5*Wsf2y 1.5*Wsf1y 1];
% Fit model to data.
[peakfit1,~] = fit(xData1,yData1,f1t,opts);
fitcoeff1 = coeffvalues(peakfit1);

peaksy = feval(peakfit1,indexy); % the fitted thermal plus sf fraction
residualy = (lineouty-peaksy); % fitting residual

amplitudey = fitcoeff1(1:4);
centery = fitcoeff1(5:7);
widthy = fitcoeff1(8:10);
offsety = fitcoeff1(11);

thermal2y = amplitudey(1)*exp(-(1/2)*((indexy-centery(1))./widthy(1)).^2)+offsety;
SFPeaks1y = amplitudey(2)*exp(-(1/2)*((indexy-(centery(3)+centery(2)))./widthy(2)).^2)+amplitudey(3)*exp(-(1/2)*((indexy-(centery(3)-centery(2)))./widthy(2)).^2);
SFPeaks3y = amplitudey(4)*exp(-(1/2)*((indexy-centery(3))./widthy(3)).^2);

if eachplot == 1
    % close all
    figure
    subplot(1,2,1);
    plot(indexy,peaksy,'r-',indexy,lineouty,'g-',indexy,thermal2y,'b-',indexy,SFPeaks1y,'k-',indexy,SFPeaks3y,'y-',indexy,residualy-0.1,'c-')
    subplot(1,2,2);
    % plot(indexy,SFPeaks1y);
    imagesc(F);
    axis tight
end


if amplitudey(1)>200*amplitudey(2)
    visibility = -2;
elseif amplitudey(1)>200*amplitudey(3)
    visibility = -3;
elseif amplitudey(1)>200*amplitudey(4)
    visibility = -4;
else visibility = ((amplitudey(2)/(amplitudey(2)+2*thermal2y(round(centery(3)+centery(2)))))+(amplitudey(3)/(amplitudey(3)+2*thermal2y(round(centery(3)-centery(2)))))+(amplitudey(4)/(amplitudey(4)+2*thermal2y(round(centery(3))))))/5;
end

ave_vis = visibility;

end






