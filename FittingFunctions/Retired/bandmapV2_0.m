function [imbalance] = bandmapV2_0(a, eachplot)

% F = double(-log(a(500:700,680:720,1)));
% % imagesc(F);
%a=-log(u);
F = a';

lmaxx = size(F(1,:),2);
indexx = 1:lmaxx;
lmaxy = size(F(:,1),1);
indexy = 1:lmaxy;

COMx = zeros(1,lmaxx);
COMy = zeros(1,lmaxy);

for n=1:lmaxy;
    COMx(n) = round(sum(F(n,:).*indexx)/sum(F(n,:),2));
end
for n=1:lmaxx;
    COMy(n) = round(sum(F(:,n).*indexy')/sum(F(:,n),1));
end

centerx = round(sum(COMx(1,(round(lmaxx/2)-3):(round(lmaxx/2)+3)))/size(COMx(1,(round(lmaxx/2)-3):(round(lmaxx/2)+3)),2));
% centery = round(sum(COMy(1,(round(lmaxy/2)-3):(round(lmaxy/2)+3)))/size(COMy(1,(round(lmaxy/2)-3):(round(lmaxy/2)+3)),2));

lineouty = sum(F(:,max(1,(centerx-3)):min(lmaxy,(centerx+3))),2)/7;
% figure;
% plot(lineouty);
diffline = diff(lineouty);
crossings = find(diffline(1:(lmaxy-2)).*diffline(2:(lmaxy-1))<0);

for it=1:length(crossings);
    if diffline(max(1,(crossings(it)-1)))<diffline(min(lmaxy,(crossings(it)+1)));
        crossings(it)=1;
    end
end

[~,peakindex] = sort(lineouty(crossings),'descend');
peakpositions = crossings(peakindex)+1;

if numel(peakpositions)<3
    imbalance = 0;
    disp('ERR:Not enough peaks to fit')
end



% Moving to fourier space for instituting some pre-fitting and cuts
fouriery = fft(lineouty); % fourier transform

%Low pass filtering for thermal fit
filter = 50; % filter order, the higher the number the sharper the cutoff is, 12 seems acceptable I wouldn't play with this number
lowpassy = 2*real(ifft(fouriery.*(1-exp(-((indexy'-lmaxy/2)./(lmaxy/2)).^filter)))); % low pass filtering the lineout

[lowpassvaluey,lowpasscentery] = max(lowpassy); % finding the amplitude and center for inital guesses for the fitter



% Now thermal peak fitting for the first guess in the y direction
[xData0, yData0] = prepareCurveData(indexy',lowpassy);

% Set up fittype and options for the thermal fraction
f0t = fittype('(a*exp(-(1/2)*((x-b)/c)^2))+d','independent','x','dependent','y');
opts = fitoptions(f0t);
opts.Display = 'Off';
opts.Lower = [0 0 0 -100];
% fitting paramaters format, alphabetical: [a b c d]
opts.StartPoint = [lowpassvaluey lowpasscentery lmaxy/4 0];
opts.Upper = [2*max(max(F)) lmaxy 10*lmaxy 100];
% Fit model to data.
[peakfit0,~] = fit(xData0,yData0,f0t,opts);
fitcoeff0 = coeffvalues(peakfit0); % coefficients used as initial guesses for the final fitter
% curve that will serve as an initial thermal component guess
peaks0 = feval(peakfit0,indexy);





% inital guesses for the fittter based on what we've done so far
Athermaly = fitcoeff0(1);
Cthermaly = fitcoeff0(2);
Wthermaly = fitcoeff0(3);
Othermaly = fitcoeff0(4);
if lineouty(peakpositions(2))>lineouty(peakpositions(3));
    Asf2y = lineouty(peakpositions(2)) - 0.8.*peaks0(peakpositions(2));
    Asf3y = lineouty(peakpositions(3)) - 0.8.*peaks0(peakpositions(3));
elseif lineouty(peakpositions(2))<lineouty(peakpositions(3));
    Asf2y = lineouty(peakpositions(3)) - 0.8.*peaks0(peakpositions(3));
    Asf3y = lineouty(peakpositions(2)) - 0.8.*peaks0(peakpositions(2));
end
Asf1y = lineouty(peakpositions(1)) - 0.8.*peaks0(peakpositions(2));
Csfdiffy = abs(sum(diff(peakpositions(1:2))));
Csf1y = peakpositions(1);
Wsf2y = 3; %sum(diff(sort(peakguess(peakindex(1:5)))))/25;
Wsf1y = 3; %sum(diff(sort(peakguess(peakindex(1:5)))))/25;

% Now full SF peak fitting
[xData1, yData1] = prepareCurveData(indexy',lineouty);

% Set up fittype and options for the sf peaks
f1t = fittype('(a*exp(-(1/2)*((x-b)/c)^2))+d+(a1*exp(-(1/2)*((x-(b3+b1))/c1)^2))+(a2*exp(-(1/2)*((x-(b3-b1))/c1)^2))+(a5*exp(-(1/2)*((x-b3)/c5)^2))','independent','x','dependent','y');
opts = fitoptions(f1t);
opts.Display = 'Off';
% format: [a a1 a2 a3 a4 a5 b b1 b3 c c1 c5 d]
opts.Lower = [0 0 0 0 3*Cthermaly/4 Csfdiffy/2 0.8*Csf1y Wthermaly/3 Wsf2y/3 Wsf1y/3 -5];
opts.StartPoint = [Athermaly Asf2y Asf3y Asf1y Cthermaly Csfdiffy Csf1y Wthermaly Wsf2y Wsf1y Othermaly];
opts.Upper = [200 200 200 200 5*Cthermaly/4 2*Csfdiffy 1.2*Csf1y 3*Wthermaly 1.5*Wsf2y 1.5*Wsf1y 5];
% Fit model to data
[peakfit1,~] = fit(xData1,yData1,f1t,opts);
fitcoeff1 = coeffvalues(peakfit1);

peaksy = feval(peakfit1,indexy); % the fitted thermal plus sf fraction
residualy = (lineouty-peaksy); % fitting residual


amplitudey = fitcoeff1(1:4);
centery = fitcoeff1(5:7);
widthy = fitcoeff1(8:10);
offsety = fitcoeff1(11);

thermal2y = amplitudey(1)*exp(-(1/2)*((indexy-centery(1))./widthy(1)).^2)+offsety;
SFPeaks1y = amplitudey(2)*exp(-(1/2)*((indexy-(centery(3)+centery(2)))./widthy(2)).^2);
SFPeaks2y = amplitudey(3)*exp(-(1/2)*((indexy-(centery(3)-centery(2)))./widthy(2)).^2);
SFPeaks3y = amplitudey(4)*exp(-(1/2)*((indexy-centery(3))./widthy(3)).^2);


initialguess = Athermaly*exp(-(1/2)*((indexy-Cthermaly)./Wthermaly).^2)+Othermaly+Asf2y*exp(-(1/2)*((indexy-(Csf1y+Csfdiffy))./Wsf2y).^2)+Asf3y*exp(-(1/2)*((indexy-(Csf1y-Csfdiffy))./Wsf2y).^2)+Asf1y*exp(-(1/2)*((indexy-Csf1y)./Wsf1y).^2);


%close all
if eachplot == 1
    figure
    plot(indexy,peaksy,indexy,lineouty,indexy,thermal2y,indexy,SFPeaks1y,indexy,SFPeaks2y,indexy,SFPeaks3y,indexy,residualy-0.1)
    axis tight
    figure;
    plot(indexy,lineouty,indexy,initialguess,indexy,peaks0)
    axis tight
end

% box=5;
% 
% line1 = sum(sum(lineouty((peakpositions(1)-box):(peakpositions(1)+box)),1));
% line2 = sum(sum(lineouty((peakpositions(2)-box):(peakpositions(2)+box)),1));
% line3 = sum(sum(lineouty((peakpositions(3)-box):(peakpositions(3)+box)),1));
% 
% imbalance = (line1-line2-line3)/(line1+line2+line3);

Ncount1 = sum(SFPeaks1y,2);
Ncount2 = sum(SFPeaks2y,2);
Ncount3 = sum(SFPeaks3y,2);
NcountTh = sum(thermal2y,2);

imbalance = (Ncount3-Ncount2-Ncount1)/(Ncount1+Ncount2+Ncount3+NcountTh);

end
