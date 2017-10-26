function fitresult = BEC5Fit(a, eachplot)
%{
fitresult details:
1: Condensate fraction
2: N Total
3: Thermal X Width
4: Thermal Y Width
5: Center in X
6: Center in Y
7: Condensate peak X width
8: Condensate peak Y width
9: Condensate fraction in side peaks
10: Side peaks X position
11: Side peaks Y position
%}

    OD=-log(a);
    roi=[0, 0, size(OD,1), size(OD,2)];
    xRegion = round(roi(1))+(1:round(roi(4)));
    yRegion = round(roi(2))+(1:round(roi(3)));
    
    therm = imgaussfilt(OD,round(roi(3)/10))./1.5;
    subtherm = imgaussfilt(OD - therm,1);
    
    X = reshape(repmat(xRegion',1,numel(yRegion))',1,[]);    
    Y = repmat(yRegion,1,numel(xRegion));
    Z = reshape(OD,1,[]);
    
    
    xProjection = sum(OD,1);
    yProjection = sum(OD,2)';
    
    xProjTherm = sum(therm,1);
    yProjTherm = sum(therm,2)';
    xProjCond = sum(subtherm,1);
    yProjCond = sum(subtherm,2)';
    
    [m,xpeak]=max(xProjection);
    xpeak=xRegion(xpeak);
    ind=find(xProjection<=m/2);
    xFWHM=max(diff(ind));
    if isempty(xFWHM)
        xFWHM=xpeak;
    end
    [m,ypeak]=max(yProjection);
    ypeak=yRegion(ypeak);
    ind=find(yProjection<=m/2);
    yFWHM=max(diff(ind));
    if isempty(yFWHM)
        yFWHM=ypeak;
    end
    
    
    [m,xThPeak]=max(xProjTherm);
    xThPeak=xRegion(xThPeak);
    ind=find(xProjTherm<=m/2);
    xThFWHM=max(diff(ind));
    if isempty(xThFWHM)
        xThFWHM=xThPeak;
    end
    [m,yThPeak]=max(yProjTherm);
    yThPeak=yRegion(yThPeak);
    ind=find(yProjTherm<=m/2);
    yThFWHM=max(diff(ind));
    if isempty(yThFWHM)
        yThFWHM=yThPeak;
    end
    
    
     [m,xCondPeak]=max(xProjCond);
    xCondPeak=xRegion(xCondPeak);
    ind=find(xProjCond<=m/2);
    xCondFWHM=max(diff(ind));
    if isempty(xCondFWHM)
        xCondFWHM=xCondPeak;
    end
    [m,yCondPeak]=max(yProjCond);
    yCondPeak=yRegion(yCondPeak);
    ind=find(yProjCond<=m/2);
    yCondFWHM=max(diff(ind));
    if isempty(yCondFWHM)
        yCondFWHM=yCondPeak;
    end
    
    
    % Find position of SF peak maxima in y
    ysubtherm = max(sum(subtherm,2),0);
    ysubdiff = diff(ysubtherm,1);
    ymaxima = find([0; ysubdiff].*[ysubdiff; 0]<0);
    if size(ymaxima,1)<3
        ymaxima = [ypeak-2*yFWHM;ypeak;ypeak+2*yFWHM];
    end
    
    [ymaxvalues,ymaxindicies] = sort(ysubtherm(ymaxima),'descend');
    y1guess = (abs(ymaxima(ymaxindicies(1))-ymaxima(ymaxindicies(2)))+abs(ymaxima(ymaxindicies(1))-ymaxima(ymaxindicies(3))))/2;
   
     % Find position of SF peak maxima in x
    xsubtherm = max(sum(subtherm,1),0);
    xsubdiff = diff(xsubtherm,1)';
    xmaxima = find([0; xsubdiff].*[xsubdiff; 0]<0);
    if size(xmaxima,1)<3
        xmaxima = [xpeak-2*xFWHM;xpeak;xpeak+2*xFWHM];
    end
    
    [xmaxvalues,xmaxindicies] = sort(xsubtherm(xmaxima),'descend');
    x1guess = (abs(xmaxima(xmaxindicies(1))-xmaxima(xmaxindicies(2)))+abs(xmaxima(xmaxindicies(1))-xmaxima(xmaxindicies(3))))/2;

    
    % Sum total number
    nTotal=abs(sum(xProjection));

    % Clean up data
    [xData, yData, zData] = prepareSurfaceData( X, Y, Z ); 
%     size(xData)
%     size(yData)
%     size(zData)
    
    % Set up fittype and options.
    % gbec is http://mathworld.wolfram.com/Polylogarithm.html
    ft = fittype( '(1-cf)*nTotal/(2*pi*sx*sy)/1.202*gbec(2,exp(-((x-x0)^2/sx^2+(y-y0)^2/sy^2)/2),3)+nTotal*(cf-cfsx-cfsy)*5/(2*pi*rx*ry)*max((1-(x-x0)^2/rx^2-(y-y0)^2/ry^2),0)^(3/2)+nTotal*cfsy*5/(4*pi*rx*ry)*max((1-(x-(x0-x1))^2/rx^2-(y-(y0-y1))^2/ry^2),0)^(3/2)+nTotal*cfsy*5/(4*pi*rx*ry)*max((1-(x-(x0+x1))^2/rx^2-(y-(y0+y1))^2/ry^2),0)^(3/2)+nTotal*cfsx*5/(4*pi*tx*ty)*max((1-(x-(x0-x2))^2/tx^2-(y-(y0-y2))^2/ty^2),0)^(3/2)+nTotal*cfsx*5/(4*pi*tx*ty)*max((1-(x-(x0+x2))^2/tx^2-(y-(y0+y2))^2/ty^2),0)^(3/2)', 'independent', {'x', 'y'}, 'dependent', 'z', 'coefficients', {'cf','nTotal','sx','sy','x0','y0','rx','ry','cfsy','x1','y1','tx','ty','cfsx','x2','y2'});
    
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );

    opts.Lower =        [0     0.7*nTotal   xThFWHM/10    yThFWHM/10    xpeak-30    ypeak-30   xCondFWHM/10     yCondFWHM/10    0     0              0              xCondFWHM/10    yCondFWHM/10    0       0           0];
    opts.StartPoint =   [0.8    nTotal      xThFWHM       yThFWHM       xpeak       ypeak      xCondFWHM        yCondFWHM       0.1   1              y1guess        xCondFWHM       yCondFWHM       0.1     x1guess     1];
    opts.Upper =        [1      1.3*nTotal  Inf           Inf           xpeak+30    ypeak+30   Inf              Inf             0.5   xThFWHM        y1guess*10     inf             Inf             0.5     x1guess*10  yThFWHM];
    
    opts.MaxFunEvals = 1000;
    opts.MaxIter = 600;
    opts.TolX = 1e-6;
    opts.TolFun = 1e-6;
    opts.DiffMinChange = 1e-8;
    opts.DiffMaxChange = 1e-2;
    
    opts.Display = 'Off';
    
    % Fit model to data.
    [result, ~] = fit( [xData, yData], zData, ft, opts );
    fitresult=coeffvalues(result);
    
    %Image plotting the fitting results versus the raw image for debugging purposes
    if eachplot
%         figure;
%         imagesc(reshape(result(xData,yData),[size(OD,1),size(OD,2)]));
%         
        residual = result(xData,yData)-zData;
        xFitP = sum(reshape(result(xData,yData),size(OD)),1);
        yFitP = sum(reshape(result(xData,yData),size(OD)),2);
        xThermP = sum(reshape(fitresult(2)*(1-fitresult(1))/(2*pi*fitresult(3)*fitresult(4))/1.202*gbec(2,exp(-((X-fitresult(5)).^2./fitresult(3)^2+(Y-fitresult(6)).^2./fitresult(4)^2)/2),3),size(OD)),1);
        yThermP = sum(reshape(fitresult(2)*(1-fitresult(1))/(2*pi*fitresult(3)*fitresult(4))/1.202*gbec(2,exp(-((X-fitresult(5)).^2./fitresult(3)^2+(Y-fitresult(6)).^2./fitresult(4)^2)/2),3),size(OD)),2);
        squareres = reshape(residual,[size(OD,1),size(OD,2)]);
        figure;
        subplot(4,4,[1 2 3 5 6 7 9 10 11]);
        imagesc(squareres);
        axis off;
        subplot(4,4,[13 14 15]);
        xfig = plot(xRegion,xProjection,'.',xRegion,xFitP,'-',xRegion,xThermP,'-');
        title('x Projection');
        axis tight;
        subplot(4,4,[4 8 12]);
        yfig = plot(yRegion,yProjection,'.',yRegion,yFitP,'-',yRegion,yThermP,'-');
        rotate(yfig,[0,0,1],-90);
        title('y Projection');
        axis tight;
    end
        
end

% Fake peaks:
% for i=1:100;
% for j=1:150;
% a(i,j) = (1-0.55*exp((-1/2)*(((i-45)^2)/10))*exp((-1/2)*(((j-55)^2)/10))-0.2*exp((-1/2)*(((i-65)^2)/10))*exp((-1/2)*(((j-55)^2)/10))-0.2*exp((-1/2)*(((i-25)^2)/10))*exp((-1/2)*(((j-55)^2)/10))-0.3*exp((-1/2)*(((i-45)^2)/120))*exp((-1/2)*(((j-55)^2)/30)));
% end
% end
% BEC fitting function:
% gbec is http://mathworld.wolfram.com/Polylogarithm.html
% ft = fittype( 'nTotal*(1-cf)/(2*pi*sx*sy)/1.202*gbec(2,exp(-((x-x0)^2/sx^2+(y-y0)^2/sy^2)/2),3)+nTotal*cf*5/(2*pi*rx*ry)*max((1-(x-x0)^2/rx^2-(y-y0)^2/ry^2),0)^(3/2)', 'independent', {'x', 'y'}, 'dependent', 'z' );
