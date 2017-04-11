function fitresult = BEC5Fit(a, eachplot)
%{
fitresult details:
1: Condensate fraction
2: N Total
3: BEC X Width
4: BEC Y Width
5: Thermal X Width scaling relative to BEC width
6: Thermal Y Width scaling relative to BEC width
7: X Peak Height
8: Y Peak Height
%}

    OD=-log(a);
    %
    roi=[0, 0, size(OD,1), size(OD,2)];
    xRegion = round(roi(1))+(1:round(roi(4)));
    yRegion = round(roi(2))+(1:round(roi(3)));
    X = reshape(repmat(xRegion',1,numel(yRegion))',1,[]);    
    Y = repmat(yRegion,1,numel(xRegion));
    Z = reshape(OD,1,[]);
    
    
%    OD=OD(yRegion,xRegion);
    xProjection = sum(OD,1);
    yProjection = sum(OD,2)';
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
%     xMoment2 = sqrt(sum((xProjection.*(xRegion - xpeak)).^2)/size(xRegion,2));
%     yMoment2 = sqrt(sum((yProjection.*(yRegion - ypeak)).^2)/size(yRegion,2));
    xMoment4 = sqrt(sqrt(sum((xProjection.*(xRegion - xpeak)).^4)/size(xRegion,2)));
    yMoment4 = sqrt(sqrt(sum((yProjection.*(yRegion - ypeak)).^4)/size(yRegion,2)));
    
   if xFWHM>xMoment4/2;
       cfStart = 0.1;
   else cfStart = 0.9;
   end
    
    % Sum total number
    nTotal=abs(sum(xProjection));

    % Clean up data
    [xData, yData, zData] = prepareSurfaceData( X, Y, Z ); 
%     size(xData)
%     size(yData)
%     size(zData)
    
    % Set up fittype and options.
    %gbec is http://mathworld.wolfram.com/Polylogarithm.html
    ft = fittype( 'nTotal*(1-cf)/(2*pi*(rx*sx)*(ry*sy)*1.2021)*gbec(2,exp(-((x-x0)^2/(rx*sx)^2+(y-y0)^2/(ry*sy)^2)/2),5)+nTotal*cf*5/(2*pi*rx*ry)*max((1-(x-x0)^2/rx^2-(y-y0)^2/ry^2),0)^(3/2)', 'independent', {'x', 'y'}, 'dependent', 'z', 'coefficients', {'cf','nTotal','rx','ry','sx','sy','x0','y0',});
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    
    opts.Lower =        [0          0.7*nTotal  xFWHM/10    yFWHM/10   1.2      1.2      xpeak-30     ypeak-30];
    opts.StartPoint =   [cfStart    nTotal      xFWHM       yFWHM      1.5      1.2   	 xpeak        ypeak];
    opts.Upper =        [1          1.3*nTotal  5*xFWHM     5*yFWHM    inf      inf      xpeak+30     ypeak+30];
    
    opts.MaxFunEvals = 800;
    opts.MaxIter = 600;
    opts.TolX = 1e-6;
    opts.TolFun = 1e-6;
    opts.DiffMinChange = 1e-8;
    opts.DiffMaxChange = 1e-2;
    
    opts.Display = 'Off';
    
    % Fit model to data.
    [result, ~] = fit( [xData, yData], zData, ft, opts );
    rawfitresult=coeffvalues(result);
    fitresult = [rawfitresult(1), rawfitresult(2), rawfitresult(3), rawfitresult(4), rawfitresult(3)*rawfitresult(5),rawfitresult(4)*rawfitresult(6), rawfitresult(7), rawfitresult(8)];
    
    %Image plotting the fitting results versus the raw image for debugging purposes
    if eachplot
        residual = result(xData,yData)-zData;
        xFitP = sum(reshape(result(xData,yData),size(OD)),1);
        yFitP = sum(reshape(result(xData,yData),size(OD)),2);
        xThermP = sum(reshape(rawfitresult(2)*(1-rawfitresult(1))/(2*pi*rawfitresult(3)*rawfitresult(5)*rawfitresult(4)*rawfitresult(6))/1.202*gbec(2,exp(-((X-rawfitresult(7)).^2./(rawfitresult(3)*rawfitresult(5))^2+(Y-rawfitresult(8)).^2./(rawfitresult(4)*rawfitresult(6))^2)/2),3),size(OD)),1);
        yThermP = sum(reshape(rawfitresult(2)*(1-rawfitresult(1))/(2*pi*rawfitresult(3)*rawfitresult(5)*rawfitresult(4)*rawfitresult(6))/1.202*gbec(2,exp(-((X-rawfitresult(7)).^2./(rawfitresult(3)*rawfitresult(5))^2+(Y-rawfitresult(8)).^2./(rawfitresult(4)*rawfitresult(6))^2)/2),3),size(OD)),2);
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


% % Fake data
% for i=1:100;
%     for j=1:150;
%         a(i,j) = exp(-(0.80*(max(1-(((i-45)/30)^2)-(((j-65)/15)^2),0)^(3/2))+(1-0.80-0.05)*exp((-1/2)*(((i-45)^2)/300))*exp((-1/2)*(((j-65)^2)/300))));
%     end
% end
% clear i j;
% imagesc(a);

 