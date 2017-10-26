conn = database('becivdatabase', 'root', 'w0lfg4ng', 'Server', 'thedarkknight', 'Vendor', 'MySQL');

imageID = 160901:160905;     %Enter imageID here
imageListLength = length(imageID);
imageNameList = zeros(imageListLength,1);
imageList = cell(imageListLength,1);
for i=1:imageListLength
    sqlquery2=['SELECT data FROM images WHERE imageID = '  int2str(imageID(i))];
    curs2=exec(conn, sqlquery2);
    curs2=fetch(curs2);
    bdata=curs2.Data;
    close(curs2);
    
    sqlquery=['SELECT name FROM images WHERE imageID =', num2str(imageID(i))];
    curs1=exec(conn, sqlquery);
    curs1=fetch(curs1);
    info = curs1.Data;
    close(curs1);
    
    blobdata=typecast(cell2mat(bdata),'int16');
    s=[1024,1024,3];            %Size of the data
    a=Blob2Matlab(blobdata,s);
    imgmode=1;      %Image mode (1-Normal, 2-Kinetics)
    framenum=1;     % Frame Number (1-Final, 2-PWA, 3-PWOA, 4-DF)
    b=data_evaluation(a,imgmode,framenum);
    %density =  imrotate(-log(double(b)),47);
    density =  imrotate(double(b),47);
    densityRot = density(730:815,280:1180);
    imageList{i} = densityRot;
    tmp = strsplit(info{1});
    imageNameList(i) = str2num(tmp{3});
end

first = zeros(1,imageListLength);
second = zeros(1,imageListLength);
third = zeros(1,imageListLength);
for i=1:imageListLength
    tmp = KapitzaDiracV1_0(imageList{i},0);
    first(i) = tmp(1);
    second(i) = tmp(2);
    third(i) = tmp(3);
end
figure(10);
%fitplt = axes('Units','pixels','Position',[1050,480,500,450]);
%fitres = uicontrol('Style','edit','String','Results','min', 0, 'max', 100,'Position',[1370,260,180,170]); %Results for fit function
%%
leastsquare=cell(imageListLength,1);
        ls=@(theta) 0;
        for i=1:imageListLength
            leastsquare{i}=@(theta)(besselj(0,theta*imageNameList(i))^2-first(i))^2+(besselj(1,theta*imageNameList(i))^2-second(i))^2+(besselj(2,theta*imageNameList(i))^2-third(i))^2;
            ls=@(theta) ls(theta)+leastsquare{i}(theta);
        end
        theta0=0.2;
        theta=fminsearch(ls, theta0);
        volts_5recoil=(0.0000125*1000*2*pi*5*2.02781)/(2*theta);
        resultx = zeros(1,imageListLength);
        resulty0 = zeros(1,imageListLength);
        resulty1 = zeros(1,imageListLength);
        resulty2 = zeros(1,imageListLength);
        for i=1:imageListLength
            resulty0(i)=first(i);
            resulty1(i)=second(i);
            resulty2(i)=third(i);
            resultx(i)=theta*imageNameList(i);
        end
        showdata1='';
        for i=1:length(resulty0)
            showdata1{i}=sprintf('%s , %s', num2str(imageNameList(i)), num2str(resulty0(i)));
        end
        showdata=sprintf('\n%s', showdata1{:});
        %axes(fitplt);
        plot(resultx,resulty0,'bo',resultx,resulty1,'ro',resultx,resulty2,'go');
        hold on
        fplot(@(x) besselj(0,x)^2, [min(resultx) max(resultx)],'b');
        fplot(@(x) besselj(1,x)^2, [min(resultx) max(resultx)],'r');
        fplot(@(x) besselj(2,x)^2, [min(resultx) max(resultx)],'g');
        hold off
        ['Theta: ' num2str(theta) ]
        ['Volts/5 Erec: ' num2str(volts_5recoil)] 
        ['Bessel 0 Data: '] 
        showdata