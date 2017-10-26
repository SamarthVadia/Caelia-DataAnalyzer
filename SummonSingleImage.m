    imageID = 161462;

    conn = database('becivdatabase', 'root', 'w0lfg4ng', 'Server', 'thedarkknight', 'Vendor', 'MySQL');

    sqlquery2=['SELECT data FROM images WHERE imageID = '  int2str(imageID)];
    curs2=exec(conn, sqlquery2);
    curs2=fetch(curs2);
    bdata=curs2.Data;
    close(curs2);
    
    sqlquery=['SELECT name FROM images WHERE imageID =', num2str(imageID)];
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
    %density =  imrotate(double(b),49);
    density =  imrotate(b,49);
density = density(565:966,685:775);
a=-log(density);
F = a;    
%%
r = size(F(:,1),1); % the y-image size
index = (1:r)'; % a column vector enumerating the pixels of the image

% Taking a sum to use as the main fitting data
projection = sum(F,2);