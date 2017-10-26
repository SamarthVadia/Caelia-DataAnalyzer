% Save entered imageID as tiff file at desired location

conn = database('becivdatabase', 'root', 'w0lfg4ng', 'Server', 'thedarkknight', 'Vendor', 'MySQL');

imageID = 154712;     %Enter imageID here

sqlquery2=['SELECT data FROM images WHERE imageID = '  int2str(imageID)];
curs2=exec(conn, sqlquery2);
curs2=fetch(curs2);
bdata=curs2.Data;
close(curs2);
blobdata=typecast(cell2mat(bdata),'int16');
s=[1024,1024,3];            %Size of the data
a=Blob2Matlab(blobdata,s);
imgmode=1;      %Image mode (1-Normal, 2-Kinetics)
framenum=1;     % Frame Number (1-Final, 2-PWA, 3-PWOA, 4-DF)
b=data_evaluation(a,imgmode,framenum);
data=cast(b,'single');

%[filename, pathname] = uiputfile('.mat','Save as');

%imwrite(data,'a.png','png');
%c=imread('a.png');
%imwrite(c,[pathname '\' filename],'tiff');