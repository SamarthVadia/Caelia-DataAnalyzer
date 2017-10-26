function CoM = CoMV2_2(a, eachplot)

OD = -log(a);

[u,v] = meshgrid(1:size(OD,2),1:size(OD,1));

[~,xcom1] = max(OD.*u);
[~,xcom2] = max(OD.*v);
[~,ycom1] = max(OD.*v,[],2);
[~,ycom2] = max(OD.*u,[],2);

centerx = (mode(xcom1)+mode(xcom2))/2;
centery = (mode(ycom1)+mode(ycom2))/2;

CoM = [centerx, centery];

end