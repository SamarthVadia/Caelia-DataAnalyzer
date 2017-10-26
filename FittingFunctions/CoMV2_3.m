function CoM = CoMV2_3(a, eachplot)

OD = -log(a);

[u,v] = meshgrid(1:size(OD,2),1:size(OD,1));

[~,xcom1] = max(OD.*u);
[~,xcom2] = max(OD.*v);
[~,ycom1] = max(OD.*u,[],2);
[~,ycom2] = max(OD.*v,[],2);

centerx1 = mode(xcom1);
centery1 = mode(ycom1);
centerx2 = mode(xcom2);
centery2 = mode(ycom2);

centerx = (centerx1 + centerx2)/2;
centery = (centery1 + centery2)/2;

CoM = [centerx, centery];

end




