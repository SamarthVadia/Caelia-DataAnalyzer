function CoM = CoMV2_0(a, eachplot)

OD = -log(a);

[u,v] = meshgrid(1:size(OD,2),1:size(OD,1));

[~,xcom] = max(OD.*u);
[~,ycom] = max(OD.*v);

centerx = mode(xcom);
centery = mode(ycom);

CoM = [centerx,centery];

end




