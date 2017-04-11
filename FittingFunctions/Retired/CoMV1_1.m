function CoM = CoMV1_1(a, eachplot)

OD = -log(a);

[r,s] = size(OD);



indexx = (1:s)';
indexy = (1:r)'; % column vectors enumerating the pixels of the lineout

COMx = zeros(1,r); % Preallocating space for xCoM
COMy = zeros(1,s); % for yCoM

for n=1:r;
    COMx(n) = ((OD(n,:)*indexx)/sum(OD(n,:),2));
end
clear n;
for n=1:s;
    COMy(n) = ((OD(:,n)'*indexy)/sum(OD(:,n),1));
end

centerx = (sum(COMx(1,(round(r/2)-5):(round(r/2)+5)))/size(COMx(1,(round(r/2)-5):(round(r/2)+5)),2));
centery = (sum(COMy(1,(round(s/2)-5):(round(s/2)+5)))/size(COMy(1,(round(s/2)-5):(round(s/2)+5)),2));
CoM = [centerx, centery];

end




