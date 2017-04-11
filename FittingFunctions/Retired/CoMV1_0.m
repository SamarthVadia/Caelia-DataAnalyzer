function [centery] = CoMV1_0(a, eachplot)

F = a';
G = a;

r = size(F(:,1),1); % the y-image size
s = size(G(:,1),1); % the x-image size

indexy = (1:r)'; % a column vector enumerating the pixels of the lineout
indexx = (1:s)';

% prefourier = fft(F)


COMx = zeros(1,s); % Preallocating space for ycom
COMy = zeros(1,r); % for xcom

for n=1:r;
    COMx(n) = round(sum(F(n,:).*indexx')/sum(F(n,:),2));
end
for n=1:s;
    COMy(n) = round(sum(F(:,n).*indexy)/sum(F(:,n),1));
end

centerx = round(sum(COMx(1,(round(s/2)-3):(round(s/2)+3)))/size(COMx(1,(round(s/2)-3):(round(s/2)+3)),2));
centery = round(sum(COMy(1,(round(r/2)-3):(round(r/2)+3)))/size(COMy(1,(round(r/2)-3):(round(r/2)+3)),2));

end






