
clear all 
close all 

samplenumber=[5,5,4,4,3];
G1=[1501,2430,1078,1691,435];
G3=[276,505,206,380,90];
G9=[51, 101, 40, 64, 17];
N1=[21, 29, 11, 18, 1];
M1=[9, 6, 0, 6, 3];
T1=[3, 1, 1, 0, 0];

G2=[216, 276, 176, 171, 55];
G5=[2, 3, 3, 1, 1];
G13=[17, 26, 12, 30, 10];
G12=[6, 16, 9, 8, 2];
G19=[2, 11, 3, 7, 1];
G10=[25, 41, 26, 47, 2];
G4=[39, 70, 36, 52, 10];

G6=[51, 134, 38, 78, 20];
G7=[3, 25, 10, 8, 4];
G8=[2, 1, 6, 0, 0];
G14=[1, 5, 2, 4, 1];
G17=[2, 6, 1, 0, 0];
G24=[0, 2, 0, 2, 0];
G11=[3, 1, 0, 0, 0];

N11=[17, 15, 11, 22, 4];
N2=[2, 5, 4, 5, 0];
N3=[5, 12, 5, 4, 3];

M4=[2, 13, 4, 9, 3];
M8=[1, 3, 0, 0, 0];
N7=[2, 0, 0, 2, 0];
N12=[1, 1, 0, 3, 0];
N4=[0, 1, 2, 0, 0];
N9=[0, 4, 0, 4, 0];
N5=[1, 2, 1, 1, 2];
N6=[1, 2, 0, 1, 0];

M9=[0, 2, 0, 2, 1];
M2=[2, 1, 1, 0, 0];
M20=[2, 4, 1, 4, 0];
M7=[4, 1, 1, 0, 0];

P8=[3, 1, 0, 2, 0];
G81=[5, 5, 2, 0, 0];




data=[G81];

for i=1:size(data,1)
    data(i,:)=data(i,:)./samplenumber;
end

newdata=data(:,[1,3,2,4,5]);


h=figure;
set(gcf, 'PaperSize', [2.5 2]); %7
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 2.5 2]);

bar(newdata)
%pause 
low=round(0.9*min(newdata));
up=ceil(max(newdata));
ylim([low,up])

mid1=round(low+0.25*(up-low));
mid2=round(low+0.50*(up-low));
mid3=round(low+0.75*(up-low));

if up>10
    set(gca,'ytick',[low,mid1,mid2,mid3,up])
end
set(gca,'xticklabel',{'DTwt','DTmt','PTwt','PTmt','DUwt'},'fontsize',6)

print(gcf, '-dpng', ['motif','_']);

