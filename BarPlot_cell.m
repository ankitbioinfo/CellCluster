
clear all 
close all 

samplenumber=[5,5,4,4,3];
G1=[827,1274,542,826,224];
G3=[142,208,90,179,45];
G9=[16, 19, 15,22,6];
N1=[5, 8, 7, 5, 2];
M1=[0, 0, 2, 1, 0];


G2=[68,94,51,58,21];
G5=[2, 2, 0, 2, 0];
G13=[7, 4, 4, 8, 3];
G12=[5, 3, 5, 7, 0];
G19=[2, 3, 0, 1, 0];
G10=[12, 18, 8, 9, 5];
G4=[15, 31, 7, 20, 5];

G6=[31, 42, 21, 17, 7];
G7=[7, 3, 0, 3, 1];
G8=[2, 3, 0, 0, 1];
G14=[1, 0, 1, 0, 3];
G17=[0,1,0,0,0];
G24=[0,0,0,1,0];
G11=[0,0,0,0,1];

N11=[4, 10, 7, 6, 1];
N2=[1, 3, 1, 1, 0];
N3=[2, 4, 0, 2, 0];

M4=[2, 2, 0, 1, 0];
M8=[1, 0, 1, 0, 0];
N7=[0,1,0,0,0];
N12=[0,0,1,0,0];
N4=[1, 1, 0, 0, 1];
N9=[0, 1, 0, 1, 0];
N5=[];
N6=[];

M9=[0,0,1,0,0];
M2=[2, 0, 0, 0, 0];
M20=[0, 0, 1, 1, 0];
M7=[0,1,0,0,0];

P8=[1, 0, 0, 0,0];
G81=[0, 1, 0, 1, 1];


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


mid1=round(low+0.25*(up-low));
mid2=round(low+0.50*(up-low));
mid3=round(low+0.75*(up-low));

if up>10
    ylim([low,up])
    set(gca,'ytick',[low,mid1,mid2,mid3,up])
end
set(gca,'xticklabel',{'DTwt','DTmt','PTwt','PTmt','DUwt'},'fontsize',6)

print(gcf, '-dpng', ['motif','_']);

