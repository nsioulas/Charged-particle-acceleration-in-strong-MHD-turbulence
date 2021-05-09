lsize = 16; % Label fontsize
nsize = 16; % Axis fontsize
%% ENERGY SPLIT
start=1;    rd1start=1;      end2=50; %number of bins
ending=10;     rd1end=length(tM);
split2=logspace(0,log10(max(W_esc_eV(:))),end2);



W_esc_half=zeros(1,end2);
p1=zeros(end2,length(tM));

Rsm=zeros(end2,length(tM));
Dw1sm=zeros(end2,length(tM));
Dw2sm=zeros(end2,length(tM));


if prosexe2==0
    prosexe2=1;
dw1=zeros(nP,length(tM));
dw2=zeros(nP,length(tM));
for i=1:nP
for j=1:length(tM)
if Wkinet_eV(i,j)>0
dw1(i,j)=(Wkinet_eV(i,j)-W0_eV(i));
dw2(i,j)=(Wkinet_eV(i,j)-W0_eV(i)).^2;
else
dw1(i,j)=0;
dw2(i,j)=0;
end
end
end


%end
rmsd_=rmsd;
dw1_=dw1;
dw2_=dw2;

for jN=1
rmsd=rmsd_;

dw1=dw1_;
dw2=dw2_;
for i=1:nP
if W_esc_eV(i)>=split2(jN)
W_esc_half(jN)=(1/2)*(split2(jN));
rmsd(i,:)=0;
dw1(i,:)=0;
dw2(i,:)=0;
end
end
for j=1:length(tM)
p1(jN,j)=nnz(rmsd(:,j));
Rsm(jN,j)=sum(rmsd(:,j))/p1(jN,j);
Dw1sm(jN,j)=sum(dw1(:,j))/p1(jN,j);
Dw2sm(jN,j)=sum(dw2(:,j))/p1(jN,j);
end
dw1=dw1_;
dw2=dw2_;
end
for jN=2:end2
disp(['percentage: ' num2str(jN/end2)])
rmsd=rmsd_;
dw1=dw1_;
dw2=dw2_;
W_esc_half(jN)=(1/2)*(split2(jN)-split2(jN-1))+split2(jN-1);
for i=1:nP
if W_esc_eV(i)<=split2(jN-1) || W_esc_eV(i)>=split2(jN)
rmsd(i,:)=0;
dw1(i,:)=0;
dw2(i,:)=0;
end
end
for j=1:length(tM)
p1(jN,j)=nnz(rmsd(:,j));
Rsm(jN,j)=sum(rmsd(:,j))/p1(jN,j);
Dw1sm(jN,j)=sum(dw1(:,j))/p1(jN,j);
Dw2sm(jN,j)=sum(dw2(:,j))/p1(jN,j);
end
rmsd=rmsd_;
dw1=dw1_;
dw2=dw2_;
%dv2=dv2_;
%vkinetsq=vkinetsq_;
end
figure(2)
for u=1:end2
loglog( tM(start:2*ending),Dw1sm(u,start:2*ending))
hold on
end
hold off
figure(3)
for u=1:end2
loglog( tM(start:2*ending),Dw2sm(u,start:2*ending))
hold on
end
hold off
%}
a1=zeros(1,end2);
a2=zeros(1,end2);
a3=zeros(1,end2);
a4=zeros(1,end2);
new=ones(end2,length(tM));
tM1=new.*tM;       %creates a new matrix so tM1 has the same dimensions as Rsm

end
%% For rmsd
for jN=1:end2
figure(8)
disp(['percentage: ' num2str(jN/end2)])
fitResults1 = polyfit(log10(tM1(jN,start:ending)),log10(Rsm(jN,start:ending)),1); %prosexe edw to /2
pol1=polyval(fitResults1,log10(tM1(jN,start:ending)));
a1(jN)=fitResults1(1);
a11=a1>-10; %non nan values of a1 and W_esc
end
scatter((W_esc_half(a11)),a1(a11),'MarkerEdgeColor',[0 .5 .5],...
'MarkerFaceColor',[0 .7 .7],...
'LineWidth',1.5)
set(gca, 'YScale', 'linear')
set(gca, 'XScale', 'log')
hold on
ans=nnz(log10(W_esc_half)<4);
a12=a1(a11);
av1=mean(a12(1:ans));
average=av1*ones(1,ans+1);
new=(W_esc_half(1:ans+1));
semilogx(new,average,'cyan','LineWidth',2)
anew=(log10(W_esc_half)>4);
%semilogx(W_esc_half(anew),a1(anew),'blue','LineWidth',2)
xlabel('Escape energy [eV]','FontSize',lsize)
ylabel('a_{r}','FontSize',lsize)
xlim([1,10^8])
ylim([min(a1)-0.1 max(a1)+0.3])
XTick = [ 10^(-1)   10^(1)  10^(3)  10^(5)  10^(7) 10^(9)];
set(gca,'xtick',XTick)
%YTick = [10^(-12) 10^(-10) 10^(-8) 10^(-6) 10^(-4) 10^(-2) 10^(0)];
%set(gca,'ytick',YTick)
set(gca, 'Fontsize', nsize)
set(gcf,'paperpositionmode','auto');
set(gcf,'windowstyle','normal');
set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))
set(gca,'fontweight','normal')
box on
opts.Colors     = get(groot,'defaultAxesColorOrder');
opts.saveFolder = 'img/';
opts.width      = 12;
opts.height     = 10;
opts.fontType   = 'Times';
% creating the zoom-in inset
% creating the zoom-in inset
ax=axes;
set(ax,'units','normalized','position',[0.18,0.7,0.25,0.25])
box(ax,'on')
for u=1:end2
loglog( tM(start:ending),Rsm(u,start:ending),'LineWidth',0.5)
hold on
end
set(gca,'fontweight','normal','FontSize', 9)
set(ax,'xlim',[tM(start),tM(ending)],'ylim',[min(Rsm(:)),10*max(Rsm(:))])
%XTick = [10^(-1) 10^(0) ];%  10^(8)];
%set(gca,'xtick',XTick)
YTick = [ 10^0 10^2 10^(4) 10^(6) 10^(8) 10^(10)  10^(12) 10^(14) 10^(16) 10^(18) 10^(20)];
set(gca,'ytick',YTick)
str = sprintf('<Är^{2}>-t');
title(str ,'FontSize', 10);
set(get(gca,'title'),'Position',[2*tM(start) 0.5*max(Rsm(:)) 0.8])
ax = gca;
ax.YAxisLocation = 'right';
hold off
if imodel == 1

    temp=['png_files_upper\f3b.png']; saveas(gca,temp);
    temp1=['eps_files_upper\f3b.eps']; saveas(gca,temp1); 
    temp2=['fig_files_upper\f3b.fig']; saveas(gca,temp2); 
    
 elseif imodel == 12
          temp=['new_png\f3b.png']; saveas(gca,temp);
    temp1=['new_eps\f3b.eps']; saveas(gca,temp1); 
    temp2=['new_fig\f3b.fig']; saveas(gca,temp2); 
 elseif imodel == 2 
     
      temp=['new_png_rec\f3b P=',num2str(P_reconnect),'.png']; saveas(gca,temp);
    temp1=['new_eps_rec\f3b P=',num2str(P_reconnect),'.eps']; saveas(gca,temp1); 
    temp2=['new_fig_rec\f3b P=',num2str(P_reconnect),'.fig']; saveas(gca,temp2); 
 end

%% for w_esc-dw
for jN=1:end2
figure(10)
disp(['percentage: ' num2str(jN/end2)])
fitResults3 = polyfit(log10(tM1(jN,start:ending)),log10(Dw1sm(jN,start:ending)),1); %prosexe edw to /2
pol3=polyval(fitResults3,log10(tM1(jN,start:ending)));
a4(jN)=fitResults3(1);
a44=a4>-10; %non nan values of a1 and W_esc
end
scatter((W_esc_half(a44)),a4(a44),'MarkerEdgeColor',[0 .5 .5],...
'MarkerFaceColor',[0 .7 .7],...
'LineWidth',1.5)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'linear')
hold on
%loglog((W_esc_half(a44)),a4(a44),'blue','LineWidth',0.1)
ans2=nnz(log10(W_esc_half)<4);
a42=a4(a44);
av2=mean(a42(1:ans2));
average2=av2*ones(1,ans2);
new2=(W_esc_half(1:ans2));
semilogx(new2,average2,'cyan','LineWidth',2)
anew=(log10(W_esc_half)>4);
%semilogx(W_esc_half(anew),a4(anew),'blue','LineWidth',2)
xlabel('Escape energy [eV]','FontSize',lsize)
ylabel('a_{w}','FontSize',lsize)
xlim([1,10^8])
%ylim([min(a42)-0.1 max(a42)+0.5])
XTick = [ 10^(-1)   10^(1)  10^(3)  10^(5)  10^(7) 10^(9)];
set(gca,'xtick',XTick)
%YTick = [10^(-12) 10^(-10) 10^(-8) 10^(-6) 10^(-4) 10^(-2) 10^(0)];
%set(gca,'ytick',YTick)
set(gca, 'Fontsize', nsize)
set(gcf,'paperpositionmode','auto');
set(gcf,'windowstyle','normal');
set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))
set(gca,'fontweight','normal')
box on
opts.Colors     = get(groot,'defaultAxesColorOrder');
opts.saveFolder = 'img/';
opts.width      = 12;
opts.height     = 10;
opts.fontType   = 'Times';
% creating the zoom-in inset
ax=axes;
set(ax,'units','normalized','position',[0.18,0.69,0.25,0.25])
box(ax,'on')
for u=1:end2
loglog( tM(start:ending),Dw1sm(u,start:ending),'LineWidth',0.5)
hold on
end
set(gca,'fontweight','normal','FontSize', 9)
set(ax,'xlim',[tM(start),tM(ending)],'ylim',[min(Dw1sm(:)),10*max(Dw1sm(:))])
%XTick = [10^(-1) 10^(0) ];%  10^(8)];
%set(gca,'xtick',XTick)
YTick = [ 10^0 10^2 10^(4) 10^(6) 10^(8) 10^(10)  10^(12) 10^(14)];
set(gca,'ytick',YTick)
str = sprintf('<Äw>-t');
title(str ,'FontSize', 10);
set(get(gca,'title'),'Position',[2*tM(start) 0.1*max(Dw1sm(:)) 0.8])
ax = gca;
ax.YAxisLocation = 'right';
hold off
if imodel == 1

    temp=['png_files_upper\f6c.png']; saveas(gca,temp);
    temp1=['eps_files_upper\f6c.eps']; saveas(gca,temp1); 
    temp2=['fig_files_upper\f6c.fig']; saveas(gca,temp2); 
    
 elseif imodel == 12
          temp=['new_png\f6c.png']; saveas(gca,temp);
    temp1=['new_eps\f6c.eps']; saveas(gca,temp1); 
    temp2=['new_fig\f6c.fig']; saveas(gca,temp2); 
 elseif imodel == 2 
     
      temp=['new_png_rec\f6c P=',num2str(P_reconnect),'.png']; saveas(gca,temp);
    temp1=['new_eps_rec\f6c P=',num2str(P_reconnect),'.eps']; saveas(gca,temp1); 
    temp2=['new_fig_rec\f6c P=',num2str(P_reconnect),'.fig']; saveas(gca,temp2); 
 end

%% for w_esc-dw^2
for jN=1:end2
figure(11)
disp(['percentage: ' num2str(jN/end2)])
fitResults3 = polyfit(log10(tM1(jN,start:ending)),log10(Dw2sm(jN,start:ending)),1); %prosexe edw to /2
pol3=polyval(fitResults3,log10(tM1(jN,start:ending)));
a4(jN)=fitResults3(1);
a44=a4>-10; %non nan values of a1 and W_esc
end
scatter((W_esc_half(a44)),a4(a44),'MarkerEdgeColor',[0 .5 .5],...
'MarkerFaceColor',[0 .7 .7],...
'LineWidth',1.5)
set(gca, 'YScale', 'linear')
set(gca, 'XScale', 'log')
hold on
%loglog((W_esc_half(a44)),a4(a44),'blue','LineWidth',0.1)
ans2=nnz(log10(W_esc_half)<4);
a42=a4(a44);
av2=mean(a4(1:ans2));
average2=av2*ones(1,ans2+1);
new2=(W_esc_half(1:ans2+1));
semilogx(new2,average2,'cyan','LineWidth',2)
anew=(log10(W_esc_half)>4);
%semilogx(W_esc_half(anew),a44(anew),'blue','LineWidth',2)
xlabel('Escape energy [eV]','FontSize',lsize)
ylabel('a_{w^{2}}','FontSize',lsize)
xlim([10^(-1),10^8])
ylim([min(a42)-0.1 ,max(a42)+1])
XTick = [ 10^(-1)   10^(1)  10^(3)  10^(5)  10^(7) 10^(9)];
set(gca,'xtick',XTick)
set(gca, 'Fontsize', nsize)
set(gcf,'paperpositionmode','auto');
set(gcf,'windowstyle','normal');
set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))
set(gca,'fontweight','normal')
box on
opts.Colors     = get(groot,'defaultAxesColorOrder');
opts.saveFolder = 'img/';
opts.width      = 12;
opts.height     = 10;
opts.fontType   = 'Times';
% creating the zoom-in inset
ax=axes;
set(ax,'units','normalized','position',[0.18,0.7,0.25,0.25])
box(ax,'on')
for u=1:end2
loglog( tM(start:ending),Dw2sm(u,start:ending),'LineWidth',0.5)
hold on
end
set(gca,'fontweight','normal','FontSize', 9)
set(ax,'xlim',[0.9*tM(start),tM(ending)],'ylim',[min(Dw2sm(start,:)),10*max(Dw2sm(:))])
%XTick = [10^(-2) 10^(-1) 10^(0)];%  10^(8)];
%set(gca,'xtick',XTick)
YTick = [ 10^(5) 10^(8) 10^(11) 10^(14)  10^(17) 10^(20)];
set(gca,'ytick',YTick)
str = sprintf('<Äw^{2}>-t');
title(str ,'FontSize', 10);
set(get(gca,'title'),'Position',[2*tM(start) 0.01*max(Dw2sm(:)) 1])
ax = gca;
ax.YAxisLocation = 'right';
if imodel == 1

    temp=['png_files_upper\f7b.png']; saveas(gca,temp);
    temp1=['eps_files_upper\f7b.eps']; saveas(gca,temp1); 
    temp2=['fig_files_upper\f7b.fig']; saveas(gca,temp2); 
    
 elseif imodel == 12
          temp=['new_png\f7b.png']; saveas(gca,temp);
    temp1=['new_eps\f7b.eps']; saveas(gca,temp1); 
    temp2=['new_fig\f7b.fig']; saveas(gca,temp2); 
 elseif imodel == 2 
     
      temp=['new_png_rec\f7b P=',num2str(P_reconnect),'.png']; saveas(gca,temp);
    temp1=['new_eps_rec\f7b P=',num2str(P_reconnect),'.eps']; saveas(gca,temp1); 
    temp2=['new_fig_rec\f7b P=',num2str(P_reconnect),'.fig']; saveas(gca,temp2); 
 end
