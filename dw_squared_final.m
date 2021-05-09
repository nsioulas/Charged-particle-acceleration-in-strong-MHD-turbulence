
lsize = 16; % Label fontsize
    nsize = 16; % Axis fontsize
    
    clear Wkinet W0 W_esc
   prosexe=prosexe+1; 

start=1;           rd1start=450;           rd1start1=800;
ending=ntM;       rd1ending=800;     rd1ending1=890;
%ntM=ending;
if prosexe ==1
    
    dw1=zeros(nP,ntM); 
  dw2=zeros(nP,ntM);                    
%W0_eV = W0 .* erg2eV;
    
for i=1:nP
   for j=1:ntM
       
       if Wkinet_eV(i,j)>0
       dw1(i,j)=(Wkinet_eV(i,j)-W0_eV(i));
       dw2(i,j)=(Wkinet_eV(i,j)-W0_eV(i))^2;
      
       else
         dw1(i,j)=0;
      dw2(i,j)=0;
       end 
    end 
end

prosexe= prosexe+1;


kinit=zeros(1,ntM);
kinit2=zeros(1,ntM);
kinit3=zeros(1,ntM);
for jN=1:ntM
nP_in=nnz(Wkinet_eV(:,jN));
fprintf('for timestep %d fraction of particles remaining in the box is %2.2f \n',jN,nP_in/nP);
kinit(jN)=(sum(Wkinet_eV(:,jN)))/nP_in;
kinit2(jN)=(sum(dw1(:,jN)))/nP_in;
kinit3(jN)=(sum(dw2(:,jN)))/nP_in;
end
clear dw1 dw2
end

%{
yval_peak=findpeaks(diff(kinit));

 [~,itM] = min(abs(yval_peak(1)-kinit));            % finds a slope change in the graph  Usefull
[~,itM1] = min(abs(yval_peak(2)-kinit)); 

rd1start=1;
rd1ending=itM;


rd1start=itM;
rd1ending=itM1;

%}

c12 = struct('rr', [0.9047, 0.1918, 0.1988], ...  %Your required color
    'bb', [0.2941, 0.5447, 0.7494], ... %Your required color
    'um', [0.0824, 0.1294, 0.4196], ... %ultra marine
    'cy', [0.0910, 0.9725, 0.9412], ... %cyan
    'or', [0.985, 0.597, 0.268] );   %orange

figure(2) % for <Wkinet>-time
fitResults1 = polyfit(log10(tM(rd1start:rd1ending)),log10(kinit(rd1start:rd1ending)),1);
pol=polyval(fitResults1,log10(tM(rd1start:rd1ending)));
a=fitResults1(1);
b=fitResults1(2);
fit3=tM(rd1start:(length(pol)+rd1start-1)).^a*10^(b);


fitResults1 = polyfit(log10(tM(rd1start1:rd1ending1)),log10(kinit(rd1start1:rd1ending1)),1);
pol1=polyval(fitResults1,log10(tM(rd1start1:rd1ending1)));
a1=fitResults1(1);
b1=fitResults1(2);
fit4=tM(rd1start1:(length(pol1)+rd1start1-1)).^a1*10^(b1);



loglog((tM(start:ending)),(kinit(start:ending)),'Color',c12.or,'LineWidth',2)
hold on
loglog((tM(rd1start:(length(pol)+rd1start-1))),fit3,'black--','LineWidth',2)
%loglog((tM(rd1start1:(length(pol1)+rd1start1-1))),fit4,'b--','LineWidth',2)
legend({'Kinetic energy ', strcat('Power-Law fit ,index:_{<w>}=',(num2str(a))), strcat('Power-Law fit ,index:_{<w>}=',(num2str(a1)))},'Location','Southeast','FontSize',14)
xlabel('time [sec]', 'Fontsize', lsize)
ylabel('W_{kinet} [eV]', 'Fontsize', lsize)
xlim([1.e-7 1.e-0])
ylim([10^2 1.01*10^3])
xticks = [1.e-7  1.e-5  1.e-3  1.e-1 ];
 set(gca,'XTick',xticks)
set(gca, 'Fontsize', nsize)
set(gcf,'paperpositionmode','auto');
    set(gcf,'windowstyle','normal');
    set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))
set(gca,'fontweight','normal')


opts.Colors     = get(groot,'defaultAxesColorOrder');
opts.saveFolder = 'img/';
opts.width      = 12;
opts.height     = 10;
opts.fontType   = 'Times';
if tfin>0.9
    if imodel == 1
    temp=['png_files_upper\f6b.png']; saveas(gca,temp);
    temp1=['eps_files_upper\f6b.eps']; saveas(gca,temp1); 
    temp2=['fig_files_upper\f6b.fig']; saveas(gca,temp2); 
    
 elseif imodel == 12
          temp=['new_png\f6b.png']; saveas(gca,temp);
    temp1=['new_eps\f6b.eps']; saveas(gca,temp1); 
    temp2=['new_fig\f6b.fig']; saveas(gca,temp2); 
 elseif imodel == 2 
     
      temp=['new_png_rec\f6b P=',num2str(P_reconnect),'.png']; saveas(gca,temp);
    temp1=['new_eps_rec\f6b P=',num2str(P_reconnect),'.eps']; saveas(gca,temp1); 
    temp2=['new_fig_rec\f6b P=',num2str(P_reconnect),'.fig']; saveas(gca,temp2); 
    end
else
    if imodel == 1
     temp=['png_files_upper\f6b_small.png']; saveas(gca,temp);
    temp1=['eps_files_upper\f6b_small.eps']; saveas(gca,temp1); 
    temp2=['fig_files_upper\f6b_small.fig']; saveas(gca,temp2); 
    
 elseif imodel == 12
          temp=['new_png\f6b_small.png']; saveas(gca,temp);
    temp1=['new_eps\f6b_small.eps']; saveas(gca,temp1); 
    temp2=['new_fig\f6b_small.fig']; saveas(gca,temp2); 
 elseif imodel == 2 
     
      temp=['new_png_rec\f6b_small P=',num2str(P_reconnect),'.png']; saveas(gca,temp);
    temp1=['new_eps_rec\f6b_small P=',num2str(P_reconnect),'.eps']; saveas(gca,temp1); 
    temp2=['new_fig_rec\f6b_small P=',num2str(P_reconnect),'.fig']; saveas(gca,temp2); 

    end
end

figure(3) % for <W>-time
fitResults2 = polyfit(log10(tM(rd1start:rd1ending)),log10(kinit2(rd1start:rd1ending)),1);
pol=polyval(fitResults2,log10(tM(rd1start:rd1ending)));
c=fitResults2(1);
d=fitResults2(2);
fit4=tM(rd1start:(length(pol)+rd1start-1)).^c*10^(d);


fitResults1 = polyfit(log10(tM(rd1start1:rd1ending1)),log10(kinit2(rd1start1:rd1ending1)),1);
pol1=polyval(fitResults1,log10(tM(rd1start1:rd1ending1)));
a1=fitResults1(1);
b1=fitResults1(2);
fit5=tM(rd1start1:(length(pol1)+rd1start1-1)).^a1*10^(b1);


loglog((tM(start:ending)),(kinit2(start:ending)),'Color',c12.or,'LineWidth',2)
%text(0.15,2*10^11,polyfit_str,'FontSize', 12, 'Color', 'g', 'FontWeight', 'normal');
hold on
loglog((tM(rd1start:(length(pol)+rd1start-1))),fit4,'black--','LineWidth',2)
%loglog((tM(rd1start1:(length(pol1)+rd1start1-1))),fit5,'b--','LineWidth',2)
legend({'Mean energy displacement ', strcat('Power-Law fit ,index:_{<w>}=',(num2str(c))), strcat('Power-Law fit ,index:_{<w>}=',(num2str(a1)))},'Location','SouthEast','FontSize',14)
xlabel('Time [sec]', 'Fontsize', lsize)
ylabel('<Äw>  [eV]', 'Fontsize', lsize)
xlim([1.e-7 1.e-0])
%ylim([10^3 10^4])
xticks = [1.e-7  1.e-5  1.e-3  1.e-1 ];
 set(gca,'XTick',xticks)
set(gca, 'Fontsize', nsize)
set(gcf,'paperpositionmode','auto');
    set(gcf,'windowstyle','normal');
    set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))
set(gca,'fontweight','normal')


opts.Colors     = get(groot,'defaultAxesColorOrder');
opts.saveFolder = 'img/';
opts.width      = 12.1;
opts.height     = 10;
opts.fontType   = 'Times';


if tfin>0.9
    if imodel == 1
    temp=['png_files_upper\f6b_2.png']; saveas(gca,temp);
    temp1=['eps_files_upper\f6b_2.eps']; saveas(gca,temp1); 
    temp2=['fig_files_upper\f6b_2.fig']; saveas(gca,temp2); 
    
 elseif imodel == 12
          temp=['new_png\f6b_2.png']; saveas(gca,temp);
    temp1=['new_eps\f6b_2.eps']; saveas(gca,temp1); 
    temp2=['new_fig\f6b_2.fig']; saveas(gca,temp2); 
 elseif imodel == 2 
     
      temp=['new_png_rec\f6b_2 P=',num2str(P_reconnect),'.png']; saveas(gca,temp);
    temp1=['new_eps_rec\f6b_2 P=',num2str(P_reconnect),'.eps']; saveas(gca,temp1); 
    temp2=['new_fig_rec\f6b_2 P=',num2str(P_reconnect),'.fig']; saveas(gca,temp2); 
    end
else
    if imodel == 1
     temp=['png_files_upper\f6b_2_small.png']; saveas(gca,temp);
    temp1=['eps_files_upper\f6b_2_small.eps']; saveas(gca,temp1); 
    temp2=['fig_files_upper\f6b_2_small.fig']; saveas(gca,temp2); 
    
 elseif imodel == 12
          temp=['new_png\f6b_2_small.png']; saveas(gca,temp);
    temp1=['new_eps\f6b_2_small.eps']; saveas(gca,temp1); 
    temp2=['new_fig\f6b_2_small.fig']; saveas(gca,temp2); 
 elseif imodel == 2 
     
      temp=['new_png_rec\f6b_2_small P=',num2str(P_reconnect),'.png']; saveas(gca,temp);
    temp1=['new_eps_rec\f6b_2_small P=',num2str(P_reconnect),'.eps']; saveas(gca,temp1); 
    temp2=['new_fig_rec\f6b_2_small P=',num2str(P_reconnect),'.fig']; saveas(gca,temp2); 

    end
end


figure(4) % for <W2>-time
fitResults2 = polyfit(log10(tM(rd1start:rd1ending)),log10(kinit3(rd1start:rd1ending)),1);
pol=polyval(fitResults2,log10(tM(rd1start:rd1ending)));
c=fitResults2(1);
d=fitResults2(2);
fit4=tM(rd1start:(length(pol)+rd1start-1)).^c*10^(d);


fitResults1 = polyfit(log10(tM(rd1start1:rd1ending1)),log10(kinit3(rd1start1:rd1ending1)),1);
pol1=polyval(fitResults1,log10(tM(rd1start1:rd1ending1)));
a1=fitResults1(1);
b1=fitResults1(2);
fit6=tM(rd1start1:(length(pol1)+rd1start1-1)).^a1*10^(b1);


loglog((tM(start:ending)),(kinit3(start:ending)),'Color',c12.or,'LineWidth',2)
%text(0.15,2*10^11,polyfit_str,'FontSize', 12, 'Color', 'g', 'FontWeight', 'normal');
hold on
loglog((tM(rd1start:(length(pol)+rd1start-1))),fit4,'black--','LineWidth',2)
%loglog((tM(rd1start1:(length(pol1)+rd1start1-1))),fit6,'c--','LineWidth',2)
legend({'Mean squared displacement in energy ', strcat('Power-Law fit ,index: a_{<w^{2}>}=',(num2str(c,'%2.2f'))), strcat('Power-Law fit ,index: a_{<w^{2}>}=',(num2str(a1','%2.2f'))),'On frac,but no DW ', strcat('Power-Law fit ,index: a_{<w^{2}>}=',(num2str(c,'%2.2f'))), strcat('Power-Law fit ,index: a_{<w^{2}>}=',(num2str(a1','%2.2f'))),'Energiz at t=0 ', strcat('Power-Law fit ,index: a_{<w^{2}>}=',(num2str(c,'%2.2f'))), strcat('Power-Law fit ,index: a_{<w^{2}>}=',(num2str(a1','%2.2f')))},'Location','SouthEast','FontSize',12)
xlabel('Time [sec]', 'Fontsize', lsize)
ylabel('<(Äw)^{2}>  [eV^{2}]', 'Fontsize', lsize)
xlim([1.e-7 1.e-0])
xticks = [1.e-7  1.e-5  1.e-3  1.e-1 ];
 set(gca,'XTick',xticks)
set(gca, 'Fontsize', nsize)
set(gcf,'paperpositionmode','auto');
    set(gcf,'windowstyle','normal');
    set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))
set(gca,'fontweight','normal')


opts.Colors     = get(groot,'defaultAxesColorOrder');
opts.saveFolder = 'img/';
opts.width      = 12;
opts.height     = 10;
opts.fontType   = 'Times';

if imodel == 1

    temp=['png_files_upper\f7a.png']; saveas(gca,temp);
    temp1=['eps_files_upper\f7a.eps']; saveas(gca,temp1); 
    temp2=['fig_files_upper\f7a.fig']; saveas(gca,temp2); 
    
 elseif imodel == 12
          temp=['new_png\f7a.png']; saveas(gca,temp);
    temp1=['new_eps\f7a.eps']; saveas(gca,temp1); 
    temp2=['new_fig\f7a.fig']; saveas(gca,temp2); 
 elseif imodel == 2 
     
      temp=['new_png_rec\f7a P=',num2str(P_reconnect),'.png']; saveas(gca,temp);
    temp1=['new_eps_rec\f7a P=',num2str(P_reconnect),'.eps']; saveas(gca,temp1); 
    temp2=['new_fig_rec\f7a P=',num2str(P_reconnect),'.fig']; saveas(gca,temp2); 
end
 


if tfin>0.9
    if imodel == 1
    temp=['png_files_upper\f7a.png']; saveas(gca,temp);
    temp1=['eps_files_upper\f7a.eps']; saveas(gca,temp1); 
    temp2=['fig_files_upper\f7a.fig']; saveas(gca,temp2); 
    
 elseif imodel == 12
          temp=['new_png\f7a.png']; saveas(gca,temp);
    temp1=['new_eps\f7a.eps']; saveas(gca,temp1); 
    temp2=['new_fig\f7a.fig']; saveas(gca,temp2); 
 elseif imodel == 2 
     
      temp=['new_png_rec\f7a P=',num2str(P_reconnect),'.png']; saveas(gca,temp);
    temp1=['new_eps_rec\f7a P=',num2str(P_reconnect),'.eps']; saveas(gca,temp1); 
    temp2=['new_fig_rec\f7a P=',num2str(P_reconnect),'.fig']; saveas(gca,temp2); 
    end
else
    if imodel == 1
     temp=['png_files_upper\f7a_small.png']; saveas(gca,temp);
    temp1=['eps_files_upper\f7a_small.eps']; saveas(gca,temp1); 
    temp2=['fig_files_upper\f7a_small.fig']; saveas(gca,temp2); 
    
 elseif imodel == 12
          temp=['new_png\f7a_small.png']; saveas(gca,temp);
    temp1=['new_eps\f7a_small.eps']; saveas(gca,temp1); 
    temp2=['new_fig\f7a_small.fig']; saveas(gca,temp2); 
 elseif imodel == 2 
     
      temp=['new_png_rec\f7a_small P=',num2str(P_reconnect),'.png']; saveas(gca,temp);
    temp1=['new_eps_rec\f7a_small P=',num2str(P_reconnect),'.eps']; saveas(gca,temp1); 
    temp2=['new_fig_rec\f7a_small P=',num2str(P_reconnect),'.fig']; saveas(gca,temp2); 

    end
end

