
 lsize = 16; % Label fontsize
    nsize = 16; % Axis fontsize
 

rd1start=1;  start=1;      start1=600;           

rd1end=1.e3; ending=780;     ending1=750;     


prosexe1=prosexe1+1;
if prosexe1 ==1
fprintf(' ---- First definition ------ \n');

rd1=zeros(1,rd1end);
for jN=1:ntM
  nP_in=nnz(rmsd(:,jN));
  %nP_in = nnz(is_in(:,jN)); %HI: the safe way of counting   
  fprintf('for timestep %d fraction of particles remaining in the box is %2.2f \n',jN,nP_in/nP);
  rd1(jN)=sum(rmsd(1:nP,jN))/nP_in;
end
prosexe1= prosexe1+1;
end
%{
yval_peak=findpeaks(diff(rd1));

 [~,itM] = min(abs(yval_peak(1)-rd1));            % finds a slope change in the graph  Usefull
[~,itM1] = min(abs(yval_peak(2)-rd1)); 
start=1;
ending=itM;

start1=itM;
ending1=itM1;
%}

figure(4)

fitResults1 = polyfit(log10(tM(start:ending)),log10(rd1(start:ending)),1);
pol=polyval(fitResults1,log10(tM(start:ending)));
a=fitResults1(1);
b=fitResults1(2);
polyfit_str = ['<(?r)^2> ~ ? ^{'  (num2str(a)) ' } '] % + ' num2str(b)]  polyfit_str will be : y = 4*x + 2

fitResults1 = polyfit(log10(tM(start1:ending1)),log10(rd1(start1:ending1)),1);
pol1=polyval(fitResults1,log10(tM(start1:ending1)));
c=fitResults1(1);
d=fitResults1(2);
polyfit_str = ['<(?r)^2> ~ ? ^{'  (num2str(c)) ' } '] % + ' num2str(b)]  polyfit_str will be : y = 4*x + 2







fit1=tM(start:ending).^a*10^(b);
fit2=tM(start1:ending1).^c*10^(d);



c12 = struct('rr', [0.9047, 0.1918, 0.1988], ...  %Your required color
    'bb', [0.2941, 0.5447, 0.7494], ... %Your required color
    'um', [0.0824, 0.1294, 0.4196], ... %ultra marine
    'cy', [0.0910, 0.9725, 0.9412], ... %cyan
    'or', [0.985, 0.197, 0.08] );   %orange

%
loglog((tM(rd1start:rd1end)),(rd1(rd1start:rd1end)),'Color',[0.2941, 0.4447, 0.6494],'LineWidth',2)
hold on
loglog((tM(start:(length(pol)+start-1))),fit1,'--','Color',[0.9 0.3 0.35],'LineWidth',2)
%loglog((tM(start1:(length(pol1)+start1-1))),fit2,'black--','LineWidth',2)
legend({ strcat('Mean squared displacement'),strcat(' Power-law fit, a_{r}: ', num2str(a,'%2.2f')),strcat(' Power-law fit, a_{r}: ', num2str(c,'%2.2f'))},'Location', 'SouthEast','FontSize',14)%,'fit for [0-20 sec]'})
legend({'Mean squared displacement ', strcat('Power-Law fit ,index: a_{<r^{2}>}=',(num2str(a,'%2.2f'))), strcat('Power-Law fit ,index: a_{<r^{2}>}=',(num2str(c','%2.2f'))),'On frac,but no DW ', strcat('Power-Law fit ,index: a_{<r^{2}>}=',(num2str(a,'%2.2f'))), strcat('Power-Law fit ,index: a_{<r^{2}>}=',(num2str(c','%2.2f'))),'Energiz at t=0 ', strcat('Power-Law fit ,index: a_{<r^{2}>}=',(num2str(a,'%2.2f'))), strcat('Power-Law fit ,index: a_{<r^{2}>}=',(num2str(c','%2.2f')))},'Location','SouthEast','FontSize',12)

 ylabel('<(Är)^2> [cm^2]', 'Fontsize', lsize)
xlabel('Time [sec]', 'Fontsize', lsize)
xlim([1.e-7 1.e-0])
YTick = [ 10^0  10^(4)  10^(8)  10^(12)  10^(16)  10^(20)];
set(gca,'ytick',YTick)

XTick = [ 1.e-7 1.e-5 10^(-3) 10^(-1) 10^(0)];% 10^(-1)  10^(0)];
set(gca,'xtick',XTick)
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
    temp=['png_files_upper\f2a.png']; saveas(gca,temp);
    temp1=['eps_files_upper\f2a.eps']; saveas(gca,temp1); 
    temp2=['fig_files_upper\f2a.fig']; saveas(gca,temp2); 
    
 elseif imodel == 12
          temp=['new_png\f2a.png']; saveas(gca,temp);
    temp1=['new_eps\f2a.eps']; saveas(gca,temp1); 
    temp2=['new_fig\f2a.fig']; saveas(gca,temp2); 
 elseif imodel == 2 
     
      temp=['new_png_rec\f2a P=',num2str(P_reconnect),'.png']; saveas(gca,temp);
    temp1=['new_eps_rec\f2a P=',num2str(P_reconnect),'.eps']; saveas(gca,temp1); 
    temp2=['new_fig_rec\f2a P=',num2str(P_reconnect),'.fig']; saveas(gca,temp2); 
    end
else
    if imodel == 1
     temp=['png_files_upper\f2a_small.png']; saveas(gca,temp);
    temp1=['eps_files_upper\f2a_small.eps']; saveas(gca,temp1); 
    temp2=['fig_files_upper\f2a_small.fig']; saveas(gca,temp2); 
    
 elseif imodel == 12
          temp=['new_png\f2a_small.png']; saveas(gca,temp);
    temp1=['new_eps\f2a_small.eps']; saveas(gca,temp1); 
    temp2=['new_fig\f2a_small.fig']; saveas(gca,temp2); 
 elseif imodel == 2 
     
      temp=['new_png_rec\f2a_small P=',num2str(P_reconnect),'.png']; saveas(gca,temp);
    temp1=['new_eps_rec\f2a_small P=',num2str(P_reconnect),'.eps']; saveas(gca,temp1); 
    temp2=['new_fig_rec\f2a_small P=',num2str(P_reconnect),'.fig']; saveas(gca,temp2); 

    end
end
