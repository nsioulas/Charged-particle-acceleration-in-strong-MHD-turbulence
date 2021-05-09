 Eeff_mean_coeff=1.e3;     % Multiplied by ED to give Eeff_mean
Eeef_sigma_coeff=(1/3)*Eeff_mean_coeff; % Multiplied by ED to give Eeff_mean

Eeff_mean=Eeff_mean_coeff*ED;
Eeef_sigma=Eeef_sigma_coeff*ED;
 Eeff = normrnd(Eeff_mean,Eeef_sigma,1.e3);
 
 

lsize=16;  
  nsize=16;
start=2;    rd1start=2;      end2=30; %number of bins
ending=end2;     rd1end=end2;  %nP=length(rmsd);
tacc=2;


dr1=Eeff(:);
dr1=dr1(dr1>0);
dr1_=dr1(:);
nP=length(dr1);
%split2=linspace((min(Wkinet_eV(:,tacc))),(max(Wkinet_eV(:,tacc))),end2);
%split2=logspace(0,log10(max(Wkinet_eV(:,tacc))),end2);
split2=logspace(log10(min(dr1)),log10(max(dr1)),end2);
dr_half=zeros(1,end2);
p1=zeros(1,end2);






for jN=1
 
  dr1=dr1_;
   dr_half(jN)=(1/2)*(split2(jN));
for i=1:nP
        if dr1(i,1)>=split2(jN)

        
dr1(i,1)=0;

       end
end

    
        p1(1,jN)=nnz(dr1(:));
 
     
    
   % dw1=dw1_;
 %dr1_=dr1;
end

dr1=dr1_;
 
for jN=2:end2
    disp(['percentage: ' num2str(jN/end2)])
   dr1=dr1_;
    dr_half(jN)=(1/2)*(split2(jN)-split2(jN-1))+split2(jN-1);
for i=1:nP
       
        if dr1(i,1)<=split2(jN-1) || dr1(i,1)>=split2(jN)

       
dr1(i,1)=0;

        end
       
end
       
p1(1,jN)=nnz(dr1(:));
        
    dr1=dr1_;
end

total=sum(p1);
relative=zeros(1,end2);
for jN=1:end2
     if jN==1
         %relative(1,jN)=p1(1,jN)/(total*(split2(jN)));
      relative(1,jN)=p1(1,jN)/(total);
     else
    %relative(1,jN)=p1(1,jN)/(total*(split2(jN)-split2(jN-1)));
    relative(1,jN)=p1(1,jN)/(total);
    
     end
end


figure(1)
fitResults1 = polyfit(log10(dr_half(start:ending)),log10(relative(start:ending)),1);
pol=polyval(fitResults1,log10(dr_half(start:ending)));
a=fitResults1(1);
b=fitResults1(2);
polyfit_str = ['<(?r)^2> ~ ? ^{'  (num2str(a)) ' } '] % + ' num2str(b)]  polyfit_str will be : y = 4*x + 2


fit1=dr_half(start:ending).^a*10^(b);




plot((dr_half(rd1start:rd1end)),(relative(rd1start:rd1end)),'Color',[0.2941, 0.4447, 0.6494],'LineWidth',2)
%loglog((dr_half(rd1start:rd1end))/(ED),(relative(rd1start:rd1end)),'w-','LineWidth',2)

hold on
loglog(dr_half(1),fit1(1),'w','LineWidth',2)
loglog(dr_half(1),fit1(1)/1.0001,'w','LineWidth',2)
legend({ strcat('Distribution of E_{eff}'),strcat('Mean :',num2str(Eeff_mean_coeff,2),' \cdot E_{D}'),strcat('Sigma :',num2str(Eeef_sigma_coeff,3),' \cdot E_{D}')},'FontSize',12,'Location','NorthWest')%,'fit for [0-20 sec]'})
xlabel('E_{eff}  [statV/cm]','FontSize',lsize)
ylabel('p(E_{eff})','FontSize',lsize)

%xlim([5.e-4,5.e-3])
%ylim([0,2.5*10^(-1)])
%XTick = [ 10^(2) 10^(4)  10^(6) 10^(8) 10^(10)];
%set(gca,'xtick',XTick)

%YTick = [ 10^(-14) 10^(-12) 10^(-10) 10^(-8) 10^(-6) 10^(-4) 10^(-2) 10^(0)];
%set(gca,'ytick',YTick)


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


 
hold off



if imodel == 1

    temp=['png_files_upper\f30.png']; saveas(gca,temp);
    temp1=['eps_files_upper\f30.eps']; saveas(gca,temp1); 
    temp2=['fig_files_upper\f30.fig']; saveas(gca,temp2); 
    
 elseif imodel == 12
          temp=['new_png\f30.png']; saveas(gca,temp);
    temp1=['new_eps\f30.eps']; saveas(gca,temp1); 
    temp2=['new_fig\f30.fig']; saveas(gca,temp2); 
 elseif imodel == 2 
     
      temp=['new_png_rec\f30 P=',num2str(P_reconnect),'.png']; saveas(gca,temp);
    temp1=['new_eps_rec\f30 P=',num2str(P_reconnect),'.eps']; saveas(gca,temp1); 
    temp2=['new_fig_rec\f30 P=',num2str(P_reconnect),'.fig']; saveas(gca,temp2); 
 end

