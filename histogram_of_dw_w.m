  
 lsize = 16; % Label fontsize
    nsize = 16; % Axis fontsize

 maen1=0;                 %1 :shows mean of Äw
                         %0 :shows median of Äw
                         %2 :both
                         %1: dw/w gia kathe swmatidio ksexwrista
 mono=2;                 %2: oxi dw/w alla to P(dw)
                         %  avrage of dw/w gia kathe swmatidio

devide_by_area=true;


start=9;    rd1start=2;      end2=100; %number of bins
ending=end2-1;     rd1end=end2;  %nP=length(rmsd);
tacc=2;
if mono==1      %  dw/w gia kathe swmatidio ksexwrista
 C1 = a11(a11>0)';
dr1=C1;
dr1_=C1;
elseif mono==2  %  oxi dw/w alla to P(dw)
    AA=cell2mat(dW_Cell)*erg2eV;
  AA = rmmissing(AA);
  AA=AA(AA>0);
    dr1=AA;
   dr1_=AA;  
elseif mono==3   % avrage of dw/w gia kathe swmatidio
  dr1=mDelW_W;
dr1_=mDelW_W;  
end
nP=length(dr1);
if (mono==1) || (mono==3)
 %split2=linspace(min(dr1),max(dr1),end2);
 split2=logspace(log10(min(dr1)),log10(max(dr1)),end2);
else
split2=logspace(log10(min(dr1)),log10(max(dr1)),end2);

end


dr_half=zeros(1,end2);
p1=zeros(1,end2);




for jN=1
 
  dr1=dr1_;
   dr_half(jN)=(1/2)*(split2(jN));
for i=1:nP
        if dr1(1,i)>=split2(jN)

        
dr1(1,i)=0;

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
       
        if dr1(1,i)<=split2(jN-1) || dr1(1,i)>=split2(jN)

       
dr1(1,i)=0;

        end
       
end
       
p1(1,jN)=nnz(dr1(:));
        
    dr1=dr1_;
end

total=sum(p1);
relative=zeros(1,end2);
if devide_by_area
for jN=1:end2
     if jN==1
         relative(1,jN)=p1(1,jN)/(total*(split2(jN)));
         %relative(1,jN)=p1(1,jN)/total;
     else
    relative(1,jN)=p1(1,jN)/(total*(split2(jN)-split2(jN-1)));
    %relative(1,jN)=p1(1,jN)/total;
    
     end
end

else
    for jN=1:end2
     if jN==1
         %relative(1,jN)=p1(1,jN)/(total*(split2(jN)));
         relative(1,jN)=p1(1,jN)/total;
     else
    %relative(1,jN)=p1(1,jN)/(total*(split2(jN)-split2(jN-1)));
    relative(1,jN)=p1(1,jN)/total;
    
     end
end
    
end


figure(1)
fitResults1 = polyfit(log10(dr_half(start:ending)),log10(relative(start:ending)),1);
pol=polyval(fitResults1,log10(dr_half(start:ending)));
a=fitResults1(1);
b=fitResults1(2);
polyfit_str = ['<(?r)^2> ~ ? ^{'  (num2str(a)) ' } '] % + ' num2str(b)]  polyfit_str will be : y = 4*x + 2


fit1=dr_half(start:ending).^a*10^(b);



loglog((dr_half(rd1start:rd1end)),(relative(rd1start:rd1end)),'Color',[0 0.75 1],'LineWidth',2)
hold on




if imodel==2 || imodel==12
    if maen1==1 
[itM2,itM] = min(abs(mean(dr1)-dr_half));
    elseif maen1==0
[itM4,itM3] = min(abs(median(dr1)-dr_half));
 elseif maen1==2
     [itM2,itM] = min(abs(mean(dr1)-dr_half));
     [itM4,itM3] = min(abs(median(dr1)-dr_half));
    end
    
    [aa1,bb1]=sort(relative);
   % bb=mean(dr1)*ones(length(relative),1);
    %cc=logspace(log10(min(relative)),log10(max(relative)),(length(relative)));
    
     loglog((dr_half(start:(length(pol)+start-1))),fit1,'black--','LineWidth',2)
     if maen1==1
     scatter(mean(dr1),relative(itM),'filled','red','LineWidth',10)
     elseif maen1==0
     scatter(median(dr1),relative(itM3),'filled','blue','LineWidth',10)
     else
     scatter(mean(dr1),relative(itM),'filled','red','LineWidth',10)
     scatter(median(dr1),relative(itM3),'filled','blue','LineWidth',10)
     end
elseif imodel==1
    [itM2,itM] = min(abs(median(dr1)-dr_half));
    [itM4,itM3] = min(abs(median(dr1)-dr_half));
    [aa1,bb1]=sort(relative);
    bb=median(dr1)*ones(length(relative),1);
    cc=logspace(log10(min(relative)),log10(max(relative)),(length(relative)));
    if maen1
      scatter(mean(dr1),relative(itM),'filled','red','LineWidth',10) 
    else
    scatter(median(dr1),relative(itM),'filled','red','LineWidth',10)
    end
    
    
end
    %loglog(bb,cc,'black','LineWidth',0.01)
hold off
if (mono==1) || (mono==3)
xlabel('Äw/w')
ylabel('P(Äw/w)')
%format shortE
if imodel==2
    if maen1==1 
      legend({ strcat('Distribution of Äw/w '),strcat(' Power-law fit, k:',num2str(a,2)),strcat('<Äw/w> :',num2str(mean(dr1),'%.2d'))},'FontSize',12,'Location','NorthEast')%,'fit for [0-20 sec]'})
    elseif maen1==0
        legend({ strcat('Distribution of Äw/w '),strcat(' Power-law fit, k:',num2str(a,2)),strcat('Median of Äw/w : ', num2str(median(dr1),'%.2d'))},'FontSize',12,'Location','SouthWest')%,'fit for [0-20 sec]'})
    elseif maen1==2
        legend({ strcat('Distribution of Äw/w '),strcat(' Power-law fit, k:',num2str(a,2)),strcat('Mean of Äw/w : ', num2str(mean(dr1,2),'%.2d')),strcat('Median of Äw/w : ', num2str(median(dr1,2),'%.2d'))},'FontSize',12,'Location','SouthWest')
    end
    
    elseif imodel==1
        if maen1
            legend({ strcat('Distribution of Äw/w '),strcat('<Äw/w> :',num2str(mean(dr1),'%.2d'))},'FontSize',12,'Location','NorthWest')%,'fit for [0-20 sec]'})
        
        else
    legend({ strcat('Distribution of Äw/w '),strcat('<Äw/w> : ',num2str(median(dr1),' %.2d'))},'FontSize',12,'Location','NorthWest')%,'fit for [0-20 sec]'})
        end
end
%xlim([1.e-5,2*dr_half(end2)])
%ylim([min(relative),1.e2])
else
    xlabel('Äw [eV]')
ylabel('P(Äw)')
legend({ strcat('Distribution of Äw'),strcat(' Power-law fit, k:',num2str(a,2)),strcat('Mean of Äw : ', num2str(mean(dr1),' %.2d'))},'FontSize',12,'Location','SouthWest')%,'fit for [0-20 sec]'})
end
%xlim([1.e-8,2.e0])
%ylim([9.e-6,1.e7])
%{
YTick = [10^(-5) 10^(-3) 10^(-1) 10^(1)  10^(3) 10^(5)  10^(7)];
set(gca,'ytick',YTick)

XTick = [ 10^(-8) 10^(-6) 10^(-4) 10^(-2) 10^(-0)];% 10^(-4) 10^(-2) 10^(0)];
set(gca,'xtick',XTick)
%}

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
box on

 
hold off

if mono==1
if imodel == 1

    temp=['new_png_rec\f18.png']; saveas(gca,temp);
    temp1=['new_eps_rec\f18.eps']; saveas(gca,temp1); 
    temp2=['new_fig_rec\f18.fig']; saveas(gca,temp2); 
    
 elseif imodel == 12
          temp=['new_png\f18 PdW=',num2str(PdW),'.png']; saveas(gca,temp);
    temp1=['new_eps\f18 PdW=',num2str(PdW),'.eps']; saveas(gca,temp1); 
    temp2=['new_fig\f18 PdW=',num2str(PdW),'.fig']; saveas(gca,temp2); 
 elseif imodel == 2 
     
      temp=['new_png_rec\f18 P=',num2str(P_reconnect),'.png']; saveas(gca,temp);
    temp1=['new_eps_rec\f18 P=',num2str(P_reconnect),'.eps']; saveas(gca,temp1); 
    temp2=['new_fig_rec\f18 P=',num2str(P_reconnect),'.fig']; saveas(gca,temp2); 
end
 
elseif mono==2
    if imodel == 1

    temp=['new_png_rec\f19.png']; saveas(gca,temp);
    temp1=['new_eps_rec\f19.eps']; saveas(gca,temp1); 
    temp2=['new_fig_rec\f19.fig']; saveas(gca,temp2); 
    
 elseif imodel == 12
          temp=['new_png\f19 PdW=',num2str(PdW),'.png']; saveas(gca,temp);
    temp1=['new_eps\f19 PdW=',num2str(PdW),'.eps']; saveas(gca,temp1); 
    temp2=['new_fig\f19 PdW=',num2str(PdW),'.fig']; saveas(gca,temp2); 
 elseif imodel == 2 
     
      temp=['new_png_rec\f19 P=',num2str(P_reconnect),'.png']; saveas(gca,temp);
    temp1=['new_eps_rec\f19 P=',num2str(P_reconnect),'.eps']; saveas(gca,temp1); 
    temp2=['new_fig_rec\f19 P=',num2str(P_reconnect),'.fig']; saveas(gca,temp2); 
end
    
else
    
    if imodel == 1

    temp=['new_png_rec\f19b.png']; saveas(gca,temp);
    temp1=['new_eps_rec\f19b.eps']; saveas(gca,temp1); 
    temp2=['new_fig_rec\f19b.fig']; saveas(gca,temp2); 
    
 elseif imodel == 12
          temp=['new_png\f19b PdW=',num2str(PdW),'.png']; saveas(gca,temp);
    temp1=['new_eps\f19b PdW=',num2str(PdW),'.eps']; saveas(gca,temp1); 
    temp2=['new_fig\f19b PdW=',num2str(PdW),'.fig']; saveas(gca,temp2); 
 elseif imodel == 2 
     
      temp=['new_png_rec\f19b P=',num2str(P_reconnect),'.png']; saveas(gca,temp);
    temp1=['new_eps_rec\f19b P=',num2str(P_reconnect),'.eps']; saveas(gca,temp1); 
    temp2=['new_fig_rec\f19b P=',num2str(P_reconnect),'.fig']; saveas(gca,temp2); 
end
    
    
    
end
    
    