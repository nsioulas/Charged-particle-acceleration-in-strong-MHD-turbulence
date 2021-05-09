 


lsize = 16; % Label fontsize
    nsize = 16; % Axis fontsize
    
    
    upper_step_limit=true; % If (F) dr can theoretically take values up to infinity 
%% Generating power-law distributed random numbers with the method of transformation

 dim=[nP,N-1];           % dimensions of dr matrix
                       
 P_law_index=1.2;      % P(dr)~dr^(P_law_index)
 
 Dr_min=10^2;          % Define minimum step size
 Dr_max=10^10;         % Define maximum step size
 
 
 if upper_step_limit 

   dr=((Dr_min^(1-P_law_index)) + (Dr_max^(1-P_law_index)-Dr_min^(1-P_law_index)).*rand(dim)).^(1/(1-P_law_index));
 else
   dr=Dr_min.*(1-rand(dim)).^(1/(1-P_law_index));
 end


%% Box size

L=10^(10);                        % length of the box size [cm]
L2=L/2;  
dmin=-L2;
dmax=L2;
 
    
    
    
    
    
    
    
start=2;    rd1start=2;      end2=30; %number of bins
ending=end2;     rd1end=end2;  %nP=length(rmsd);
tacc=2;


dr1=dr(:);
dr1_=dr(:);
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
         relative(1,jN)=p1(1,jN)/(total*(split2(jN)));
         %relative(1,jN)=p1(1,jN)/(total);
     else
    relative(1,jN)=p1(1,jN)/(total*(split2(jN)-split2(jN-1)));
    %relative(1,jN)=p1(1,jN)/(total);
    
     end
end


figure(1)
fitResults1 = polyfit(log10(dr_half(start:ending)),log10(relative(start:ending)),1);
pol=polyval(fitResults1,log10(dr_half(start:ending)));
a=fitResults1(1);
b=fitResults1(2);
polyfit_str = ['<(?r)^2> ~ ? ^{'  (num2str(a)) ' } '] % + ' num2str(b)]  polyfit_str will be : y = 4*x + 2


fit1=dr_half(start:ending).^a*10^(b);


loglog((dr_half(rd1start:rd1end)),(relative(rd1start:rd1end)),'blue-','LineWidth',2)
hold on
loglog((dr_half(start:(length(pol)+start-1))),fit1,'cyan--','LineWidth',2)
legend({ strcat('Distribution of Är '),strcat(' Power-law fit, k:',num2str(a))},'FontSize',12)%,'fit for [0-20 sec]'})
xlabel('Är  [cm]','FontSize',lsize)
ylabel('p(Är)','FontSize',lsize)

xlim([Dr_min,Dr_max])
ylim([1.e-13,10^(-2)])
XTick = [ 10^(2) 10^(4)  10^(6) 10^(8) 10^(10)];
set(gca,'xtick',XTick)

YTick = [ 10^(-14) 10^(-12) 10^(-10) 10^(-8) 10^(-6) 10^(-4) 10^(-2) 10^(0)];
set(gca,'ytick',YTick)


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




if upper_step_limit
      if imodel == 1
    temp=['png_files_upper\f1a.png']; saveas(gca,temp);
    temp1=['eps_files_upper\f1a.eps']; saveas(gca,temp1); 
    temp2=['fig_files_upper\f1a.fig']; saveas(gca,temp2); 
      else
          
          temp=['png_files_upper\dr_hist_new_noEnerg,Lmin=',num2str(Dr_min),'MLim=',num2str(Dr_max),'nP=',num2str(nP),'.png']; saveas(gca,temp);
    temp1=['eps_files_upper\dr_hist_new_noEnerg,Lmin=',num2str(Dr_min),'MLim=',num2str(Dr_max),'nP=',num2str(nP),'.eps']; saveas(gca,temp1); 
    temp2=['fig_files_upper\dr_hist_new_noEnerg,Lmin=',num2str(Dr_min),'MLim=',num2str(Dr_max),'nP=',num2str(nP),'.fig']; saveas(gca,temp2);
      
      end
    else
     if imodel == 1
 temp3=['png_files\dr_hist_new,Lmin=',num2str(Dr_min),'nP=',num2str(nP),'.png']; saveas(gca,temp3);
 
    temp4=['eps_files\dr_hist_new,Lmin=',num2str(Dr_min),'nP=',num2str(nP),'.eps']; saveas(gca,temp4); 
    temp5=['fig_files\dr_hist_new,Lmin=',num2str(Dr_min),'nP=',num2str(nP),'.fig']; saveas(gca,temp5); 
    
     else
         
         temp3=['png_files\dr_hist_new_noEnerg,Lmin=',num2str(Dr_min),'nP=',num2str(nP),'.png']; saveas(gca,temp3);
    temp4=['eps_files\dr_hist_new_noEnerg,Lmin=',num2str(Dr_min),'nP=',num2str(nP),'.eps']; saveas(gca,temp4); 
    temp5=['fig_files\dr_hist_new_noEnerg,Lmin=',num2str(Dr_min),'nP=',num2str(nP),'.fig']; saveas(gca,temp5); 
         
         
     end
end
