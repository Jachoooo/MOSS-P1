clear all; close all; clc; format long;

out=LTspice2Matlab('projekt4.raw'); 
freq_vect_spice=out.freq_vect(:); 
przebieg_spice=out.variable_mat(:);
przebieg_spice_db=20*log10(abs(przebieg_spice));
l=length(freq_vect_spice);
%%
is=1e-15;
bf=100;
br=10;
nf=1;
nr=1;
%rb=50;
tf=0.1e-9;
tr=10e-9;
cjc=2e-12;
cje=0.2e-12;
vjc=0.5;
vje=0.6;
mjc=0.5;
mje=0.5;
ut=0.026;

%punkt pracy wyliczony w LTSpice
ube1=3-2.50752;
ube2=2.50752-1.89566;
ubc1=3-12;
ubc2=2.50752-12;

[JB1,GBE1,GBC1,JC1,GCE1,GCC1,CBE1,CBC1]=BJT(ube1,ubc1,is,nf,nr,ut,bf,br,tf,tr,cje,cjc,vje,mje,1);
[JB2,GBE2,GBC2,JC2,GCE2,GCC2,CBE2,CBC2]=BJT(ube2,ubc2,is,nf,nr,ut,bf,br,tf,tr,cje,cjc,vje,mje,1);

rb=1e3;
RE=100e3;

%analiza AC
fmin = 10;
fmax = 1e9;
npoints = 100;


prompt = 'Jaka skala wykresow? 1 - log, 2 - lin ';
inp=input(prompt);

if inp==1
    skala = 'log';
    deltaf = 10^(1 / npoints);
else
    if inp==2
    	skala = 'lin';
        deltaf = (fmax - fmin) / (npoints - 1);
    else 
         msg = 'Zly numer skali';
         error(msg)
    end
end


f = fmin;
i = 1;
%%
while(1)
    
    omega = 2*pi*f;
    YBE1= 1i*omega*CBE1;
    YBC1= 1i*omega*CBC1;
    YBE2= 1i*omega*CBE2;
    YBC2= 1i*omega*CBC2;
    
    Y=zeros(9,9);
    
    Y(1,8)=1;
    Y(1,6)=-1;
    
    Y(2,2)=YBC1+YBC2-GCC1-GCC2;
    Y(2,3)=-YBC2-GCE1+GCE2+GCC2;
    Y(2,4)=-GCE2;
    Y(2,5)=-YBC1+GCE1+GCC1;
    Y(2,7)=1;
    
    Y(3,2)=-YBC2+GBC1+GCC1-GBC2;
    Y(3,3)=YBE1+YBE2+YBC2+GBE1+GCE1+GBE2+GBC2;
    Y(3,4)=-YBE2-GBE2;
    Y(3,5)=-YBE1-GBE1-GCC1-GCE1-GBC1;
    
    Y(4,2)=GBC2+GCC2;
    Y(4,3)=-YBE2-GBE2-GBC2-GCE2-GCC2;
    Y(4,4)=YBE2+GBE2+GCE2;
    Y(4,9)=1;
    
    Y(5,2)=-YBC1-GBC1;
    Y(5,3)=-YBE1-GBE1;
    Y(5,5)=YBE1+YBC1+GBE1+GBC1;
    Y(5,8)=-1;
    
    Y(6,1)=1;
    
    Y(7,2)=1;
    
    Y(8,1)=1;
    Y(8,5)=-1;
    Y(8,8)=-rb;
    
    Y(9,4)=1;
    Y(9,9)=-RE;
    
     B = [0; 0; 0; 0; 0; 1; 0; 0; 0];
     
     X=Y\B;
     syg = X(4);
     ff(i) = f;
  
     ch_amp(i) = abs(syg);
     ch_faz(i) = (angle(syg) * 180 / pi); 
     i=i+1;

     if (skala == 'lin') 
         f=f+deltaf; 
     else
         f=f*deltaf;
     end
     
     if (f > fmax) % Koniec analizy
      if(skala == 'lin')
          figure(1); 
          plot(ff, mag2db(ch_amp), 'r-'); %matlab
          hold on;
          plot(freq_vect_spice,przebieg_spice_db, 'g-'); %spice
          title('Charakterystyka amplitudowa');
          xlabel('Czestotliwosc[Hz]'); ylabel('Wzmocnienie [db]');
          legend('matlab','spice');
          grid on
          hold off;

          figure(2); 
          plot(ff, ch_faz, 'r-');
          hold on
          plot(freq_vect_spice, angle(przebieg_spice)*180/pi, 'g-');
          title('Charakterystyka fazowa');
          xlabel('Czestotliwosc [Hz]'); ylabel('Faza [stopnie]');
          legend('matlab','spice');
          grid on
          hold off
          
      blad_wzm=abs(przebieg_spice_db'-mag2db(ch_amp));
      blad_fazy=abs((angle(przebieg_spice) * 180 / pi)'-ch_faz); 
      blad_wzm_wzg=blad_wzm./przebieg_spice_db';
      blad_fazy_wzg=blad_fazy./((angle(przebieg_spice) * 180 / pi)');
      
      figure()      
      subplot(2,2,1)
      plot(freq_vect_spice, blad_wzm);
      grid on
      title('Blad ch-ki amplitudowej');
      xlabel('Czestotliwosc [Hz]'); ylabel('Blad'); 
      subplot(2,2,2)
      plot(freq_vect_spice, blad_fazy);
      grid on
    
      title('Blad ch-ki fazowej');
      xlabel('Czestotliwosc [Hz]'); ylabel('Blad'); 

      subplot(2,2,3)
      plot(freq_vect_spice, blad_wzm_wzg);
      grid on
      title('Blad wzgledny ch-ki amplitudowej');
      xlabel('Czestotliwosc [Hz]'); ylabel('Blad'); 
      subplot(2,2,4)
      plot(freq_vect_spice, blad_fazy_wzg);
      grid on
      title('Blad wzgledny ch-ki fazowej');
      xlabel('Czestotliwosc [Hz]'); ylabel('Blad'); 
      
      else
          figure(1); 
          semilogx(ff, mag2db(ch_amp), 'r-');
          hold on;
          semilogx(freq_vect_spice, przebieg_spice_db, 'g-');
          legend('matlab','spice');
          xlabel('Czestotliwosc[Hz]'); ylabel('Wzmocnienie [db]');
          title('Charakterystyka amplitudowa');
          grid on
          hold off;

          figure(2); 
          semilogx(ff, ch_faz, 'r-');
          hold on;
          semilogx(freq_vect_spice, angle(przebieg_spice) * 180 / pi, 'g-');
          legend('matlab','spice');
          xlabel('Czestotliwosc [Hz]'); ylabel('Faza [stopnie]');
          title('Charakterystyka fazowa');
          grid on
          hold off;
          
      blad_wzm=abs(przebieg_spice_db'-(20*log10(ch_amp)));
      blad_fazy=abs((angle(przebieg_spice) * 180 / pi)'-ch_faz); 
      blad_wzm_wzg=blad_wzm./abs(przebieg_spice_db');
      blad_fazy_wzg=blad_fazy./abs((angle(przebieg_spice) * 180 / pi)');
      
      figure()      
      subplot(2,2,1)
      semilogx(freq_vect_spice, blad_wzm);
      grid on
      title('Blad ch-ki amplitudowej');
      xlabel('Czestotliwosc [Hz]'); ylabel('Blad'); 
      subplot(2,2,2)
      semilogx(freq_vect_spice, blad_fazy);
      grid on
    
      title('Blad ch-ki fazowej');
      xlabel('Czestotliwosc [Hz]'); ylabel('Blad'); 

      subplot(2,2,3)
      semilogx(freq_vect_spice, blad_wzm_wzg);
      grid on
      title('Blad wzgledny ch-ki amplitudowej');
      xlabel('Czestotliwosc [Hz]'); ylabel('Blad'); 
      subplot(2,2,4)
      semilogx(freq_vect_spice, blad_fazy_wzg);
      grid on
      title('Blad wzgledny ch-ki fazowej');
      xlabel('Czestotliwosc [Hz]'); ylabel('Blad'); 
      end
      
      
      
     
    return;
    end
end
