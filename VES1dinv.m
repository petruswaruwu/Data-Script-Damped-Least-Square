%VES1dinv
%1D Inversion of Schlumberger Sounding Data
%By Yunus Levent Ekinci and Alper Demirci
close all;
clear all;
clc;
format long;
%load atau input data observasi 
pp = load ('C:\Users\user\Documents\Tugas Akhir\Matlab\TQ.dat'),1; 
x=pp(:,1);           %AB/2
roa=pp(:,2);         %Resistivitas dobs
r=[100 70 50 20];    %Model nilai resistivity awal sembarang
t=[1 3 10];          %Model ketebalan awal sembarang
m=[r t];             %Matriks m 
rinitial=r;          %Resistivitas inisial
tinitial=t;          %Ketebalan inisial
lr=length(r);        %Jumlah data resistivitas
lt=length(t);        %Jumlah data ketebalan
kr=1e-3;             %Faktor damping
iteration=1;         %Jumlah iterasi
maxiteration=10;     %Jumlah max iterasi
dfit=1;              %Inisial perbandingan misfit

%Proses Iterasi 
while iteration<maxiteration    %Ketentuan
    r=m(1:lr);                  %Menguraikan data r
    t=m(1+lr:lr+lt);            %Menguraikan data t
    for i=1:length(x)           %Komposisi i
        s=pp(i);                %s menjadi data load komposisi i
        [g]=VES1dmod(r,t,s);    %Pemodelan kedepan
        roa1(i,:)=g;            %roa1 = Resitivitas Kalkulasi Pertama
    end

    
%Fungsi Objektif1    
%Perhitungan Misfit 1    
    e1=log(roa)-log(roa1);      %e1 = Perbandingan resistivtas dobs degan dcal1
    dd=e1;                      %dd = hasil perbandingan
    misfit1=e1'*e1;             %Missfit pertama

    
%Perbdaninagn Misfit 1 dengan Faktor Damping 
    if misfit1<kr                       %ketentuan
        loglog(x,roa,'r.-',x,roa1,'k'); %Menampilkan grafik log pada kedua sumbu 
        axis([1 250 0 100])             %Batas Sb-X dan Sb-Y
        xlabel('AB/2(m)');              %Label Sb-X
        ylabel('Apparent Resistivity(Ohm-m)'); %Label Sb-Y
        h=legend('obs','clc');                 %Legenda
       break
    end
  
    
%Proses Modifikasi Parameter Model
%Metode Levenberg-Marquardt (Damped Least Sqaure)
%Dekomposisi matriks (SVD)
    [A]=jacobian(pp,x,r,t,lr,lt,roa,roa1);  %[A] = Matriks Jacobian
    [U S V]=svd(A,0);                       %Persamaan SVD di Matriks Jacobi
    ss=length(S);                           %ss = Panjang Matriks Diagonal S
    say=1;                                  
    k=0;
    while say<ss            
        diagS=diag(S);
        beta=S(say)*(dfit^(1/say));
        if beta<10^-5
            beta=0.001*say;
        end
        for i4=1:ss
            SS(i4,i4)=S(i4,i4)/(S(i4,i4)^2+beta);
        end
%U = Matriks data dan parameter
%V = Matriks parameter dan parameter
%S = Matriks diagonal parameter dan parameter
%ss = Jumlah data S


%Solusi persamaan Damped Least Square (DLS) dengan SVD
        dmg=V*SS*U'*dd; %Persamaan solusi DLS
        mg=exp(log(m)+dmg');
        r=mg(1:lr);
        t=mg(1+lr:lr+lt);
        for i5=1:length(x)
            s=pp(i5);
            [g]=VES1dmod(r,t,s);
            roa4(i5,:)=g;
        end
%Keterangan:
%dmg = delta resistivitas
%V = Vdiag
%SS = Matriks Diagonal 
%U' = Matrik U transpose
%dd = delta data
%mg = Resistivitas Kalkulasi
%s = data load
%roa4 = Hasil Proses DLS dalam SVD

        
%Fungsi Objektif2
%Perbandingann misfit 1 dan misfit 2
        e2=log(roa)-log(roa4); %Perbandingan Resistivitas dobs dengan dcal2
        misfit2=e2'*e2;        %Misfit kedua
        if misfit2>misfit1     %ketentuan
            ('Beta Control')   
            say=say+1;          
            k=k+1;
            if k==s-1
                iteration=maxiteration;
                say=ss+1;
            end

 

%Mencari Nilai Misfit Yang Diperoleh
        else
            say=ss+1;                       %jika Jumlah data Matriks S+1
            m=mg;                           %m = resistivitas hasil kalkulasi
            dfit=(misfit1-misfit2)/misfit1; %dfit = delta Xr (misfit total)
            iteration=iteration+1;
            a=iteration;
            if dfit<kr                      %ketentuan
             iteration=maxiteration;
             say=say+1;
            end
        end
    end
%Proses Iterasi Selesai   
 

%Plot Kurva Hubungan Antara Data Observasi Dengan Data Kalkulasi
    subplot(1,2,1)                          %Plot Gambar pada bagian kiri
    loglog(x,roa,'r.-.',x,roa4,'k');        %Menampilkan grafik log kedua sumbu
    axis([1 150 0 100])                     %Batas Sb-X dan Sb-Y
    xlabel('AB/2(m)');                      %label X
    ylabel('Apparent Resistivity(Ohm-m)');  %Label Y
    h=legend('observed','calculated');      %Legenda
    pause(0.001)                            %Pembentukan Kurva Terdelay 
end


%Pemodelan 1-D di sisi bagian kanan
observed=roa;    %data Observasi
calculated=roa4; %data kalkulasi
format bank;     %menampilkan nilai matriks dua angka desimal
m;               %resistivitas kalkulasi       
rr=[0,r];        %matriks resistivitas kalkulasi
tt=[0,cumsum(t),max(t)*10]; %Ketebalan Kalkulasi
subplot(1,2,2),             %Plot grafik pada gambar kanan
stairs(rr,tt,'r-');         %Plot data obs dan data cal dalam grafik tangga
rrr=[0,rinitial];           %Resistivitas awal
ttt=[0,cumsum(tinitial),max(t)*10]; %Ketebalan awal
hold on;
subplot(1,2,2),               %Plot grafik pada gambar kanan
stairs(rrr,ttt,'k--');        %Plot grafik tangga
set(gca,'Ydir','normal');     %Mengatur Nilai Sb-Y nilai terkecil dari bawah
set(gca,'Xscale','log');      %Nilai X menjadi persamaan log
ylim([0 sum(tt)-max(tt)]);    %Batas Sb-Y
xlim([1 2*max(rr)]);          %Batas Sb-X 
xlabel('Resistivity(Ohm.m)'); %Label Sb-X
ylabel('Depth(m)');           %Label Sb-Y
rms=norm(roa4-roa)/sqrt(length(roa)) %Persamaan RMS Error
title(['Iteration',num2str(a),'rms',num2str(rms), '%']) %Judul
h=legend('calculated model','starting model');          %Legenda
%Pemodelan 1-D Gambar 1 selesai


%Gabungan Antara Pemodelan 1-D dan Persebaran Data obs serta Kurva Kalkulasi 
figure(2)
loglog(x,roa,'r.',x,roa4,'k');  %Menampilkan grafik log pada kedua sumbu 
ylim([0 100]);                  %Batas Sb-Y
xlim([0.5 250]);                %batas Sb-X
xlabel('AB/2 (m)/ kedalaman (m)');  %Judul X
ylabel('True Resistivity (Ohm-m)'); %Judul Y
hold on
disp('Nilai Kedalaman (tt) kecuali 0.1 dan angka terakhir =');   %Tambahan
tt(1)=0.1
disp('Ketebalan =');
H=[cummax(t)]                 %H = Ketebalan Kalkulasi                                        
disp('Nilai True Resistivity (Ohm.m) kecuali angka terakhir =');
rn=[r,r(length(r))]           %rn = Resistivitas Kalkulasi
stairs(tt,rn,'-b');           %Plot Grafik Tangga
    set(gca,'Ydir','normal'); %Grafik Sb-Y, nilai kecil dari bawah
    set(gca,'Xscale','log');  %Data Sb-X membentuk kurva log
    set(gca,'Yscale','log');  %Data Sb-Y membentuk kurva log
title(['Inverse Damped Least Square'])
legend('observed data','calculated data','final model')
hold off
%Selesai

