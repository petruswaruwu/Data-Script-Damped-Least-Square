%Function VES1dmod
function[g]=VES1dmod(r,t,s)
q=13;       %Panjang nilai koefisien filter 
f=10;       %Besar alfa pada koefisien filter
m=4.438;    %Frekuensi Sampling
x=0;
e=exp(0.5*log(10)/m);   %Nilai eksponensial
h=2*q-2;                %koefisien filter linear terjauh
u=s*exp(-f*log(10)/m-x);%Nilai koefisien filter linear 
l=length(r);            %Jumlah data resistivitas
n=1;                    %koefisien linear awal


% Fungsi Transformasi Resistivitas
for i=1:n+h
    w=l;
    v=r(l);
    while w>1
        w=w-1;
        aa=tanh(t(w)/u);
        v=(v+r(w)*aa)/(1+v*aa/r(w));%Persamaan Transformasi
    end
    a(i)=v;
    u=u*e;
end
i=1;
%Keterangan:
%v = transformasi resistivitas
%r(w) = Resistivitas setiap lapisan
%aa = Nilai fungsi tanh


%Metode Filter Linear
g=105*a(i)-262*a(i+2)+416*a(i+4)-746*a(i+6)+1065*a(i+8);
g=g-4390*a(i+10)+13396*a(i+12)-27841*a(i+14);
g=g+16448*a(i+16)+8183*a(i+18)+2525*a(i+20);
g=(g+336*a(i+22)+225*a(i+24))/10000;
return
%g = Resistivitas Semu
%i = lapisan ke-k
