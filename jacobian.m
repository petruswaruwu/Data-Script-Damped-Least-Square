%Function Jacobian
function[A]=jacobian(pp,x,r,t,lr,lt,roa,roa1)
par=0.1;


%Turunan resistivitas terhadap (A1)
r2=r;
for i2=1:lr;
    r2(i2)=(r(i2)*par)+r(i2);
    for ii=1:length(x)
        s=pp(ii);
        [g]=VES1dmod(r2,t,s);
        roa2(ii,:)=g;
    end
    A1(:,i2)=[(roa2-roa1)/(r(i2)*par)]*r(i2)./roa;
    r2=r;
end
%Keterangan:
%roa1 = resistivitas dcal1
%roa2 = resistivitas dcal2


%Turunan resistivitas terhadap ketebalan (A2)
t2=t;
for i3=1:lt
    t2(i3)=(t(i3)*par)+t(i3);
    for ii=1:length(x)
        s=pp(ii);
        [g]=VES1dmod(r,t2,s);
        roa3(ii,:)=g;
    end
    A2(:,i3)=[(roa3-roa1)/(t(i3)*par)]*t(i3)./roa;
    t2=t;
end
%Keterangan:
%roa1 = resistivitas dcal1
%roa3 = resistivitas dcal3


A=[A1 A2];
return