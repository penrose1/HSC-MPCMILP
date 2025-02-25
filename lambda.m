Pv=[1;2;3;4;5];
Pw=[11;12;13;14;15];
lambda = ones(10,1);
Aug_pv=zeros(length(Pv));
Aug_pw=zeros(length(Pw));
for i=1:length(Pv)
    Aug_pv(i,i) = Pv(i);
end
for i=1:length(Pw)
    Aug_pw(i,i) = Pw(i);
end

Aug_P = [Aug_pv Aug_pw]
res = Aug_P*lambda