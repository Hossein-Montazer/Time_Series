%% ARMAX 
%% LOAD MATERIAL 89

clc
clear

load fjarrvarme89.dat

ObsNbr = fjarrvarme89(:,1);
P = fjarrvarme89(:,2);
AAT = fjarrvarme89(:,3);
SWT = fjarrvarme89(:,4);
Y = fjarrvarme89(:,5);
M = fjarrvarme89(:,6);
D = fjarrvarme89(:,7);
H = fjarrvarme89(:,8);

subplot(211);
plot(P)
subplot(212);
plot(AAT);

%% LOAD MATERIAL 90

clc
clear

load fjarrvarme90.dat

ObsNbr = fjarrvarme90(:,1);
P = fjarrvarme90(:,2);
AAT = fjarrvarme90(:,3);
SWT = fjarrvarme90(:,4);
Y = fjarrvarme90(:,5);
M = fjarrvarme90(:,6);
D = fjarrvarme90(:,7);
H = fjarrvarme90(:,8);

subplot(211);
plot(P)
subplot(212);
plot(AAT);

%% REMOVE OUTLIERS (Not very effective method)

subplot(221);
plot(P);
title('P');
subplot(222);
plot(AAT);
title('AAT');

percntiles = prctile(AAT,[0.5 99.5]); 
outlierIndex = AAT < percntiles(1) | AAT > percntiles(2);
AAT_o = AAT;
AAT_o(outlierIndex) = [];

percntiles = prctile(P,[0.5 99.5]); 
outlierIndex = P < percntiles(1) | P > percntiles(2);
P_o = P;
P_o(outlierIndex) = [];

subplot(223);
plot(P_o);
title('P - Outliers');
subplot(224);
plot(AAT_o);
title('AAT - Outlier');

%% TRANSFORM

nbr = 0;
delay = 0;
AAT_s = AAT(1680*nbr + 1 - delay :(nbr + 1)*1680 - delay);
%bcNormPlot(AAT_s)
subplot(221);
plot(AAT_s);
title('AAT');
AAT_s_t = log(AAT_s);
ma = mean(AAT_s_t);
AAT_s_t_m = AAT_s_t-ma;
subplot(222);
plot(AAT_s_t_m);
title('AAT transform');

P_s = P(1680*nbr + 1 - delay:(nbr + 1)*1680 - delay);
P_s([290 291 292 293 478 479 480 481 482 781 782 783 784 785 941 942 943 1264 1265 1266]) = NaN;
%P_s([290 291 292 293 478 479 480 481 482 ]) = NaN;
P_s = misdata(iddata(P_s)); 
P_s = P_s.OutputData;
%bcNormPlot(P_s)
subplot(223);
plot(P_s);
title('P');
P_s_t = log(P_s);
mp = mean(P_s_t);
P_s_t_m = P_s_t-mp;
subplot(224);

plot(P_s_t_m);
title('P transform');

%% ACF/PACF/NORM P 

ehat = P_s_t_m;
figure; 

%subplot(1,3,1);
acf(ehat,round(length(ehat)/8), 2/sqrt(length(ehat)), 1, 0 ,0);
%title('MA-C');
%subplot(1,3,2);
%pacf(ehat,round(length(ehat)/8),2/sqrt(length(ehat)),  1, 0);
%title('AR-A');
%subplot(1,3,3);
%normplot(ehat);
%figure;
%plot(tacf(ehat,round(length(ehat)/8), 2/sqrt(length(ehat)), 1, 0 ,0));
%% ACF/PACF/NORM AAT

ehat = AAT_s_t_m;
figure; 

subplot(3,1,1);
acf(ehat,round(length(ehat)/4), 2/sqrt(length(ehat)), 1, 0 ,0);
title('MA-C');
subplot(3,1,2);
pacf(ehat,round(length(ehat)/4),2/sqrt(length(ehat)),  1, 0);
title('AR-A');
subplot(3,1,3);
plotNTdist(ehat);
figure;
plot(tacf(ehat,round(length(ehat)/8), 2/sqrt(length(ehat)), 1, 0 ,0));
%% DETREND

subplot(211);
plot(P_s_t_m);
[P_s_t_m_dt,DT] = deSeason(P_s_t_m,[24]);
subplot(212);
plot(P_s_t_m_dt);

ehat = P_s_t_m_dt;

figure; 

subplot(3,1,1);
acf(ehat,round(length(ehat)/4), 2/sqrt(length(ehat)), 1, 0 ,0);
title('MA-C');
subplot(3,1,2);
pacf(ehat,round(length(ehat)/4),2/sqrt(length(ehat)),  1, 0);
title('AR-A');
subplot(3,1,3);
plotNTdist(ehat);

%% DESEASON 

subplot(211);
plot(P_s_t_m_dt);
[P_s_t_m_dt_ds,DS] = deSeason(P_s_t_m_dt,[24]);
subplot(212);
plot(P_s_t_m_dt_ds);

ehat = P_s_t_m_dt_ds;

figure; 

subplot(3,1,1);
acf(ehat,round(length(ehat)/4), 2/sqrt(length(ehat)), 1, 0 ,0);
title('MA-C');
subplot(3,1,2);
pacf(ehat,round(length(ehat)/4),2/sqrt(length(ehat)),  1, 0);
title('AR-A');
subplot(3,1,3);
plotNTdist(ehat);

%% DESEASON P NO TREND

subplot(121);
plot(P_s_t_m);
[P_s_t_m_ds,DS] = deSeason(P_s_t_m,[24]);
%DS = [1 zeros(1, 23), -1 -1];
%P_s_t_m_ds = filter([1 zeros(1, 23), -1 -1], 1,P_s_t_m);
subplot(122);
plot(P_s_t_m_ds);

ehat = P_s_t_m_ds;

figure; 

subplot(1,3,1);
acf(ehat,round(length(ehat)/4), 2/sqrt(length(ehat)), 1, 0 ,0);
title('MA-C');
subplot(1,3,2);
pacf(ehat,round(length(ehat)/4),2/sqrt(length(ehat)),  1, 0);
title('AR-A');
subplot(1,3,3);
plotNTdist(ehat);



%% DESEASON AAT NO TREND

subplot(211);
plot(AAT_s_t_m);
[AAT_s_t_m_ds,DS] = deSeason(AAT_s_t_m,[24]);
subplot(212);
plot(AAT_s_t_m_ds);

%% CROSSCORR P AAT

M = 100;
stem(-M:M, crosscorr(AAT_s_t_m,P_s_t_m,M));
title('Cross correlation function');
xlabel('Lag');
hold on
plot(-M:M, 2/sqrt(length(AAT)) * ones(1,2*M+1),'--');
plot(-M:M, -2/sqrt(length(AAT)) * ones(1,2*M+1),'--');
hold off;

%% DESCIDE MODEL
[P_s_t_m_ds,DS] = deSeason(P_s_t_m,[24]);
%[AAT_s_t_m_ds,DS] = deSeason(AAT_s_t_m,[24]);
%[SWT_s_t_m_ds,DS] = deSeason(SWT_s_t_m,[24]);

I = P_s_t_m_ds;
%U = AAT_s_t_m_ds;
U = AAT_s_t_m;
Z = iddata(I,U);

delay = 11;

B = [zeros(1,delay) 1 1];
A = [1 1 1 1 0 0 0 1 1];
C = [1 zeros(1,23) 0 0];

model = idpoly(A,B,C);
model.Structure.c.Free = [1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1];
model.Structure.b.Free = B;
model.Structure.a.Free = A;

thx = armax(Z,model);

present(thx)

ehat = resid(thx,iddata(I,U));
ehat = ehat.OutputData;

figure;

subplot(1,3,1);
acf(ehat,round(length(ehat)/4), 2/sqrt(length(ehat)), 1, 0 ,0);
title('MA-C');
subplot(1,3,2);
pacf(ehat,round(length(ehat)/4),2/sqrt(length(ehat)),  1, 0);
title('AR-A');
subplot(1,3,3);
plotNTdist(ehat);
figure;
stem(tacf(ehat,round(length(ehat)/8), 2/sqrt(length(ehat)), 1, 0 ,0));
hold on
plot( (2/sqrt(length(ehat))) * ones(1,(length(ehat)/8)),'--');
plot( -(2/sqrt(length(ehat))) * ones(1,(length(ehat)/8)),'--');
figure; 

whitenessTest(ehat);

thx.a = conv(thx.a,DS);
%thx.b = conv(thx.b,DS);


%% PREDICTION

C = thx.c;
A = thx.a;
B = thx.b;

K = [1];

for i = 1:length(K);
    
k = K(i);

[Cs,As] = equalLength(C,A);

[Fk, Gk] = deconv(conv([1 zeros(1,k-1)], Cs),As);

BF = conv(B,Fk);

[Cs, BF] = equalLength(Cs,BF);

[Fhat,Ghat] = deconv(conv([1 zeros(1,k-1)],BF),Cs);

% Test other data

load fjarrvarme89.dat;
load fjarrvarme90.dat;

P = fjarrvarme89(:,2);
AAT = fjarrvarme89(:,3);

nbr = 0;
delay = 0;

tP_s = P(1680*nbr + 1 + delay:(nbr + 1)*1680 + delay);
tP_s_log = log(tP_s);
tmp = mean(tP_s_log);
tP_test = tP_s_log-tmp;

tP_test([290 291 292 293 478 479 480 481 482 781 782 783 784 785 941 942 943 1264 1265 1266]) = NaN;
%tP_test([290 291 292 293 478 479 480 481 482 ]) = NaN;
tP_test = misdata(iddata(tP_test)); 
tP_test = tP_test.OutputData;

tAAT_s = AAT(1680*nbr + 1 + delay:(nbr + 1)*1680 + delay);
tAAT_s_log = log(tAAT_s);
tma = mean(tAAT_s_log);
tAAT_test = tAAT_s_log-tma;

yhat_k = filter(Ghat,Cs,tAAT_test) + filter(Fhat,1,tAAT_test) + filter(Gk,Cs,tP_test);

figure(1);
subplot(1,length(K),i);
yhat_k_p = exp(yhat_k) + tmp;
P_p = exp(tP_test) + tmp;
plot(yhat_k_p(K(i):end));
hold on 
plot(P_p);
legend('yhat' ,'y')
title(strcat('K = ',int2str(K(i))));
% EVALUATION

ehat = yhat_k_p-P_p;
ehat = ehat(length(Gk):end);

disp('Variance');
disp(var(ehat));

figure; 

subplot(1,3,1);
acf(ehat,round(length(ehat)/4), 2/sqrt(length(ehat)), 1, 0 ,0);
title('MA-C');
subplot(1,3,2);
pacf(ehat,round(length(ehat)/4),2/sqrt(length(ehat)),  1, 0);
title('AR-A');
subplot(1,3,3);
plotNTdist(ehat);

figure;

whitenessTest(ehat);
cup = mean(ehat) + 1.96 * sqrt(var(ehat));
clo = -(mean(ehat) + 1.96 * sqrt(var(ehat))); 

disp('Procent out');
Nbrout = sum(ehat>cup) + sum(ehat<clo);  
disp( Nbrout / length(ehat) );

end

%%

K = 72;

ttP_s = P(1680:1680 + K-1);
ttP_s_log = log(ttP_s);
ttmp = mean(ttP_s_log);
ttP_test =  ttP_s_log-ttmp;
yhat_k = forecast(thx,iddata(tP_test,tAAT_test),K);

yhat_k_p = exp(yhat_k.OutputData) + tmp;
P_p = exp(ttP_test) + ttmp;

plot(yhat_k_p);
hold on 
plot(P_p);
legend('yhat' , 'y')
var(yhat_k_p-P_p)

%%

%% HEMTENTA
% PROBLEM 2

clc
clear
size = 10;
a = 0.2;

x = randn(1,size);
y_def = filter(1,[1 a], x);

for i = 1:size-1
   
   predict1(i) = y_def(i)*(-a);
   predict2(i) = y_def(i);
   
   error1(i) = y_def(i + 1) - predict1(i);
   error2(i) = y_def(i + 1) - predict2(i);
      
end

vest = var(predict1)/var(predict2);
vteo = (a^2);

disp('Variance y_k')
disp(var(predict2));
disp('Variance -a*y_k')
disp(var(predict1));
disp('Estimated variance');
disp(vest);
disp('Theoretical variance');
disp(vteo);
