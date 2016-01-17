%% ARMAX 
%
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

subplot(131);
plot(P(1:1680))
title('P');
ylabel('Amplitude');
xlabel('Time');
subplot(132);
plot(AAT(1:1680));
title('AAT');
ylabel('Amplitude');
xlabel('Time');
subplot(133);
plot(SWT(1:1680));
title('SWT');
ylabel('Amplitude');
xlabel('Time');
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

subplot(311);
plot(P)
subplot(312);
plot(AAT);
subplot(313);
plot(SWT);

%% REMOVE OUTLIERS (Not very effective method)

subplot(321);
plot(P);
title('P');
subplot(322);
plot(AAT);
title('AAT');
subplot(323);
plot(SWT);
title('SWT');

%%
percntiles = prctile(AAT,[0.5 99.5]); 
outlierIndex = AAT < percntiles(1) | AAT > percntiles(2);
AAT(outlierIndex) = [];
%%
percntiles = prctile(P,[0.5 99.5]); 
outlierIndex = P < percntiles(1) | P > percntiles(2);
P(outlierIndex) = [];
%%
percntiles = prctile(SWT,[1 99]); 
outlierIndex = SWT < percntiles(1) | SWT > percntiles(2);
SWT(outlierIndex) = NaN;
SWT = misdata(iddata(SWT)); 
SWT = SWT.OutputData;
%%

subplot(324);
plot(P);
title('P - Outliers');
subplot(325);
plot(AAT);
title('AAT - Outlier');
subplot(326)
plot(SWT);
title('SWT - Outlier');

%% TRANSFORM

nbr = 0;
delay = 0;
AAT_s = AAT(1680*nbr + 1 + delay :(nbr + 1)*1680 + delay);
%AAT_s = AAT(delay*nbr + 1 : (nbr + 1)*delay);
%bcNormPlot(AAT_s)
subplot(321);
plot(AAT_s);
title('AAT');
AAT_s_t = log(AAT_s);
ma = mean(AAT_s_t);
AAT_s_t_m = AAT_s_t-ma;
subplot(322);
plot(AAT_s_t_m);
title('AAT transform');

P_s = P(1680*nbr + 1 + delay:(nbr + 1)*1680 + delay);
%P_s = P(delay*nbr + 1 : (nbr + 1)*delay);
%bcNormPlot(P_s)
subplot(323);
plot(P_s);
title('P');
P_s_t = log(P_s);
mp = mean(P_s_t);
P_s_t_m = P_s_t-mp;
subplot(324);
plot(P_s_t_m);
title('P transform');


SWT_s = SWT(1680*nbr + 1 + delay:(nbr + 1)*1680 + delay);
%SWT_s = SWT(delay*nbr + 1 : (nbr + 1)*delay);
%bcNormPlot(SWT_s)
subplot(325);
plot(SWT_s);
title('SWT');
SWT_s_t = log(SWT_s);
mp = mean(SWT_s_t);
SWT_s_t_m = SWT_s_t-mp;
subplot(326);
plot(SWT_s_t_m);
title('SWT transform');

%% REMOVE OUTLIERS (Effective method)

P_s_t_m([ 478 479 480 481 482]) = NaN;
P_s_t_m = misdata(iddata(P_s_t_m)); 
P_s_t_m = P_s_t_m.OutputData;

%% ACF/PACF/NORM P 

ehat = P_s_t_m;
figure; 

subplot(1,3,1);
acf(ehat,round(length(ehat)/4), 2/sqrt(length(ehat)), 1, 0 ,0);
title('ACF');
ylabel('Amplitude');
xlabel('Time');
subplot(1,3,2);
pacf(ehat,round(length(ehat)/4),2/sqrt(length(ehat)),  1, 0);
title('PACF');
ylabel('Amplitude');
xlabel('Time');
subplot(1,3,3);
plotNTdist(ehat);

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

%% ACF/PACF/NORM SWT

ehat = SWT_s_t_m;
figure; 

subplot(3,1,1);
acf(ehat,round(length(ehat)/4), 2/sqrt(length(ehat)), 1, 0 ,0);
title('MA-C');
subplot(3,1,2);
pacf(ehat,round(length(ehat)/4),2/sqrt(length(ehat)),  1, 0);
title('AR-A');
subplot(3,1,3);
plotNTdist(ehat);

%% DETREND

subplot(211);
plot(P_s_t_m);
[P_s_t_m_dt,DT] = deSeason(P_s_t_m,[24]);
subplot(212);
plot(P_s_t_m_dt);

%% DESEASON

subplot(211);
plot(P_s_t_m_dt);
[P_s_t_m_dt_ds,DS] = deSeason(P_s_t_m_dt,[24]);
subplot(212);
plot(P_s_t_m_dt_ds);

%% DESEASON P NO TREND

subplot(211);
plot(P_s_t_m);
[P_s_t_m_ds,DS] = deSeason(P_s_t_m,[24]);
subplot(212);
plot(P_s_t_m_ds);

%% DESEASON AAT NO TREND

subplot(211);
plot(AAT_s_t_m);
[AAT_s_t_m_ds,DS] = deSeason(AAT_s_t_m,[24]);
subplot(212);
plot(AAT_s_t_m_ds);

%% DESEASON SWT NO TREND

subplot(211);
plot(SWT_s_t_m);
[SWT_s_t_m_ds,DS] = deSeason(SWT_s_t_m,[24]);
subplot(212);
plot(SWT_s_t_m_ds);

%% CROSSCORR SWT P 

figure;

subplot(121);
M = 100;
stem(-M:M, crosscorr(SWT_s_t_m,P_s_t_m,M));
title('Cross correlation function SWT P');
xlabel('Lag');
ylabel('Amplitude');
hold on
plot(-M:M, 2/sqrt(length(AAT)) * ones(1,2*M+1),'--');
plot(-M:M, -2/sqrt(length(AAT)) * ones(1,2*M+1),'--');
hold off;

 %CROSSCORR AAT P 
 
subplot(122);

M = 100; 
stem(-M:M, crosscorr(AAT_s_t_m,P_s_t_m,M));
title('Cross correlation function AAT P');
xlabel('Lag');
ylabel('Amplitude');

hold on
plot(-M:M, 2/sqrt(length(AAT)) * ones(1,2*M+1),'--');
plot(-M:M, -2/sqrt(length(AAT)) * ones(1,2*M+1),'--');
hold off;


%% DESCIDE MODEL

subplot(121);
plot(P_s_t_m - SWT_s_t_m*1.5)
title('P-SWT');
subplot(122);
plot(P_s_t_m);
title('P');
%%

[P_s_t_m_ds,DS] = deSeason(P_s_t_m,[24]);
%P_s_t_m_ds = P_s_t_m_ds(24:end);
[AAT_s_t_m_ds,DS] = deSeason(AAT_s_t_m,[24]);
%[SWT_s_t_m_ds,DS] = deSeason(SWT_s_t_m,[24]);

U = [SWT_s_t_m AAT_s_t_m];
I = P_s_t_m_ds;
Z = iddata(I,U);

B = cell(1,2);
delay1 = 0;
delay2 = 11;
B{1,1} = [1 1 1 1];
B{1,2} = [zeros(1,delay2) 1];
A = [1 1 1 1 0 0 0 1 1];
C = [1 zeros(1,23) 0 0];

model = idpoly(A,B,C);
model.Structure.c.Free = [1 1 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1];
model.Structure.b(1).Free = B{1,1};
model.Structure.a.Free = A;

thx = pem(Z,model);

present(thx);

ehat = resid(thx,Z);
ehat = ehat.OutputData;

% Evalueate residual

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

whitenessTest(ehat);

thx.a = conv(thx.a,DS);

%% PREDICTION

C = thx.c;
A = thx.a;

B1 = cell2mat(thx.b(1));
B2 = cell2mat(thx.b(2));

K = [1];

for i = 1:length(K);
    
k = K(i);

[Cs,As] = equalLength(C,A);
[Fk, Gk] = deconv(conv([1 zeros(1,k-1)], Cs),As);

BF1 = conv(B1,Fk);
BF2 = conv(B2,Fk);

[Cs1, BF1] = equalLength(Cs,BF1);
[Cs2, BF2] = equalLength(Cs,BF2);

[Fhat1,Ghat1] = deconv(conv([1 zeros(1,k-1)],BF1),Cs);
[Fhat2,Ghat2] = deconv(conv([1 zeros(1,k-1)],BF2),Cs);

nbr = 0;
delay = 0;

load fjarrvarme90.dat;
load fjarrvarme89.dat;

P = fjarrvarme89(:,2);
AAT = fjarrvarme89(:,3);
SWT = fjarrvarme89(:,4);

tP_s = P(1680*nbr + 1 + delay:(nbr + 1)*1680 + delay);
%tP_s = P(delay*nbr + 1 : (nbr + 1)*delay);

tP_s_log = log(tP_s);
tmp = mean(tP_s_log);
tP_test = tP_s_log-tmp;

%tP_test([290 291 292 293 478 479 480 481 482]) = NaN;
%tP_test = misdata(iddata(tP_test)); 
%tP_test = tP_test.OutputData;

tAAT_s = AAT(1680*nbr + 1 + delay:(nbr + 1)*1680 + delay);
%tAAT_s = P(delay*nbr + 1 : (nbr + 1)*delay);
tAAT_s_log = log(tAAT_s);
tma = mean(tAAT_s_log);
tAAT_test = tAAT_s_log-tma;

tSWT_s = SWT(1680*nbr + 1 + delay:(nbr + 1)*1680 + delay);
%tSWT_s = P(delay*nbr + 1 : (nbr + 1)*delay);
tSWT_s_log = log(tSWT_s);
tms = mean(tSWT_s_log);
tSWT_test = tSWT_s_log-tms;

yhat_k = filter(Ghat1,C,tSWT_test) + filter(Fhat1,1,tSWT_test) + filter(Ghat2,C,tAAT_test)  + filter(Gk,C,tP_test) ;
%yhat_k = filter(Ghat1,Cs,SWT_s_t_m) + filter(Ghat2,Cs,AAT_s_t_m) + filter(Gk,Cs,P_s_t_m);

figure;

figure(1);
subplot(1,length(K),i);
yhat_k_p = exp(yhat_k) + tmp;
P_p = exp(tP_test) + tmp;
plot(yhat_k_p);
hold on 
plot(P_p);
legend('yhat' ,'y')
title(strcat('K = ',int2str(K(i))));

% EVALUATION

figure;

ehat = yhat_k_p-P_p;
ehat = ehat(length(Gk):end);

disp('Variance');
disp(var(ehat));

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
yhat_k = forecast(thx,iddata(tP_test,[tSWT_test,tAAT_test]),K);

yhat_k_p = exp(yhat_k.OutputData) + tmp;
P_p = exp(ttP_test) + ttmp;

plot(yhat_k_p);
hold on 
plot(P_p);
legend('yhat' , 'y')
var(yhat_k_p-P_p)