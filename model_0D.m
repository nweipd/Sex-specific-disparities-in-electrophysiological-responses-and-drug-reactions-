function [output]=model_0D(t,X,flag_gender,flag_cell,flag_ode,flag_hormone,ko,fkatp,flag_drug)

constants_0D

N=1; 

%endo = 0, epi = 1, M=2

%extracellular ionic concentrations
nao=140.0;
cao=1.8;


SFNa=0.75;  % change based on ischemia 
SFCaL=0.75; % change based on ischemia

%%% the following is the start for each type 

celltype_vector=flag_cell;

%give names to the state vector values
v=X(1);
nai=X(2);
nass=X(3);
ki=X(4);
kss=X(5);
cai=X(6);
cass=X(7);
cansr=X(8);
cajsr=X(9);
m=X(10);
hf=X(11);
hs=X(12);
j=X(13);
hsp=X(14);
jp=X(15);
mL=X(16);
hL=X(17);
hLp=X(18);
a=X(19);
iF=X(20);
iS=X(21);
ap=X(22);
iFp=X(23);
iSp=X(24);
d=X(25);
ff=X(26);
fs=X(27);
fcaf=X(28);
fcas=X(29);
jca=X(30);
nca=X(31);
ffp=X(32);
fcafp=X(33);
xrf=X(34);
xrs=X(35);
xs1=X(36);
xs2=X(37);
xk1=X(38);
Jrelnp=X(39);
Jrelp=X(40);
CaMKt=X(41);


%%%%%
GKs_rt = 0;
GKr_rt = 0;
GK1_rt = 0;
Gto_rt = 0;
pCa_rt = 0;
NaK_rt = 0;
Gup_rt = 0;
CaM_rt = 0;
    

%%%%%%%%%%%%%%%%%%% drug effects %%%%%%%%%%%%%%%%%%%%%
dose=4; 
if flag_drug==1 % Bepridil CiPA
    Cmax=33;
    IC50_IKr=50;
    hill_IKr=0.9;
    IC50_INaL=1813.9;
    hill_INaL=1.4;
    IC50_ICaL=2808.1;
    hill_ICaL=0.6;
    IC50_INa=2929.3;
    hill_INa=1.2;
    IC50_Ito=8594;
    hill_Ito=3.5;
    IC50_IKs=28628.3;
    hill_IKs=0.7;
    GKrfc=drug_effects(dose,Cmax,IC50_IKr,hill_IKr);
    GNaLfc=drug_effects(dose,Cmax,IC50_INaL,hill_INaL);
    GCaLfc=drug_effects(dose,Cmax,IC50_ICaL,hill_ICaL);
    GNafc=drug_effects(dose,Cmax,IC50_INa,hill_INa);
    Gtofc=drug_effects(dose,Cmax,IC50_Ito,hill_Ito);
    GK1fc=1;
    GKsfc=drug_effects(dose,Cmax,IC50_IKs,hill_IKs);
elseif flag_drug==2 % Sotalol CiPA
    Cmax=14690;
    IC50_IKr=110600;
    hill_IKr=0.8;
    IC50_ICaL=7061527;
    hill_ICaL=0.9;
    IC50_INa=1140000000;
    hill_INa=0.5;
    IC50_Ito=43143455;
    hill_Ito=0.7;
    IC50_IK1=3050260;
    hill_IK1=1.2;
    IC50_IKs=4221856;
    hill_IKs=1.2;
    GKrfc=drug_effects(dose,Cmax,IC50_IKr,hill_IKr);
    GNaLfc=1;
    GCaLfc=drug_effects(dose,Cmax,IC50_ICaL,hill_ICaL);
    GNafc=drug_effects(dose,Cmax,IC50_INa,hill_INa);
    Gtofc=drug_effects(dose,Cmax,IC50_Ito,hill_Ito);
    GK1fc=drug_effects(dose,Cmax,IC50_IK1,hill_IK1);
    GKsfc=drug_effects(dose,Cmax,IC50_IKs,hill_IKs);
elseif flag_drug==3 % Diazepam
    Cmax=29;
    IC50_IKr=53200;
    hill_IKr=1;
    IC50_ICaL=30500;
    hill_ICaL=1;
    IC50_INa=306400;
    hill_INa=1;
    GKrfc=drug_effects(dose,Cmax,IC50_IKr,hill_IKr);
    GNaLfc=1;
    GCaLfc=drug_effects(dose,Cmax,IC50_ICaL,hill_ICaL);
    GNafc=drug_effects(dose,Cmax,IC50_INa,hill_INa);
    Gtofc=1;
    GK1fc=1;
    GKsfc=1;
elseif flag_drug==4 % Diltiazem CiPA
    Cmax=122;
    IC50_IKr=13150;
    hill_IKr=0.9;
    IC50_INaL=21868.5;
    hill_INaL=0.7;
    IC50_ICaL=112.1;
    hill_ICaL=0.7;
    IC50_INa=110859;
    hill_INa=0.7;
    IC50_Ito=2820000000;
    hill_Ito=0.2;
    GKrfc=drug_effects(dose,Cmax,IC50_IKr,hill_IKr);
    GNaLfc=drug_effects(dose,Cmax,IC50_INaL,hill_INaL);
    GCaLfc=drug_effects(dose,Cmax,IC50_ICaL,hill_ICaL);
    GNafc=drug_effects(dose,Cmax,IC50_INa,hill_INa);
    Gtofc=drug_effects(dose,Cmax,IC50_Ito,hill_Ito);
    GK1fc=1;
    GKsfc=1;

elseif flag_drug==5 % Mibefradil1
    Cmax=12;
    IC50_IKr=1700;
    hill_IKr=1;
    IC50_ICaL=510;
    hill_ICaL=1;
    IC50_INa=5600;
    hill_INa=1;
    GKrfc=drug_effects(dose,Cmax,IC50_IKr,hill_IKr);
    GNaLfc=1;
    GCaLfc=drug_effects(dose,Cmax,IC50_ICaL,hill_ICaL);
    GNafc=drug_effects(dose,Cmax,IC50_INa,hill_INa);
    Gtofc=1;
    GK1fc=1;
    GKsfc=1;
elseif flag_drug==6 % Mibefradil2
    Cmax=12;
    IC50_IKr=1800;
    hill_IKr=1;
    IC50_ICaL=156;
    hill_ICaL=1;
    IC50_INa=980;
    hill_INa=1;
    GKrfc=drug_effects(dose,Cmax,IC50_IKr,hill_IKr);
    GNaLfc=1;
    GCaLfc=drug_effects(dose,Cmax,IC50_ICaL,hill_ICaL);
    GNafc=drug_effects(dose,Cmax,IC50_INa,hill_INa);
    Gtofc=1;
    GK1fc=1;
    GKsfc=1;
elseif flag_drug==7 % Nifedipine1
    Cmax=8;
    IC50_IKr=44000;
    hill_IKr=1;
    IC50_ICaL=12;
    hill_ICaL=1;
    IC50_INa=88500;
    hill_INa=1;
    GKrfc=drug_effects(dose,Cmax,IC50_IKr,hill_IKr);
    GNaLfc=1;
    GCaLfc=drug_effects(dose,Cmax,IC50_ICaL,hill_ICaL);
    GNafc=drug_effects(dose,Cmax,IC50_INa,hill_INa);
    Gtofc=1;
    GK1fc=1;
    GKsfc=1;
elseif flag_drug==8  % Nifedipine2
    Cmax=7.7;
    IC50_IKr=275000;
    hill_IKr=1;
    IC50_ICaL=60;
    hill_ICaL=1;
    IC50_INa=37000;
    hill_INa=1;
    GKrfc=drug_effects(dose,Cmax,IC50_IKr,hill_IKr);
    GNaLfc=1;
    GCaLfc=drug_effects(dose,Cmax,IC50_ICaL,hill_ICaL);
    GNafc=drug_effects(dose,Cmax,IC50_INa,hill_INa);
    Gtofc=1;
    GK1fc=1;
    GKsfc=1;
elseif flag_drug==9 % Nitrendipine1
    Cmax=3;
    IC50_IKr=24600;
    hill_IKr=1;
    IC50_ICaL=25;
    hill_ICaL=1;
    IC50_INa=21600;
    hill_INa=1;
    GKrfc=drug_effects(dose,Cmax,IC50_IKr,hill_IKr);
    GNaLfc=1;
    GCaLfc=drug_effects(dose,Cmax,IC50_ICaL,hill_ICaL);
    GNafc=drug_effects(dose,Cmax,IC50_INa,hill_INa);
    Gtofc=1;
    GK1fc=1;
    GKsfc=1;
elseif flag_drug==10 % Nitrendipine2
    Cmax=3.02;
    IC50_IKr=10000 ;
    hill_IKr=1;
    IC50_ICaL=0.35 ;
    hill_ICaL=1;
    IC50_INa=36000;
    hill_INa=1;
    GKrfc=drug_effects(dose,Cmax,IC50_IKr,hill_IKr);
    GNaLfc=1;
    GCaLfc=drug_effects(dose,Cmax,IC50_ICaL,hill_ICaL);
    GNafc=drug_effects(dose,Cmax,IC50_INa,hill_INa);
    Gtofc=1;
    GK1fc=1;
    GKsfc=1;
elseif flag_drug==11 % Prenylamine
    Cmax=17;
    IC50_IKr=65 ;
    hill_IKr=1;
    IC50_ICaL=1240 ;
    hill_ICaL=1;
    IC50_INa=2520;
    hill_INa=1;
    GKrfc=drug_effects(dose,Cmax,IC50_IKr,hill_IKr);
    GNaLfc=1;
    GCaLfc=drug_effects(dose,Cmax,IC50_ICaL,hill_ICaL);
    GNafc=drug_effects(dose,Cmax,IC50_INa,hill_INa);
    Gtofc=1;
    GK1fc=1;
    GKsfc=1;
elseif flag_drug==12 % Propranolol
    Cmax=26;
    IC50_IKr=2828 ;
    hill_IKr=1;
    IC50_ICaL=18000 ;
    hill_ICaL=1;
    IC50_INa=2100;
    hill_INa=1;
    GKrfc=drug_effects(dose,Cmax,IC50_IKr,hill_IKr);
    GNaLfc=1;
    GCaLfc=drug_effects(dose,Cmax,IC50_ICaL,hill_ICaL);
    GNafc=drug_effects(dose,Cmax,IC50_INa,hill_INa);
    Gtofc=1;
    GK1fc=1;
    GKsfc=1;
elseif flag_drug==13  % Verapamil CiPA
    Cmax=81;
    IC50_IKr=288;
    hill_IKr=1;
    IC50_INaL=7028;
    hill_INaL=1;
    IC50_ICaL=201.8;
    hill_ICaL=1.1;
    IC50_Ito=13429.2;
    hill_Ito=0.8;
    IC50_IK1=349000000;
    hill_IK1=0.3;
    GKrfc=drug_effects(dose,Cmax,IC50_IKr,hill_IKr);
    GNaLfc=drug_effects(dose,Cmax,IC50_INaL,hill_INaL);
    GCaLfc=drug_effects(dose,Cmax,IC50_ICaL,hill_ICaL);
    GNafc=1;
    Gtofc=drug_effects(dose,Cmax,IC50_Ito,hill_Ito);
    GK1fc=drug_effects(dose,Cmax,IC50_IK1,hill_IK1);
    GKsfc=1;
elseif flag_drug==14  % Ranolazine CiPA
    Cmax=1948.2;
    IC50_IKr=8270;
    hill_IKr=0.9;
    IC50_INaL=7884.5;
    hill_INaL=0.9;
    IC50_INa=68774;
    hill_INa=1.4;
    IC50_IKs=36155020;
    hill_IKs=0.5;
    GKrfc=drug_effects(dose,Cmax,IC50_IKr,hill_IKr);
    GNaLfc=drug_effects(dose,Cmax,IC50_INaL,hill_INaL);
    GCaLfc=1;
    GNafc=drug_effects(dose,Cmax,IC50_INa,hill_INa);
    Gtofc=1;
    GK1fc=1;
    GKsfc=drug_effects(dose,Cmax,IC50_IKs,hill_IKs);
elseif flag_drug==0
    GKrfc=1;
    GNaLfc=1;
    GCaLfc=1;
    GNafc=1;
    Gtofc=1;
    GK1fc=1;
    GKsfc=1;
end





%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CaMK constants
KmCaMK=0.15;

aCaMK=0.05;
bCaMK=0.00068;
CaMKo=0.05;
KmCaM=0.0015;
%update CaMK
CaMKb=CaMKo.*(1.0-CaMKt)./(1.0+KmCaM./cass);
CaMKa=CaMKb+CaMKt;
dCaMKt=aCaMK.*CaMKb.*(CaMKb+CaMKt)-bCaMK.*CaMKt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%reversal potentials
ENa=(R.*T./F).*log(nao./nai);
EK=(R.*T./F).*log(ko./ki);
PKNa=0.01833;
EKs=(R.*T./F).*log((ko+PKNa.*nao)./(ki+PKNa.*nai));
%convenient shorthand calculations
vffrt=v.*F.*F./(R.*T);
vfrt=v.*F./(R.*T);

%%%%%%%%revise based on Ele's code%%%%%%%%%%
if flag_hormone==0
    factor_Kr = 1.0;
	factor_Ks = 1.0;
	factor_CaL = 1.0;
elseif flag_hormone==1
    factor_Kr = 1.0;
	factor_Ks = 1.38;
	factor_CaL = 0.94;
elseif flag_hormone==2
    factor_Kr = 1.0;
	factor_Ks = 1.4;
	factor_CaL = 0.8;
elseif flag_hormone==3
    factor_Kr = 0.98;
	factor_Ks = 1.19;
	factor_CaL = 1.0;
elseif flag_hormone==4
    factor_Kr = 0.86;
	factor_Ks = 1.19;
	factor_CaL = 1.0;
elseif flag_hormone==5
    factor_Kr = 0.9;
	factor_Ks = 1.4;
	factor_CaL = 1.0;
end
%%%%%%




%calculate INa
mss=1.0./(1.0+exp((-(v+39.57))./9.871));
tm=1.0./(6.765.*exp((v+11.64)./34.77)+8.552.*exp(-(v+77.42)./5.955));
dm=(mss-m)./tm; 
hss=1.0/(1+exp((v+78.5)/6.22)); %%%%%% revised based on modified ORd 
thf=1.0/(3.6860e-6*exp(-(v+3.8875)/7.8579)+16*exp((v-0.4963)/9.1843)); %%%%%%% revised based on modified ORd 
ths=1.0./(0.009794.*exp(-(v+17.95)./28.05)+0.3343.*exp((v+5.730)./56.66));
Ahf=0.99;
Ahs=1.0-Ahf;
dhf=(hss-hf)./thf;
dhs=(hss-hs)./ths;
h=Ahf.*hf+Ahs.*hs;
jss=hss;
tj=4.8590+1.0/(0.8628*exp(-(v+116.7258)/7.6005)+1.1096*exp((v+6.2719)/9.0358)); %%%%%%%%% revised based on modified ORd 
dj=(jss-j)./tj; 
hssp=1.0/(1+exp((v+84.7)/6.22)); %%%%%% revised based on modified ORd
thsp=3.0.*ths;
dhsp=(hssp-hsp)./thsp; 
hp=Ahf.*hf+Ahs.*hsp;
tjp=1.46.*tj;
djp=(jss-jp)./tjp;
GNa=75; % mS/μF
GNa=GNa/GNafc; % added for drug 
fINap=(1.0./(1.0+KmCaMK./CaMKa));
INa=SFNa.*GNa.*(v-ENa).*m.^3.0.*((1.0-fINap).*h.*j+fINap.*hp.*jp);

%calculate INaL
mLss=1.0./(1.0+exp((-(v+42.85))./5.264));
tmL=tm;
dmL=(mLss-mL)./tmL; 
hLss=1.0./(1.0+exp((v+87.61)./7.488));
thL=200.0;
dhL=(hLss-hL)./thL;
hLssp=1.0./(1.0+exp((v+93.81)./7.488));
thLp=3.0.*thL;
dhLp=(hLssp-hLp)./thLp;

GNaL0=0.0075;
GNaL0=GNaL0/GNaLfc; % added for drug 
% change based on celltype
%%%%%%%%
GNaL=GNaL0*ones(N,1);
% GNaL(celltype_vector==1)=GNaL0.*0.6;
%%%%%%%%
fINaLp=(1.0./(1.0+KmCaMK./CaMKa));
INaL=GNaL.*(v-ENa).*mL.*((1.0-fINaLp).*hL+fINaLp.*hLp);


%calculate Ito
ass=1.0./(1.0+exp((-(v-14.34))./14.82));
ta=1.0515./(1.0./(1.2089.*(1.0+exp(-(v-18.4099)./29.3814)))+3.5./(1.0+exp((v+100.0)./29.3814)));
da=(ass-a)./ta;
iss=1.0./(1.0+exp((v+43.94)./5.711));


% change based on celltype
%%%%
delta_epi=ones(N,1);
% delta_epi(celltype_vector==1)=1.0-(0.95./(1.0+exp((v(celltype_vector==1)+70.0)./5.0)));
%%%%%


tiF=4.562+1./(0.3933.*exp((-(v+100.0))./100.0)+0.08004.*exp((v+50.0)./16.59));
tiS=23.62+1./(0.001416.*exp((-(v+96.52))./59.05)+1.780e-8.*exp((v+114.1)./8.079));
tiF=tiF.*delta_epi;
tiS=tiS.*delta_epi;
AiF=1.0./(1.0+exp((v-213.6)./151.2));
AiS=1.0-AiF;
diF=(iss-iF)./tiF;
diS=(iss-iS)./tiS;


scaleItos=zeros(N,1); % Change based on Clancy code
if flag_gender==1
    scaleItos(celltype_vector==0)=1*(1+ Gto_rt);
    scaleItos(celltype_vector==1)=1*(0.6 + Gto_rt);
elseif flag_gender==2
    scaleItos(celltype_vector==0)=1*(0.64 + Gto_rt);
    scaleItos(celltype_vector==1)=1*(0.26 + Gto_rt);
end

i=AiF.*iF+AiS.*iS.*scaleItos; %%% Change based on Clancy code


assp=1.0./(1.0+exp((-(v-24.34))./14.82));
dap=(assp-ap)./ta;
dti_develop=1.354+1.0e-4./(exp((v-167.4)./15.89)+exp(-(v-12.23)./0.2154));
dti_recover=1.0-0.5./(1.0+exp((v+70.0)./20.0));
tiFp=dti_develop.*dti_recover.*tiF;
tiSp=dti_develop.*dti_recover.*tiS;
diFp=(iss-iFp)./tiFp;
diSp=(iss-iSp)./tiSp;

ip=AiF.*iFp+AiS.*iSp.*scaleItos; %%% change based on Clancy code
Gto0=0.02;
Gto0=Gto0/Gtofc; % added for drug 

% change based on celltype
%%%%%%%%
Gto=Gto0*ones(N,1);
% Gto(celltype_vector==1)=Gto0.*4.0;
% Gto(celltype_vector==2)=Gto0.*4.0;
%%%%%%%%%


fItop=(1.0./(1.0+KmCaMK./CaMKa));
Ito=Gto.*(v-EK).*((1.0-fItop).*a.*i+fItop.*ap.*ip);

%calculate ICaL, ICaNa, ICaK
dss=1.0./(1.0+exp((-(v+3.940))./4.230));
td=0.6+1.0./(exp(-0.05.*(v+6.0))+exp(0.09.*(v+14.0)));
dd=(dss-d)./td;


fss=1.0./(1.0+exp((v+19.58)./3.696));
tff=7.0+1.0./(0.0045.*exp(-(v+20.0)./10.0)+0.0045.*exp((v+20.0)./10.0));
tfs=1000.0+1.0./(0.000035.*exp(-(v+5.0)./4.0)+0.000035.*exp((v+5.0)./6.0));
Aff=0.6;
Afs=1.0-Aff;
dff=(fss-ff)./tff;
dfs=(fss-fs)./tfs;
f=Aff.*ff+Afs.*fs;

fcass=fss;
tfcaf=7.0+1.0./(0.04.*exp(-(v-4.0)./7.0)+0.04.*exp((v-4.0)./7.0));
tfcas=100.0+1.0./(0.00012.*exp(-v./3.0)+0.00012.*exp(v./7.0));
Afcaf=0.3+0.6./(1.0+exp((v-10.0)./10.0));
Afcas=1.0-Afcaf;
dfcaf=(fcass-fcaf)./tfcaf;
dfcas=(fcass-fcas)./tfcas;
fca=Afcaf.*fcaf+Afcas.*fcas;
tjca=75.0;
djca=(fcass-jca)./tjca;

tffp=2.5.*tff;
dffp=(fss-ffp)./tffp;

fp=Aff.*ffp+Afs.*fs;
tfcafp=2.5.*tfcaf;
dfcafp=(fcass-fcafp)./tfcafp; 
fcap=Afcaf.*fcafp+Afcas.*fcas;

Kmn=0.002;
k2n=1000.0;
km2n=jca.*1.0;
anca=1.0./(k2n./km2n+(1.0+Kmn./cass).^4.0);
dnca=anca.*k2n-nca.*km2n;

PhiCaL=4.0.*vffrt.*(cass.*exp(2.0.*vfrt)-0.341.*cao)./(exp(2.0.*vfrt)-1.0);
PhiCaNa=1.0.*vffrt.*(0.75.*nass.*exp(1.0.*vfrt)-0.75.*nao)./(exp(1.0.*vfrt)-1.0);
PhiCaK=1.0.*vffrt.*(0.75.*kss.*exp(1.0.*vfrt)-0.75.*ko)./(exp(1.0.*vfrt)-1.0);
zca=2.0;
PCa0=0.0001;
PCa0=PCa0/GCaLfc; % added for drug 

% change based on celltype
%%%%%%%
PCa=PCa0*ones(N,1);
% PCa(celltype_vector==1)=PCa0*1.2;
% PCa(celltype_vector==2)=PCa0*2.5;
%%%%%%%%


PCap=1.1.*PCa;
PCaNa=0.00125.*PCa;
PCaK=3.574e-4.*PCa;
PCaNap=0.00125.*PCap;
PCaKp=3.574e-4.*PCap;
fICaLp=(1.0./(1.0+KmCaMK./CaMKa));
ICaL=factor_CaL*(1.0-fICaLp).*PCa.*PhiCaL.*d.*(f.*(1.0-nca)+jca.*fca.*nca)+fICaLp.*PCap.*PhiCaL.*d.*(fp.*(1.0-nca)+jca.*fcap.*nca); % factor_CaL is used for hormone based on Ele and Clancy code 
ICaL=SFCaL.*ICaL; 
ICaNa=(1.0-fICaLp).*PCaNa.*PhiCaNa.*d.*(f.*(1.0-nca)+jca.*fca.*nca)+fICaLp.*PCaNap.*PhiCaNa.*d.*(fp.*(1.0-nca)+jca.*fcap.*nca);
ICaK=(1.0-fICaLp).*PCaK.*PhiCaK.*d.*(f.*(1.0-nca)+jca.*fca.*nca)+fICaLp.*PCaKp.*PhiCaK.*d.*(fp.*(1.0-nca)+jca.*fcap.*nca);

%calculate IKr
xrss=1.0./(1.0+exp((-(v+8.337))./6.789));
txrf=12.98+1.0./(0.3652.*exp((v-31.66)./3.869)+4.123e-5.*exp((-(v-47.78))./20.38));
txrs=1.865+1.0./(0.06629.*exp((v-34.70)./7.355)+1.128e-5.*exp((-(v-29.74))./25.94));
Axrf=1.0./(1.0+exp((v+54.81)./38.21));
Axrs=1.0-Axrf;
dxrf=(xrss-xrf)./txrf;
dxrs=(xrss-xrs)./txrs;
xr=Axrf.*xrf+Axrs.*xrs;
rkr=1.0./(1.0+exp((v+55.0)./75.0)).*1.0./(1.0+exp((v-10.0)./30.0));
GKr0=factor_Kr*0.046; % factor_Kr is used for hormone based on Ele and Clancy code 
GKr0=GKr0/GKrfc;  % added for drug 


% change based on gendertype and celltype 
%%%%%%%%%
GKr=GKr0*ones(N,1);
if flag_gender==1
    GKr(celltype_vector==1)=GKr0.*(1.09+GKr_rt);
    GKr(celltype_vector==0)=GKr0.*(1+GKr_rt);    
elseif flag_gender==2
    GKr(celltype_vector==1)=GKr0.*(0.875+GKr_rt);
    GKr(celltype_vector==0)=GKr0.*(0.79+GKr_rt);  
else
    GKr(celltype_vector==1)=GKr0.*1.3;
    GKr(celltype_vector==2)=GKr0.*0.8;
end
%%%%%%%%%

IKr=GKr.*sqrt(ko./5.4).*xr.*rkr.*(v-EK);



%calculate IKs
GKs0=factor_Ks*0.0034; % factor_Ks is used for hormone based on Ele and Clancy code 
GKs0=GKs0/GKsfc; % added for drug 

% change based on celltype and gendertype
%%%%%%
GKs=GKs0*ones(N,1);
xs1ss=zeros(N,1);
txs1=zeros(N,1);
txs2=zeros(N,1);
if flag_gender==1
    GKs(celltype_vector==0)=GKs0.*(1+ GKs_rt);
    xs1ss(celltype_vector==0)=1.0./(1.0+exp((-(v(celltype_vector==0)+11.60))./8.932).*(1 + GKs_rt) );
    txs1(celltype_vector==0)=817.3+1.0./(2.326e-4.*exp((v(celltype_vector==0)+48.28)./17.80).*(1 + GKs_rt) + 0.001292.*exp((-(v(celltype_vector==0)+210.0))./230.0).*(1+ GKs_rt) );
    txs2(celltype_vector==0)=1.0./(0.01.*exp((v(celltype_vector==0)-50.0)./20.0).*(1+ GKs_rt) + 0.0193.*exp((-(v(celltype_vector==0)+66.54))./31.0).*(1 + GKs_rt) );
    
    
    GKs(celltype_vector==1)=GKs0.*(1.04 + GKs_rt);
	xs1ss(celltype_vector==1)=1.0./(1.0+exp((-(v(celltype_vector==1)+11.60))./8.932) .* (1.04 + GKs_rt) );
	txs1(celltype_vector==1)=817.3+1.0./(2.326e-4.*exp((v(celltype_vector==1)+48.28)./17.80).*(1.04 + GKs_rt) + 0.001292.*exp((-(v(celltype_vector==1)+210.0))./230.0).*(1.04 + GKs_rt) );
	txs2(celltype_vector==1)=1.0./(0.01.*exp((v(celltype_vector==1)-50.0)./20.0).*(1.04 + GKs_rt) + 0.0193.*exp((-(v(celltype_vector==1)+66.54))./31.0).*(1.04 + GKs_rt) );
elseif flag_gender==2
    GKs(celltype_vector==0)=GKs0.*(0.83 + GKs_rt);
	xs1ss(celltype_vector==0)=1.0./(1.0+exp((-(v(celltype_vector==0)+11.60))./8.932).*(0.83 + GKs_rt) );
	txs1(celltype_vector==0)=817.3+1.0./(2.326e-4.*exp((v(celltype_vector==0)+48.28)./17.80).*(0.83 + GKs_rt) + 0.001292.*exp((-(v(celltype_vector==0)+210.0))./230.0).*(0.83 + GKs_rt) );
    txs2(celltype_vector==0)=1.0./(0.01.*exp((v(celltype_vector==0)-50.0)./20.0).*(0.83 + GKs_rt) + 0.0193.*exp((-(v(celltype_vector==0)+66.54))./31.0).*(0.83 + GKs_rt) );
    
    GKs(celltype_vector==1)=GKs0.*(0.87 + GKs_rt);
	xs1ss(celltype_vector==1)=1.0./(1.0+exp((-(v(celltype_vector==1)+11.60))./8.932).*(0.87 + GKs_rt) );
	txs1(celltype_vector==1)=817.3+1.0./(2.326e-4.*exp((v(celltype_vector==1)+48.28)./17.80).*(0.87 + GKs_rt) + 0.001292.*exp((-(v(celltype_vector==1)+210.0))./230.0).*(0.87 + GKs_rt) );
	txs2(celltype_vector==1)=1.0./(0.01.*exp((v(celltype_vector==1)-50.0)./20.0).*(0.87 + GKs_rt) + 0.0193.*exp((-(v(celltype_vector==1)+66.54))./31.0).*(0.87 + GKs_rt) );
else
    GKs(celltype_vector==1)=GKs0.*1.4;
    xs1ss=1.0./(1.0+exp((-(v+11.60))./8.932));
    txs1=817.3+1.0./(2.326e-4.*exp((v+48.28)./17.80)+0.001292.*exp((-(v+210.0))./230.0));
    txs2=1.0./(0.01.*exp((v-50.0)./20.0)+0.0193.*exp((-(v+66.54))./31.0));
end

xs2ss=xs1ss;
dxs1=(xs1ss-xs1)./txs1;
dxs2=(xs2ss-xs2)./txs2;
KsCa=1.0+0.6./(1.0+(3.8e-5./cai).^1.4);
IKs=GKs.*KsCa.*xs1.*xs2.*(v-EKs);



% calculate IK1 
xk1ss=1.0./(1.0+exp(-(v+2.5538.*ko+144.59)./(1.5692.*ko+3.8115)));
txk1=122.2./(exp((-(v+127.2))./20.36)+exp((v+236.8)./69.33));
dxk1=(xk1ss-xk1)./txk1;



rk1=1.0./(1.0+exp((v+105.8-2.6.*ko)./9.493));
GK10=0.1908;
GK10=GK10/GK1fc; % added for drug 
% change based on celltype and gendertype 
%%%%%%
GK1=GK10*ones(N,1);
if flag_gender==1
    GK1(celltype_vector==0)=GK10.*(1+ GK1_rt);
    GK1(celltype_vector==1)=GK10.*(0.98 + GK1_rt);
    
elseif flag_gender==2
    GK1(celltype_vector==0)=GK10.*(0.86+ GK1_rt);
    GK1(celltype_vector==1)=GK10.*(0.74 + GK1_rt);
else
    GK1(celltype_vector==1)=GK10.*1.2;
    GK1(celltype_vector==2)=GK10.*1.3;
end
%%%%%%%
IK1=GK1.*sqrt(ko).*rk1.*xk1.*(v-EK);



%%%%%% revise based on Ele's code %%%%%%
if flag_gender==1
    scaleNaCas=1;
elseif flag_gender==2
    scaleNaCas=1.15;
end
%%%%%%%


%calculate INaCa_i
kna1=15.0;
kna2=5.0;
kna3=88.12;
kasymm=12.5;
wna=6.0e4;
wca=6.0e4;
wnaca=5.0e3;
kcaon=1.5e6;
kcaoff=5.0e3;
qna=0.5224;
qca=0.1670;
hca=exp((qca.*v.*F)./(R.*T));
hna=exp((qna.*v.*F)./(R.*T));
h1=1+nai./kna3.*(1+hna);
h2=(nai.*hna)./(kna3.*h1);
h3=1.0./h1;
h4=1.0+nai./kna1.*(1+nai./kna2);
h5=nai.*nai./(h4.*kna1.*kna2);
h6=1.0./h4;
h7=1.0+nao./kna3.*(1.0+1.0./hna);
h8=nao./(kna3.*hna.*h7);
h9=1.0./h7;
h10=kasymm+1.0+nao./kna1.*(1.0+nao./kna2);
h11=nao.*nao./(h10.*kna1.*kna2);
h12=1.0./h10;
k1=h12.*cao.*kcaon;
k2=kcaoff;
k3p=h9.*wca;
k3pp=h8.*wnaca;
k3=k3p+k3pp;
k4p=h3.*wca./hca;
k4pp=h2.*wnaca;
k4=k4p+k4pp;
k5=kcaoff;
k6=h6.*cai.*kcaon;
k7=h5.*h2.*wna;
k8=h8.*h11.*wna;
x1=k2.*k4.*(k7+k6)+k5.*k7.*(k2+k3);
x2=k1.*k7.*(k4+k5)+k4.*k6.*(k1+k8);
x3=k1.*k3.*(k7+k6)+k8.*k6.*(k2+k3);
x4=k2.*k8.*(k4+k5)+k3.*k5.*(k1+k8);
E1=x1./(x1+x2+x3+x4);
E2=x2./(x1+x2+x3+x4);
E3=x3./(x1+x2+x3+x4);
E4=x4./(x1+x2+x3+x4);
KmCaAct=150.0e-6;
allo=1.0./(1.0+(KmCaAct./cai).^2.0);
zna=1.0;
JncxNa=3.0.*(E4.*k7-E1.*k8)+E3.*k4pp-E2.*k3pp;
JncxCa=E2.*k2-E1.*k1;
Gncx0=0.0008;

% change based on celltype
%%%%%%%%
Gncx=Gncx0*ones(N,1);
% Gncx(celltype_vector==1)=Gncx0.*1.1;
% Gncx(celltype_vector==2)=Gncx0.*1.4;
%%%%%%%

INaCa_i=0.8.*Gncx.*allo.*(zna.*JncxNa+zca.*JncxCa).*scaleNaCas; %% revise based on Ele's code 



%calculate INaCa_ss
h1=1+nass./kna3.*(1+hna);
h2=(nass.*hna)./(kna3.*h1);
h3=1.0./h1;
h4=1.0+nass./kna1.*(1+nass./kna2);
h5=nass.*nass./(h4.*kna1.*kna2);
h6=1.0./h4;
h7=1.0+nao./kna3.*(1.0+1.0./hna);
h8=nao./(kna3.*hna.*h7);
h9=1.0./h7;
h10=kasymm+1.0+nao./kna1.*(1+nao./kna2);
h11=nao.*nao./(h10.*kna1.*kna2);
h12=1.0./h10;
k1=h12.*cao.*kcaon;
k2=kcaoff;
k3p=h9.*wca;
k3pp=h8.*wnaca;
k3=k3p+k3pp;
k4p=h3.*wca./hca;
k4pp=h2.*wnaca;
k4=k4p+k4pp;
k5=kcaoff;
k6=h6.*cass.*kcaon;
k7=h5.*h2.*wna;
k8=h8.*h11.*wna;
x1=k2.*k4.*(k7+k6)+k5.*k7.*(k2+k3);
x2=k1.*k7.*(k4+k5)+k4.*k6.*(k1+k8);
x3=k1.*k3.*(k7+k6)+k8.*k6.*(k2+k3);
x4=k2.*k8.*(k4+k5)+k3.*k5.*(k1+k8);
E1=x1./(x1+x2+x3+x4);
E2=x2./(x1+x2+x3+x4);
E3=x3./(x1+x2+x3+x4);
E4=x4./(x1+x2+x3+x4);
KmCaAct=150.0e-6;
allo=1.0./(1.0+(KmCaAct./cass).^2.0);
JncxNa=3.0.*(E4.*k7-E1.*k8)+E3.*k4pp-E2.*k3pp;
JncxCa=E2.*k2-E1.*k1;
INaCa_ss=0.2.*Gncx.*allo.*(zna.*JncxNa+zca.*JncxCa).*scaleNaCas; %% revise based on Ele's code 

%calculate INaK
k1p=949.5;
k1m=182.4;
k2p=687.2;
k2m=39.4;
k3p=1899.0;
k3m=79300.0;
k4p=639.0;
k4m=40.0;
Knai0=9.073;
Knao0=27.78;
delta=-0.1550;
Knai=Knai0.*exp((delta.*v.*F)./(3.0.*R.*T));
Knao=Knao0.*exp(((1.0-delta).*v.*F)./(3.0.*R.*T));
Kki=0.5;
Kko=0.3582;
MgADP=0.05;
MgATP=9.8;
Kmgatp=1.698e-7;
H=1.0e-7;
eP=4.2;
Khp=1.698e-7;
Knap=224.0;
Kxkur=292.0;
P=eP./(1.0+H./Khp+nai./Knap+ki./Kxkur);
a1=(k1p.*(nai./Knai).^3.0)./((1.0+nai./Knai).^3.0+(1.0+ki./Kki).^2.0-1.0);
b1=k1m.*MgADP;
a2=k2p;
b2=(k2m.*(nao./Knao).^3.0)./((1.0+nao./Knao).^3.0+(1.0+ko./Kko).^2.0-1.0);
a3=(k3p.*(ko./Kko).^2.0)./((1.0+nao./Knao).^3.0+(1.0+ko./Kko).^2.0-1.0);
b3=(k3m.*P.*H)./(1.0+MgATP./Kmgatp);
a4=(k4p.*MgATP./Kmgatp)./(1.0+MgATP./Kmgatp);
b4=(k4m.*(ki./Kki).^2.0)./((1.0+nai./Knai).^3.0+(1.0+ki./Kki).^2.0-1.0);
x1=a4.*a1.*a2+b2.*b4.*b3+a2.*b4.*b3+b3.*a1.*a2;
x2=b2.*b1.*b4+a1.*a2.*a3+a3.*b1.*b4+a2.*a3.*b4;
x3=a2.*a3.*a4+b3.*b2.*b1+b2.*b1.*a4+a3.*a4.*b1;
x4=b4.*b3.*b2+a3.*a4.*a1+b2.*a4.*a1+b3.*b2.*a1;
E1=x1./(x1+x2+x3+x4);
E2=x2./(x1+x2+x3+x4);
E3=x3./(x1+x2+x3+x4);
E4=x4./(x1+x2+x3+x4);
zk=1.0;
JnakNa=3.0.*(E1.*a3-E2.*b3);
JnakK=2.0.*(E4.*b1-E3.*a1);


Pnak0=30;
%%%%%% revise based on Ele's code (based on celltype, remove gendertype) %%%%%
Pnak=Pnak0*ones(N,1);
Pnak(celltype_vector==0)=Pnak0.*(1.0 + NaK_rt);
Pnak(celltype_vector==1)=Pnak0.*(0.94 + NaK_rt);
%%%%%%%%
INaK=Pnak.*(zna.*JnakNa+zk.*JnakK);



%calculate IKb
xkb=1.0./(1.0+exp(-(v-14.48)./18.34));
GKb0=0.003;

% change based on celltype
%%%%%%
GKb=GKb0*ones(N,1);
% GKb(celltype_vector==1)=GKb0.*0.6;
%%%%%%
IKb=GKb.*xkb.*(v-EK);


%calculate INab
PNab=3.75e-10;
INab=PNab.*vffrt.*(nai.*exp(vfrt)-nao)./(exp(vfrt)-1.0);


%calculate ICab
PCab=2.5e-8;
ICab=PCab.*4.0.*vffrt.*(cai.*exp(2.0.*vfrt)-0.341.*cao)./(exp(2.0.*vfrt)-1.0);

%calculate IpCa
%change based on celltype and gendertype
%%%%
GpCa0=0.0005;
GpCa=GpCa0*ones(N,1);
if flag_gender==1
    GpCa(celltype_vector==0)=GpCa0.*(1+pCa_rt);
    GpCa(celltype_vector==1)=GpCa0.*(0.88+pCa_rt);
elseif flag_gender==2
    GpCa(celltype_vector==0)=GpCa0.*(1.6+pCa_rt);
    GpCa(celltype_vector==1)=GpCa0.*(1.6+pCa_rt);
end
    
IpCa=GpCa.*cai./(0.0005+cai);

%calculate Ikatp
K_o_n = 5.4; % normal extracellular potassium concentration
akik = (ko/K_o_n)^0.24;
A_atp=2.0;
K_atp=0.25;
bkik = 1.0/(1.0+(A_atp/K_atp)^2.0);
gkatp = 4.4;
Ikatp = fkatp*gkatp*akik*bkik*(v-EK);





%calculate the stimulus current, Istim
amp=-80.0;
duration=0.5;
if t<=duration
    Istim=amp;
else
    Istim=0.0;
end



%calculate diffusion fluxes
JdiffNa=(nass-nai)./2.0;
JdiffK=(kss-ki)./2.0;
Jdiff=(cass-cai)./0.2;

%calculate ryanodione receptor calcium induced calcium release from the jsr
bt=4.75;
a_rel=0.5.*bt;

Jrel_inf0=a_rel.*(-ICaL)./(1.0+(1.5./cajsr).^8.0);

% change based on celltype
%%%%
Jrel_inf=Jrel_inf0;
Jrel_inf(celltype_vector==2)=Jrel_inf0(celltype_vector==2).*1.7;
%%%%%

tau_rel=bt./(1.0+0.0123./cajsr);
tau_rel(tau_rel<0.001)=0.001; 

dJrelnp=(Jrel_inf-Jrelnp)./tau_rel;

btp=1.25.*bt;
a_relp=0.5.*btp;
Jrel_infp0=a_relp.*(-ICaL)./(1.0+(1.5./cajsr).^8.0);


% change based on celltype
%%%%
Jrel_infp=Jrel_infp0;
Jrel_infp(celltype_vector==2)=Jrel_infp0(celltype_vector==2).*1.7; 
%%%%%%
tau_relp=btp./(1.0+0.0123./cajsr);
tau_relp(tau_relp<0.001)=0.001; 
dJrelp=(Jrel_infp-Jrelp)./tau_relp;
fJrelp=(1.0./(1.0+KmCaMK./CaMKa));
Jrel=(1.0-fJrelp).*Jrelnp+fJrelp.*Jrelp;

%calculate serca pump, ca uptake flux
Jupnp0=0.004375.*cai./(cai+0.00092);
Jupp0=2.75.*0.004375.*cai./(cai+0.00092-0.00017);



 
%%% revise based on Ele's code based on celltype, and remove gendertype %%%%%
Jupnp=Jupnp0; 
Jupp=Jupp0; 
Jupnp(celltype_vector==0)=Jupnp0(celltype_vector==0).*(1+ Gup_rt);
Jupp(celltype_vector==0)=Jupp0(celltype_vector==0).*(1 + Gup_rt);
Jupnp(celltype_vector==1)=Jupnp0(celltype_vector==1).*(1.42 + Gup_rt);
Jupp(celltype_vector==1)=Jupp0(celltype_vector==1).*(1.42 + Gup_rt);
%%%

fJupp=(1.0./(1.0+KmCaMK./CaMKa));
Jleak=0.0039375.*cansr./15.0;
Jup=(1.0-fJupp).*Jupnp+fJupp.*Jupp-Jleak;

%calculate tranlocation flux
Jtr=(cansr-cajsr)./100.0;



%calcium buffer constants

cmdnmax0=0.05;

% change based on celltype and gendertype 
%%%%%%
cmdnmax=cmdnmax0*ones(N,1);

if flag_gender==1
    cmdnmax(celltype_vector==0)=cmdnmax0.*(1 + CaM_rt);
    cmdnmax(celltype_vector==1)=cmdnmax0.*(1.07 + CaM_rt);
elseif flag_gender==2
    cmdnmax(celltype_vector==0)=cmdnmax0.*(1.21 + CaM_rt);
    cmdnmax(celltype_vector==1)=cmdnmax0.*(1.41 + CaM_rt);
else
    cmdnmax(celltype_vector==1)=cmdnmax0.*1.3;
end
%%%%%%


kmcmdn=0.00238;
trpnmax=0.07;
kmtrpn=0.0005;
BSRmax=0.047;
KmBSR=0.00087;
BSLmax=1.124;
KmBSL=0.0087;
csqnmax=10.0;
kmcsqn=0.8;



%update intracellular concentrations, using buffers for cai, cass, cajsr
dnai=-(INa+INaL+3.0.*INaCa_i+3.0.*INaK+INab).*Acap./(F.*vmyo)+JdiffNa.*vss./vmyo;
dnass=-(ICaNa+3.0.*INaCa_ss).*Acap./(F.*vss)-JdiffNa;

dki=-(Ito+IKr+IKs+IK1+IKb+Istim-2.0.*INaK).*Acap./(F.*vmyo)+JdiffK.*vss./vmyo;
dkss=-(ICaK).*Acap./(F.*vss)-JdiffK;

Bcai=1.0./(1.0+cmdnmax.*kmcmdn./(kmcmdn+cai).^2.0+trpnmax.*kmtrpn./(kmtrpn+cai).^2.0);
dcai=Bcai.*(-(IpCa+ICab-2.0.*INaCa_i).*Acap./(2.0.*F.*vmyo)-Jup.*vnsr./vmyo+Jdiff.*vss./vmyo);

Bcass=1.0./(1.0+BSRmax.*KmBSR./(KmBSR+cass).^2.0+BSLmax.*KmBSL./(KmBSL+cass).^2.0);
dcass=Bcass.*(-(ICaL-2.0.*INaCa_ss).*Acap./(2.0.*F.*vss)+Jrel.*vjsr./vss-Jdiff);

dcansr=Jup-Jtr.*vjsr./vnsr;

Bcajsr=1.0./(1.0+csqnmax.*kmcsqn./(kmcsqn+cajsr).^2.0);
dcajsr=Bcajsr.*(Jtr-Jrel);

Ion=INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa_i+INaCa_ss+INaK+INab+IKb+IpCa+ICab+Istim+Ikatp;
dv=-Ion; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
if flag_ode==1
    output=[dv dnai dnass dki dkss dcai dcass dcansr dcajsr dm dhf dhs dj... 
        dhsp djp dmL dhL dhLp da diF diS dap diFp diSp dd dff dfs dfcaf...
        dfcas djca dnca dffp dfcafp dxrf dxrs dxs1 dxs2 dxk1 dJrelnp dJrelp dCaMKt]';
   
 
else
    output=[INa INaL Ito ICaL IKr IKs IK1 INaCa_i INaCa_ss INaK  IKb INab ICab IpCa Jdiff JdiffNa JdiffK Jup Jleak Jtr Jrel CaMKa Istim Ikatp]';
end











    



