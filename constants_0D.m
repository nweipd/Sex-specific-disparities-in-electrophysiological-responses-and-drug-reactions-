%physical constants 
R=8314.0;
T=310.0;
F=96485.0;

%cell geometry
L=0.01; % cm 
rad=0.0011;
vcell=1000*pi*rad*rad*L;
Ageo=2*pi*rad*rad+2*pi*rad*L;
Acap=2*Ageo;
vmyo=0.68*vcell;
vnsr=0.0552*vcell;
vjsr=0.0048*vcell;
vss=0.02*vcell;


% action potential propagation
v_start=-60;
v_end=v_start;

alternans_thr=3; 
iter_max=20;





