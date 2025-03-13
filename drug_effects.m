function[Gfc]=drug_effects(DOSE,CMAX,IC50,HILL)
Gfc=1+(DOSE*CMAX/IC50)^HILL;
