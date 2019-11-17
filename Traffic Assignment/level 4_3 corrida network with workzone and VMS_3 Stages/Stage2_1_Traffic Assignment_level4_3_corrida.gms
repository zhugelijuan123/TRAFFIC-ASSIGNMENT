$TITLE OD DEMAND ASSIGNMENT PROBLEM

OPTIONS  ITERLIM=1000, RESLIM = 1000000, SYSOUT = OFF, SOLPRINT = OFF, NLP = MINOS5,OPTCR= 0.1, LIMROW = 0, LIMCOL = 0;
SET I ZONES/1*2/;
SET K LINKS/1*55/;
SET P PATHS/1*4/;
SET GP GPATHS/1*2/;

ALIAS (I,J);
ALIAS (I,I1);
ALIAS (J,J1);
ALIAS (P,P1);
ALIAS (K,L);
ALIAS (K,K1);


$include C:\Users\zhuge\Desktop\trans committee\learning documents\Traffic Assignment\dataset for gams\INPUT_DEMAND.inc
$include C:\Users\zhuge\Desktop\trans committee\learning documents\Traffic Assignment\dataset for gams\PARAMETER PATH_LINK.inc
$include C:\Users\zhuge\Desktop\trans committee\learning documents\Traffic Assignment\dataset for gams\PARAMETER LINK_CAPACITY_WITH WORKZONE.inc
$include C:\Users\zhuge\Desktop\trans committee\learning documents\Traffic Assignment\dataset for gams\PARAMETER LINK_FFTT.inc


*redifined the path p between OD(ij) with the VMS link
PARAMETER VMS_PATH_LINK(I,J,P,K)/
   1.   2.   1.   2   1
   1.   2.   2.   2   1
/;

*link id of the VMS setting link
PARAMETER VMS(K)/
   2   1
/;

PARAMETER PATH_PRO(I,J,P)/
   1.   2.   3   0.25
   1.   2.   4   0.25
/;


PARAMETER ALPHA;
ALPHA=0.15;

PARAMETER BETA;
BETA=4;

*positive VARIABLES
POSITIVE VARIABLES
ESPATHFLOW(I,J,P)   estimated pathflow on path p of OD pair (ij) ;

POSITIVE VARIABLES
ESVMSPATHLINKPRO(I,J,P) estimated pathflow(with VMS) proportion on OD pair(ij);

POSITIVE VARIABLES
ESVMSPATHFLOW(I,J,P,K) estimated pathflow(with VMS) on OD pair(ij);

POSITIVE VARIABLES
ESLINKFLOW(K)   estimated link total flow with all the going through path;

POSITIVE VARIABLES
ESVMSLINKFLOW(K) estimated link total flow with VMS;

POSITIVE VARIABLES
LINK_TT(K)    link real travel time with BPR;

POSITIVE VARIABLES
PATH_TT(I,J,P)    path travel time of the sum of link travel time;

POSITIVE VARIABLES
PIE(I,J)   the shortest path for OD pair(ij);

POSITIVE VARIABLES
LINK_PROPORTION(I,J,K) link flow proportion for OD(ij);

VARIABLE
SO;


EQUATIONS
VMSF(I,J,P,K)       *define the path(with VMS) proportion of the corresponding link(with VMS)
ESPATHFLOWF_1(I,J)      *the sum of ESPATHFLOW(IJP)=DEMAND(IJ)
ESPATHFLOWF_2(I,J,P,K)    *ESPATHFLOW(IJP)>=0
PATH_PROF(I,J,P)
ESVMSPATHFLOWF(I,J,P,K)    *get the ESVMSPATHFLOW from the ESPATHFLOW
ESVMSPATHFLOWFF(I,J)       *make sure the ESVMSPATHFLOW(with VMS) proportion of the corresponding OD(ij) reflect the situation when there is no workzone
ESLINKFLOWF_1(K)        *calculate eslinkflow(K) from espathflow(IJP)
ESLINKFLOWF_2(I,J,P,K)        *eslinkflow(K)<=capacity(K)
ESVMSLINKFLOWF(K)       *get the ESVMSLINKFLOWF from the ESLINKFLOWF
LINK_TTF(K)             *BPR counting for link real travel time
PATH_TTF(I,J,P)             *path travel time
PIEF(I,J,P)             *shortest path for OD(IJ)
LINK_PROPORTIONF(I,J,K) *calculate link proportion for the use of ODME
SOF;

VMSF(I,J,P,K)$(VMS_PATH_LINK(I,J,P,K) AND VMS(K))..ESVMSPATHLINKPRO(I,J,P)=E=(ESVMSPATHFLOW(I,J,P,K)/(ESVMSLINKFLOW(K)+0.0001));
ESPATHFLOWF_1(I,J).. SUM(P$(DEMAND(I,J)),ESPATHFLOW(I,J,P))=E=DEMAND(I,J);
ESPATHFLOWF_2(I,J,P,K)$(DEMAND(I,J) AND PATH_LINK(I,J,P,K)).. ESPATHFLOW(I,J,P)=G=0;
PATH_PROF(I,J,P)$PATH_PRO(I,J,P)..ESPATHFLOW(I,J,P)/DEMAND(I,J)=E=PATH_PRO(I,J,P);
ESVMSPATHFLOWF(I,J,P,K)$(VMS_PATH_LINK(I,J,P,K))..ESVMSPATHFLOW(I,J,P,K)=E=ESPATHFLOW(I,J,P);
ESVMSPATHFLOWFF(I,J)$(DEMAND(I,J))..SUM((P,K)$(VMS_PATH_LINK(I,J,P,K) AND VMS(K) AND DEMAND(I,J)),ESVMSPATHFLOW(I,J,P,K)/(DEMAND(I,J)+0.0001))=E=0.50;
ESLINKFLOWF_1(K).. ESLINKFLOW(K)=E=SUM((I,J,P)$(DEMAND(I,J) AND PATH_LINK(I,J,P,K)),ESPATHFLOW(I,J,P)*PATH_LINK(I,J,P,K));
ESLINKFLOWF_2(I,J,P,K)$(DEMAND(I,J) AND PATH_LINK(I,J,P,K)).. ESLINKFLOW(K)=L=LINK_CAPACITY(K);
ESVMSLINKFLOWF(K)$VMS(K)..ESVMSLINKFLOW(K)=E=ESLINKFLOW(K);
LINK_TTF(K)..LINK_TT(K)=E=LINK_FFTT(K)*(1+ALPHA*(ESLINKFLOW(K)/LINK_CAPACITY(K)));
PATH_TTF(I,J,P).. PATH_TT(I,J,P)=E=SUM(K$(DEMAND(I,J) AND PATH_LINK(I,J,P,K)),PATH_LINK(I,J,P,K)*LINK_TT(K));
PIEF(I,J,P)$(DEMAND(I,J)).. PIE(I,J)=L=PATH_TT(I,J,P);
LINK_PROPORTIONF(I,J,K).. LINK_PROPORTION(I,J,K)=E=SUM(P$(DEMAND(I,J) AND PATH_LINK(I,J,P,K)),(ESPATHFLOW(I,J,P)*PATH_LINK(I,J,P,K))/DEMAND(I,J));
SOF..SO=E=12;


MODEL ODE /ALL/;
SOLVE ODE USING NLP MINIMIZING SO;


File TrafficAssignment /TrafficAssignment.dat/;
put TrafficAssignment;


DISPLAY SO.L;
DISPLAY ESPATHFLOW.L;
DISPLAY ESVMSPATHFLOW.L;
DISPLAY ESVMSPATHLINKPRO.L;
DISPLAY ESLINKFLOW.L;
DISPLAY ESVMSLINKFLOW.L;
DISPLAY LINK_TT.L;
DISPLAY PATH_TT.L;
DISPLAY PIE.L;
DISPLAY LINK_PROPORTION.L;
DISPLAY LINK_FFTT;




