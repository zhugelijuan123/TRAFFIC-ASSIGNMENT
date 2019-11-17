$TITLE OD DEMAND ASSIGNMENT PROBLEM

OPTIONS  ITERLIM=1000, RESLIM = 1000000, SYSOUT = OFF, SOLPRINT = OFF, NLP = MINOS5,OPTCR= 0.1, LIMROW = 0, LIMCOL = 0;
SET I ZONES/1*2/;
SET K LINKS/1*55/;
SET P PATHS/1*4/;
SET GP GPATHS/1*2/;

ALIAS (I,J);
ALIAS (K,L);

$include C:\Users\zhuge\Desktop\trans committee\learning documents\Traffic Assignment\dataset for gams\INPUT_DEMAND.inc
$include C:\Users\zhuge\Desktop\trans committee\learning documents\Traffic Assignment\dataset for gams\PARAMETER PATH_LINK.inc
$include C:\Users\zhuge\Desktop\trans committee\learning documents\Traffic Assignment\dataset for gams\PARAMETER LINK_CAPACITY.inc
$include C:\Users\zhuge\Desktop\trans committee\learning documents\Traffic Assignment\dataset for gams\PARAMETER LINK_FFTT.inc

PARAMETER ALPHA;
ALPHA=0.15;

PARAMETER BETA;
BETA=4;


*positive VARIABLES
POSITIVE VARIABLES
ESPATHFLOW(I,J,P)   estimated pathflow on path p of OD pair (ij) ;

POSITIVE VARIABLES
ESLINKFLOW(K)   estimated link total flow with all the going through path;

POSITIVE VARIABLES
LINK_TT(K)    link real travel time with BPR;

POSITIVE VARIABLES
PATH_TT(I,J,P)    path travel time of the sum of link travel time;

POSITIVE VARIABLES
PIE(I,J)   the shortest path for OD pair(ij);

POSITIVE VARIABLES
LINK_PROPORTION(I,J,K) link flow proportion for OD(ij);


EQUATIONS
ESPATHFLOWF_1(I,J)      *the sum of ESPATHFLOW(IJP)=DEMAND(IJ)
ESPATHFLOWF_2(I,J,P,K)    *ESPATHFLOW(IJP)>=0
ESLINKFLOWF_1(K)        *calculate eslinkflow(K) from espathflow(IJP)
ESLINKFLOWF_2(I,J,P,K)        *eslinkflow(K)<=capacity(K)
LINK_TTF(K)             *BPR counting for link real travel time
PATH_TTF(I,J,P)             *path travel time
PIEF(I,J,P)             *shortest path for OD(IJ)
LINK_PROPORTIONF(I,J,K) *calculate link proportion for the use of ODME;


ESPATHFLOWF_1(I,J).. SUM(P$(DEMAND(I,J)),ESPATHFLOW(I,J,P))=e=DEMAND(I,J);
ESPATHFLOWF_2(I,J,P,K)$(DEMAND(I,J) AND PATH_LINK(I,J,P,K)).. ESPATHFLOW(I,J,P)=G=0;
ESLINKFLOWF_1(K).. ESLINKFLOW(K)=E=SUM((I,J,P)$(DEMAND(I,J) AND PATH_LINK(I,J,P,K)),ESPATHFLOW(I,J,P)*PATH_LINK(I,J,P,K));
ESLINKFLOWF_2(I,J,P,K)$(DEMAND(I,J) AND PATH_LINK(I,J,P,K)).. ESLINKFLOW(K)=L=LINK_CAPACITY(K);
LINK_TTF(K)..LINK_TT(K)=E=LINK_FFTT(K)*(1+ALPHA*(ESLINKFLOW(K)/LINK_CAPACITY(K)));
PATH_TTF(I,J,P).. PATH_TT(I,J,P)=E=SUM(K$(DEMAND(I,J) AND PATH_LINK(I,J,P,K)),PATH_LINK(I,J,P,K)*LINK_TT(K));
PIEF(I,J,P)$(DEMAND(I,J)).. PIE(I,J)=L=PATH_TT(I,J,P);
LINK_PROPORTIONF(I,J,K).. LINK_PROPORTION(I,J,K)=E=SUM(P$(DEMAND(I,J) AND PATH_LINK(I,J,P,K)),(ESPATHFLOW(I,J,P)*PATH_LINK(I,J,P,K))/DEMAND(I,J));



VARIABLE
OBJ     gap of path choice;

EQUATIONS
OBJF    define objective function;

OBJF..  OBJ=E=SUM(K,ESLINKFLOW(K)*LINK_TT(K));  



MODEL ODE /ALL/;
SOLVE ODE USING NLP MINIMIZING OBJ;


File TrafficAssignment /TrafficAssignment.dat/;
put TrafficAssignment;



DISPLAY OBJ.L;
DISPLAY ESPATHFLOW.L;
DISPLAY ESLINKFLOW.L;
DISPLAY LINK_TT.L;
DISPLAY PATH_TT.L;
DISPLAY PIE.L;
DISPLAY LINK_PROPORTION.L;
DISPLAY LINK_FFTT;




