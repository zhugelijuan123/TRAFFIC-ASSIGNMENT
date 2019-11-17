$TITLE OD ESTIMATION PROBLEM

OPTIONS  ITERLIM=1000, RESLIM = 1000000, SYSOUT = OFF, SOLPRINT = OFF, NLP = MINOS5,OPTCR= 0.1, LIMROW = 0, LIMCOL = 0;
SET I NODES/1*4/;
SET K LINKS/1*5/;
SET P PATHS/1*2/;

ALIAS (I,J);
ALIAS (K,L);

PARAMETER DEMAND(I,J) /
    1.   2     600
    1.   4     2300
/;

PARAMETER PATH_LINK(I,J,P,K) /
   1.   2.   1.   1   1
   1.   2.   2.   2   1
   1.   2.   2.   3   1
   1.   4.   1.   2   1
   1.   4.   1.   4   1
   1.   4.   2.   5   1
/;

PARAMETER LINK_CAPACITY(K) /
    1   600
    2   1800
    3   500
    4   1500
    5   1500
/;

PARAMETER LINK_FFTT(K) /
    1   20
    2   18
    3   3
    4   6
    5   26
/;

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
ESLINKFLOWF_2(K)        *eslinkflow(K)<=capacity(K)
LINK_TTF(I,J,P,K)             *BPR counting for link real travel time
PATH_TTF(I,J,P)             *path travel time
PIEF(I,J,P)             *shortest path for OD(IJ)
LINK_PROPORTIONF(I,J,K) *calculate link proportion for the use of ODME;

ESPATHFLOWF_1(I,J).. SUM(P$DEMAND(I,J),ESPATHFLOW(I,J,P))=e=DEMAND(I,J);
ESPATHFLOWF_2(I,J,P,K)$PATH_LINK(I,J,P,K).. ESPATHFLOW(I,J,P)=G=0;
ESLINKFLOWF_1(K).. ESLINKFLOW(K)=E=SUM((I,J,P)$(PATH_LINK(I,J,P,K)),ESPATHFLOW(I,J,P)*PATH_LINK(I,J,P,K));
ESLINKFLOWF_2(K).. ESLINKFLOW(K)=L=LINK_CAPACITY(K);
LINK_TTF(I,J,P,K).. LINK_TT(K)=E=LINK_FFTT(K)*(1+ALPHA*((ESLINKFLOW(K)/LINK_CAPACITY(K))));
PATH_TTF(I,J,P)$(DEMAND(I,J)).. PATH_TT(I,J,P)=E=SUM(K,PATH_LINK(I,J,P,K)*LINK_TT(K));
PIEF(I,J,P)$(DEMAND(I,J)).. PIE(I,J)=L=PATH_TT(I,J,P);
LINK_PROPORTIONF(I,J,K).. LINK_PROPORTION(I,J,K)=E=SUM(P$PATH_LINK(I,J,P,K),(ESPATHFLOW(I,J,P)*PATH_LINK(I,J,P,K))/DEMAND(I,J));


VARIABLE
OBJ     gap of path choice;

EQUATIONS
OBJF    define objective function;

OBJF..  OBJ=E=sum((i,j,p), ESPATHFLOW(I,J,P)*( PATH_TT(I,J,P)-PIE(I,J) )*( PATH_TT(I,J,P)-PIE(I,J)) );



MODEL ODE /ALL/;
SOLVE ODE USING NLP MINIMIZING OBJ;

VARIABLES
ESLINKFLOW(K)      ESTIMATED link FLOW;


File estdemand /estdemand.dat/;
put estdemand;


loop((K)$(ESLINKFLOW.L(K)),put @1, K.tl, @3, ESLINKFLOW.L(K)/);
loop((I,J,P)$(ESPATHFLOW.L(I,J,P)),put @1,I.tl,@3,J.tl,@5,P.tl,@7,ESPATHFLOW.L(I,J,P),@20,PIE.L(I,J)/);

DISPLAY OBJ.L;
DISPLAY ESLINKFLOW.L;
DISPLAY ESPATHFLOW.L;
DISPLAY PIE.L;
DISPLAY LINK_TT.L;
DISPLAY PATH_TT.L;


