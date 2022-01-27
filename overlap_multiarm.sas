/*************************************************************************************
This macro condcuts a three arm overlap weight cox survival regression.
The code may be used as a base and extended for more arms as desired.

It takes as intput:
A)dataset- contains the cohort groups (GROUP1 GROUP2 GROUP3), unique id (id),
	covariates, outcome date, outcome history, 
	end of followup date (endoffolow), and the start of followup date (t0)
B)history- An indicator of prior history of the disease being examined
	Used to exclude those with the disease. May be set to 1 or a date
C)outcome- This contains the date of the outcome
D)pssrez- storage of weights for covariate balance evaluation
E)hrrz-storage of hazard ratios
F)surrez-storage of survial probability estimates
G)hdps- dataset with high dimeisonal covariate names in long form 
H)pre- dataset with predefined covariates names in long form
*************************************************************************************/


%macro disease(history, outcome, dataset, hdps, pre, psrez, hrrez, surrez);
/*Initializes a dataset for surival probalibty estimation*/
data baseline;
length group $16.;
event=0; time=180;
group="GROUP1";output;
group="GROUP2";output;
group="GROUP3";output;
run;

/*High dimeisonal covariates*/
data eee;
set &hdps;
length cov $32767.;
retain cov;
cov=catx(' ',cov,_name_); 
call symput('hdvs',cov);
run;

/*predefind covariates*/
data eee;
set &pre;
length cov $32767.;
retain cov;
cov=catx(' ',cov,_name_); 
call symput('pre',cov);
run;

/*The HD selection is only done once, in an outside program;*/
data model;
set &dataset;
by scrssn;
if a;
if &history<=0;
run;


proc logistic data=model noprint;
class sex race smoke   ;
model group(event="GROUP1")=
&hdvs &pre;
output out=ps1 (keep=id  ps_1 ) pred=ps_1 ;
run;

proc logistic data=model noprint;
class sex race smoke   ;
model group(event="GROUP2")=
&hdvs &pre;
output out=ps2 (keep=id  ps_2) pred=ps_2 ;
run;

proc logistic data=model noprint;
class sex race smoke   ;
model group(event="GROUP3")=
&hdvs &pre;
output out=ps3 (keep=id  ps_3  ) pred=ps_3 ;
run;

proc sort data=ps1; by id; run;
proc sort data=ps2; by id; run;
proc sort data=ps3; by id; run;

data &psrez ;
merge ps1 ps2 ps3 ;
by id;
base=(1/PS_1 + 1/PS_2 + 1/PS_3)
if group="GROUP1" then ps_&outcome=(1/PS_1)/base;
if group="GROUP2" then ps_&outcome=(1/PS_2)/base;
if group="GROUP3" then ps_&outcome=(1/PS_3)/base;
drop base ps_1 ps_2 PS_3;
run;


data outcome;
set &psrez (where=(ps_&outcome^=.));
if &outcome ^=. and &outcome<=endoffollow then do;
	event=1;
	time=&outcome -t0;
end;
else do;
	event=0;
	time=endoffollow-(t0+30);
end;
run;

/*Survival model*/
ods output hazardratios=hr_out;
proc phreg data=outcome covs(aggregate);
class group (ref="REF");
model time*event(0)=group;
hazardratio group/diff=pairwise;
weight ps_&outcome;
id id;
baseline out=sur_out covariates=baseline survival=sur upper=suru lower=surl;
run;


data hr;
set hr_out;
outcome="&outcome";
run;

data sur;
set sur_out;
outcome="&outcome";
run;

data &hrrez;
set &hrrez  hr;
run;

data &surrez;
set &surrez  sur;
run;
%mend;


*Storage for hazard ratios;
data  HR;
length outcome $50;
set _null_;
run;
*Storage for survival probability estiamtes;
data  SUR;
length outcome $50;
set _null_;
run;
*storage of weights;
data  PS;
set cohort (keep=id);
run;

%disease (history,outcome,dataset, hdps, pre, PS, HR,SUR);
