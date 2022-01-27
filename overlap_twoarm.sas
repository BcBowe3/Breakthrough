*************************************************************************************
This macro condcuts a two arm overlap weight cox surival regression.
It takes as intput:
A)dataset- contains the cohort groups (BREAK and REF), unique id (id),
	covariates, outcome date, outcome history, 
	end of followup date (endoffolow), and the start of followup date (t0)
B)history- An indicator of prior history of the disease being examined
	Used to exclude those with the disease. May be set to 1 or a date
C)outcome- This contains the date of the outcome
D)pssrez- storage of weights for covariate balance evaluation
E)hrrz-storage of hazard ratios
F)surrez-storage of survial probability estimates
G)hdps- dataset with high dimeisonal covariates in long form 
H)pre- dataset with predefined covariates in long form



;


%macro disease(history, outcome, dataset, hdps, pre, psrez, hrrez, surrez);
/*Initializes a dataset for surival probalibty estimation*/
data baseline;
length group $16.;
event=0; time=180;
group="BREAK";output;
group="REF";output;
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
model group(event="BREAK")=
&hdvs &pre;
output out=ps (keep=id group ps_score endoffollow &outcome ) pred=ps_score ;
run;

data &psrez ;
set ;
if group="REF" then ps_&outcome=ps;
if group="BREAK" then  ps_&outcome=(1-ps);
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
baseline out=sur covariates=baseline survival=sur upper=suru lower=surl;
run;


data hr;
set hr;
outcome="&outcome";
run;

data sur;
set sur;
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
