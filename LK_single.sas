/* LK.sas 1.1 10/26/15 */

%macro lk(var,k=10,alpha=0.05,data=_last_,data2=,where=1,id=,idxout=);

%let laglast=&syslast;
%if &data eq _last_ %then %let data=&syslast;

%if %scan(&var,2)^= %then %do;
   %put ERROR: macro %LK can only handle one variable at a time;
   %goto errsys;
   %end;

/* Subset the data to non-missing and where conditions. */
data; 
   set &data;
   where(&var^=. and (&where));
   keep &id &var;
   run;
%let temp=&syslast;

/* Get descriptive statistics. */
data; run; %let mn=&syslast;
proc means noprint data=&temp;
  var &var;
  output out=&mn n=__n sum=__sum;
run;

proc sort data=&temp;
  by &var;
  run;

/* Get k+1 observations from the upper tail of the distribution.    */
/* &temp contains x_(n-k) to calculate LK_k, while &temp2 does not. */
data;
  set &temp nobs=__nobs;
  keep &id &var __obs;
  if _n_>=(__nobs-&k) then do;
    __obs=_n_+&k-__nobs;
    output;
  end;
run;
%let temp2=&syslast;
data &temp(keep=&var); 
  set &temp2;
run;
data &temp2(keep=&id __obs);
  set &temp2;
  if __obs>0;
run;

/* Calculate the test statistics for each protential outlier. */
data &temp(keep=&var __test __sum __n __obs);  
  set &temp end=__fin;
  array __a[%eval(&k+1)];
  retain __a1-__a%eval(&k+1);
  __obs=_n_;
  __a[__obs]=&var;
  if __fin then do; 
    set &mn;
    do __i=1 to (&k); /* find as many as k outliers */
      __obs=&k+1-__i;
      __test=(__a[__obs+1]-__a[__obs])/__sum;
      &var=__a[__obs+1];
      output;
      __n=__n-1;
      __sum=__sum-__a[__obs+1];
    end;
  end;
run;

/* Estimate the number of outliers by removing one point each time. */
proc sort data=&temp;
  by __obs;
run;
data &temp;
  set &temp;
  length __flag $3.;
  __ts=1-(&alpha)**(1/(__n-1));
  retain __flag "  ";
  if __test>__ts or __flag="***" then __flag="***";
run;

proc sort data=&temp;
  by __obs;
run;
proc sort data=&temp2;
  by __obs;
run;
data &temp;
  merge &temp2 &temp;
  by __obs;
  label __n="n" __test="test" __ts="critical value" __flag="flag";
run;
proc sort data=&temp;
  by descending __obs;
run;

/* The last observation contains the number of true outliers in k protential outliers. */ 
data &temp;
  set &temp;
  retain count 0;
  if __flag="***" then count=count+1;
run;

/* Add a new observation to &data2, which contains #outliers. */
data new(keep=counts);
  set &temp end=last;
  if last then counts=count;
data &data2; set &data2 new;
data &data2(keep=counts); set &data2; where counts^=.;
run;

/* title "LK Maltiple Upper-Outlier Procedure";      */
/* title2 "Reference: S. Lalitha and N. Kumar, 2012, ""Multiple outlier test for upper outliers "; */
/* title3 "in an exponential sample"", Journal of Applied Statistics,39:6:1323-1330"; */
/*  */
/* proc print data=&temp split=" "; */
/*   var __n &var __test __ts __flag &id; */
/*   format __test __ts 6.3 &var 6.2; */
/* run; */
%if %length(&idxout)>0 %then %do;
   %if %length(&id)>0 %then %do;
   data &idxout(keep=&id);
      set &temp;
      if __flag="***";
      run;
   proc sort;
      by &id;
      run;
   %end;
   %else %do;
       %put WARNING(esd): no ID variables, IDXOUT database not created;
       %let syslast=&laglast;
       %end;
%end;
%else %do;
   %let syslast=&laglast;
%end;

proc datasets nolist;
   delete 
      %scan(&temp,2,.)
      %scan(&mn,2,.)
      %scan(&temp2,2,.);;
   quit;


title;
%goto out;
%errsys: %do;
         %put ERROR: The macro LK discontinued execution;
         %goto out;
         %end;
%out: 
%mend lk;

/* Estimation of the proportion of #outlier */
/* It will return a table containing the proportion of #outlier */
%macro lk_est(var=iteration,k=10,iteration=,alpha=0.05,data=,id=);

data proportion;
  label counts="Proportion of #outliers";
run;

%do i=1 %to &iteration;
  %lk(&var%eval(&i),k=&k,alpha=&alpha,data=&data,data2=proportion,id=&id);
%end;

/* Calculate the size of sample. */
proc means noprint data=&data;
  var &id;
  output out=temp n=__n;
run;
data _null_; set temp; call symput('numb',__n);
run;

title "LK Maltiple Upper-Outlier Procedure";
title2 "#observations=&numb, #outliers= &k, #iterations= &iteration. ";
proc freq data=proportion; tables counts; run;

%mend;
