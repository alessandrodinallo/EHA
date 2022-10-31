* Collective exercise

global dir /Users/aledinal/Dropbox/PostDoc/Teaching/EHA/data

cd   "$dir"

/* sysuse cancer, clear */ // You may also decide to use the dataset "cancer" that is provided by default by STATA
sysuse cancer, clear 
su
de
stset studytim , failure(died)

// look at new variables created by -stset- (prefixed with "_")
de _*
su _*
st
stdes
stsum


ge id = _n  /* create unique person identifier */
lab var id "subject identifier"

* Look at the data format for the first four persons
list id studytim died drug age in 1/4

* Now expand the data set so that one obs per month at risk of death

expand studytim   /* see -help expand- */

bysort id: ge seqvar = _n   // -bysort- is like -sort- followed by -by-

lab var seqvar "spell month identifier, by subject"

	// NB generation of variable using a logical condition
bysort id: ge dead = died == 1 & _n==_N
lab var dead "binary depvar for discrete hazard model"


stset seqvar, failure(dead) id(id)
stsum
de _*
su _*


ge logd = ln(seqvar)
lab var logd "ln(t)"


lab var id "subject identifier"
recode drug 1=0 2/3=1
lab var drug "receives drug?"
lab def drug 0 "placebo" 1 "drug"
lab val drug drug

sts list

sts graph, title("Survivor function, Cancer data (sts)") ///
	xtick(1(1)39)  
sts graph, title("Cumulative hazard function, Cancer data (sts)") ///
	na xtick(1(1)39)  

* smoothed hazard, default bandwidth (and default kernel)
sts graph, hazard title("Smoothed hazard function, Cancer data (sts)") ///
	xtick(1(1)39)  

****************************
**** stratified estimates ***

sts graph, title("Survivor functions, by drug, Cancer data (sts)") ///
	by(drug) lost //saving(surv5, replace)

sts list, by(drug) 
sts list, by(drug) compare
sts test drug
sts test drug, wilcoxon
	
	
**************** WEIBULL *****************************
* Weibull model, PH, coeffients
streg drug age, dist(weibull) nolog nohr

	
* Elasticity of hazard w.r.t. age (age covariate in logs) = b_age
ge lage = log(age)
streg drug lage, dist(weibull) nolog nohr
	

streg drug age, dist(weibull) nolog nohr

* Ratio of hazards for two persons with different X, same t
* h(t,X1)/h(t,X2) = exp[ (X1-X2)b ].  

di "h(t;age=y+10,drug=x)/h(t;age=y,drug=x) = "  exp(_b[age]*10)

di "h(t;age=y,drug=0)/h(t;age=y,drug=1) = "  exp(_b[drug])

di "h(t;age=y+10,drug=1)/h(t;age=y,drug=0) = "  exp(_b[age]*10 + _b[drug])
	
	
	

	
	
	
	
	

************* Cancer data **********
sysuse cancer, clear

ge id = _n  
lab var id "subject identifier"

* drug = 1 (placebo); drug =2,3 (receives drug)
ta drug died
recode drug 1=0 2/3=1
lab var drug "receives drug?"
lab def drug 0 "placebo" 1 "drug"
lab val drug drug

ta drug

************************************
* Episode-splitting --> data in person-month format

expand studytim   
bysort id: ge j = _n   
* spell month identifier, by subject
lab var j "spell month"
bysort id: ge dead = died==1 & _n==_N
lab var dead "binary depvar for discrete hazard model"

sort id j
order id j studytime died dead drug

* We don't have to -stset- the data for estimation, but might as
*    well -- it emphasises parallels with continuous time models
*    esp. when there are TVCs.
stset j, failure(dead) id(id)

* Note correspondences between generated and stset vbles
su j _t dead _d

* Create duration-specific dummy variables, one for each spell month at risk
* (the maximum number is 39 here). There are 2 methods (at least)

* Perhaps most straightforward way (also gives var label automatically)
ta j, ge(d)

ds d*

* alternative using -forvalues-
forvalues x = 1(1)39 {
	ge byte durat`x' = (j == `x')
}

ds durat*

* su durat* d*
drop durat*


ge e1 = j < 9
ge e2 = j >= 9 & j <= 17
ge e3 = j >= 18 & j <.


ge Late = j > 17 & j < .

ge dur1 = d1+d2+d3+d4+d5+d6
ge dur2 = d7+d8+d9+d10+d11+d12
ge dur3 = d13+d14+d15+d16+d17+d18
ge dur4 = d19+d20+d21+d22+d23+d24
ge dur5 = d25+d26+d27+d28+d29+d30
ge dur6 = d31+d32+d33+d34+d35+d36+d37+d38+d39


ge lnj = ln(j)
ge j2 = j^2
ge j3 = j^3


order id  studytime lnj j j2 j3 dur*  d1-d39 e1 e2 e3 died dead drug
br 
compress

	
logit dead drug age lnj,  
logit dead drug age lnj, nolog
logit dead drug age lnj, or
logit, or

logistic dead drug age lnj, 
logistic dead drug age lnj, nolog
logistic dead drug age lnj, or
logistic dead drug age lnj, log

predict h, p

bysort id (j): ge s = exp(sum(ln(1-h)))

ge h0 = h if age == 55 & drug == 0
ge h1 = h if age == 55 & drug == 1
lab var h0 "drug = 0"
lab var h1 "drug = 1"
twoway (connect h0 j , sort  msymbol(t) ) (connect h1 j, sort msymbol(o) ) ///
	, title("Time = log(t)") name(lnj, replace) 

ge s0 = s if age == 55 & drug == 0
ge s1 = s if age == 55 & drug == 1
lab var s0 "drug = 0"
lab var s1 "drug = 1" /*
twoway (connect s0 j if age == 55, sort  msymbol(t) ) (connect s1 j, sort msymbol(o) ) ///
	, title("S(j),c(t)=(q-1)ln(j),age=55") name(dlogit2, replace) */


drop h h0 h1 s s0 s1	



* cubic polynomial
logit dead drug age j j2 j3, nolog

predict h, p


bysort id (j): ge s = exp(sum(ln(1-h)))

ge h0 = h if age == 55 & drug == 0
ge h1 = h if age == 55 & drug == 1
lab var h0 "drug = 0"
lab var h1 "drug = 1"
/*
twoway (connect h0 j, sort  msymbol(t) ) (connect h1 j, sort msymbol(o) ) ///
	, title("h(j),c(j) cubic polyn.,age=55") name(dlogit3, replace) */

twoway (connect h0 j, sort  msymbol(t) ) (connect h1 j, sort msymbol(o) ) ///
	, title("Time = cubic polyn.") name(cub, replace) 

	
ge s0 = s if age == 55 & drug == 0
ge s1 = s if age == 55 & drug == 1
lab var s0 "drug = 0"
lab var s1 "drug = 1"
twoway (connect s0 j, sort  msymbol(t) ) (connect s1 j, sort msymbol(o) ) ///
	, title("S(j),c(j) cubic polyn.,age=55") name(dlogit4, replace) 


drop h h0 h1  s s0 s1

* Combine log and cubic approximation of hazard
graph combine lnj cub


* piece-wise constant baseline
logit dead  drug age dur1-dur6, nocons nolog
logit, or
logit dead drug age dur2-dur6 , nolog
logit, or

logit dead drug age e2 e3, nolog
logit, or


logit dead drug age Late, nolog
logit, or

predict h, p


bysort id (j): ge s = exp(sum(ln(1-h)))

ge h0 = h if age == 55 & drug == 0
ge h1 = h if age == 55 & drug == 1
lab var h0 "drug = 0"
lab var h1 "drug = 1"
/*
twoway (connect h0 j, sort  msymbol(t) ) (connect h1 j, sort msymbol(o) ) ///
	, title("h(j),c(j) piecewise constant,age=55") name(dlogit5, replace) */

twoway (connect h0 j, sort  msymbol(t) ) (connect h1 j, sort msymbol(o) ) ///
	, title("Time = piece-wise constant") name(piecewise, replace) 

	
ge s0 = s if age == 55 & drug == 0
ge s1 = s if age == 55 & drug == 1
lab var s0 "drug = 0"
lab var s1 "drug = 1"
twoway (connect s0 j, sort  msymbol(t) ) (connect s1 j, sort msymbol(o) ) ///
	, title("S(j),c(j) piecewise constant,age=55") name(dlogit6, replace) 

drop h h0 h1 s s0 s1




* fully non-parametric baseline

set matsize 100

logit dead drug age d1-d39, nocons nolog
logit dead drug age d1-d8 d10-d13 d15-d17 d22-d25 d28 d33 /*
      */ if (j>=1 & j<=8) | (j>=10 & j<=13) /*
      */ | (j>=15 & j<=17) | (j>=22 & j<=25) /*
      */ | j==28 | j==33 , nocons nolog
logit, or

* use -predict- to derive estimate of predicted hazard and survivor function
* and thence median duration

** (1) set the hazard for months in which no events equal to 0 (implausible!!)

predict h, p
replace h = 0 if j==9|j==14|(j >=18&j<= 21)|j==26|j==27|(j>=29&j<=32)|j>=34
bysort id (j): ge s = exp(sum(ln(1-h)))

ge h0 = h if age == 55 & drug == 0
ge h1 = h if age == 55 & drug == 1
lab var h0 "drug = 0"
lab var h1 "drug = 1"

twoway (connect h0 j, sort  msymbol(t) ) (connect h1 j, sort msymbol(o) ) ///
	, xtick(1(1)39) title("Time = fully non-parametric") name(nonpar, replace) 

/*
twoway (connect h0 j, sort  msymbol(t) ) (connect h1 j, sort msymbol(o) ) ///
	, xtick(1(1)39) title("p(j),c(j) fully non-par.,age=55") saving(dlogit7a, replace) 
*/

* Combine log and cubic approximation of hazard
graph combine piecewise nonpar
		
/*
** (2) redo model assuming piece-wise constant to 'cover' unidentified bits
**	but still has no events for one group at some times

drop h h0 h1 s

ge byte d89 = d8 + d9
ge byte d1314 = d13+d14
ge byte d1721 = d17 + d18 + d19 + d20 + d21
ge byte d2527 = d25 + d26 + d27
ge byte d2832 = d28 + d29 + d30 + d31 + d32
ge byte d3339 = d33 + d34 + d35 + d36 + d37 + d38 + d39
logit dead drug age d1-d7 d89 d10-d12 d1314 d15 d16 d1721 d22-d24 d2527 d2832 d3339, nocons nolog
*/

** (3) redo model assuming more groupig (piece-wise constant) to 'cover' unidentified bits
drop h h0 h1 s


logit dead drug age dur1 dur2 dur3 dur4 dur5 dur6, nocons nolog

predict h, p
bysort id (j): ge s = exp(sum(ln(1-h)))


ge h0 = h if age == 55 & drug == 0
ge h1 = h if age == 55 & drug == 1
lab var h0 "drug = 0"
lab var h1 "drug = 1"
	
twoway (connect h0 j, sort  msymbol(t) ) (connect h1 j, sort msymbol(o) ) ///
	, xtick(1(1)39) title("Time = piece-wise constant (more groups)") name(piecewise2, replace) 
	
* Combine log and cubic approximation of hazard
graph combine piecewise piecewise2	

ge s0 = s if age == 55 & drug == 0
ge s1 = s if age == 55 & drug == 1
lab var s0 "drug = 0"
lab var s1 "drug = 1"
twoway (connect s0 j, sort  msymbol(t) ) (connect s1 j, sort msymbol(o) ) ///
	, xtick(1(1)39) title("S(j),c(j) piecewise-constant,age=55") 


drop h h0 h1 s s0 s1 



*--------------------
*Dataset on separation
 use    sep, clear
 
 * STATA needs to be given an id (coupleid, in this case) and a time counter variable (wave)
 xtset coupleid wave

  order coupleid pidp_? wave sep union_duration union_duration_ln union hiqual_f age_xt_f income income_miss age_homo background hiqual_? education_homo background
 br
 
 global dur  union_duration_ln 
 global demo c.age_xt_f##c.age_xt_f /* ib3.age_xt_cat_f */ i.age_homo i.ethnicity i.background // i.background_? 
 global edu  i.hiqual_?	
 
    logit sep union_duration_ln c.age_xt_f##c.age_xt_f c.income##c.income i.income_miss i.hiqual_f ib1.education_homo i.union
	
	margins, at(union_duration_ln=(0(0.5)4.5))
	marginsplot

	logit sep union_duration_ln c.age_xt_f##c.age_xt_f c.income##c.income i.income_miss i.hiqual_f ib1.education_homo i.union
	margins, at(age_xt_f=(15(5)60))
	marginsplot

	recode hiqual_f (99 = 7)
    logit sep union_duration_ln c.age_xt_f##c.age_xt_f c.income##c.income i.income_miss i.hiqual_f ib1.education_homo i.union, or
	margins hiqual_f
	marginsplot, horizontal
	
	* Quadratic time function
logit sep c.union_duration##c.union_duration c.age_xt_f##c.age_xt_f c.income##c.income i.income_miss i.hiqual_f ib1.education_homo i.union

	margins, at(union_duration=(1(5)30))
	marginsplot

	
	
* Recurrent events
xtlogit sep union_duration_ln , re  // But you can choose other shapes for time
* You can add more variables to reduce omitted variable bias



* Unobserved heterogeneity ('frailty') examples

************** continuous time ********************************


use bc.dta, clear

stset t, f(dead) 

** weibull **

streg age smoking dietfat, d(weib) nohr nolog

* omit dietfat as a way of introducing unobserved heterogeneity
* basic model w/out dietfat and no frailty

streg age smoking , d(weib) nohr nolog

streg age smoking , d(weib) nolog frailty(gamma) nohr 

* now reintroduce dietfat -- does frailty disappear?

streg age smoking dietfat, d(weib) nolog frailty(gamma) nohr
streg age smoking dietfat, d(weib) nolog frailty(gamma)


*** cancer data set, Weibull model ***

use cancer, clear

recode drug 1 = 0 2/3=1
stset studytim, f(died)

streg age drug, d(weib) nolog 
streg age drug, d(weib) nolog nohr

* now see if can introduce frailty
* base model - no covariates, and then frailty models

streg , d(weib) nolog nohr
streg , d(weib) nolog nohr frailty(gamma)




* Independent competing risks example


************** discrete time *********************************
global dir /Users/aledinal/Dropbox/PostDoc/Teaching/EHA/data
cd   "$dir"

use unemp, clear
ta status exit

ge status2 = 0
replace status2 = 1 if status == 2
replace status2 = 2 if status == 3
lab def status2 0 "censored" 1 "Exit-job" 2 "Exit-other"
lab val status2 status2
ta status status2





* treat survival times as discrete
* prepare the data (expand; new censoring vbles)

expand conmths


bysort newid: ge t = _n
bysort newid (t): ge picklast = _n ==_N
lab var picklast "1=selects last month observed"

* single destination censoring vble
bysort newid (t): ge leftui = exit == 1 & _n==_N
lab var leftui "1=Exit UI"

* multiple destination censoring vbles
bysort newid (t): ge y_noUA = status == 0 & _n == _N if status < .
lab var y_noUA "1=Exhaust UI, no UA"
bysort newid (t): ge y_UA = status == 1 & _n == _N if status < .
lab var y_UA "1=Exhaust UI, UA"
bysort newid (t): ge job = status == 2 & _n == _N if status < .
lab var job "1=Exit UI to job"
bysort newid (t): ge other = status == 3 & _n == _N if status < .
lab var other "1=Exit UI to other dest."


* create some covariates
* age, famresp, tyentry, and ...
ta groupreg, ge(reg)
ge logt = ln(t)

order newid t leftui y_* status status2

************ logit versions of above eqns; cf with mlogit

* estimate a single destination PH model with log(time) hazard
logit leftui age famresp tyentry reg1-reg4 logt, nolog

* estimate the components of the competing risk model 

logit job age famresp tyentry reg1-reg4 logt, nolog
logit, or

logit other age famresp tyentry reg1-reg4 logt, nolog
logit, or

sort newid t

* MNL estimates based on revised censoring vble 

ge depvar = 0
bysort newid (t): replace depvar = 1 if status2==1 & _n==_N
bysort newid (t): replace depvar = 2 if status2==2 & _n==_N
lab def depvar 0 "unemployed" 1 "Exit-job" 2 "Exit-other"
lab val depvar depvar
ta depvar status2

order newid t status2 depvar leftui job other age famresp groupreg status 


mlogit depvar age famresp tyentry reg1-reg4 logt, 
mlogit depvar age famresp tyentry reg1-reg4 logt, nolog baseoutcome(0) 
mlogit, rrr
