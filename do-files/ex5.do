
************** Estimation: (i) parametric and Cox models *********


************* Cancer data **********
sysuse cancer, clear
de
stset studytim , failure(died)
ge id = _n  
lab var id "subject identifier"

* drug = 1 (placebo); drug =2,3 (receives drug)
ta drug died
recode drug 1=0 2/3=1
lab var drug "receives drug?"
lab def drug 0 "placebo" 1 "drug"
lab val drug drug

ta drug

**************** WEIBULL *****************************
* Weibull model, PH, coeffients

streg drug age, dist(weibull) nolog nohr

* Elasticity of hazard w.r.t. age (age covariate in levels) = b_age * age
ge elas_age = _b[age]*age
su elas_age, detail

* Elasticity of hazard w.r.t. age (age covariate in logs) = b_age
ge lage = log(age)
streg drug lage, dist(weibull) nolog nohr

* return to original regression
streg drug age, dist(weibull) nolog nohr

di "Check that -b/p from the PH spec. equals b* from the AFT spec (below):"
di " "
di "-b[drug]/p  = "  -_b[drug]/e(aux_p)
di "-b[age]/p  = "  -_b[age]/e(aux_p)
di "-b[_cons]/p  = "  -_b[_cons]/e(aux_p)


* ratios of hazards at different points in time, given same X
* For Weibull model h(t)/h(u) = (t/u)^(alpha-1)

di "h(10,X)/h(5,X) =  "  (10/5)^(e(aux_p)-1)
di "h(20,X)/h(10,X) =  "  (20/10)^(e(aux_p)-1)
di "h(20,X)/h(5,X) =  "  (20/5)^(e(aux_p)-1)
di "h(30,X)/h(20,X) =  "  (30/20)^(e(aux_p)-1)
di "h(30,X)/h(5,X) =  "  (30/5)^(e(aux_p)-1)

* Ratio of hazards for two persons with different X, same t
* h(t,X1)/h(t,X2) = exp[ (X1-X2)b ].  

di "h(t;age=y+10,drug=x)/h(t;age=y,drug=x) = "  exp(_b[age]*10)

di "h(t;age=y,drug=0)/h(t;age=y,drug=1) = "  exp(_b[drug])

di "h(t;age=y+10,drug=1)/h(t;age=y,drug=0) = "  exp(_b[age]*10 + _b[drug])

* use -predict- to derive estimate of predicted median duration

* median is [(1/L)*{ln(2)}]^(1/a)  where L = exp(xb)
* mean is (1/L)^(1/a)*G(1+(1/a)) where G(x) is the Gamma function

predict xb, xb
su xb
di "WEIBULL MODEL"
di "Pred. Median [at sample mean X] = "  (ln(2)*exp(-r(mean)))^(1/e(aux_p)) 
di "Pred. Mean [at sample mean X] = " exp(-r(mean)/e(aux_p))*exp(lngamma(1+1/e(aux_p)))

* median duration for each person in sample
ge mediand = (ln(2)*exp(-xb))^(1/e(aux_p))
* expected (mean) duration for each person in sample
ge meand = exp(-xb/e(aux_p))*exp(lngamma(1+1/e(aux_p)))

* NB Stata allows you to generate these directly:
* e.g. 	predict mediand, median time
* 	predict meand, meand time

predict mediandS, median time
predict meandS, mean time

su mediand mediandS meand meandS 

* listing predicted means and medians for different groups

sort age drug
list id age drug mediand meand if age >50 & age <60 & drug ==0 , noobs
list id age drug mediand meand if age >50 & age <60 & drug ==1 , noobs

drop xb mediand meand mediandS meandS

su age
local meana = r(mean)
local xb0 = _b[_cons] + _b[age]*`meana' + _b[drug]*0

di "Mean age = " `meana' " ,_b[_cons] + _b[age]*(mean age) + _b[drug]*0 = " `xb0'
di "Pred. Median [mean(age), drug=0] = "  (ln(2)*exp(-`xb0'))^(1/e(aux_p))
di "Pred. Mean [mean(age),drug=0] = " exp(-`xb0'/e(aux_p))*exp(lngamma(1+1/e(aux_p)))
local xb0 = _b[_cons] + _b[age]*`meana' + _b[drug]*1
di "Mean age = " `meana' " ,_b[_cons] + _b[age]*(mean age) + _b[drug]*1 = " `xb0'
di "Pred. Median [mean(age), drug=1] = "  (ln(2)*exp(-`xb0'))^(1/e(aux_p))
di "Pred. Mean [mean(age),drug=1] = " exp(-`xb0'/e(aux_p))*exp(lngamma(1+1/e(aux_p)))

* use -stcurv- to look at estimated survivor and hazard functions
*  for person with sample mean values of covariates

stcurv, survival title("Cancer data, at sample means") ///
	saving(streg1,replace) 
stcurv, hazard title("Cancer data, at sample means") ///
	saving(streg2,replace)

* out of sample predictions using the range option

stcurv, survival title("Survival, out-of-sample prediction, at sample means") ///
	saving(streg1a,replace) range(0 50)
stcurv, hazard title("Hazard, out-of-sample prediction, at sample means") ///
	saving(streg2a,replace) range(0 50)


*  comparing person with drug 0/1
stcurv, hazard at(drug=0) title("Cancer data: drug=0") ///
	saving(streg3,replace) 
stcurv, hazard title("Cancer data:drug=1") ///
	saving(streg4,replace) 
stcurv, survival title("Cancer data:drug=0") at(drug=0) ///
	saving(streg5,replace) 
stcurv, survival title("Cancer data:drug=1") at(drug=1) ///
	saving(streg6,replace) 

* put both curves on one graph

stcurv, survival title("Survival, Cancer data: drug=0,1") at(drug=0) at(drug=1) ///
	saving(streg5a,replace) 


* Weibull model, PH, exp(coefficients), via replay of earlier
streg, hr

* Weibull model, AFT, coefficients

streg drug age, dist(weibull) nolog time

* Weibull model, AFT, exp(coefficients), via replay of earlier
streg, tr

**************** LOG-LOGISTIC *****************************

****** log-logistic model, coeffs

streg drug age, dist(logl) nolog

* use -predict- to derive estimate of predicted median duration

* predicted median is 1/L where L = exp(-xb).  So 1/L = exp(xb)
* predicted mean is (1/L)*(_pi*gamma)/sin(_pi*gamma) if gamma < 1
*
* ... and gamma is less than 1 here

predict xb, xb
su xb
di "LOG-LOGISTIC MODEL"
di "Pred. Median [at sample mean X] = "  exp(r(mean))
di "Pred. Mean [at sample mean X] = "  exp(r(mean))*(_pi*e(gamma))/sin(_pi*e(gamma))

* median duration for each person in sample
* NB one can also generate these directly using -predict-
ge mediand = exp(xb)
ge meand = exp(xb)*(_pi*e(gamma))/sin(_pi*e(gamma))


predict mediandS, median time
predict meandS, mean time
su mediand mediandS meand meandS


sort age drug
list id age drug mediand meand if age >50 & age <60 & drug ==0 , noobs
list id age drug mediand meand if age >50 & age <60 & drug ==1, noobs

drop xb mediand meand mediandS meandS

* mean age is already in local macro meana 

local xb0 = _b[_cons] + _b[age]*`meana' + _b[drug]*0
di "Mean age = " `meana' " ,_b[_cons] + _b[age]*(mean age) + _b[drug]*0 = " `xb0'
di "Pred. Median [mean(age), drug=0] = "  exp(`xb0')
di "Pred. Mean [mean(age), drug=0] = "  exp(`xb0')*(_pi*e(gamma))/sin(_pi*e(gamma))
local xb0 = _b[_cons] + _b[age]*`meana' + _b[drug]*1
di "Mean age = " `meana' " ,_b[_cons] + _b[age]*(mean age) + _b[drug]*1 = " `xb0'
di "Pred. Median [mean(age), drug=1] = "  exp(`xb0')
di "Pred. Mean [mean(age), drug=1] = "  exp(`xb0')*(_pi*e(gamma))/sin(_pi*e(gamma))

* use stcurve to look at estimated survivor and hazard functions
*  for person with sample mean values of covariates

stcurv, survival title("Cancer data, at sample means") ///
	saving(streg7,replace) 
stcurv, hazard  title("Cancer data, at sample means") ///
	saving(streg8,replace)

*  comparing person with drug 0/1
stcurv, hazard at(drug=0) title("Cancer data:drug=0") ///
	saving(streg9,replace) 
stcurv, hazard at(drug=1) title("Cancer data:drug=1") ///
	saving(streg10,replace) 
stcurv, survival title("Cancer data:drug=0") at(drug=0) ///
	saving(streg11,replace) 
stcurv, survival title("Cancer data:drug=1") at(drug=1) ///
	saving(streg12,replace) 

stcurv, hazard at(drug=0) at1(drug=1) title("Cancer data:drug=0,1") ///
	saving(streg9a,replace) 
stcurv, survival at(drug=0) at1(drug=1) title("Cancer data:drug=0,1") ///
	saving(streg11a,replace) 


************* OTHER MODELS ********************

****** lognormal model, coeffs

streg drug age, dist(logn) nolog


* NB one gets exactly the same results from a censored normal
*	regression on log(survival time)!

ge cens = 1-died
ge lntime = ln(studytim)
cnreg lntime drug age, censored(cens)

* How badly does OLS do?

reg lntime drug age 		/* ignoring censoring */
reg lntime drug age if died==1 	/* dropping censored cases */

****** generalised gamma model, coeffs

streg drug age, dist(gamma) nolog

* Gamma model has ancillary parameters kappa and sigma. Special cases:
* kappa = 1: Weibull model
* kappa = 1, sigma = 1: Exponential 
* kappa = 0: Lognormal model

* Wald test of null hypothesis kappa = 0 (Lognormal model)
*	Shown in output above: see "z and P>|z|" on /kappa row

* Wald test of null hypothesis kappa = 1 (Weibull model)
* test stat = [(kappa-1)/se(kappa)]^2.  Cf. with Chi2(1).

* Code Wald test something like the following:
mat v = e(V)
mat list v
local vk = v[5,5]
di "`vk'"
di "Reject null if test statistic greater than critical value = Chi2(1)"
di "Test statistic is: " (e(kappa)-1)^2/`vk' ". Critical value: " chi2tail(1,.05)


* use stcurve to look at estimated survivor and hazard functions
*  for person with sample mean values of covariates

stcurv, survival title("Cancer data, at sample means") ///
	saving(streg13,replace) 
stcurv, hazard title("Cancer data, at sample means") 	///
	saving(streg14,replace)

*  comparing person with drug 0/1 (and mean age)
stcurv, hazard at(drug=0) title("Cancer data:drug=0") ///
	saving(streg15,replace) 
stcurv, hazard at(drug=1) title("Cancer data:drug=1") ///
	saving(streg16,replace) 
stcurv, survival at(drug=0) title("Cancer data:drug=0")  ///
	saving(streg17,replace) 
stcurv, survival at(drug=1) title("Cancer data:drug=1")  ///
	saving(streg18,replace) 


stcurv, hazard at(drug=0) at1(drug=1) title("Cancer data:drug=0,1") ///
	saving(streg15a,replace) 
stcurv, survival at(drug=0) at1(drug=1) title("Cancer data:drug=0,1")  ///
	saving(streg17a,replace) 

******************* Cox Model ***************

stcox drug age, nohr
stcox, hr

stcox drug age, nohr bases(s0) basech(ch0) 

/* Stata 7 graphics
gr7 s0 _t, c(J) s(i) sort xlab ylab b1("Baseline S(t),age=0,Cox model") /*
	*/ l1(Survival) saving(stcox1,replace)
gr7 ch0 _t, c(J) s(i) sort xlab ylab b1("Baseline Cum.Haz.,age=0,Cox model") /* 
	*/ l1(Cumulative Hazard) saving(stcox2,replace)
*/

twoway line s0 _t, sort connect(J) title("Baseline S(t),age=0,Cox model") ///
	saving(stcox1, replace) 
twoway line ch0 _t, sort connect(J) title("Baseline Cum.Haz.,age=0,Cox model") ///
	saving(stcox2, replace) 


* Generally better practice to generate baseline functions at 
* values of x within the range of the data.  
* Above example generated them at age=0!  (NB mean age=55.9)
* So redo at age = 55
ge age55 = age-55
stcox drug age55, nohr bases(s1) basech(ch1)

twoway line s1 _t, sort connect(J) title("Baseline S(t),age=55,Cox model") ///
	saving(stcox3, replace) 
twoway line ch1 _t, sort connect(J) title("Baseline Cum.Haz.,age=55,Cox model") ///
	saving(stcox4, replace) 


* How do the results compare?
graph combine stcox1.gph stcox2.gph stcox3.gph stcox4.gph, saving(stcox5, replace)

stcox age, nohr strata(drug)


*** Piece-wise constant exponential model as per Lesson 3 ***
*** split at times 8, 17

use cancer, clear
ge id = _n
recode drug 1=0 2/3=1
lab var drug "receives drug?"
lab def drug 0 "placebo" 1 "drug"
lab val drug drug

stset studytim, f(died) id(id)

* split at survival times 8 and 17 
* ... key intervals are (0,8], (8,17], (17,39]

stsplit ehaz, at(8 17)
ta ehaz, ge(e)

* report hazard ratios
streg drug age e2 e3, dist(exp) nolog
streg, nohr
* report coefficients
streg drug age e2 e3, dist(exp) nolog nohr
streg drug age e1 e2 e3, dist(exp) nolog nocons nohr
streg drug age e1 e2 e3, dist(exp) nolog  


*****************************************
*************** Marriage data **********
*
* married variable distinguishes between legal marriages and cohabiting partnerships

use duration, clear
stset time , failure(status)
ge id = _n  
lab var id "subject identifier"
 
* create dummy variable for legal marital status
ge single = married == 1
lab define single 1 "Cohabiting" 0 "Legally married"
lab val single single
lab var single "Legal status of partnership"

streg single sex, dist(weibull) nolog nohr
streg, hr

* use -predict- to derive estimate of predicted median duration

* median is [(1/L)*{ln(2)}]^(1/a)  where L = exp(xb)
* mean is (1/L)^(1/a)*G(1+(1/a)) where G(x) is the Gamma function

predict xb, xb
su xb
di "WEIBULL MODEL"
di "Pred. Median [at sample mean X] = "  (ln(2)*exp(-r(mean)))^(1/e(aux_p)) 
di "Pred. Mean [at sample mean X] = " exp(-r(mean)/e(aux_p))*exp(lngamma(1+1/e(aux_p)))

* median duration for each person in sample
ge mediand = (ln(2)*exp(-xb))^(1/e(aux_p))
* expected (mean) duration for each person in sample
ge meand = exp(-xb/e(aux_p))*exp(lngamma(1+1/e(aux_p)))

sort mediand
list id married sex mediand meand , noobs

drop xb mediand meand

local xb0 = _b[_cons] + _b[single]*1 + _b[sex]*0
di "Pred. Median [single,woman] = "  (ln(2)*exp(-`xb0'))^(1/e(aux_p))
di "Pred. Mean [single,woman] = " exp(-`xb0'/e(aux_p))*exp(lngamma(1+1/e(aux_p)))

local xb0 = _b[_cons] + _b[single]*0 + _b[sex]*0
di "Pred. Median [married,woman] = "  (ln(2)*exp(-`xb0'))^(1/e(aux_p))
di "Pred. Mean [married,woman] = " exp(-`xb0'/e(aux_p))*exp(lngamma(1+1/e(aux_p)))

local xb0 = _b[_cons] + _b[single]*1 + _b[sex]*1
di "Pred. Median [single,man] = "  (ln(2)*exp(-`xb0'))^(1/e(aux_p))
di "Pred. Mean [single,man] = " exp(-`xb0'/e(aux_p))*exp(lngamma(1+1/e(aux_p)))

local xb0 = _b[_cons] + _b[single]*0 + _b[sex]*1
di "Pred. Median [married,man] = "  (ln(2)*exp(-`xb0'))^(1/e(aux_p))
di "Pred. Mean [married,man] = " exp(-`xb0'/e(aux_p))*exp(lngamma(1+1/e(aux_p)))

* use -stcurv- to look at estimated survivor and hazard functions

*  comparing 4 groups as above
stcurv, hazard at(single=1 sex=0) title("Cohabiting,woman") ///
	saving(streg19,replace) 
stcurv, hazard at(single=0 sex=0) title("Married,woman") ///
	saving(streg20,replace) 
stcurv, hazard at(single=1 sex=1) title("Cohabiting,man") ///
	saving(streg21,replace) 
stcurv, hazard at(single=0 sex=1) title("Married,man") ///
	saving(streg22,replace) 

stcurv, hazard at(single=0 sex=0) at1(single=0 sex=1) ///
	at2(single=1 sex=0) at3(single=1 sex=1) ///
	title("Hazard, by type") ///
	saving(streg22a,replace) 

*sample means
stcurv, hazard title(All sample) saving(streg23,replace) 

*  comparing 4 groups as above
stcurv, survival at(single=1 sex=0) title("Cohabiting,woman") ///
	saving(streg24,replace) 
stcurv, survival at(single=0 sex=0) title("Married,woman") ///
	saving(streg25,replace) 
stcurv, survival at(single=1 sex=1) title("Cohabiting,man") ///
	saving(streg26,replace) 
stcurv, survival at(single=0 sex=1) title("Married,man") ///
	saving(streg27,replace) 

stcurv, survival at(single=0 sex=0) at1(single=0 sex=1) ///
	at2(single=1 sex=0) at3(single=1 sex=1) ///
	title("Survivor function, by type") ///
	saving(streg27a,replace) 


*sample means
stcurv, survival title(All sample) saving(streg28,replace) 


* Weibull model, PH, exp(coefficients), via replay of earlier
streg, hr

* Weibull model, AFT, coefficients

streg single sex, dist(weibull) nolog time

* Weibull model, AFT, exp(coefficients), via replay of earlier
streg, tr

****** log-logistic model, coeffs

streg single sex, dist(logl) nolog

* use -predict- to derive estimate of predicted median duration

* predicted median is 1/L where L = exp(-xb).  So 1/L = exp(xb)
* predicted mean is (1/L)*(_pi*gamma)/sin(_pi*gamma) if gamma < 1

predict xb, xb
su xb
di "LOG-LOGISTIC MODEL"
di "Pred. Median [at sample mean X] = "  exp(r(mean))
di "Pred. Mean [at sample mean X] = "  exp(r(mean))*_pi*e(gamma)/sin(_pi*e(gamma))

* median and mean duration for each person in sample
ge mediand = exp(xb)
ge meand = exp(xb)*_pi*e(gamma)/sin(_pi*e(gamma))

sort mediand
list id married sex mediand  meand , noobs

drop xb mediand meand

local xb0 = _b[_cons] + _b[single]*1 + _b[sex]*0
di "Pred. Median [single,woman] = "   exp(`xb0')
di "Pred. Mean [single,woman] = "   exp(`xb0')*_pi*e(gamma)/sin(_pi*e(gamma))

local xb0 = _b[_cons] + _b[single]*0 + _b[sex]*0
di "Pred. Median [married,woman] = "  exp(`xb0')
di "Pred. Mean [married,woman] = "  exp(`xb0')*_pi*e(gamma)/sin(_pi*e(gamma))

local xb0 = _b[_cons] + _b[single]*1 + _b[sex]*1
di "Pred. Median [single,man] = "   exp(`xb0')
di "Pred. Mean [single,man] = "   exp(`xb0')*_pi*e(gamma)/sin(_pi*e(gamma))

local xb0 = _b[_cons] + _b[single]*0 + _b[sex]*1
di "Pred. Median [married,man] = "  exp(`xb0')
di "Pred. Mean [married,man] = "  exp(`xb0')*_pi*e(gamma)/sin(_pi*e(gamma))

* use -stcurv- to look at estimated survivor and hazard functions

*  comparing 4 groups as above
stcurv, hazard at(single=1 sex=0) title("Cohabiting,woman") ///
	saving(streg29,replace) 
stcurv, hazard at(single=0 sex=0) title("Married,woman") ///
	saving(streg30,replace) 
stcurv, hazard at(single=1 sex=1) title("Cohabiting,man") ///
	saving(streg31,replace) 
stcurv, hazard at(single=0 sex=1) title("Married,man") ///
	saving(streg32,replace) 

stcurv, hazard at(single=0 sex=0) at1(single=0 sex=1) ///
	at2(single=1 sex=0) at3(single=1 sex=1) ///
	title("Hazard, by type") ///
	saving(streg32a,replace) 

*sample means
stcurv, hazard  title(All sample) saving(streg33,replace) 

*  comparing 4 groups as above
stcurv, survival at(single=1 sex=0) title("Cohabiting,woman") ///
	saving(streg34,replace) 
stcurv, survival at(single=0 sex=0) title("Married,woman") ///
	saving(streg35,replace) 
stcurv, survival at(single=1 sex=1) title("Cohabiting,man") ///
	saving(streg36,replace) 
stcurv, survival at(single=0 sex=1) title("Married,man") ///
	saving(streg37,replace) 

stcurv, survival at(single=0 sex=0) at1(single=0 sex=1) ///
	at2(single=1 sex=0) at3(single=1 sex=1) ///
	title("Survivor function, by type") ///
	saving(streg37a,replace) 


*sample means
stcurv, survival  title(All sample) saving(streg38,replace) 


* Cox 
stcox single sex, nohr
stcox, hr

stcox single sex, nohr bases(s0) basech(ch0) 

twoway line s0 _t, sort connect(J) title("Baseline S(t),age=0,Cox model") ///
	saving(stcox6, replace) 
twoway line ch0 _t, sort connect(J) title("Baseline Cum.Haz.,age=0,Cox model") ///
	saving(stcox7, replace) 


