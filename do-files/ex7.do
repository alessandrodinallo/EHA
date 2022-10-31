
* Unobserved heterogeneity ('frailty') examples

************** continuous time ********************************



use bc.dta, clear

stset t, f(dead) 

** weibull **

streg age smoking dietfat, d(weib) nohr nolog

* omit dietfat as a way of introducing unobserved heterogeneity
* basic model w/out dietfat and no frailty

streg age smoking , d(weib) nohr nolog

predict xb, xb
su xb
di "Pred. Median [at sample mean X] = "  (ln(2)*exp(-r(mean)))^(1/e(aux_p)) 
* median duration for each person in sample

* NB Stata allows you to generate these directly:
predict mediand, time
su mediand
su mediand, de
drop xb mediand

streg age smoking , d(weib) nolog

* basic model w/out dietfat but now allowing for frailty

streg age smoking , d(weib) nohr nolog frailty(gamma)

predict xb, xb
su xb
di "Pred. Median [at sample mean X] = "  ((2^e(theta) - 1)/ (e(theta)*exp(r(mean))))^(1/e(aux_p))

* median duration for each person in sample
* NB Stata onwards allows you to generate these directly:
predict mediand, time
su mediand
su mediand, de


* repeat but get predicted median at mean frailty value
predict mediand2, time alpha1
su mediand2
su mediand2, de

drop mediand mediand2 xb


streg age smoking , d(weib) nohr nolog frailty(invgauss)
predict xb, xb
su xb
di "Pred. Median [at sample mean X] = "  ( ((1+e(theta)*ln(2))^2 - 1)/(2*e(theta)*exp(r(mean))) )^(1/e(aux_p))

* median duration for each person in sample
* NB Stata allows you to generate these directly:
predict mediand, time
su mediand
su mediand, de
drop xb mediand

streg age smoking , d(weib) nolog frailty(gamma) nohr 
streg age smoking , d(weib) nolog frailty(invgauss) nohr 


* now reintroduce dietfat -- does frailty disappear?

streg age smoking dietfat, d(weib) nolog frailty(gamma) nohr
streg age smoking dietfat, d(weib) nolog frailty(invg) nohr


streg age smoking dietfat, d(weib) nolog frailty(gamma)
streg age smoking dietfat, d(weib) nolog frailty(invg)


**** Repeat above using log-logistic model (AFT)
**** [But note that data were simulated/created assuming a Weibull model!]

streg age smoking dietfat, d(logl) nolog

* omit dietfat as a way of introducing unobserved heterogeneity
* basic model w/out dietfat and no frailty

streg age smoking , d(logl) nolog

* basic model w/out dietfat but now allowing for frailty
streg age smoking , d(logl) nolog frailty(gamma)
* streg age smoking , d(logl) nolog frailty(invgauss)


* now reintroduce dietfat -- does frailty disappear?

* models don't converge (so limit # iterations)
streg age smoking dietfat, d(logl) frailty(gamma) iter(20)
streg age smoking dietfat, d(logl) frailty(gamma) iter(30)  difficult tech(dfp 5 nr 5)

streg age smoking dietfat, d(logl) frailty(invg) iter(20)
streg age smoking dietfat, d(logl) frailty(invg) iter(30) difficult tech(dfp 5 nr 5)

* Repeat Weibull example and show that non-informative episode-splitting
*  does not affect results

use bc, clear
ge id = _n
stset t, f(dead) id(id)
streg age smoking , d(weib) nohr nolog frailty(gamma)

stsplit time, every(1)
stset
streg age smoking , d(weib) nohr nolog frailty(gamma)

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
streg , d(weib) nolog nohr frailty(invg)

********************************************************
** marriage data set, Weibull model ***
********************************************************

use duration, clear
de
ge nopart = married == 1
stset time, f(status)

* full model

streg nopart sex, d(weib) nolog nohr 
streg nopart sex, d(weib) nolog 


* model with no covariates, with and without frailty
streg  , d(weib) nolog nohr
streg  , d(weib) nolog nohr frailty(gamma)
streg  , d(weib) nolog nohr frailty(invg)




**************************************************************
************** discrete time (bc data) ***********************
**************************************************************



* beware t is a non-integer vble!

use bc, clear
su t, de

* convert survival times to discrete integers (rounding up)
* Variable td indexes intervals

ge td = ceil(t)  // has same effect as: ge td = round(t+.49,1)

su t td

sort t

ge id = _n
expand td 	// expand on td (not t) since time intervals indexed on td 
bysort id: ge newt = _n
bysort id: ge died = dead==1 & _n==_N
ge logt = ln(newt)

tab newt died

** log(time) specification for hazard, PH model ***


cloglog died logt age smoking dietfat, nolog
cloglog died logt age smoking, nolog

* (log)normal error
xtcloglog died logt age smoking , nolog i(id)

* Gamma distribution
ssc install pgmhaz8
pgmhaz8 logt age smoking , id(id) seq(newt) dead(died) nolog


pgmhaz8 logt age smoking , id(id) seq(newt) dead(died) 

/* to generate error messages

pgmhaz8 logt age smoking , id(id) seq(newt) dead(died) trace

pgm_ll does not compute a continuous nonconstant function
could not calculate numerical derivatives
r(430);

*/


* Now add dietfat back in :

pgmhaz8 logt age smoking dietfat, id(id) seq(newt) dead(died) iter(25)

** repeat for logistic hazard ***

logit died logt age smoking dietfat, nolog
logit died logt age smoking , nolog
xtlogit died logt age smoking,  nolog i(id)
xtlogit died logt age smoking dietfat,  nolog i(id)



**** Predictions from -pgmhaz-


/* Complementary log-log ('cloglog') hazard function:

     	p(t)   = 1 - exp[-L*exp(c(t))]

     	p(t)   = 1 - exp[-exp(z)]

where z = Xb + c(t).  E.g. Xb + c.log(t)

The discrete time survival function S(t)

   	S(t) = (1-p(1))*(1-p(2))*(1-p(3))*...*(1-p(t))

To calculate S(t), it is easier to re-write it as

		    s=t
	S(t) = exp{ SUM [ln(1-p(s))] }.
		    s=1

*/

* First re-do earlier model
pgmhaz8 logt age smoking , id(id) seq(newt) dead(died) nolog

* Now generate predicted probabilities
mat b = e(b)
mat list b
scalar n = colsof(b) - 1 
scalar list n

mat b = b[1,2..n]  /* exclude coeffs on dur dep var(s), and ln_varg */
mat list b
mat score xb = b
sum xb if newt == 1

ge z0 = r(mean) + _b[logt]*logt
ge z1 = _b[_cons] + _b[age]*45 + _b[smoking]*1 + _b[logt]*logt
ge z2 = _b[_cons] + _b[age]*30 + _b[smoking]*0 + _b[logt]*logt

sort id newt
by id: gen p0 = 1 - exp(-exp(z0))
lab var p0 "Predicted h(t) at mean of covariates"
by id: gen p1 = 1 - exp(-exp(z1))
lab var p1 "Predicted h(t),age=45,smoking"
by id: gen p2 = 1 - exp(-exp(z2))
lab var p2 "Predicted h(t),age=30,non-smoking"

by id: ge s0 =  exp(sum(ln(1-p0)))
lab var s0 "Predicted S(t) at mean of covariates"
by id: ge s1 =  exp(sum(ln(1-p1)))
lab var s1 "Predicted S(t),age=45,smoking"
by id: ge s2 =  exp(sum(ln(1-p2)))
lab var s2 "Predicted S(t),age=30,non-smoking"

twoway (connect p0 newt , sort  msymbol(t) ) ///
	(connect p1 newt, sort msymbol(o) ) ///
   	(connect p2 newt, sort msymbol(x) )  ///
 	, title("Predicted discrete hazard rates from -pgmhaz8-") ///
	saving(pgmh1, replace) ytitle("p(t)")

twoway (connect s0 newt , sort  msymbol(t) )   ///
	(connect s1 newt, sort msymbol(o) )   ///
	(connect s2 newt , sort  msymbol(x) )  ///
	, title("Predicted survivor functions from -pgmhaz8-") ///
	saving(pgms1, replace) ytitle("S(t)")

drop z0 z1 z2 xb p0 p1 p2 s0 s1 s2


*----------------------------------------------------
* now repeat with normal heterogeneity

xtcloglog died logt age smoking , i(id)  nolog
predict p, pu0

* Now generate predicted probabilities
mat b = e(b)
mat list b
scalar n = colsof(b) - 1 
scalar list n

mat b = b[1,2..n]  /* exclude coeffs on dur dep var(s), and lnsig2u */
		/* code relies on vbles being in correct order */
mat list b
mat score xb = b
sum xb if newt == 1
ge z0 = r(mean) + _b[logt]*logt
ge z1 = _b[_cons] + _b[age]*45 + _b[smoking]*1 + _b[logt]*logt
ge z2 = _b[_cons] + _b[age]*30 + _b[smoking]*0 + _b[logt]*logt

sort id newt
by id: gen p0 = 1 - exp(-exp(z0))
lab var p0 "Predicted h(t) at mean of covariates"
by id: gen p1 = 1 - exp(-exp(z1))
lab var p1 "Predicted h(t),age=45,smoking"
by id: gen p2 = 1 - exp(-exp(z2))
lab var p2 "Predicted h(t),age=30,non-smoking"


by id: ge s0 =  exp(sum(ln(1-p0)))
lab var s0 "Predicted S(t) at mean of covariates"
by id: ge s1 =  exp(sum(ln(1-p1)))
lab var s1 "Predicted S(t),age=45,smoking"
by id: ge s2 =  exp(sum(ln(1-p2)))
lab var s2 "Predicted S(t),age=30,non-smoking"

twoway (connect p0 newt , sort  msymbol(t) ) ///
	(connect p1 newt, sort msymbol(o) ) ///
	 (connect p2 newt, sort msymbol(x) ) ///
	, title("Predicted discrete hazard rates from -pgmhaz8-") ///
	saving(pgmh2, replace) ytitle("p(t)")

twoway (connect s0 newt , sort  msymbol(t) ) ///
	(connect s1 newt, sort msymbol(o) ) ///
	(connect s2 newt , sort  msymbol(x) ) 		///
	, title("Predicted survivor functions from -pgmhaz8-") ///
	saving(pgms2, replace) 	ytitle("S(t)")


drop z0 z1 z2 xb p0 p1 p2 s0 s1 s2

sort id newt
by id: ge s =  exp(sum(ln(1-p)))
gsort -s
list id age smoking newt s if abs(s-.5) < .05
drop s p


*----------------------------------------------------
* now repeat with H-S mass point heterogeneity
*----------------------------------------------------
hshaz logt age smoking , id(id) seq(newt) dead(died) nolog

* NB if add dietfat back in, model unable to be fitted 
* 	(initial values not feasible)
* 	r(1400);

/*
hshaz logt age smoking dietfat, id(id) seq(newt) dead(died)
*/


* Now generate predicted probabilities
mat b = e(b)
mat list b
scalar n = colsof(b) - 2 	// 2 parameters (m2, logitp2 ) 
scalar list n

mat b = b[1,2..n]  /* exclude coeffs on dur dep var(s), and m2, logitp2 */
mat list b
mat score xb = b
sum xb if newt == 1
* Type 1 person (intercept is b0 + m1 = b0)
ge z01 =  r(mean) + _b[logt]*logt
ge z11 = _b[_cons] + _b[age]*45 + _b[smoking]*1 + _b[logt]*logt
ge z21 = _b[_cons] + _b[age]*30 + _b[smoking]*0 + _b[logt]*logt

* Type 2 person (intercept is b0 + m2)
ge z02 = e(m2) +  r(mean) + _b[logt]*logt
ge z12 = e(m2) + _b[_cons] + _b[age]*45 + _b[smoking]*1 + _b[logt]*logt
ge z22 = e(m2) + _b[_cons] + _b[age]*30 + _b[smoking]*0 + _b[logt]*logt

sort id newt
by id: gen p01 = 1 - exp(-exp(z01))
lab var p01 "Type 1: Predicted h(t) at mean of covariates"
by id: gen p11 = 1 - exp(-exp(z11))
lab var p11 "Type 1: Predicted h(t),age=45,smoking"
by id: gen p21 = 1 - exp(-exp(z21))
lab var p21 "Type 1: Predicted h(t),age=30,non-smoking"

by id: gen p02 = 1 - exp(-exp(z02))
lab var p02 "Type 2: Predicted h(t) at mean of covariates"
by id: gen p12 = 1 - exp(-exp(z12))
lab var p12 "Type 2: Predicted h(t),age=45,smoking"
by id: gen p22 = 1 - exp(-exp(z22))
lab var p22 "Type 2: Predicted h(t),age=30,non-smoking"

by id: ge s01 =  exp(sum(ln(1-p01)))
lab var s01 "Type 1: Predicted S(t) at mean of covariates"
by id: ge s11 =  exp(sum(ln(1-p11)))
lab var s11 "Type 1: Predicted S(t),age=45,smoking"
by id: ge s21 =  exp(sum(ln(1-p21)))
lab var s21 "Type 1: Predicted S(t),age=30,non-smoking"

by id: ge s02 =  exp(sum(ln(1-p02)))
lab var s02 "Type 2: Predicted S(t) at mean of covariates"
by id: ge s12 =  exp(sum(ln(1-p12)))
lab var s12 "Type 2: Predicted S(t),age=45,smoking"
by id: ge s22 =  exp(sum(ln(1-p22)))
lab var s22 "Type 2: Predicted S(t),age=30,non-smoking"

* compare Type 1 and Type 2 at mean of covariates

twoway (connect p01 newt , sort  msymbol(t) ) ///
	(connect p02 newt, sort msymbol(o) ) ///
	, title("Predicted hazard rates from -hshaz-") ///
	saving(hshaz_h1, replace) ytitle("p(t)")

twoway (connect s01 newt , sort  msymbol(t) ) ///
	(connect s02 newt, sort msymbol(o) ) ///
	, title("Predicted survivor functions from -hshaz-") ///
	saving(hshaz_s1, replace) ytitle("S(t)")


capture drop z01 z11 z21 xb p01 p11 p21 ///
	s01 s11 s21 z02 z12 z22 p02 p12 p22 s02 s12 s22



************************************************
** cancer data ***


use cancer, clear
ge id = _n


recode drug 1=0 2/3=1
lab var drug "receives drug?"
lab def drug 0 "placebo" 1 "drug"
lab val drug drug

su
expand studytim
bysort id: ge t = _n
bysort id: ge dead = died == 1 & _n ==_N
ge logt = ln(t)
compress
su

tab t dead

* (log)normal PH model 

cloglog dead drug age logt, nolog
xtcloglog  dead drug age logt, nolog i(id)


* Gamma PH 

pgmhaz8  drug age logt, id(id) s(t) d(dead) nolog
pgmhaz8, eform


* logistic

logit dead drug age logt, nolog
xtlogit  dead drug age logt, nolog i(id)



