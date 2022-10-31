
************** Estimation: (i) discrete time models *********


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



/*
(NB) Fully non-parametric baseline hazard case.
 
-- essential to check whether events occur at each value of 
   t (here 't' or _t).
   The hazard cannot be estimated for values of t with no events
   (cf. non-parametric baseline hazard in Cox model).

For the duration intervals with no events, then either
one must refine the grouping on the duration dimension (a piece-wise
constant baseline) is an example of this or one must drop the relevant
person months from the estimation. (Cf. the discussion of identification
of the logit model in Reference Manuals under "perfect predictors")

   Check like this:
*/

ta j dead

/*
There are no deaths during months 9, 14, 18-21, 26, 27, 29-32, 34-39,
and so a month-specific hazard rate cannot be estimated for these
intervals.
*/


/**********************************************************
****************** CLOGLOG HAZARD MODELS *****************
* Compare model estimated with different baseline hazard specifications.
* Use -predict- to derive estimate of predicted hazard and survivor function
* and thence median duration.  First use within-sample info.


* log(j) baseline [and 'or' option; logit versus logistic]

* cloglog = glm, f(b) l(c). Can also use -glm- .
* See help -glm- and note glm Deviance = -2*LogL from cloglog

glm dead drug age lnj, f(b) l(c)  
glm, eform

cloglog dead drug age lnj, nolog
predict h, p

cloglog, eform	// replay results, but this time with hazard ratio



bysort id (j): ge s = exp(sum(ln(1-h)))

ge h0 = h if age == 55 & drug == 0
ge h1 = h if age == 55 & drug == 1
lab var h0 "drug = 0"
lab var h1 "drug = 1"
twoway (connect h0 j, sort  msymbol(t) ) (connect h1 j, sort msymbol(o) ) ///
	, xlabel(0(10)50) title("p(j),c(j)=(q-1)ln(j),age=55") saving(dclog1, replace) 

ge s0 = s if age == 55 & drug == 0
ge s1 = s if age == 55 & drug == 1
lab var s0 "drug = 0"
lab var s1 "drug = 1"
twoway (connect s0 j, sort  msymbol(t) ) (connect s1 j, sort msymbol(o) ) ///
	, xlabel(0(10)50) title("S(j),c(j)=(q-1)ln(j),age=55") saving(dclog2, replace) 

drop h h0 h1 s s0 s1

* cubic polynomial
cloglog dead drug age j j2 j3,  nolog
cloglog, eform

cloglog dead drug age lnj, nolog
predict h, p
cloglog, eform

bysort id (j): ge s = exp(sum(ln(1-h)))

ge h0 = h if age == 55 & drug == 0
ge h1 = h if age == 55 & drug == 1
lab var h0 "drug = 0"
lab var h1 "drug = 1"
twoway (connect h0 j, sort  msymbol(t) ) (connect h1 j, sort msymbol(o) ) ///
	, xlabel(0(10)50) title("p(j),c(j) cubic polyn.,age=55") saving(dclog3, replace) 


ge s0 = s if age == 55 & drug == 0
ge s1 = s if age == 55 & drug == 1
lab var s0 "drug = 0"
lab var s1 "drug = 1"
twoway (connect s0 j, sort  msymbol(t) ) (connect s1 j, sort msymbol(o) ) ///
	, xlabel(0(10)50) title("S(j),c(j) cubic polyn.,age=55") saving(dclog4, replace) 

drop h h0 h1 s s0 s1


* piece-wise constant baseline
cloglog dead dur1-dur6 drug age, nocons  nolog
cloglog dead dur2-dur6 drug age, nolog

cloglog dead drug age e2 e3, nolog
cloglog, eform

cloglog dead drug age Late, nolog
predict h, p
cloglog, eform



bysort id (j): ge s = exp(sum(ln(1-h)))

ge h0 = h if age == 55 & drug == 0
ge h1 = h if age == 55 & drug == 1
lab var h0 "drug = 0"
lab var h1 "drug = 1"
twoway (connect h0 j, sort  msymbol(t) ) (connect h1 j, sort msymbol(o) ) ///
	, xlabel(0(10)50) title("p(j),c(j) piece-wise cnst.,age=55") saving(dclog5, replace) 

ge s0 = s if age == 55 & drug == 0
ge s1 = s if age == 55 & drug == 1
lab var s0 "drug = 0"
lab var s1 "drug = 1"
twoway (connect s0 j, sort  msymbol(t) ) (connect s1 j, sort msymbol(o) ) ///
	, xlabel(0(10)50) title("S(j),c(j) piece-wise cnst.,age=55") saving(dclog6, replace) 


drop h h0 h1 s s0 s1


* fully non-parametric baseline
cloglog dead drug age d1-d39, nocons nolog
cloglog dead drug age d1-d8 d10-d13 d15-d17 d22-d25 d28 d33 /*
      */ if (j>=1 & j<=8) | (j>=10 & j<=13) /*
      */ | (j>=15 & j<=17) | (j>=22 & j<=25) /*
      */ | j==28 | j==33 , nocons nolog

* use -predict- to derive estimate of predicted hazard and survivor function
* and thence median duration

predict h, p
replace h = 0 if j==9|j==14|(j>=18 & j <=21)|j==26|j==27|(j>=29&j<=32)|j>=34


bysort id (j): ge s = exp(sum(ln(1-h)))

ge h0 = h if age == 55 & drug == 0
ge h1 = h if age == 55 & drug == 1
lab var h0 "drug = 0"
lab var h1 "drug = 1"
twoway (connect h0 j, sort  msymbol(t) ) (connect h1 j, sort msymbol(o) ) ///
	, xtick(1(1)39) xlabel(0(10)50) title("p(j),c(j) fully non-par.,age=55") saving(dclog7a, replace) 

/*
* (2) redo model assuming piece-wise constant to 'cover' unidentified bits
drop h h0 h1 s

cloglog dead drug age d1-d7 d89 d10-d12 d1314 d15 d16 d1721 d22-d24 d2527 d2832 d3339, nocons

predict h, p
bysort id (j): ge s = exp(sum(ln(1-h)))
*/

** (3) redo model assuming more grouping (piece-wise constant) to 'cover' unidentified bits
drop h h0 h1 s

cloglog dead drug age dur1 dur2 dur3 dur4 dur5 dur6, nocons nolog

predict h, p
bysort id (j): ge s = exp(sum(ln(1-h)))


ge h0 = h if age == 55 & drug == 0
ge h1 = h if age == 55 & drug == 1
lab var h0 "drug = 0"
lab var h1 "drug = 1"
twoway (connect h0 j, sort  msymbol(t)  connect(J)) 	///
	(connect h1 j, sort msymbol(o)  connect(J)) 	///
	, xtick(1(1)39) title("p(j),c(j) piecewise-constant,age=55") ///
	saving(dclog7, replace) 

ge s0 = s if age == 55 & drug == 0
ge s1 = s if age == 55 & drug == 1
lab var s0 "drug = 0"
lab var s1 "drug = 1"
twoway (connect s0 j, sort  msymbol(t) ) (connect s1 j, sort msymbol(o) ) ///
	, xtick(1(1)39) title("S(j),c(j) piecewise-constant.,age=55") saving(dclog8, replace) 


drop h h0 h1 s s0 s1 

*************************************************************
* Redo using articial observations and out-of-sample prediction
*************************************************************


*** drug == 0, age == mean
local newn = _N + 50
su id
local idmax = r(max)

set obs `newn'

cloglog dead drug age lnj, nolog

replace id = `idmax' + 1 if id==.
sort id

bysort id (j): replace j = _n if id==(`idmax' + 1)
replace lnj = ln(j) if  id==(`idmax' + 1)
replace drug = 0 if  id==(`idmax' + 1)


su age if j==1  /* one obs per person */
local meana = r(mean)
di "Average age =  "  `meana'

replace age = `meana' if id==(`idmax' + 1)  

predict h0 if id==(`idmax' + 1), p
lab var h0 "drug = 0"


bysort id (j): ge s0 = exp(sum(ln(1-h0))) if id==(`idmax' + 1)
lab var s0 "drug = 0"
twoway connect h0 j, sort  msymbol(t)  saving(dclog9, replace) 
twoway connect s0 j, sort  msymbol(t)  saving(dclog10, replace) 


* sample size right?
su id drug age j lnj studytim

*** now repeat for drug == 1, age == mean

local newn = _N + 50
su id
local idmax = r(max)

set obs `newn'


replace id = `idmax' + 1 if id==.
sort id
bysort id (j): replace j = _n if id==(`idmax' + 1)
replace lnj = ln(j) if id==(`idmax' + 1) 
replace drug = 1 if id==(`idmax' + 1) 
replace age = `meana' if id==(`idmax' + 1) 

predict h1 if  id==(`idmax' + 1), p
lab var h1 "drug = 1"


bysort id (j): ge s1 = exp(sum(ln(1-h1))) if id==(`idmax' + 1) 
lab var s1 "drug = 1"

twoway connect h0 j, sort  msymbol(t)  saving(dclog11, replace) 
twoway connect s0 j, sort  msymbol(t)  saving(dclog12, replace) 


* sample size right?
su id drug age j lnj studytim
twoway (connect h0 j, sort  msymbol(t) ) (connect h1 j, sort msymbol(o) ) ///
	, xlabel(0(10)50) title("p(j),c(j)=(q-1)ln(j),age=mean") saving(dclog13, replace) 
twoway (connect s0 j, sort  msymbol(t) ) (connect s1 j, sort msymbol(o) ) ///
	, xlabel(0(10)50) title("S(j),c(j)=(q-1)ln(j),age=mean") saving(dclog14, replace) 



sort j
list id j s0 s1 if ((abs(s0-.5) < .1)|(abs(s1-.5) < .1)) & id > `idmax'-1

drop h0 h1 s0 s1
drop if id > `idmax'-1
graph combine dlogit13.gph dclog13.gph dlogit14.gph  dclog14.gph  ///
	, title("logistic (LHS) versus cloglog(RHS)") saving(dclog15,replace)



****************** LOGISTIC HAZARD MODELS *****************/
* Compare model estimated with different baseline hazard specifications.
* Use -predict- to derive estimate of predicted hazard and survivor function
* and thence median duration.  First use within-sample info.

* log(j) baseline [and 'or' option; logit versus logistic]

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
	, title("Time = log(t)") name(lnj, replace) saving(dlogit1, replace) 

ge s0 = s if age == 55 & drug == 0
ge s1 = s if age == 55 & drug == 1
lab var s0 "drug = 0"
lab var s1 "drug = 1"
twoway (connect s0 j if age == 55, sort  msymbol(t) ) (connect s1 j, sort msymbol(o) ) ///
	, title("S(j),c(t)=(q-1)ln(j),age=55") saving(dlogit2, replace) 


drop h h0 h1 s s0 s1

* cubic polynomial
logit dead drug age j j2 j3, nolog

predict h, p


bysort id (j): ge s = exp(sum(ln(1-h)))

ge h0 = h if age == 55 & drug == 0
ge h1 = h if age == 55 & drug == 1
lab var h0 "drug = 0"
lab var h1 "drug = 1"
twoway (connect h0 j, sort  msymbol(t) ) (connect h1 j, sort msymbol(o) ) ///
	, title("h(j),c(j) cubic polyn.,age=55") saving(dlogit3, replace) 

twoway (connect h0 j, sort  msymbol(t) ) (connect h1 j, sort msymbol(o) ) ///
	, title("Time = cubic polyn.") name(cub, replace) saving(cub, replace) 

	
ge s0 = s if age == 55 & drug == 0
ge s1 = s if age == 55 & drug == 1
lab var s0 "drug = 0"
lab var s1 "drug = 1"
twoway (connect s0 j, sort  msymbol(t) ) (connect s1 j, sort msymbol(o) ) ///
	, title("S(j),c(j) cubic polyn.,age=55") saving(dlogit4, replace) 


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
twoway (connect h0 j, sort  msymbol(t) ) (connect h1 j, sort msymbol(o) ) ///
	, title("h(j),c(j) piecewise constant,age=55") saving(dlogit5, replace) 

twoway (connect h0 j, sort  msymbol(t) ) (connect h1 j, sort msymbol(o) ) ///
	, title("Time = piece-wise constant") name(piecewise, replace) saving(dlogit5, replace) 

	
ge s0 = s if age == 55 & drug == 0
ge s1 = s if age == 55 & drug == 1
lab var s0 "drug = 0"
lab var s1 "drug = 1"
twoway (connect s0 j, sort  msymbol(t) ) (connect s1 j, sort msymbol(o) ) ///
	, title("S(j),c(j) piecewise constant,age=55") saving(dlogit6, replace) 

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
	, xtick(1(1)39) title("Time = fully non-parametric") name(nonpar, replace) saving(dlogit7a, replace) 


twoway (connect h0 j, sort  msymbol(t) ) (connect h1 j, sort msymbol(o) ) ///
	, xtick(1(1)39) title("p(j),c(j) fully non-par.,age=55") saving(dlogit7a, replace) 


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
twoway (connect h0 j, sort  msymbol(t) connect(J) ) ///
	(connect h1 j, sort msymbol(o) connect(J) ) ///
	, xtick(1(1)39) title("p(j),c(j) piecewise-constant,age=55") ///
	saving(dlogit7, replace) 
	
twoway (connect h0 j, sort  msymbol(t) ) (connect h1 j, sort msymbol(o) ) ///
	, xtick(1(1)39) title("Time = piece-wise constant (more groups)") name(piecewise2, replace) saving(dlogit7a, replace) 
	
* Combine log and cubic approximation of hazard
graph combine piecewise piecewise2	

ge s0 = s if age == 55 & drug == 0
ge s1 = s if age == 55 & drug == 1
lab var s0 "drug = 0"
lab var s1 "drug = 1"
twoway (connect s0 j, sort  msymbol(t) ) (connect s1 j, sort msymbol(o) ) ///
	, xtick(1(1)39) title("S(j),c(j) piecewise-constant,age=55") saving(dlogit8, replace) 


drop h h0 h1 s s0 s1 



*************************************************************
* Redo using articial observations and out-of-sample prediction
*************************************************************

*** drug == 0, age == mean
local newn = _N + 50
su id
local idmax = r(max)

set obs `newn'

logit dead drug age lnj

replace id = `idmax' + 1 if id==.
sort id

bysort id (j): replace j = _n if id==(`idmax' + 1)
replace lnj = ln(j) if  id==(`idmax' + 1)
replace drug = 0 if  id==(`idmax' + 1)


su age if j==1  /* one obs per person */
local meana = r(mean)
di "Average age =  "  `meana'

replace age = `meana' if id==(`idmax' + 1)  

predict h0 if id==(`idmax' + 1), p
lab var h0 "drug = 0"


bysort id (j): ge s0 = exp(sum(ln(1-h0))) if id==(`idmax' + 1)
lab var s0 "drug = 0"

twoway connect h0 j, sort  msymbol(t)  saving(dlogit9, replace) 
twoway connect s0 j, sort  msymbol(t)  saving(dlogit10, replace) 

* sample size right?
su id drug age j lnj studytim

*** now repeat for drug == 1, age == 55

local newn = _N + 50
su id
local idmax = r(max)

set obs `newn'

replace id = `idmax' + 1 if id==.
sort id
bysort id (j): replace j = _n if id==(`idmax' + 1)
replace lnj = ln(j) if id==(`idmax' + 1) 
replace drug = 1 if id==(`idmax' + 1) 
replace age = `meana' if id==(`idmax' + 1) 

predict h1 if  id==(`idmax' + 1), p
lab var h1 "drug = 1"


bysort id (j): ge s1 = exp(sum(ln(1-h1))) if id==(`idmax' + 1) 
lab var s1 "drug = 1"
twoway connect h0 j, sort  msymbol(t)  saving(dlogit11, replace) 
twoway connect s0 j, sort  msymbol(t)  saving(dlogit12, replace) 

* sample size right?
su id drug age j lnj studytim

twoway (connect h0 j, sort  msymbol(t) ) (connect h1 j, sort msymbol(o) ) ///
	, xlabel(0(10)50) title("p(j),c(j)=(q-1)ln(j),age=mean") saving(dlogit13, replace) 

twoway (connect s0 j, sort  msymbol(t) ) (connect s1 j, sort msymbol(o) ) ///
	, xlabel(0(10)50) title("S(j),c(j)=(q-1)ln(j),age=mean") saving(dlogit14, replace) 

sort j
list id j s0 s1 if ((abs(s0-.5) < .1)|(abs(s1-.5) < .1)) & id > `idmax'-1

drop h0 h1 s0 s1
drop if id > `idmax'-1




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
