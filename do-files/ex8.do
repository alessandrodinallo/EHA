
* Independent competing risks example


************** discrete time *********************************

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

* treat single risk as 'leaving UI register'
* treat competing risks as 'leaving to a job', 'leaving to other'


* estimate a single destination PH model with log(time) hazard
cloglog leftui age famresp tyentry reg1-reg4 logt, nolog
cloglog, eform


* estimate the components of the competing risk model 

cloglog job age famresp tyentry reg1-reg4 logt, nolog
cloglog, eform

cloglog other age famresp tyentry reg1-reg4 logt, nolog
cloglog, eform


************ logit versions of above eqns; cf with mlogit

* estimate a single destination PH model with log(time) hazard
logit leftui age famresp tyentry reg1-reg4 logt, nolog

* estimate the components of the competing risk model 

logit job age famresp tyentry reg1-reg4 logt, nolog
logit, or

logit other age famresp tyentry reg1-reg4 logt, nolog
logit, or

sort newid t

* 
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

********************************************************
********** Continuous time model ***********************
********************************************************


* pretend survival times are continuous (not discrete intervals)
use unemp, clear
ta status exit

* create some covariates
* age, famresp, tyentry, and ...
ta groupreg, ge(reg)

**** Weibull model ****

* treat single risk as 'leaving UI register'
* treat competing risks as 'leaving to a job', 'leaving to other'

* single risk model

stset conmths, f(exit)
streg age famresp tyentry reg1-reg4, dist(weib) nolog nohr
local llleft = e(ll)
local n_rest = e(k)  	// # restrictions for test: see below 
di "Number restrictions  = " `n_rest'
local n_exit = e(N_fail)
di "Number exits, all kinds = " `n_exit'
streg, hr

* competing risks

* (1) exit to a job
stset conmths, fail(status==2)  
streg age famresp tyentry reg1-reg4, dist(weib) nolog nohr
local lljob = e(ll)
local n_job = e(N_fail)
di "Number exits to job = " `n_job'
streg, hr

* (1) exit to other destination
stset conmths, fail(status==3)  
streg age famresp tyentry reg1-reg4, dist(weib) nolog nohr
local lloth = e(ll)
local n_oth = e(N_fail)
di "Number exits, other kind = " `n_oth'
streg, hr

* Testing proportionality of risks of exit to job versus other exits:
*
* Narendranathan & Stewart (1991) 'Testing the proportionality of 
*  cause-specific hazards in competing risk models', Oxford Bulletin
*  of Economics and Statistics 53, 331-340
* 
* Are the exits to different states behaviourally distinct or is the
* state exited to incidental?
* Set of restrictions implied is equality of all parameters except intercepts
* [In Weibull model, # restrictions = #covariates -1 (intercept) + 1 (shape)]
*
* test statistic = 2[ln(L_CR) - ln(L_SR) - {sum_j (n_j)*ln(p_j)} ]
*			where p_j = n_j/sum_k{n_k}
*			n_j = # exiting to state j	
*
*		d.of f. = # restrictions
*
* NB their particular test was developed in terms of a continuous time PH model 
*    (so can't use for discrete!)

local test = 2*(`lljob' + `lloth' - `llleft' ///
	- (`n_job')*ln(`n_job'/`n_exit') - (`n_oth')*ln(`n_oth'/`n_exit') )

* test statistic is Chi-sq with dof = # restrictions 

di "The test statistic is " `test' 
di "with p-value = " chiprob(`n_rest', `test')

