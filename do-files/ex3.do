

************** Cancer data ********************

 /* sysuse cancer, clear */ // You may also decide to use the dataset "cancer" that is provided by default by STATA
use cancer, clear 
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

* Compare this data format with the earlier one, taking the same four persons
list id studytim seqvar died dead age if id <= 4

de  /* There should be 744 obs in the data set now */
su

stset seqvar, failure(dead) id(id)
stsum
de _*
su _*


ge logd = ln(seqvar)
lab var logd "ln(t)"


*** using stsplit to accomplish same results

sysuse cancer, clear
ge id = _n
* Use of -stsplit- later requires that we set id(.)
stset studytim, f(died) id(id) 
ge T = studytim

stsplit time, every(1) 
stset
sort id time
ge seqvar = time + 1
by id: list id died T studytim time seqvar , noobs
stsum

*** Now treat the survival times as as continuous
*** to illustrate time-varying covariate creation
*** -- episode split for piece-wise exponential model

sysuse cancer, clear
ge id = _n
stset studytim, f(died) id(id) 
de _*
su _*

* Vble to pick out 3 subjects with illustrative times
ge pick = id == 3 | id == 17 | id == 48 
sort id
list id studytim _t died _d _t0 if pick

* split at survival times 8 and 17 
* ... key intervals are (0,8], (8,17], (17,infinity]

stsplit ehaz, at(8 17)
ta ehaz, ge(e)
sort id ehaz
list id studytim _t died _d _t0 e* if pick


*** Episode-splitting for the Cox model for discrete TVCs
*** -- split at failure times

use cancer, clear
ge id = _n
stset studytim, f(died)  id(id) 
de _*
su _*

stsplit , at(failures)
de _*
su _*

sort id studytim
list id studytim _t died _d _t0  if (id==3|id==17|id==48), noobs

* Now you would manually generate the relevant TVCs


***************** Strike data **********************
use kennan, clear
de
stset time, f(status)
stsum
di r(risk)
di r(N_sub)
local nrisk = r(risk)  /* see below (some fancy stuff using local macros) */
local nsub = r(N_sub)

ge id = _n
lab var id "strike identifier"

* Look at the data format for the first four obs
list id time status in 1/4

* Now expand the data set so that one obs per month at risk of 'failure'
expand time   
bysort id: ge seqvar = _n  
lab var seqvar "spell month identifier, by strike"
quietly by id: ge dead = status & _n==_N
lab var dead "binary depvar for discrete hazard model"

* Compare this data format with the earlier one, same four obs
list id time seqvar status dead if id <= 4
di "Should now be `nrisk' spell months in total, for `nsub' subjects"
de  
stset seqvar, f(dead) id(id)
stsum
ge logd = ln(seqvar)
lab var logd "ln(t)"
su

********* Marriage data ******************************
use duration, clear
de
stset time, f(status)
stsum
di r(risk)
di r(N_sub)
local nrisk = r(risk)  /* see below (some fancy stuff using local macros) */
local nsub = r(N_sub)

ge id = _n
lab var id "person identifier"
* Look at the data format for the first four obs
list id time status in 1/4
* Now expand the data set so that one obs per month at risk of 'failure'
expand time   
bysort id: ge seqvar = _n  
lab var seqvar "spell month identifier, by person"
bysort id: ge dead = status & _n==_N
lab var dead "binary depvar for discrete hazard model"
* Compare this data format with the earlier one, same four obs
list id time seqvar status dead if id <= 4
di "Should now be `nrisk' spell months in total, for `nsub' subjects"
de  
stset seqvar dead, id(id)
stsum
ge logd = ln(seqvar)
lab var logd "ln(t)"
su

******************* HIV data - spell date exercise ******

use hmohiv, clear
de, de
ge start = date(entdate, "DM19Y")
ge end = date(enddate, "DM19Y")
format start %td
format end %td
ge T = end - start
lab var T "Survival time (days)"
su
list in 1/5, noobs
* survival time measured in months
stset time, f(censor)
stsum
* survival time measured in days
stset T, f(censor)
stsum


**** Unemployment data and TVC ****

* episode splitting using -expand-
use unemp, clear
stset conmths, fail(exit) id(newid)
stsum
expand conmths
bysort newid: ge t = _n
bysort newid: ge ended = exit==1 & _n==_N
stset t, f(ended) id(newid)
stsum
* generate TVC
ge netrr = .
replace netrr = rn1 if t < 7
replace netrr = rn2 if t >= 7 & t < 13
replace netrr = rn3 if t >=13 & t ~=.
su
de, de /* look at data set size */

* now use -stplit- to do the organisation instead
use unemp, clear
stset conmths, fail(exit) id(newid)
stsplit time, at(1(1)24)
stset
stsum
* generate TVC
ge netrr = .
replace netrr = rn1 if _t < 7
replace netrr = rn2 if _t >= 7 & _t < 13
replace netrr = rn3 if _t >=13 & _t ~=.
su
de, de // look at data set size 

* list the data, seeing how the new TVC (netrr) varies with survival time values

list newid time _t netrr in 1/100, sepby(newid) noobs

* now pretend that conmths is a continuous time vble
* and repeat exercise creating TVC.
* Values of replacement rate change at 7 and 13
* so split there

use unemp, clear
stset conmths, fail(exit) id(newid)
stsplit rrtime, at(7 13)
stset
stsum
* generate TVC
ta rrtime, ge(period)  // creates 3 new dummy variables named period1, period2, period3
de period*
ge netrr = .
replace netrr = rn1 if period1==1
replace netrr = rn2 if period2==1
replace netrr = rn3 if period3==1
su
de, de // look at data set size 

list newid rrtime _t netrr period* in 1/100, sepby(newid) noobs


