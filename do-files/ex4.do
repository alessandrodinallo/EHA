
************** Survivor & hazard functions ************


************* Cancer data **********
sysuse cancer, clear
su
de
stset studytim , failure(died)
ge id = _n  
lab var id "subject identifier"
recode drug 1=0 2/3=1
lab var drug "receives drug?"
lab def drug 0 "placebo" 1 "drug"
lab val drug drug

sts list

sts graph, title("Survivor function, Cancer data (sts)") ///
	xtick(1(1)39) saving(surv1, replace)
sts graph, title("Cumulative hazard function, Cancer data (sts)") ///
	na xtick(1(1)39) saving(integh1, replace)

* smoothed hazard, default bandwidth (and default kernel)
sts graph, hazard title("Smoothed hazard function, Cancer data (sts)") ///
	xtick(1(1)39) saving(smoothhaz1, replace)

* smaller width
sts graph, hazard width(2) title("Smoothed hazard function, Cancer data (sts)") xtick(1(1)39)

sts gen s = s
lab var s "KM survivor function, from -sts gen-"
sts gen cumh = na
lab var cumh "NA cumulative hazard function, from -sts gen-"
sts gen deltach = h 
lab var deltach "Hazard contribution, from -sts gen-"

twoway line s _t, sort connect(J) title("Survivor function, Cancer data (sts gen)") ///
	ytitle("S(t)") xtitle("survival time, t") xtick(1(1)39) saving(surv2, replace) 

twoway line cumh _t, sort connect(J) title("Cumulative hazard function, Cancer data (sts gen)") ///
	ytitle("H(t)") xtitle("survival time, t") xtick(1(1)39) saving(chaz1, replace) 

twoway line deltach _t, sort connect(J) title("Hazard contribution, Cancer data (sts gen)") ///
	ytitle("deltaH(t)") xtitle("survival time, t") xtick(1(1)39) saving(deltach1, replace) 

* Inspecting the numbers themselves:
* the next 2 lines provide way of selecting just one value of t for the -list-ing
* (remember that _t is a synonym for studytim now that data are -stset- )


egen tagt = tag(_t) 
sort _t
list _t deltach cumh s if tagt == 1


****************************
**** stratified estimates ***

sts graph, title("Survivor functions, by drug, Cancer data (sts)") ///
	by(drug) lost saving(surv5, replace)

sts list, by(drug) 
sts list, by(drug) compare
sts test drug
sts test drug, wilcoxon


**** check plausibility of PH, Weibull, log-logistic models, using graphs
**** "parallel lines"?

* PH, by drug receipt?

	
sts gen cumhg = na, by(drug)
ge lch = log(cumhg)
separate lch, by(drug) 	// creates lch0 lch1
de lch*


twoway connect lch0 lch1 _t, ytitle("log(cumulative hazard)") ///
	title("Proportional hazard check - parallel lines?") ///
	xtitle("survival time") xtick(1(1)39) sort saving(PHtest1, replace)

* Weibull? How straight & parallel are the lines?

ge logt = log(_t)
twoway connect lch0 lch1 logt, ytitle("log(cumulative hazard)") ///
	xtitle("log(survival time)") ///
	title("Weibull model check - parallel straight lines?") ///
	sort saving(Weibtest1, replace)

* Log-logistic?  How straight & parallel are the lines?


sts gen sg = s, by(drug) 
ge sls = log(sg/(1-sg))
separate sls, by(drug)   // creates sls0 sls1
de sls*


twoway connect sls0 sls1 logt, ytitle("log[S/(1-S)]") xtitle("log(survival time)")  ///
	title("Log-logistic model check - parallel straight lines?") sort saving(LogLtest1, replace)


**************************************
******** lifetable estimates *********
**************************************

use cancer, clear
su
de
stset studytim , failure(died)
ge id = _n  
lab var id "subject identifier"
recode drug 1=0 2/3=1
lab var drug "receives drug?"
lab def drug 0 "placebo" 1 "drug"
lab val drug drug


ltable studytim died


* To save the survivor function graph produced by -ltable-, use
* 	 graph save [graphname] filename [, asis replace ]
* after running the command

ltable studytim died, graph title("Survivor function, Cancer data (ltable)")
graph save surv3, replace

* In Stata 8, -ltable- will not graph the hazard 
* but you can get it via the -saving- option which saves the life table data to a file, 
*	or you can use -sts- instead (as we did earlier)
* Here's how the former method works


ltable studytim died, hazard saving(surv3a, replace) 

* graphs of the hazard + CIs
preserve
use surv3a, clear
de // see the variables corresponding to the table
* simple graph
twoway connect lhazard hazard uhazard t0
* fancier:
twoway (rspike uhazard lhazard t0) (connect hazard t0), ytitle(hazard) ///
	title("Hazard estimate from -ltable-, with 95% CI") ///
	xtitle("survival time") xtick(1(1)39) legend(off) saving(ltable-haz, replace)
restore



* Alternatively use version 7
version 7: ltable studytim died, haz graph saving(haz3a, replace)

*** repeat using noadjust option ********************

ltable studytim died, noadjust 
ltable studytim died, hazard noadjust 

ltable studytim died, graph title("Survivor function, Cancer data (ltable)")  noadjust 
graph save surv4, replace

	

ltable studytim died, hazard saving(surv3a-noa, replace) 

* graphs of the hazard + CIs

preserve
use surv3a-noa, clear
de
twoway (rspike uhazard lhazard t0) (connect hazard t0), ytitle(hazard) ///
	title("Hazard estimate from -ltable-, with 95% CI") ///
	xtitle("survival time") xtick(1(1)39) legend(off) saving(ltable-haz-noa, replace)
restore


***********************************
* now confirm that you get the same results from data in
* person-month format

expand studytim   
bysort id: ge seqvar = _n   
lab var seqvar "spell month identifier, by subject"
bysort id: ge dead = (died==1)  & _n==_N
lab var dead "binary depvar for discrete hazard model"


stset seqvar dead, id(id)
stsum
sts
sts list
ltable studytim died, gr tvid(id)	
ltable studytim died, hazard tvid(id)



*******&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
************* Strike data **********

use kennan, clear
su
de
stset time , failure(status)
ge id = _n  
lab var id "subject identifier"

sts list
sts graph, title("Survivor function, Strike data (sts)") saving(surv6, replace)
sts graph, title("Cumulative hazard function, Strike data (sts)") ///
	na saving(integh2, replace)

sts gen s = s
lab var s "KM survivor function from -sts gen-"
sts gen cumh = na
lab var cumh "NA cumulative hazard function, from -sts gen-"
sts gen deltach = h 
lab var deltach "Hazard contribution, from -sts gen-"

twoway line s _t, sort connect(J) title("Survivor function, Strike data (sts gen)") ///
	saving(surv7, replace) 

twoway line cumh _t, sort connect(J) title("Cumulative hazard function, Strike data (sts gen)") ///
	saving(chaz4, replace) 

twoway line deltach _t, sort connect(J) title("Hazard contribution, Strike data (sts gen)") ///
	saving(deltach2, replace) 

// listing the data to see the estimates at each value of t

list id _t status deltach cumh s

// ... but notice that this command produces, for each t, repeated values for individuals
//	-- what we really want is one estimate listed for each t value. Here's how:

sort _t
egen first = tag(_t)
list id _t status deltach cumh s if first == 1


ltable time status
ltable time status, graph title("Survivor function, Strike data (ltable)") noadjust
graph save surv8, replace 
ltable time status, hazard   noadjust


************* Marriage data **********

use duration, clear
su
de
stset time , failure(status)
ge id = _n  
lab var id "subject identifier"

sts list
sts graph, b1("Survivor function, Marriage data (sts)") saving(surv9, replace)
sts graph, b1("Cumulative hazard function, Marriage data (sts)") ///
	na saving(integh3, replace)

sts gen s = s
lab var s "KM survivor function from -sts gen-"
sts gen cumh = na
lab var cumh "NA cumulative hazard function, from -sts gen-"
sts gen deltach = h 
lab var deltach "Hazard contribution, from -sts gen-"


twoway line s _t, sort connect(J) title("Survivor function, Marriage data (sts gen)") ///
	saving(surv10, replace) 

twoway line cumh _t, sort connect(J) title("Cumulative hazard function, Marriage data (sts gen)") ///
	saving(chaz5, replace) 

twoway line deltach _t, sort connect(J) title("Hazard contribution, Marriage data (sts gen)") ///
	saving(deltach3, replace) 


list id _t status deltach cumh s

ltable time status
ltable time status, graph b1("Survivor function, Marriage data (ltable)") 
graph save surv11, replace
ltable time status, hazard 

sts graph, title("Survivor functions, by sex, Marriage data (sts)") ///
	by(sex) lost saving(surv12, replace)
sts list, by(sex) 
sts list, by(sex) compare

sts test sex
sts test sex, wilcoxon


