

** (A) CONTINUOUS TIME MODELS
** ==========================
clear all
/* 
Weibull model
-------------
hazard function: h(t) = a*L*t^(a-1) where L = exp(Xb)
survivor function: S(t) = exp(-L*t^a)
failure function: F(t) = 1 - S(t)
integrated hazard function: H(t) = -ln[S(t)] = L*t^a

Parameters alpha (a) and lambda (L) = exp(Xb).

NB Exponential model is the case when a=1
*/

* Weibull hazard: lambda=1, alpha varying

* lambda = 1, alpha = .5
twoway  function y = .5*x^(-.5), range(0 5)

* lambda = 1, alpha = 1.5
twoway  function y = 1.5*x^(.5), range(0 5)

* lambda = 1, alpha = 1
twoway  function y = 1*x^(0), range(0 5)

* lambda = 1, alpha = 2
twoway  function y = 2*x, range(0 5)

* Produce graph in which overlay all on top of each other

twoway  (function y = .5*x^(-.5), range(0 5) yvarlab("a=.5") )   ///
     	( function y = 1.5*x^(.5), range(0 5) yvarlab("a=1.5") ) ///
     	( function y = 1*x^(0), range(0 5) yvarlab("a=1") )	///
     	( function y = 2*x, range(0 2) yvarlab("a=2") )			///
	, saving(weib1, replace)				///
	title("Weibull hazard: lambda=1, alpha varying")	///
	ytitle(hazard) xtitle(t)	


* Weibull hazard: lambda varying, alpha =.5

* lambda = 2, alpha = .5
twoway  function y = x^(-.5), range(0 5)

* lambda = 1, alpha = .5
twoway  function y = .5*x^(-.5), range(0 5)

* lambda = .5, alpha = .5
twoway  function y = .25*x^(-.5), range(0 5)	

* Produce graph in which overlay all on top of each other

twoway  ( function y = x^(-.5), range(0 5) yvarlab("L=2") )	 ///
	( function y = .5*x^(-.5), range(0 5) yvarlab("L=1") ) ///
     	( function y = 0.25*x^(-.5), range(0 5) yvarlab("L=.5") )	///
	, saving(weib2, replace)				///
	title("Weibull hazard: lambda varying, alpha = 0.5")	///
	ytitle(hazard) xtitle(t)	
	



* Weibull hazard: lambda varying, alpha =1.5

* lambda=2,alpha=1.5
twoway  function y = 2*1.5*x^(.5), range(0 5)

* lambda=1,alpha=1.5
twoway  function y = 1.5*x^(.5), range(0 5)

* lambda=.5,alpha=1.5
twoway  function y = .5*1.5*x^(.5), range(0 5)

twoway  ( function y = 2*1.5*x^(.5), range(0 5) yvarlab("L=2") ) ///
	( function y = 1.5*x^(.5), range(0 5) yvarlab("L=1.5") )  ///
     	( function y = .5*1.5*x^(.5), range(0 5) yvarlab("L=.5")) ///
	, saving(weib3, replace)				///
	title("Weibull hazard: lambda varying, alpha = 1.5")	///
	ytitle(hazard) xtitle(t)	

* Weibull survivor function: lambda=1, alpha varying


* lambda=1,alpha=.5
twoway function y = exp(-x^.5), range(0 5)

* lambda=1,alpha=1.5
twoway function y = exp(-x^1.5), range(0 5)

* lambda=1,alpha=1
twoway function y = exp(-x), range(0 5)

* lambda=1,alpha=2
twoway function y = exp(-x^2), range(0 5)


twoway  ( function y = exp(-x^.5), range(0 5) yvarlab("a=.5") )		///
    	( function y = exp(-x), range(0 5) yvarlab("a=1") )			///
	( function y = exp(-x^1.5), range(0 5) yvarlab("a=1.5") )		///
     	( function y = exp(-x^2), range(0 5) yvarlab("a=2") )			///
	, saving(weib4, replace)				///
	title("Weibull survivor function: lambda = 1, alpha varying")	///
	ytitle(survival probability) xtitle(t)	


* Weibull integrated hazard function: lambda=1, alpha varying (so H(t) = t^a in this case)
twoway  ( function y = x^.5, range(0 5) yvarlab("a=.5") )		///
    	( function y = x, range(0 5) yvarlab("a=1") )			///
	( function y = x^1.5, range(0 5) yvarlab("a=1.5") )		///
     	( function y = x^2, range(0 5) yvarlab("a=2") )			///
	, saving(weib5, replace)				///
	title("Weibull integrated hazard function: lambda = 1, alpha varying")	///
	ytitle("Integrated hazard") xtitle(t)	



/* 
Weibull model: median and mean duration
---------------------------------------

Median duration t' satisfies S(t') = 0.5

2 approximate derivations: 
(i) read off survivor function at S(t) = .5
(ii) use -list- and interpolate
Else derive exactly using -display- and inverting the survivor function:
     median  t' = [(1/L)*{ln(2)}]^(1/a)

NB for other percentiles in addition to median = p50, substitute relevant
   values for 0.5 in formula above

Mean duration t* = G(1+(1/a))*(1/L)^(1/a) 
			where G(x) is the Gamma function
			i.e. exp(lngamma(x)) in Stata


First create an artificial data set to represent
'time', t: 100 equally-spaced points between t=0 and t=5 
*/
set obs 100
ge t = (_n-1)/20
su t
compress
* Weibull survivor function: lambda=1, alpha varying
ge sw1 = exp(-t^.5)
lab var sw1 "lambda=1,alpha=.5"
ge sw2 = exp(-t^1.5)
lab var sw2 "lambda=1,alpha=1.5"
ge sw3 = exp(-t)
lab var sw3 "lambda=1,alpha=1"
ge sw4 = exp(-t^2)
lab var sw4 "lambda=1,alpha=2"

* Median: derive by interpolation
di "lambda=1,alpha=.5 (sw1)"
list t sw1 if abs(sw1-.5) < .1
di "lambda=1,alpha=1.5 (sw2)"
list t sw2 if abs(sw2-.5) < .1
di "lambda=1,alpha=1 (sw3)"
list t sw3 if abs(sw3-.5) < .1
di "lambda=1,alpha=2 (sw4)"
list t sw4 if abs(sw4-.5) < .1

*  Median: derive using exact formula
di "median duration for lambda=1, alpha=.5 is "  (ln(2))^(1/.5)
di "median duration for lambda=1, alpha= 1 is "  (ln(2))^(1/1)
di "median duration for lambda=1, alpha= 1.5 is "  (ln(2))^(1/1.5)
di "median duration for lambda=1, alpha= 2 is "  (ln(2))^(1/2)
di "median duration for lambda=1, alpha= 4 is "  (ln(2))^(1/4)

di "median duration for lambda=.5, alpha=.5 is "  (2*ln(2))^(1/.5)
di "median duration for lambda=.5, alpha= 1 is "  (2*ln(2))^(1/1)
di "median duration for lambda=.5, alpha= 1.5 is "  (2*ln(2))^(1/1.5)
di "median duration for lambda=.5, alpha= 2 is "  (2*ln(2))^(1/2)
di "median duration for lambda=.5, alpha= 4 is "  (2*ln(2))^(1/4)

di "median duration for lambda=2, alpha=.5 is "  (.5*ln(2))^(1/.5)
di "median duration for lambda=2, alpha= 1 is "  (.5*ln(2))^(1/1)
di "median duration for lambda=2, alpha= 1.5 is "  (.5*ln(2))^(1/1.5)
di "median duration for lambda=2, alpha= 2 is "  (.5*ln(2))^(1/2)
di "median duration for lambda=2, alpha= 4 is "  (.5*ln(2))^(1/4)

*  Mean: derive using formula

di "mean duration for lambda=1, alpha=.5 is "  exp(lngamma(3))
di "mean duration for lambda=1, alpha= 1 is "  exp(lngamma(2))
di "mean duration for lambda=1, alpha= 1.5 is " exp(lngamma(1+(1/1.5)))
di "mean duration for lambda=1, alpha= 2 is "  exp(lngamma(1.5))
di "mean duration for lambda=1, alpha= 4 is "  exp(lngamma(1.25))

di "mean duration for lambda=.5, alpha=.5 is " exp(lngamma(3))*(1/.5)^(1/.5)
di "mean duration for lambda=.5, alpha= 1 is "  exp(lngamma(2))*(1/.5)^(1/1)
di "mean duration for lambda=.5, alpha= 1.5 is "  exp(lngamma(1+(1/1.5)))*(1/.5)^(1/1.5)
di "mean duration for lambda=.5, alpha= 2 is "  exp(lngamma(1.5))*(1/.5)^(1/2)
di "mean duration for lambda=.5, alpha= 4 is "  exp(lngamma(1.25))*(1/.5)^(1/4)

di "mean duration for lambda=2, alpha=.5 is "  exp(lngamma(3))*(1/2)^(1/.5)
di "mean duration for lambda=2, alpha= 1 is "  exp(lngamma(2))*(1/2)^(1/1)
di "mean duration for lambda=2, alpha= 1.5 is "  exp(lngamma(1+(1/1.5)))*(1/2)^(1/1.5)
di "mean duration for lambda=2, alpha= 2 is "  exp(lngamma(1.5))*(1/2)^(1/2)
di "mean duration for lambda=2, alpha= 4 is "  exp(lngamma(1.25))*(1/2)^(1/4)


/*
Log-logistic model
------------------
hazard function h(t) = L^(1/g) * t^[(1/g)-1]
                       ----------------------
                       g*[ 1 + (L*t)^(1/g) ]

survivor function S(t) =         1
                          -----------------
   			   1 + (L*t)^(1/g) 

Parameters gamma (g) and lambda (L) = exp(-Xb).

*/

* Log-logistic hazard: lambda varying, gamma varying

ge hll1 = t/ (.5 * (1 + t^2) )
lab var hll1 "phi=1,gamma=.5"
ge hll2 = t^(-.5) / (2 * (1 + t^.5) )
lab var hll2 "phi=1,gamma=2"
ge hll3 = (4*t) / (.5 * (1 + (2*t)^2) ) 
lab var hll3 "phi=2,gamma=.5"
ge hll4 = ( 2^.5 * t^(-.5) ) / (2 * (1 + (2*t)^.5 ) )
lab var hll4 "phi=2,gamma=2"

/*
* Stata 7 graphics
gr7 hll1 hll2 hll3 hll4 t, xlabel(0,1,2,3,4,5) ylabel /*
   */ b1("Log-logistic hazard") l1("h(t)") c(llll) s(.pod) /*
   */ saving(logl1,replace)
*/

twoway  ( function y = x/ (.5 * (1 + x^2) ), range(0 5) yvarlab("phi=1,g=.5") )		///
    	( function y = x^(-.5) / (2 * (1 + x^.5) ), range(0 5) yvarlab("phi=1,g=2") )		///
	( function y = (4*x) / (.5 * (1 + (2*x)^2) ), range(0 5) yvarlab("phi=2,g=.5") )		///
     	( function y = ( 2^.5 * x^(-.5) ) / (2 * (1 + (2*x)^.5 ) ), range(0 5) yvarlab("phi=2,g=2") )	///
	, saving(logl1, replace)				///
	title("Log-logistic hazard function: phi=exp(-Xb) and gamma varying")	///
	ytitle(hazard) xtitle(t)



* Log-logistic survivor function: lambda varying, gamma varying

ge sll1 = 1 /  (1 + t^2)
lab var sll1 "lambda=1,gamma=.5"
ge sll2 = 1 / (1 + t^.5) 
lab var sll2 "lambda=1,gamma=2"
ge sll3 = 1 / (1 + (2*t)^2) 
lab var sll3 "lambda=2,gamma=.5"
ge sll4 = 1 / (1 + (2*t)^.5 ) 
lab var sll4 "lambda=2,gamma=2"

/*
* Stata 7 graphics

gr7 sll1 sll2 sll3 sll4 t, xlabel(0,1,2,3,4,5) ylabel  /*
   */ b1("Log-logistic survivor function") l1("S(t)") c(llll) s(.pod) /*
   */ saving(logl2,replace)
*/

twoway  ( function y = 1 /  (1 + x^2), range(0 5) yvarlab("phi=1,g=.5") )		///
    	( function y = 1 / (1 + x^.5), range(0 5) yvarlab("phi=1,g=2") )		///
	( function y = 1 / (1 + (2*x)^2) , range(0 5) yvarlab("phi=2,g=.5") )		///
     	( function y = 1 / (1 + (2*x)^.5 ) , range(0 5) yvarlab("phi=2,g=2") )	///
	, saving(logl2, replace)				///
	title("Log-logistic survivor function: phi=exp(-Xb) and gamma varying")	///
	ytitle(survival probability) xtitle(t)



/*
Log-logistic model: median duration
---------------------------------------

Median duration t' satisfies S(t') = 0.5

2 approximate derivations: 
(i) read off survivor function at S(t) = .5
(ii) use -list- and interpolate
Else derive exactly using -display- and inverting the survivor function:
     median  t' = (1/L)*[(1/.5) - 1]^g = (1/L)

NB for other percentiles in addition to median = p50, substitute relevant
   values for 0.5 in formula above.

Mean duration t* = (1/L)*(g*_pi)/sin(g*_pi) if g < 1.

*/


* Median: derive by interpolation
di "lambda=1,gamma=.5 (sll1 case)"
list t sll1 if abs(sll1-.5) < .1
di "lambda=1,gamma=2 (sll2 case)"
list t sll2 if abs(sll2-.5) < .1
di "lambda=2,gamma=.5 (sll3 case)"
list t sll3 if abs(sll3-.5) < .1
di "lambda=2,gamma=2 (sll4 case)"
list t sll4 if abs(sll4-.5) < .1

* Median: derive using formula
di "median duration for lambda=1, gamma=.5 is "  1
di "median duration for lambda=1, gamma= 2 is "  1
di "median duration for lambda=2, gamma= .5 is " .5
di "median duration for lambda=2, gamma= 2 is "  .5

* Mean: derive using formula
di "mean duration for lambda=1, gamma=.5 is "  (1/1)*(.5*_pi)/sin(.5*_pi)
di "median duration for lambda=2, gamma= .5 is " (1/2)*(.5*_pi)/sin(.5*_pi)



/*
(B) DISCRETE TIME MODELS
========================

Let z(t) = c(t) + Xb for a representative person in month t,  
			where c(t) is the baseline hazard function
			and Xb includes an intercept term.

Logistic discrete time hazard function:

	log[p(t)/(1-p(t))] = z(t)

      =>     p(t)   = 1/[1 + exp( -z(t) )] 

Complementary log-log ('cloglog') discrete time hazard function:

	log[-log( 1-p(t) )] = z(t)

      =>     p(t)   = 1 - exp[-exp( z(t) )] 

*/

* Look at the shapes of each hazard function against z

clear
set obs 101
ge z = ((_n-1)-50)/10
ge logitp  = 1/(1 + exp(-z) )
ge clogp  = 1-exp(-exp(z) )

twoway line logitp clogp z
drop z


/* 
We can re-write the hazard functions using the components of z:

Logistic hazard function:

	p(t)   	= 1/[1 + exp( -c(t) - Xb )] 
		
		= 1/[1 + (1/L)*exp(-c(t))]  where L = exp(Xb)

Complementary log-log ('cloglog') hazard function:

      =>     p(t)   = 1 - exp[-L*exp(c(t))] 


The discrete time survival function S(t)

   	S(t) = (1-p(1))*(1-p(2))*(1-p(3))*...*(1-p(t))

NB if p(t) does not depend on t (constant hazard), then
		S(t) = (1-p)^t
This yields a Geometric duration distribution
(discrete time analogue of Exponential distribution
 in continuous time model with constant hazard)

To calculate S(t), it is easier to re-write it as

		    s=t
	S(t) = exp{ SUM [ln(1-p(s))] }.
		    s=1

It is much harder to derive closed form solutions for the mean and 
median durations in discrete time models compared to continuous time
ones.

The median duration t' is defined implicitly by S(t') = 0.5

		      max(t)	
The mean duration t* = SUM  [ t*p(t) ].
		       t=0

Moreover we have not yet specified the shape of the baseline 
hazard function.

Consider the case c(t) = (q-1)*log(t)

[This is the discrete time analogue to the Weibull model--see why below.]

In this case, the logistic hazard simplifies to:

	p(t) = 1/[1 + (1/L)*t^(1-q)]  

and the cloglog hazard to:

	p(t)   = 1 - exp[-L*t^(q-1)]

Now let us look at the hazard and survival functions for specific values
of q and L:

First create an artificial data set to represent
'time'. Index the intervals ('months') by integers 1,...,10.

*/

clear
set obs 10
ge t = _n
compress


* Logistic hazard model: q = 0.5, 1, 2; L = 0.5, 1, 2
* ===================================================
*	p(t) = 1/[1 + (1/L)*t^(1-q)]   

ge hlog1 =  1/(1 + t^(-.5))
lab var hlog1 "lambda=1,q=1.5"
ge hlog2 =  1/2
lab var hlog2 "lambda=1,q=1"
ge hlog3 =  1/(1 + t^.5)
lab var hlog3 "lambda=1,q=0.5"

/* Stata 7 graphics
gr7 hlog1 hlog2 hlog3 t,xlabel ylabel xtick(0,1,2,3,4,5,6,7,8,9,10) /*
	*/ b1("Logistic hazard,c(t)=(q-1)*log(t)") /*
	*/ l1("h(t)") c(lll) s(.po) saving(logist1,replace)

*/

twoway connect hlog1 hlog2 hlog3 t ///
	, title("Logistic hazard,c(t)=(q-1)*log(t)")  ///
	ytitle("h(t)") xtick(0(1)10) saving(logist1,replace)


ge hlog4 =  1/(1 + 2*(t^(-.5)))
lab var hlog4 "lambda=.5,q=1.5"
ge hlog5 =  1/(1 + .5*(t^(-.5)))
lab var hlog5 "lambda=2,q=1.5"

/* Stata 7 graphics
gr7 hlog1 hlog4 hlog5 t,xlabel ylabel xtick(0,1,2,3,4,5,6,7,8,9,10) /*
	*/ b1("Logistic hazard,c(t)=(q-1)*log(t)") /*
	*/ l1("h(t)") c(lll) s(.po) saving(logist2,replace)
*/

twoway connect hlog1 hlog4 hlog5 t ///
	, title("Logistic hazard,c(t)=(q-1)*log(t)")  ///
	ytitle("h(t)") xtick(0(1)10) saving(logist2,replace)


ge hlog6 =  1/(1 + 2*(t^(.5)))
lab var hlog6 "lambda=.5,q=0.5"
ge hlog7 =  1/(1 + .5*(t^(.5)))
lab var hlog7 "lambda=2,q=0.5"

/* Stata 7 graphics
gr7 hlog3 hlog6 hlog7 t,xlabel ylabel xtick(0,1,2,3,4,5,6,7,8,9,10) /*
	*/ b1("Logistic hazard,c(t)=(q-1)*log(t)") /*
	*/ l1("h(t)") c(lll) s(.po) saving(logist3,replace)
*/

twoway connect hlog3 hlog6 hlog7 t ///
	, title("Logistic hazard,c(t)=(q-1)*log(t)")  ///
	ytitle("h(t)") xtick(0(1)10) saving(logist3,replace)



sort t
ge slog1 =  exp(sum(ln(1-hlog1)))
lab var slog1 "lambda=1,q=1.5"
ge slog2 =  exp(sum(ln(1-hlog2)))
lab var slog2 "lambda=1,q=1"
ge slog3 =  exp(sum(ln(1-hlog3)))
lab var slog3 "lambda=1,q=0.5"
ge slog4 =  exp(sum(ln(1-hlog4)))
lab var slog4 "lambda=.5,q=1.5"
ge slog5 =  exp(sum(ln(1-hlog5)))
lab var slog5 "lambda=2,q=1.5"

/* Stata 7 graphics
gr7 slog1 slog2 slog3 t,xlabel ylabel xtick(0,1,2,3,4,5,6,7,8,9,10) /*
	*/ b1("Logistic survivor function,c(t)=(q-1)*log(t)") /*
	*/ l1("S(t)") c(lll) s(.po) saving(logist4,replace)

gr7 slog1 slog4 slog5 t, xlabel ylabel xtick(0,1,2,3,4,5,6,7,8,9,10) /*
	*/ b1("Logistic survivor function,c(t)=(q-1)*log(t)") /*
	*/ l1("S(t)") c(lll) s(.po) saving(logist5,replace)

gr7 using logist1 logist2 logist4 logist5, saving(logist6,replace)
*/

twoway connect slog1 slog2 slog3 t ///
	, title("Logistic survivor function,c(t)=(q-1)*log(t)")  ///
	ytitle("S(t)") xtick(0(1)10) saving(logist4,replace)

twoway connect slog1 slog4 slog5 t ///
	, title("Logistic survivor function,c(t)=(q-1)*log(t)")  ///
	ytitle("S(t)") xtick(0(1)10) saving(logist5,replace)

graph combine logist1.gph logist2.gph logist4.gph logist5.gph, saving(logist6,replace)


* Cloglog hazard model: q = 0.5, 1, 2; L = 0.5, 1, 2
* ===================================================
*	p(t)   = 1 - exp[-L*t^(q-1)]

ge hclog1 =  1-exp(-1*t^.5)
lab var hclog1 "lambda=1,q=1.5"
ge hclog2 =  1-exp(-1)
lab var hclog2 "lambda=1,q=1"
ge hclog3 = 1-exp(-1*t^(-.5))
lab var hclog3 "lambda=1,q=0.5"

/* Stata 7 graphics
gr7 hclog1 hclog2 hclog3 t,xlabel ylabel xtick(0,1,2,3,4,5,6,7,8,9,10) /*
	*/ b1("Cloglog hazard,c(t)=(q-1)*log(t)") /*
	*/ l1("h(t)") c(lll) s(.po) saving(clog1,replace)
*/

twoway connect hclog1 hclog2 hclog3 t ///
	, title("cloglog hazard,c(t)=(q-1)*log(t)")  ///
	ytitle("h(t)") xtick(0(1)10) saving(clog1,replace)


ge hclog4 =  1-exp(-.5*t^.5)
lab var hclog4 "lambda=.5,q=1.5"
ge hclog5 =  1-exp(-2*t^.5)
lab var hclog5 "lambda=2,q=1.5"

/* Stata 7 graphics
gr7 hclog1 hclog4 hclog5 t,xlabel ylabel xtick(0,1,2,3,4,5,6,7,8,9,10) /*
	*/ b1("cloglog hazard,c(t)=(q-1)*log(t)") /*
	*/ l1("h(t)") c(lll) s(.po) saving(clog2,replace)
*/

twoway connect hclog1 hclog4 hclog5 t ///
	, title("cloglog hazard,c(t)=(q-1)*log(t)")  ///
	ytitle("h(t)") xtick(0(1)10) saving(clog2,replace)



ge hclog6 =  1-exp(-.5*t^(-.5))
lab var hclog6 "lambda=.5,q=0.5"
ge hclog7 =  1-exp(-2*t^(-.5))
lab var hclog7 "lambda=2,q=0.5"


/* Stata 7 graphics
gr7 hclog3 hclog6 hclog7 t,xlabel ylabel xtick(0,1,2,3,4,5,6,7,8,9,10) /*
	*/ b1("Cloglog hazard,c(t)=(q-1)*log(t)") /*
	*/ l1("h(t)") c(lll) s(.po) saving(clog3,replace)
*/

twoway connect hclog3 hclog6 hclog7 t ///
	, title("cloglog hazard,c(t)=(q-1)*log(t)")  ///
	ytitle("h(t)") xtick(0(1)10) saving(clog3,replace)


sort t
ge sclog1 =  exp(sum(ln(1-hclog1)))
lab var sclog1 "lambda=1,q=1.5"
ge sclog2 =  exp(sum(ln(1-hclog2)))
lab var sclog2 "lambda=1,q=1"
ge sclog3 =  exp(sum(ln(1-hclog3)))
lab var sclog3 "lambda=1,q=0.5"
ge sclog4 =  exp(sum(ln(1-hclog4)))
lab var sclog4 "lambda=.5,q=1.5"
ge sclog5 =  exp(sum(ln(1-hclog5)))
lab var sclog5 "lambda=2,q=1.5"

/* Stata 7 graphics
gr7 sclog1 sclog2 sclog3 t,xlabel ylabel xtick(0,1,2,3,4,5,6,7,8,9,10) /*
	*/ b1("Cloglog survivor function,c(t)=(q-1)*log(t)") /*
	*/ l1("S(t)") c(lll) s(.po) saving(clog4,replace)

gr7 sclog1 sclog4 sclog5 t, xlabel ylabel xtick(0,1,2,3,4,5,6,7,8,9,10) /*
	*/ b1("Logistic survivor function,c(t)=(q-1)*log(t)") /*
	*/ l1("S(t)") c(lll) s(.po) saving(clog5,replace)

gr7 using clog1 clog2 clog4 clog5, saving(clog6,replace)
*/


twoway connect sclog1 sclog2 sclog3 t ///
	, title("cloglog survivor function,c(t)=(q-1)*log(t)")  ///
	ytitle("S(t)") xtick(0(1)10) saving(clog4,replace)

twoway connect sclog1 sclog4 sclog5 t ///
	, title("cloglog survivor function,c(t)=(q-1)*log(t)")  ///
	ytitle("S(t)") xtick(0(1)10) saving(clog5,replace)

graph combine clog1.gph clog2.gph clog4.gph clog5.gph, saving(clog6,replace)




/*
* Compare logistic and cloglog hazard and survivor functions directly
gr7 using logist1 logist4 clog1 clog4, saving(clog7,replace)
*/

graph combine logist1.gph logist4.gph clog1.gph clog4.gph, saving(clog7,replace)

/* 
Given closed form for survivor function not generally available 
in discrete time case, typically one interpolates from the estimated 
survivor function
*/


* Median (logistic hazard): derive by interpolation
di "lambda=1,q=1.5 (slog1 case)"
list t slog1 if abs(slog1-.5) < .1
di "lambda=1,q=1 (slog2 case)"
list t slog2 if abs(slog2-.5) < .1
di "lambda=1,q=0.5 (slog3 case)"
list t slog3 if abs(slog3-.5) < .1
di "lambda=0.5,q=1.5 (slog4 case)"
list t slog4 if abs(slog4-.5) < .1


* Median (cloglog hazard): derive by interpolation
di "lambda=1,q=1.5 (sclog1 case)"
list t sclog1 if abs(sclog1-.5) < .1
di "lambda=1,q=1 (sclog2 case)"
list t sclog2 if abs(sclog2-.5) < .1
di "lambda=1,q=0.5 (sclog3 case)"
list t sclog3 if abs(sclog3-.5) < .1
di "lambda=0.5,q=1.5 (sclog4 case)"
list t sclog4 if abs(sclog4-.5) < .1

/* 
Closed forms for the survivor functions are
available in the case when the hazard does not
depend on time: p(t) = p
e.g. q=1 in c(t) = (q-1)*log(t) case.

	S(t) = (1-p)^t
=>	t = log[S(t)]/log(1-p)

=>	median t' = log(.5)/log(1-p)

In the logistic hazard case

	p = 1/[1+(1/L)] = L/(1+L) => 1-p = 1/(1+L)
=> 	log(1-p) = -log(1+L)

Hence	t' = -log(.5)/log(1+L) = log(2)/log(1+L)

In the cloglog hazard case

	p = 1-exp(-L) => 1-p = exp(-L)
=> 	log(1-p) = -L

Hence	t' = -log(.5)/L = log(2)/L

*/

di "Median (logistic hazard, lambda=.5,q=1) = " ln(2)/ln(1.5)
di "Median (logistic hazard, lambda=1,q=1) = " ln(2)/ln(2)
di "Median (logistic hazard, lambda=2,q=1) = " ln(2)/ln(3)

di "Median (cloglog hazard, lambda=.5,q=1) = " ln(2)/.5
di "Median (cloglog hazard, lambda=1,q=1) = " ln(2)
di "Median (cloglog hazard, lambda=2,q=1) = " ln(2)/2

// Now think about how you might draw density and integrated hazard functions!

/* Weibull density and integrated hazard

f(t) = a*L*t^(a-1)*exp(-L*t^a)

H(t) = L*t^a

*/

* Weibull density function: lambda=1, alpha varying
twoway  (function y = .5*x^(-.5)*exp(-x^.5), range(0 5) yvarlab("a=.5") )   ///
     	( function y = 1.5*x^(.5)*exp(-x^1.5), range(0 5) yvarlab("a=1.5") ) ///
     	( function y = 1*x^(0)*exp(-x^1), range(0 5) yvarlab("a=1") )	///
     	( function y = 2*x*exp(-x^2), range(0 5) yvarlab("a=2") )			///
	, saving(weib6, replace)				///
	title("Weibull density function: lambda=1, alpha varying")	///
	ytitle(density) xtitle(t)	

* Weibull density function: lambda varying, alpha = .5
twoway  (function y = .5*x^(-.5)*exp(-.5*x^.5), range(0 5) yvarlab("L=.5") )   ///
     	( function y = .5*x^(-.5)*exp(-x^.5), range(0 5) yvarlab("L=1") ) ///
     	( function y = .5*x^(-.5)*exp(-2*x^.5), range(0 5) yvarlab("L=2") )	///
	, saving(weib7, replace)				///
	title("Weibull density function: lambda varying, alpha = .5")	///
	ytitle(density) xtitle(t)	

* Weibull density function: lambda varying, alpha = 2
twoway  (function y = 2*x*exp(-.5*x^2), range(0 5) yvarlab("L=.5") )   ///
        ( function y = 2*x*exp(-x^2), range(0 5) yvarlab("L=1") )	///
     	( function y = 2*x*exp(-2*x^2), range(0 5) yvarlab("L=2") )	///
	, saving(weib8, replace)				///
	title("Weibull density function: lambda varying, alpha = 2")	///
	ytitle(density) xtitle(t)

* Weibull integrated hazard function: lambda=1, alpha varying 
twoway  ( function y = x^.5, range(0 5) yvarlab("a=.5") )		///
    	( function y = x, range(0 5) yvarlab("a=1") )			///
	( function y = x^1.5, range(0 5) yvarlab("a=1.5") )		///
     	( function y = x^2, range(0 5) yvarlab("a=2") )			///
	, saving(weib5, replace)				///
	title("Weibull integrated hazard function: lambda = 1, alpha varying")	///
	ytitle("Integrated hazard") xtitle(t)	




/* Log-logistic density and integrated hazard

hazard function h(t) = L^(1/g) * t^[(1/g)-1]
                       ----------------------
                       g*[ 1 + (L*t)^(1/g) ]

survivor function S(t) =         1
                          -----------------
   			   1 + (L*t)^(1/g) 

Parameters gamma (g) and lambda (L) = exp(-Xb).

f(t) = h(t)*S(t)

H(t) = -logS(t) = log( 1 + (L*t)^(1/g) )

*/

twoway  ( function y = (x/ (.5 * (1 + x^2) ))/(1 + x^2), range(0 5) yvarlab("phi=1,g=.5") )		///
    	( function y = (x^(-.5) / (2 * (1 + x^.5) ))/(1 + x^.5), range(0 5) yvarlab("phi=1,g=2") )		///
	( function y = ((4*x) / (.5 * (1 + (2*x)^2) ))/(1 + (2*x)^2), range(0 5) yvarlab("phi=2,g=.5") )		///
     	( function y = (( 2^.5 * x^(-.5) ) / (2 * (1 + (2*x)^.5 ) ))/(1 + (2*x)^.5 ), range(0 5) yvarlab("phi=2,g=2") )	///
	, saving(logl3, replace)				///
	title("Log-logistic density function: phi=exp(-Xb); gamma varying")	///
	ytitle(density) xtitle(t)



twoway  ( function y = log(1 + x^2), range(0 5) yvarlab("phi=1,g=.5") )		///
    	( function y = log(1 + x^.5), range(0 5) yvarlab("phi=1,g=2") )		///
	( function y = log(1 + (2*x)^2) , range(0 5) yvarlab("phi=2,g=.5") )		///
     	( function y = log(1 + (2*x)^.5 ) , range(0 5) yvarlab("phi=2,g=2") )	///
	, saving(logl4, replace)				///
	title("Log-logistic integrated hazard: phi=exp(-Xb); gamma varying")	///
	ytitle(integrated hazard) xtitle(t)




* Now you can repeat for other cases such as lognormal case and Gompertz, etc.

