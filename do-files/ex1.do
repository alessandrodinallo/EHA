

cd "  " // [Insert your path in the quotation marks]


use duration
* 1. What are the variable names, and how many observations?
describe

* 2. What are the value labels for the variables? [Hint: -help label-]
label list

* 3. How many men are there in the data? (# and % of all persons)
ta sex

* 4. How many women are married? (# and % of all women)
ta sex married, row col
ta married if sex == 0

* 5. What does the distribution of marriage durations look like?
*   -- e.g. min, max, median, mean; general shape etc
list time
ta time
ge dur = time
recode dur min/5=1 6/10=2 11/15=3 16/20=4 21/25=5 26/30=6 31/max=7
label var dur "=time,recoded"
label def dur 1 "min-5" 2 "6-10" 3 "11-15" 4 "16-20" 5 "21-25" /*
	*/  6 "26-30" 7 "31-max"
lab val dur dur
ta dur
inspect time
su time
su time, de
 
histogram time, bin(10)  saving(junk1.gph,replace)

* 6. How do mean marriage durations differ between
*     persons with censored and uncensored spells?
su time if status==0
su time if status==1
sort status
by status: su time

* Let's create a new (fictional) variable 'age'.
* For each person, we'll set it equal to their observation number
*  (for observations sorted in ascending order by 'time') plus 30
* Also, for practice sake, we'll assume any ages greater than 50  
*  are 'missing'.
sort time
ge age = _n + 30
replace age = . if age > 50

* 7. How many cases are missing a value for age?
count if age == .

* 8. What is the mean, median, min and max of age? 
su age, de

* Showing off a little, using -display- and the
* stored results from -summarise-  (See Ref Man 4, summarize):
di "Mean(age) = " r(mean) " Median(age) = " r(p50)
di "Min(age) = " r(min) " Max(age) = " r(max)

* 9. Create a new age variable with values equal
*    to the square each person's age?
ge age2 = age^2 

* 10. Create a new age variable with values equal
*    to the natural logarithm of each person's age?
ge lnage = ln(age)

* 11. Create a new categorical variable summarising age 
*     in 5-year bands (30-34, 35-39, 40-44, 45-49, 50-).
*     Give the new variable a useful variable name
*     and label its values.  Hint: -help label-* 
*     NB there are several ways of creating the variable.  
*     Hint: methods include using generate and replace commands; 
*     generate command with its recode function; recode command.
*     Beware of how you treat the missing values -- 
*     ensure that the categorical age variable is missing 
*     whenever age is missing.

ge agecat1 = .
replace agecat1 = 1 if age < 35
replace agecat1 = 2 if age >= 35 & age < 40
replace agecat1 = 3 if age >= 40 & age < 45
replace agecat1 = 4 if age >= 45 & age < 50
replace agecat1 = 5 if age >= 50 & age < .
   /* Note the "& age < ." qualifier in the last line. Without
   this those persons with missing values on age would have agecat1=3
   rather than missing, as required.  */
ta agecat1
ta agecat1, missing

ge agecat2 = age if age ~=.  	// NB use of 'if'
recode agecat2 min/34=1 35/39=2 40/44=3 45/49=4 50/max=5 
ta agecat2

ge agecat3 = recode(age,34,39,44,49)
    /* see -help functions- under Special Functions 
       This has the advantage of being succinct and
       automatically treating missing values OK --
       but I never remember the syntax! Especially when 
       I want to add value labels */

lab var agecat1 "Age,5-year bands"
lab var agecat2 "Age,5-year bands"
lab var agecat3 "Age,5-year bands"
lab def agecat 1 "-34" 2 "35-39" 3 "40-44" 4 "45-49" 5 "50-"
lab val agecat1 agecat
lab val agecat2 agecat

ta agecat1
ta agecat2
ta agecat3

ta agecat1, nolabel
ta agecat3, nolabel

* 12. Create a set of dummy variables based on the categorical variable
*     derived in Q11.  
*     NB there are several ways of doing this.  You need to decide how to
*     treat the missing values.  Should they equal zero for each of the
*     dummies?  Or equal to missing for all the dummies?
*     Use the second convention in this exercise.
*     Hint: use generate command and logical expressions, or the tabulate
*     command

tab agecat1 if age < ., gen(agect1)
tab1 agect1*

ge agect21 = (agecat1 == 1)  if age < .
ge agect22 = (agecat1 == 2)  if age < .
ge agect23 = (agecat1 == 3)  if age < .
ge agect24 = (agecat1 == 4)  if age < .
ge agect25 = (agecat1 == 5)  if age < .
tab1 agect2*

* 13.  -describe- your data set, then -compress- it, and finally
*      -describe- it again.  Look at what happened to the 
*      size of the file in Kb, and the storage types of variables.
de
compress
de

* 14.  Now drop the log(age) variable from Q11 from the data set, and 
*      save the data in a file--call it "mydurat.dta" 
*      (We use a new name so as not to over-write the original data set.)
drop lnage
save mydurat, replace  
    /* ",replace" required if want to over-write earlier versions with the
        same filename  */


