/*==============================================================================
MORT3-13a: MULTIPLE CAUSE-OF-DEATH ANALYSIS (underlying vs any mention)

PURPOSE
	Among the 5,922 deaths in the BCCSS five-year survivor cohort, use the
	underlying cause (ucause) AND the contributory causes (mcause_1..mcause_15)
	on the death certificate to describe, for each ICD chapter group:
	  (i)   deaths with the group as underlying cause
	  (ii)  deaths with the group mentioned anywhere on the certificate
	  (iii) ratio (ii)/(i)
	  (iv)  deaths where the group is contributory only
	overall and by attained age at death (<40 vs 40+)          [Table A]
	and, among deaths from recurrence, SPN and circulatory disease, the % with
	infection / respiratory / circulatory / genitourinary (renal) mentioned
	as a contributory cause                                     [Table B]
	plus a horizontal bar chart underlying vs any mention.

	A classifier (program -codgrp-) maps any ICD code + ICD revision (6-10)
	to the chapter groups used in the paper; it is validated against codcat
	(gold standard) for the underlying cause.

	NOTE: contributory causes are only present on certificates from 1993
	onwards (ONS multiple-cause coding; ICD-9 1993-2000 and ICD-10 2001-2022).
	ICD-6/7/8 certificates and ICD-9 certificates before 1993 have NO
	contributory causes. The any-mention columns are therefore restricted to
	deaths in 1993-2022 (restricted denominators are reported in the table).

INPUTS
	$temp/x-mort3-prepforstset.dta
OUTPUTS
	$output/table_multiple_cause.docx   (Tables A and B)
	$output/table_multiple_cause.csv    (numbers behind Table A, long)
	$output/table_multiple_cause_B.csv  (numbers behind Table B)
	$output/fig_multiple_cause.png

Date: 23 Oct 2026    Author: RCR (Claude-assisted)
==============================================================================*/

do "/Users/raoulreulen/Documents/GitHub/mort3/mort3-00c-rr-setwd-mac.do"

*-------------------------------------------------------------------------------
* 0. CLASSIFIER: ICD code (str) + icdver -> chapter group
*    groups: 1 infection 2 neoplasm 3 blood 4 endocrine 5 mental 6 nervous
*            7 circulatory 8 respiratory 9 digestive 10 musculoskeletal/skin
*            11 genitourinary 12 congenital/perinatal 13 other (symptoms,
*            pregnancy, ill-defined) 14 external 15 COVID-19 (ICD-10 U07)
*            0 = blank / unclassifiable
*    circulatory sub-group: 1 cardiac 2 cerebrovascular 3 other circulatory
*-------------------------------------------------------------------------------
cap program drop codgrp
program define codgrp
	// syntax: codgrp newvar codevar icdvervar
	args new code ver
	tempvar c l n3 sub
	gen str4 `c' = upper(trim(`code'))
	gen str1 `l' = substr(`c',1,1)
	// numeric 3-digit part (ICD-6/7/8/9 codes stored as 3 or 4 chars without
	// decimal point, e.g. "1919"=191.9, "586"=586; E-codes "E912"=E912;
	// ICD-10 stored as letter + 2/3 digits, e.g. "C719"=C71.9, "B24"=B24)
	gen `n3' = real(substr(`c',1,3)) if !inrange(`l',"A","Z")
	replace `n3' = real(substr(`c',2,3)) if inlist(`l',"E","N") & `ver'<=9
	gen `sub' = real(substr(`c',2,2)) if inrange(`l',"A","Z") & `ver'==10
	gen `new' = 0
	gen `new'_circ = .

	*--- ICD-8 and ICD-9 (same chapter boundaries) ---
	replace `new' = 1  if inrange(`ver',8,9) & !inlist(`l',"E","V") & inrange(`n3',0,139)
	replace `new' = 2  if inrange(`ver',8,9) & !inlist(`l',"E","V") & inrange(`n3',140,239)
	replace `new' = 4  if inrange(`ver',8,9) & !inlist(`l',"E","V") & inrange(`n3',240,279)
	replace `new' = 3  if inrange(`ver',8,9) & !inlist(`l',"E","V") & inrange(`n3',280,289)
	replace `new' = 5  if inrange(`ver',8,9) & !inlist(`l',"E","V") & inrange(`n3',290,319)
	replace `new' = 6  if inrange(`ver',8,9) & !inlist(`l',"E","V") & inrange(`n3',320,389)
	replace `new' = 7  if inrange(`ver',8,9) & !inlist(`l',"E","V") & inrange(`n3',390,459)
	replace `new' = 8  if inrange(`ver',8,9) & !inlist(`l',"E","V") & inrange(`n3',460,519)
	replace `new' = 9  if inrange(`ver',8,9) & !inlist(`l',"E","V") & inrange(`n3',520,579)
	replace `new' = 11 if inrange(`ver',8,9) & !inlist(`l',"E","V") & inrange(`n3',580,629)
	replace `new' = 13 if inrange(`ver',8,9) & !inlist(`l',"E","V") & inrange(`n3',630,679)  // pregnancy
	replace `new' = 10 if inrange(`ver',8,9) & !inlist(`l',"E","V") & inrange(`n3',680,739)
	replace `new' = 12 if inrange(`ver',8,9) & !inlist(`l',"E","V") & inrange(`n3',740,779)  // congenital + perinatal
	replace `new' = 13 if inrange(`ver',8,9) & !inlist(`l',"E","V") & inrange(`n3',780,799)  // symptoms/ill-defined
	replace `new' = 14 if inrange(`ver',8,9) & !inlist(`l',"E","V") & inrange(`n3',800,999)  // nature of injury (N codes)
	replace `new' = 14 if inrange(`ver',8,9) & `l'=="E" & inrange(`n3',800,999)              // external cause E codes
	replace `new' = 13 if inrange(`ver',8,9) & `l'=="V"                                      // supplementary V codes
	replace `new'_circ = 1 if inrange(`ver',8,9) & `new'==7 & inrange(`n3',390,429)
	replace `new'_circ = 2 if inrange(`ver',8,9) & `new'==7 & inrange(`n3',430,438)
	replace `new'_circ = 3 if inrange(`ver',8,9) & `new'==7 & mi(`new'_circ)

	*--- ICD-6 and ICD-7 (1949-1967; chapter boundaries differ) ---
	replace `new' = 1  if inrange(`ver',6,7) & `l'!="E" & inrange(`n3',1,138)
	replace `new' = 2  if inrange(`ver',6,7) & `l'!="E" & inrange(`n3',140,239)
	replace `new' = 4  if inrange(`ver',6,7) & `l'!="E" & inrange(`n3',240,289)
	replace `new' = 3  if inrange(`ver',6,7) & `l'!="E" & inrange(`n3',290,299)
	replace `new' = 5  if inrange(`ver',6,7) & `l'!="E" & inrange(`n3',300,329)
	replace `new' = 6  if inrange(`ver',6,7) & `l'!="E" & inrange(`n3',330,399)
	replace `new' = 7  if inrange(`ver',6,7) & `l'!="E" & inrange(`n3',400,469)
	replace `new' = 8  if inrange(`ver',6,7) & `l'!="E" & inrange(`n3',470,529)
	replace `new' = 9  if inrange(`ver',6,7) & `l'!="E" & inrange(`n3',530,589)
	replace `new' = 11 if inrange(`ver',6,7) & `l'!="E" & inrange(`n3',590,639)
	replace `new' = 13 if inrange(`ver',6,7) & `l'!="E" & inrange(`n3',640,689)  // pregnancy
	replace `new' = 10 if inrange(`ver',6,7) & `l'!="E" & inrange(`n3',690,749)
	replace `new' = 12 if inrange(`ver',6,7) & `l'!="E" & inrange(`n3',750,779)
	replace `new' = 13 if inrange(`ver',6,7) & `l'!="E" & inrange(`n3',780,799)
	replace `new' = 14 if inrange(`ver',6,7) & `l'!="E" & inrange(`n3',800,999)
	replace `new' = 14 if inrange(`ver',6,7) & `l'=="E"
	replace `new'_circ = 1 if inrange(`ver',6,7) & `new'==7 & inrange(`n3',400,447)   // heart + hypertension; vascular lesions of CNS (330-334) sit in the nervous chapter in ICD-6/7
	replace `new'_circ = 3 if inrange(`ver',6,7) & `new'==7 & mi(`new'_circ)

	*--- ICD-10 ---
	replace `new' = 1  if `ver'==10 & inlist(`l',"A","B")
	replace `new' = 2  if `ver'==10 & (`l'=="C" | (`l'=="D" & `sub'<=48))
	replace `new' = 3  if `ver'==10 & `l'=="D" & inrange(`sub',50,89)
	replace `new' = 4  if `ver'==10 & `l'=="E"
	replace `new' = 5  if `ver'==10 & `l'=="F"
	replace `new' = 6  if `ver'==10 & inlist(`l',"G","H")
	replace `new' = 7  if `ver'==10 & `l'=="I"
	replace `new' = 8  if `ver'==10 & `l'=="J"
	replace `new' = 9  if `ver'==10 & `l'=="K"
	replace `new' = 10 if `ver'==10 & inlist(`l',"L","M")
	replace `new' = 11 if `ver'==10 & `l'=="N"
	replace `new' = 13 if `ver'==10 & `l'=="O"            // pregnancy
	replace `new' = 12 if `ver'==10 & inlist(`l',"P","Q") // perinatal, congenital
	replace `new' = 13 if `ver'==10 & `l'=="R"            // symptoms/ill-defined
	replace `new' = 14 if `ver'==10 & inlist(`l',"S","T","V","W","X","Y")
	replace `new' = 15 if `ver'==10 & `l'=="U" & `sub'==7  // U07.1/U07.2 COVID-19
	replace `new' = 14 if `ver'==10 & `l'=="U" & `sub'==50 // U50.9 (ONS: accident awaiting inquest) -> external
	replace `new' = 13 if `ver'==10 & `l'=="U" & `new'==0  // other U codes
	replace `new'_circ = 1 if `ver'==10 & `new'==7 & inrange(`sub',0,52)
	replace `new'_circ = 2 if `ver'==10 & `new'==7 & inrange(`sub',60,69)
	replace `new'_circ = 3 if `ver'==10 & `new'==7 & mi(`new'_circ)

	replace `new' = 0 if `c'==""
	label define lgrp 0 "blank/unclassified" 1 "infection" 2 "neoplasm" 3 "blood" 		///
		4 "endocrine" 5 "mental" 6 "nervous" 7 "circulatory" 8 "respiratory" 			///
		9 "digestive" 10 "musculoskeletal/skin" 11 "genitourinary" 					///
		12 "congenital/perinatal" 13 "other" 14 "external" 15 "COVID-19", modify
	label values `new' lgrp
	label define lcirc 1 "cardiac" 2 "cerebrovascular" 3 "other circulatory", modify
	label values `new'_circ lcirc
end

*-------------------------------------------------------------------------------
* 1. LOAD DEATHS, INSPECT CODE FORMATS
*-------------------------------------------------------------------------------
use "$temp/x-mort3-prepforstset", clear
keep if allcauses==1
count
assert r(N)==5922

tab icdver, m
gen yod = year(dod)
gen ageatdeath = (dod - dob)/365.25
gen age40 = ageatdeath>=40
label define l40 0 "<40" 1 "40+"
label values age40 l40

* number of contributory causes
gen ncontrib = 0
forval i=1/15 {
	replace ncontrib = ncontrib + (trim(mcause_`i')!="")
}
gen anycontrib = ncontrib>0
tab icdver anycontrib, m row
tab yod anycontrib if yod>=1985, m
bysort icdver: sum ncontrib
sum yod if anycontrib==1
local firstyr = r(min)
di as res "First year of death with any contributory cause: `firstyr'"

* sample of code formats by ICD revision
forval v=6/10 {
	di as res "---- ICD-`v': sample of underlying codes by codcat"
	preserve
	keep if icdver==`v'
	bysort codcat (ucause): gen k=_n
	list codcat ucause mcause_1 mcause_2 if k<=3, noobs clean
	restore
}
tab ucause if codcat==0

*-------------------------------------------------------------------------------
* 2. APPLY CLASSIFIER TO UNDERLYING CAUSE AND VALIDATE AGAINST codcat
*-------------------------------------------------------------------------------
codgrp ugrp ucause icdver
tab codcat ugrp, m
tab circsubcat ugrp_circ if codcat==7, m

* agreement (codcat groups 12 perinatal & 13 other are both 'other' in paper;
* codcat 0 unknown = U07 COVID-19 in ICD-10 + blanks)
gen agree = (ugrp==codcat)
replace agree = 1 if codcat==0 & ugrp==15      // COVID-19: codcat has no group
tab agree, m
di as res "Disagreements (excluding COVID-19 / codcat unknown):"
list icdver ucause codcat ugrp mcause_1 if agree==0, noobs clean

* RECONCILE: use codcat (gold standard) for the underlying cause group, with
* COVID-19 (codcat unknown, U07) kept as its own group 15; the classifier is
* used for the contributory causes.
gen ugrp_final = codcat
replace ugrp_final = 15 if codcat==0 & ugrp==15
label values ugrp_final lgrp
tab ugrp_final, m

*-------------------------------------------------------------------------------
* 3. APPLY CLASSIFIER TO CONTRIBUTORY CAUSES; any-mention flags
*-------------------------------------------------------------------------------
forval i=1/15 {
	codgrp mgrp`i' mcause_`i' icdver
}
* check: any non-blank contributory code that could not be classified
gen nunclass = 0
forval i=1/15 {
	replace nunclass = nunclass + (trim(mcause_`i')!="" & mgrp`i'==0)
}
tab nunclass
forval i=1/15 {
	list icdver mcause_`i' if trim(mcause_`i')!="" & mgrp`i'==0, noobs clean
}

* per-group flags: contrib_g = mentioned in any contributory position
* (ONS mentions include the underlying cause itself in the mcause list, so
*  'contributory' is defined as mentioned in a group other than the underlying
*  group, i.e. contrib_g=1 only if g != underlying group)
forval g=1/15 {
	gen ment`g' = 0
	forval i=1/15 {
		replace ment`g' = 1 if mgrp`i'==`g'
	}
	gen under`g'   = ugrp_final==`g'
	gen contrib`g' = ment`g'==1 & under`g'==0
	gen any`g'     = under`g'==1 | contrib`g'==1
}
* circulatory sub-groups (cardiac / cerebrovascular) ; underlying uses circsubcat
forval s=1/2 {
	gen umentc`s' = 0
	forval i=1/15 {
		replace umentc`s' = 1 if mgrp`i'==7 & mgrp`i'_circ==`s'
	}
	gen underc`s'   = codcat==7 & circsubcat==`s'
	gen contribc`s' = umentc`s'==1 & underc`s'==0
	gen anyc`s'     = underc`s'==1 | contribc`s'==1
}

* restriction: deaths with multiple-cause coding (1993 onwards)
gen mcc = yod>=`firstyr'
tab mcc icdver, m
tab mcc anycontrib, row

save "$temp/x-mort3-multicause", replace

*-------------------------------------------------------------------------------
* 4. TABLE A: underlying / any mention / ratio / contributory-only
*    by cause group, overall and by attained age (<40, 40+), restricted to
*    deaths from `firstyr' onwards
*-------------------------------------------------------------------------------
tempfile tabA
local first = 1
local groups 1 2 3 4 5 6 7 71 72 8 9 10 11 12 13 14 15
foreach g of local groups {
	local gname : label lgrp `g'
	if `g'==71 local gname "  cardiac"
	if `g'==72 local gname "  cerebrovascular"
	local v `g'
	if `g'==71 local v c1
	if `g'==72 local v c2
	foreach s in all 0 1 {
		preserve
		if "`s'"=="all" keep if mcc==1
		else keep if mcc==1 & age40==`s'
		qui count
		local N = r(N)
		qui count if under`v'==1
		local nu = r(N)
		qui count if any`v'==1
		local na = r(N)
		qui count if contrib`v'==1
		local nc = r(N)
		clear
		set obs 1
		gen grp = `g'
		gen str40 group = "`gname'"
		gen str5 agestr = "`s'"
		gen N = `N'
		gen underlying = `nu'
		gen anymention = `na'
		gen ratio = `na'/`nu'
		gen contribonly = `nc'
		if `first'==0 append using `tabA'
		save `tabA', replace
		local first = 0
		restore
	}
}
* all underlying deaths (all years) for reference column
use `tabA', clear
replace agestr = "overall" if agestr=="all"
replace agestr = "<40" if agestr=="0"
replace agestr = "40+" if agestr=="1"
gen period = "`firstyr'-2022"
order grp group agestr period N underlying anymention ratio contribonly
sort grp agestr
export delimited using "$output/table_multiple_cause.csv", replace

* wide for docx
gen ustr = string(underlying)
gen astr = string(anymention)
gen rstr = cond(underlying>0, string(ratio,"%4.1f"), "-")
gen cstr = string(contribonly)
gen cpct = cond(anymention>0, string(100*contribonly/anymention,"%3.0f")+"%", "-")
gen cstr2 = cstr + " (" + cpct + ")"
gen s = cond(agestr=="overall",1,cond(agestr=="<40",2,3))
keep grp group s ustr astr rstr cstr2 N
reshape wide ustr astr rstr cstr2 N, i(grp group) j(s)
sort grp
local N1 = N1[1]
local N2 = N2[1]
local N3 = N3[1]
drop N1 N2 N3 grp
order group ustr1 astr1 rstr1 cstr21 ustr2 astr2 rstr2 cstr22 ustr3 astr3 rstr3 cstr23
tempfile tabAwide
save `tabAwide'

*-------------------------------------------------------------------------------
* 5. TABLE B: among deaths from recurrence, SPN, circulatory: % with
*    infection / respiratory / circulatory / genitourinary as contributory
*-------------------------------------------------------------------------------
use "$temp/x-mort3-multicause", clear
keep if mcc==1
gen ucat3 = .
replace ucat3 = 1 if recur==1
replace ucat3 = 2 if spn==1
replace ucat3 = 3 if circulation==1
label define lu3 1 "recurrence/progression" 2 "SPN" 3 "circulatory"
label values ucat3 lu3
tempfile tabB
local first = 1
forval u=1/3 {
	local uname : label lu3 `u'
	preserve
	keep if ucat3==`u'
	qui count
	local N = r(N)
	foreach g in 1 8 7 11 {
		qui count if contrib`g'==1
		local n`g' = r(N)
	}
	clear
	set obs 1
	gen str30 underlying = "`uname'"
	gen N = `N'
	foreach g in 1 8 7 11 {
		local gname : label lgrp `g'
		gen n_`g' = `n`g''
		gen pct_`g' = 100*`n`g''/`N'
		gen str20 s_`g' = string(`n`g'') + " (" + string(100*`n`g''/`N',"%4.1f") + "%)"
	}
	if `first'==0 append using `tabB'
	save `tabB', replace
	local first = 0
	restore
}
use `tabB', clear
gen ord = cond(underlying=="recurrence/progression",1,cond(underlying=="SPN",2,3))
sort ord
drop ord
rename (n_1 pct_1 n_8 pct_8 n_7 pct_7 n_11 pct_11) ///
	(n_infection pct_infection n_respiratory pct_respiratory n_circulatory pct_circulatory n_genitourinary pct_genitourinary)
export delimited underlying N n_* pct_* using "$output/table_multiple_cause_B.csv", replace
* for circulatory underlying, circulatory contributory is by definition 0 -> show n/a
replace s_7 = "n/a" if underlying=="circulatory"
keep underlying N s_1 s_8 s_7 s_11
gen Nstr = string(N)
order underlying Nstr s_1 s_8 s_7 s_11
drop N
tempfile tabBwide
save `tabBwide'

*-------------------------------------------------------------------------------
* 6. PUTDOCX
*-------------------------------------------------------------------------------
cap putdocx clear
putdocx begin, font("arial narrow", 10) landscape
putdocx paragraph
putdocx text ("Table A. Deaths with each cause group recorded as the underlying cause, mentioned anywhere on the death certificate (underlying or contributory), and as a contributory cause only, among deaths with multiple-cause coding (`firstyr'-2022; N=`N1' overall, `N2' at attained age <40 years, `N3' at 40+ years), BCCSS."), bold
use `tabAwide', clear
putdocx table tA = data("group ustr1 astr1 rstr1 cstr21 ustr2 astr2 rstr2 cstr22 ustr3 astr3 rstr3 cstr23"), ///
	border(all, nil) border(top) border(bottom)
putdocx table tA(1,.), addrows(2, before)
putdocx table tA(1,1) = ("Cause group")
putdocx table tA(1,2) = ("All ages (N=`N1')"), colspan(4) halign(center)
putdocx table tA(1,3) = ("Attained age <40 (N=`N2')"), colspan(4) halign(center)
putdocx table tA(1,4) = ("Attained age 40+ (N=`N3')"), colspan(4) halign(center)
forval b=0/2 {
	local c = 2 + 4*`b'
	putdocx table tA(2,`c')   = ("Underlying")
	putdocx table tA(2,`=`c'+1') = ("Any mention")
	putdocx table tA(2,`=`c'+2') = ("Ratio any/underlying")
	putdocx table tA(2,`=`c'+3') = ("Contributory only, n (% of any mention)")
}
putdocx table tA(1,.), bold border(bottom)
putdocx table tA(2,.), bold border(bottom)
putdocx table tA(.,1), halign(left)
forval c=2/13 {
	putdocx table tA(.,`c'), halign(center)
}
putdocx paragraph
putdocx text ("Contributory causes are available only on certificates coded with multiple causes (deaths from `firstyr'; ICD-9 1993-2000, ICD-10 2001-2022). 'Contributory only' = mentioned on the certificate but in a different group from the underlying cause. Underlying cause groups follow codcat; COVID-19 (ICD-10 U07) shown separately. Cardiac = ICD-10 I00-I52 / ICD-9 390-429; cerebrovascular = I60-I69 / 430-438.")

putdocx paragraph
putdocx text ("Table B. Among deaths from recurrence/progression, SPN and circulatory disease (underlying cause; deaths `firstyr'-2022), number (%) with infection, respiratory, circulatory and genitourinary (renal) disease recorded as a contributory cause."), bold
use `tabBwide', clear
putdocx table tB = data("underlying Nstr s_1 s_8 s_7 s_11"), border(all, nil) border(top) border(bottom)
putdocx table tB(1,.), addrows(1, before)
putdocx table tB(1,1) = ("Underlying cause")
putdocx table tB(1,2) = ("Deaths")
putdocx table tB(1,3) = ("Infection, n (%)")
putdocx table tB(1,4) = ("Respiratory, n (%)")
putdocx table tB(1,5) = ("Circulatory, n (%)")
putdocx table tB(1,6) = ("Genitourinary, n (%)")
putdocx table tB(1,.), bold border(bottom)
putdocx table tB(.,1), halign(left)
forval c=2/6 {
	putdocx table tB(.,`c'), halign(center)
}
putdocx save "$output/table_multiple_cause.docx", replace

*-------------------------------------------------------------------------------
* 7. FIGURE: horizontal bars, underlying vs any mention, main groups
*-------------------------------------------------------------------------------
import delimited "$output/table_multiple_cause.csv", clear varnames(1)
keep if agestr=="overall" & inlist(grp,1,2,3,4,5,6,7,8,9,10,11,14,15)
gsort -anymention
gen ord = _n
labmask ord, values(group)
graph hbar (asis) anymention underlying, over(ord, sort(ord) label(labsize(small))) ///
	bar(1, color(gs10)) bar(2, color(navy)) ///
	legend(order(2 "Underlying cause" 1 "Any mention on certificate") rows(1) pos(6) size(small)) ///
	ytitle("Number of deaths (`firstyr'-2022)", size(small)) ylabel(, labsize(small)) ///
	blabel(bar, size(vsmall) format(%4.0f)) ///
	graphregion(color(white)) plotregion(margin(small))
graph export "$output/fig_multiple_cause.png", replace width(2400)

exit
