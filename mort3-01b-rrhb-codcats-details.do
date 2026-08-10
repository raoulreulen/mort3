*-----------------------------------------------------------------
* Cuases of death in more details within each category*
*-----------------------------------------------------------------
// Cause of death categories sorted from highest observed number to lowest observed numbers// 
* 1st change the name of the label of codcat variable to "Cause of death category" instead of the current label "SPN Category"
label var codcat "Cause of death category"
tab codcat, sort 
tab codcat, nolabel sort

*---------------
* Circulation 
*---------------
tab ucause if codcat==7, sort
gen code_icd = ucause + " (ICD-" + string(icdver) + ")"
tab code_icd if codcat==7, sort
drop code_icd

*---------------
* Respiratory
*---------------

// This section looks into respiratory deaths (codcat==8) to identify the most common specific causes of death and their ICD version within this category// 

* count how many respiratory deaths (codcat==8) fall under each raw ucause code
tab ucause if codcat==8

* same as above, but sorted from most to least frequent 
tab ucause if codcat==8, sort

* check the split of ICD versions among respiratory deaths
tab icdver if codcat==8

* This combines ucause code abd ICD version into one variable so they can be tabulated together, then creates the table with frequencies sorted most to least common, then removed the temorary combined variable to keep data unchanged
gen code_icd = ucause + " (ICD-" + string(icdver) + ")"
tab code_icd if codcat==8, sort
drop code_icd

*-------------
* External
*-------------

// This section looks into external deaths (codcat==14) to identify the most common specific causes of death and their ICD version within this category//

tab ucause if codcat==14, sort

gen code_icd = ucause + " (ICD-" + string(icdver) + ")"
tab code_icd if codcat==14, sort
drop code_icd

*--------------
* Nervous 
*--------------
// This section looks into nervous causes of deaths (codcat==6) to identify the most common specific causes of death and their ICD version within this category//

tab ucause if codcat==6, sort
gen code_icd = ucause + " (ICD-" + string(icdver) + ")"
tab code_icd if codcat==6, sort
drop code_icd

*--------------
* Digestive 
*--------------

// This section looks into digestive causes of deaths (codcat==9) to identify the most common specific causes of death and their ICD version within this category//

tab ucause if codcat==9, sort

gen code_icd = ucause + " (ICD-" + string(icdver) + ")"
tab code_icd if codcat==9, sort
drop code_icd

*-------------
* Infection
*-------------
// This section looks into infection causes of deaths (codcat==1) to identify the most common specific causes of death and their ICD version within this category//

tab ucause if codcat==1, sort 
gen code_icd = ucause + " (ICD-" + string(icdver) + ")"
tab code_icd if codcat==1, sort
drop code_icd

*-------------
* Perinatal 
*-------------
 
 // This section looks into perinatal deaths (codcat==12) to identify the most common specific causes of death and their ICD version within this category//
tab ucause if codcat==12, sort
gen code_icd = ucause + " (ICD-" + string(icdver) + ")"
tab code_icd if codcat==12, sort
drop code_icd

*-------------
* Endocrine
*-------------
 
// This section looks into endocrine deaths (codcat==4) to identify the most common specific causes of death and their ICD version within this category//
 
tab ucause if codcat==4, sort
gen code_icd = ucause + " (ICD-" + string(icdver) + ")"
tab code_icd if codcat==4, sort
drop code_icd
 
*---------------
* Genitourinary 
*---------------
 
// This section looks into genitourinary deaths (codcat==11) to identify the most common specific causes of death and their ICD version within this category//
tab ucause if codcat==11, sort
gen code_icd = ucause + " (ICD-" + string(icdver) + ")"
tab code_icd if codcat==11, sort
drop code_icd

*-------------
* Mental 
*-------------
 
// This section looks into mental deaths (codcat==5) to identify the most common specific causes of death and their ICD version within this category//
 
tab ucause if codcat==5, sort
gen code_icd = ucause + " (ICD-" + string(icdver) + ")"
tab code_icd if codcat==5, sort
drop code_icd

*------------------
* Musculoskeletal
*------------------
 
// This section looks into musculoskeletal/skin deaths (codcat==10) to identify the most common specific causes of death and their ICD version within this category//
 
tab ucause if codcat==10, sort
gen code_icd = ucause + " (ICD-" + string(icdver) + ")"
tab code_icd if codcat==10, sort
drop code_icd
 
*-----------------
* Other
*-----------------
 
// This section looks into other/unspecified deaths (codcat==13) to identify the most common specific causes of death and their ICD version within this category//
tab ucause if codcat==13, sort
gen code_icd = ucause + " (ICD-" + string(icdver) + ")"
tab code_icd if codcat==13, sort
drop code_icd

*----------------
* Blood
*----------------
 
// This section looks into blood-related deaths (codcat==3) to identify the most common specific causes of death and their ICD version within this category//
tab ucause if codcat==3, sort
gen code_icd = ucause + " (ICD-" + string(icdver) + ")"
tab code_icd if codcat==3, sort
drop code_icd

*----------------
* Unknown
*----------------

tab ucause if codcat==0, sort
gen code_icd = ucause + " (ICD-" + string(icdver) + ")"
tab code_icd if codcat==0, sort
drop code_icd


*-----------------------------------------------------------------------------
* EXPORT ALL CATEGORY DRILLDOWN TABLES TO ONE WORD DOCUMENT
*-----------------------------------------------------------------------------
putdocx clear
putdocx begin
putdocx paragraph, style(Title)
putdocx text ("Cause of Death Detail by Category")

levelsof codcat, local(catlist)
foreach c of local catlist {
    local label : label CAT `c'
    
    * count total deaths in this category, store in a local macro
    quietly count if codcat==`c'
    local n = r(N)
    
    putdocx paragraph, style(Heading2)
    putdocx text ("`label' (codcat==`c')")
    
    * add a line showing the total, as a sanity check
    putdocx paragraph
    putdocx text ("Total deaths in this category: `n'")
    
    frame put ucause icdver if codcat==`c', into(tempcat)
    frame tempcat: gen code_icd = ucause + " (ICD-" + string(icdver) + ")"
    frame tempcat: contract code_icd, freq(count)
    frame tempcat: gsort -count
    frame tempcat: putdocx table tbl_`c' = data(code_icd count), varnames
    frame drop tempcat
}

putdocx save "$temp/all_cod_drilldown.docx", replace

*-----------------------------------------------------------------------------
* EXPORT ALL CATEGORY DRILLDOWN TABLES TO ONE EXCEL FILE (separate sheet each)
*-----------------------------------------------------------------------------
levelsof codcat, local(catlist)
foreach c of local catlist {
    local label : label CAT `c'
    
    frame put ucause icdver if codcat==`c', into(tempcat)
    frame tempcat: gen code_icd = ucause + " (ICD-" + string(icdver) + ")"
    frame tempcat: contract code_icd, freq(count)
    frame tempcat: gsort -count
    frame tempcat: export excel code_icd count using "$temp/all_cod_drilldown.xlsx", ///
        sheet("`label'") sheetreplace firstrow(variables)
    frame drop tempcat
}
