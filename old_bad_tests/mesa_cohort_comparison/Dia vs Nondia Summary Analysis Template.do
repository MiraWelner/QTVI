clear all
set more off

*------------------------------------------------------------------
* SET FOLDER PATHS  // Update as per your location
*------------------------------------------------------------------

* Do not include a trailing "\" at the end of these paths
local diabetic_dir    "F:\MESA Data analysis\v2.4\Output_diabetic\qtvi_marker_path"
local nondiabetic_dir "F:\MESA Data analysis\v2.4\Output_nondiabetic\qtvi_marker_path"
local output_dir      "F:\MESA Data analysis\v2.4\STATA Output"


local raw_pattern "*_template_markings.csv"

capture mkdir "`output_dir'"

* Temporary master datasets
tempfile diabetic_master nondiabetic_master


*------------------------------------------------------------------
* 1. APPEND ALL DIABETIC CSV FILES
*------------------------------------------------------------------

cd "`diabetic_dir'"

local diabetic_files : dir "." files "`raw_pattern'"

* Stop if no matching files are found
if trim(`"`diabetic_files'"') == "" {
    display as error "No diabetic CSV files matching `raw_pattern' were found in:"
    display as error "`diabetic_dir'"
    exit 601
}

display as result "Diabetic files found:"
display `"`diabetic_files'"'

local n_diabetic : word count `diabetic_files'
display as result "Number of diabetic files: `n_diabetic'"

* Import the first diabetic file
local first : word 1 of `diabetic_files'

display as text "Importing diabetic file 1 of `n_diabetic': `first'"

import delimited "`first'", clear varnames(1)

* Remove these variables if they already exist in an input file
capture drop group
capture drop sourcefile

gen str20 group = "Diabetic"
gen strL sourcefile = "`first'"

save `diabetic_master', replace


* Import and append remaining diabetic files
if `n_diabetic' > 1 {

    forvalues i = 2/`n_diabetic' {

        local f : word `i' of `diabetic_files'

        display as text ///
            "Importing diabetic file `i' of `n_diabetic': `f'"

        import delimited "`f'", clear varnames(1)

        capture drop group
        capture drop sourcefile

        gen str20 group = "Diabetic"
        gen strL sourcefile = "`f'"

        append using `diabetic_master'
        save `diabetic_master', replace
    }
}

* Save diabetic master files
use `diabetic_master', clear

save "`output_dir'\diabetic_master.dta", replace

export delimited using ///
    "`output_dir'\diabetic_master.csv", replace

display as result ///
    "Diabetic master completed: " _N " observations."


*------------------------------------------------------------------
* 2. APPEND ALL NONDIABETIC CSV FILES
*------------------------------------------------------------------

cd "`nondiabetic_dir'"

local nondiabetic_files : dir "." files "`raw_pattern'"

* Stop if no matching files are found
if trim(`"`nondiabetic_files'"') == "" {
    display as error ///
        "No non-diabetic CSV files matching `raw_pattern' were found in:"
    display as error "`nondiabetic_dir'"
    exit 601
}

display as result "Non-diabetic files found:"
display `"`nondiabetic_files'"'

local n_nondiabetic : word count `nondiabetic_files'
display as result ///
    "Number of non-diabetic files: `n_nondiabetic'"

* Import the first non-diabetic file
local first : word 1 of `nondiabetic_files'

display as text ///
    "Importing non-diabetic file 1 of `n_nondiabetic': `first'"

import delimited "`first'", clear varnames(1)

* Remove these variables if they already exist in an input file
capture drop group
capture drop sourcefile

gen str20 group = "NonDiabetic"
gen strL sourcefile = "`first'"

save `nondiabetic_master', replace


* Import and append remaining non-diabetic files
if `n_nondiabetic' > 1 {

    forvalues i = 2/`n_nondiabetic' {

        local f : word `i' of `nondiabetic_files'

        display as text ///
            "Importing non-diabetic file `i' of `n_nondiabetic': `f'"

        import delimited "`f'", clear varnames(1)

        capture drop group
        capture drop sourcefile

        gen str20 group = "NonDiabetic"
        gen strL sourcefile = "`f'"

        append using `nondiabetic_master'
        save `nondiabetic_master', replace
    }
}

* Save non-diabetic master files
use `nondiabetic_master', clear

save "`output_dir'\nondiabetic_master.dta", replace

export delimited using ///
    "`output_dir'\nondiabetic_master.csv", replace

display as result ///
    "Non-diabetic master completed: " _N " observations."


*------------------------------------------------------------------
* 3. COMBINE DIABETIC AND NONDIABETIC GROUPS
*------------------------------------------------------------------

use `diabetic_master', clear
append using `nondiabetic_master'

* Confirm that the grouping variable was created correctly
tabulate group

* Make QRS numeric if it was imported as a string
capture confirm numeric variable qrs_ch1_ms

if _rc {
    destring qrs_ch1_ms, replace force
}

* Drop observations with missing QRS duration
drop if missing(qrs_ch1_ms)

* Optional: sort the combined dataset
sort group file_id bin_index

* Save final combined dataset
save "`output_dir'\diabetes_qrs_master.dta", replace

export delimited using ///
    "`output_dir'\diabetes_qrs_master.csv", replace

	
display as result "Combined master dataset completed."
display as result "Final number of observations: " _N	
	
*------------------------------------------------------------------
* 4. STAT ANALSIS OUTPUT
*------------------------------------------------------------------


// tabulate group
// summarize qrs_ch1_ms, detail
// bysort group: summarize qt_ch1_ms, detail
//
// ttest qrs_ch1_ms, by(group)

// tabstat qrs_ch1_ms, ///
//     by(group) ///
//     statistics(mean median sd p25 p75 min max n)
//
// bysort group: tabstat qrs_ch1_ms, ///
//     by(file_id) ///
//     statistics(mean median sd p25 p75 min max n)
//	
//
// tabstat qt_ch1_ms, ///
//     by(group) ///
//     statistics(mean median sd p25 p75 min max n)
//
// tabstat qt_ch1_ms, ///
//     by(file_id) ///
//     statistics(mean median sd p25 p75 min max n)

// bysort group: tabstat qrs_ch1_ms, ///
//     by(file_id) ///
// 	statistics(mean median sd p25 p75 min max n)


// PPG Var Gen  //
gen ppg_amp_y_mv_user = ppg_peak_y_mv_user - ppg_onset_y_mv_user
gen ppg_width_x_ms_user = ppg_end_x_ms_user - ppg_onset_x_ms_user

****************************************************

local vars ///
    qrs_ch1_ms ///
    qt_ch1_ms ///
    q_peak_ch1_y_mv ///
    r_peak_ch1_y_mv ///
    s_peak_ch1_y_mv ///
    t_peak_ch1_y_mv ///
	ppg_amp_y_mv_user ///
	ppg_width_x_ms_user

foreach v of local vars {

    di "======================================================"
    di "Variable: `v'"
    di "======================================================"

    bysort group: tabstat `v', ///
        by(file_id) ///
        statistics(mean median sd p25 p75 min max n)

    *----------------------------------
    * Two-sample t-test
    *----------------------------------

    ttest `v', by(group)

    di ""
    di ""

}




// ranksum qrs_ch1_ms, by(group)
	
*------------------------------------------------------------------
* 5. OVERLAID HISTOGRAM COMPARISON
*------------------------------------------------------------------

twoway ///
    (histogram qrs_ch1_ms if group == "Diabetic", ///
        percent ///
        color(cranberry%35) ///
        lcolor(cranberry%80) ///
        start(60) width(5)) ///
    (histogram qrs_ch1_ms if group == "NonDiabetic", ///
        percent ///
        color(navy%35) ///
        lcolor(navy%80) ///
        start(60) width(5)), ///
    legend(order(1 "Diabetic (5)" 2 "Non-diabetic (3)") ///
        position(1) ring(0)) ///
    xtitle("QRS duration (ms)") ///
    ytitle("Percent") ///
    title("QRS Duration Distribution") ///
    ylabel(, nogrid) ///
	xlabel(, nogrid) ///
    graphregion(color(white)) ///
    name(qrs_hist_compare, replace)

// graph export "`output_dir'\QRS_histogram_diabetic_vs_nondiabetic.png", ///
//     replace width(2400)
	
	
twoway ///
    (histogram qt_ch1_ms if group == "Diabetic", ///
        percent ///
        color(cranberry%35) ///
        lcolor(cranberry%80) ///
        start(100) width(5)) ///
    (histogram qt_ch1_ms if group == "NonDiabetic", ///
        percent ///
        color(navy%35) ///
        lcolor(navy%80) ///
        start(100) width(5)), ///
    legend(order(1 "Diabetic (5)" 2 "Non-diabetic (3)") ///
        position(1) ring(0)) ///
    xtitle("QT duration (ms)") ///
    ytitle("Percent") ///
    title("QT Duration Distribution") ///
    ylabel(, nogrid) ///
	xlabel(, nogrid) ///
    graphregion(color(white)) ///
    name(qt_hist_compare, replace)

// graph export "`output_dir'\QT_histogram_diabetic_vs_nondiabetic.png", ///
//     replace width(2400)


twoway ///
    (histogram q_peak_ch1_y_mv if group == "Diabetic", ///
        percent ///
        color(cranberry%35) ///
        lcolor(cranberry%80) ///
        start() width(0.2)) ///
    (histogram q_peak_ch1_y_mv if group == "NonDiabetic", ///
        percent ///
        color(navy%35) ///
        lcolor(navy%80) ///
        start() width(0.2)), ///
    legend(order(1 "Diabetic (5)" 2 "Non-diabetic (3)") ///
        position(1) ring(0)) ///
    xtitle("Q_peak height (mv)") ///
    ytitle("Percent") ///
    title("Q_peak Height Distribution") ///
    ylabel(, nogrid) ///
	xlabel(, nogrid) ///
    graphregion(color(white)) ///
    name(q_peak_hist_compare, replace)
	
	

twoway ///
    (histogram r_peak_ch1_y_mv if group == "Diabetic", ///
        percent ///
        color(cranberry%35) ///
        lcolor(cranberry%80) ///
        start() width(0.2)) ///
    (histogram r_peak_ch1_y_mv if group == "NonDiabetic", ///
        percent ///
        color(navy%35) ///
        lcolor(navy%80) ///
        start() width(0.2)), ///
    legend(order(1 "Diabetic (5)" 2 "Non-diabetic (3)") ///
        position(1) ring(0)) ///
    xtitle("R_peak height (mv)") ///
    ytitle("Percent") ///
    title("R_peak Height Distribution") ///
    ylabel(, nogrid) ///
	xlabel(, nogrid) ///
    graphregion(color(white)) ///
    name(r_peak_hist_compare, replace)

	
twoway ///
    (histogram s_peak_ch1_y_mv if group == "Diabetic", ///
        percent ///
        color(cranberry%35) ///
        lcolor(cranberry%80) ///
        start() width(0.2)) ///
    (histogram s_peak_ch1_y_mv if group == "NonDiabetic", ///
        percent ///
        color(navy%35) ///
        lcolor(navy%80) ///
        start() width(0.2)), ///
    legend(order(1 "Diabetic (5)" 2 "Non-diabetic (3)") ///
        position(1) ring(0)) ///
    xtitle("S_peak height (mv)") ///
    ytitle("Percent") ///
    title("S_peak Height Distribution") ///
    ylabel(, nogrid) ///
	xlabel(, nogrid) ///
    graphregion(color(white)) ///
    name(s_peak_hist_compare, replace)
	

twoway ///
    (histogram t_peak_ch1_y_mv if group == "Diabetic", ///
        percent ///
        color(cranberry%35) ///
        lcolor(cranberry%80) ///
        start() width(0.2)) ///
    (histogram t_peak_ch1_y_mv if group == "NonDiabetic", ///
        percent ///
        color(navy%35) ///
        lcolor(navy%80) ///
        start() width(0.2)), ///
    legend(order(1 "Diabetic (5)" 2 "Non-diabetic (3)") ///
        position(1) ring(0)) ///
    xtitle("T_peak height (mv)") ///
    ytitle("Percent") ///
    title("T_peak Height Distribution") ///
    ylabel(, nogrid) ///
	xlabel(, nogrid) ///
    graphregion(color(white)) ///
    name(t_peak_hist_compare, replace)


twoway ///
    (histogram ppg_amp_y_mv_user if group == "Diabetic", ///
        percent ///
        color(cranberry%35) ///
        lcolor(cranberry%80) ///
        start() width(0.05)) ///
    (histogram ppg_amp_y_mv_user if group == "NonDiabetic", ///
        percent ///
        color(navy%35) ///
        lcolor(navy%80) ///
        start() width(0.05)), ///
    legend(order(1 "Diabetic (5)" 2 "Non-diabetic (3)") ///
        position(1) ring(0)) ///
    xtitle("PPG_amp height (mv)") ///
    ytitle("Percent") ///
    title("PPG Amplitude Distribution") ///
    ylabel(, nogrid) ///
	xlabel(, nogrid) ///
    graphregion(color(white)) ///
    name(PPG_amp_hist_compare, replace)	
	

twoway ///
    (histogram ppg_width_x_ms_user if group == "Diabetic", ///
        percent ///
        color(cranberry%35) ///
        lcolor(cranberry%80) ///
        start() width(5)) ///
    (histogram ppg_width_x_ms_user if group == "NonDiabetic", ///
        percent ///
        color(navy%35) ///
        lcolor(navy%80) ///
        start() width(5)), ///
    legend(order(1 "Diabetic (5)" 2 "Non-diabetic (3)") ///
        position(1) ring(0)) ///
    xtitle("PPG width duration (ms)") ///
    ytitle("Percent") ///
    title("PPG Width Duration Distribution") ///
    ylabel(, nogrid) ///
	xlabel(, nogrid) ///
    graphregion(color(white)) ///
    name(PPG_width_hist_compare, replace)	
	
///////////////////////////////////////////


