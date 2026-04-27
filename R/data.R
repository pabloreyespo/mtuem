#' enut-i
#'
#' Processed dataset from the first National Time-Use Survey (ENUT I), applied by the
#' Instituto Nacional de Estadísticas de Chile in 2015. Contains both the original 11
#' aggregated time-use activity categories and the new 10-category time allocation
#' structure, along with aggregated household expenditures imputed from the VIII Encuesta
#' de Presupuestos Familiares (EPF). Used as the primary input for structural time-use
#' models via \code{get_data()} and \code{get_data_tc()}.
#'
#' All income and expenditure variables are expressed in weekly thousands of Chilean pesos.
#' Time variables are expressed in weekly hours, normalized to sum to 168.
#'
#' \describe{
#'   \item{id_persona}{Individual identifier}
#'   \item{id_hogar}{Household identifier}
#'
#'   \item{es_trabajador}{1 if employed, CAE = "Ocupada(o)", and positive paid work time}
#'   \item{es_familia}{1 if es_trabajador, lives with partner, and has children in household}
#'
#'   \item{dia_semana}{Weekday of diary (1 = Monday to 5 = Friday)}
#'   \item{dia_fin_semana}{Weekend day of diary (6 = Saturday, 7 = Sunday)}
#'
#'   \item{parentesco}{Relationship to household head}
#'   \item{n_menores_0_5}{Number of household members under age 6 (ages 0-5)}
#'   \item{n_menores_6_11}{Number of household members aged 6-11}
#'   \item{n_menores_0_14}{Number of household members under age 15 (ages 0-14)}
#'   \item{n_menores_12_17}{Number of household members aged 12-17}
#'   \item{n_menores}{Number of household members under 18}
#'   \item{n_mayores}{Number of adult household members (18+)}
#'   \item{n_tiempo}{Number of household members who reported time use}
#'   \item{n_trabajadores}{Number of employed workers in household}
#'   \item{n_profesionales}{Number of household members with completed university or higher}
#'   \item{n_tercera_edad}{Number of household members aged 60+}
#'   \item{hay_tercera_edad}{1 if household contains at least one elderly member aged 60+ who is not the respondent}
#'   \item{n_personas}{Total household size}
#'   \item{edad_promedio}{Mean age of household members}
#'   \item{tiene_hijos}{1 if respondent has children living in the household}
#'   \item{en_pareja}{1 if respondent is in a couple relationship (married or cohabiting)}
#'   \item{vive_pareja}{1 if respondent lives with their partner}
#'
#'   \item{servicio_domestico}{1 if household receives paid domestic service}
#'   \item{ayuda_cercanos}{1 if household receives unpaid help from relatives or neighbours}
#'   \item{fuentes_externas}{1 if household receives any external care support (servicio_domestico or ayuda_cercanos)}
#'
#'   \item{sexo}{Female indicator (1 = female, 0 = male; recoded from raw variable)}
#'   \item{edad_anios}{Age in years}
#'   \item{tramo_edad}{Age bracket: "12-24", "25-44", "45-65", or "66+"}
#'   \item{nivel_escolaridad}{Highest completed education level: 1 = ninguna, 2 = básica,
#'     3 = media, 4 = técnica, 5 = universitaria o más}
#'   \item{estudia}{1 if currently enrolled in education}
#'   \item{trabaja}{1 if currently employed (from k11_1_1)}
#'   \item{horas_trabajo_contratadas}{Horas semanales contratadas segun contrato de trabajo (from k31_1_3)}
#'   \item{horas_trabajo_habituales}{Horas habituales de trabajo en la ocupacion principal (from k31_1_1);
#'     compare with \code{t_paid_work} (diary-measured) and \code{horas_trabajo_contratadas} (contracted)}
#'   \item{dias_trabajo_semana}{Dias trabajados por semana en la ocupacion principal (from k31_1_2);
#'     used internally to scale weekday diary hours during weekend imputation}
#'   \item{quintil}{Household income quintile (1-5)}
#'   \item{macrozona}{Geographic macro-zone: "norte", "metropolitana", "centro", or "sur"}
#'   \item{region}{Region code (numeric)}
#'   \item{prop_ing_hogar}{Respondent's share of total personal income in the household; used to allocate household-level expenditures to individuals}
#'
#'   \item{cae}{Economic activity category: "Ocupada(o)", "Desocupada(o)", "Inactiva(o)", "Menor de 15 años", or "Sin clasificar"}
#'   \item{cise}{Employment status per CISE classification}
#'   \item{ciuo_agrupada}{Occupation group per CIUO-08 grouped classification}
#'
#'   \item{ing_ocuppal}{Income from main occupation (weekly, thousands CLP)}
#'   \item{ing_trab}{Total labor income (weekly, thousands CLP)}
#'   \item{ing_jub_aps}{Pension and AFP income (weekly, thousands CLP)}
#'   \item{ing_g}{Imputed group income component from household total (weekly, thousands CLP)}
#'   \item{ing_mon}{Monetary transfers and other non-labor income (weekly, thousands CLP)}
#'   \item{ing_mon_pc}{Per capita monetary income (weekly, thousands CLP)}
#'   \item{ing_gpp}{Individual share of group income (weekly, thousands CLP)}
#'   \item{ing_personal}{Personal income: ing_trab + ing_jub_aps + ing_gpp (weekly, thousands CLP)}
#'   \item{ingreso_hogar}{Total household disposable income (weekly, thousands CLP)}
#'   \item{income_person_week}{Household income divided by number of members (weekly, thousands CLP)}
#'
#'   \item{t_paid_work}{Trabajo remunerado; equivale a t_to, horas semanales}
#'   \item{t_job_search}{Busqueda de trabajo; equivale a t_to_js, horas semanales}
#'   \item{t_domestic_work}{Trabajo domestico no remunerado; suma de t_tdnr_psc + t_tdnr_lv +
#'     t_tdnr_lrc + t_tdnr_mrm + t_tdnr_admnhog + t_tdnr_comphog + t_tdnr_cmp, horas semanales}
#'   \item{t_care_work}{Trabajo de cuidado no remunerado; suma de t_tcnr_ce + t_tcnr_re +
#'     t_tcnr_oac, horas semanales}
#'   \item{t_unpaid_voluntary}{Trabajo voluntario y ayuda a otros hogares; suma de t_tvaoh_tv +
#'     t_tvaoh_oh, horas semanales}
#'   \item{t_education}{Educacion y formacion; equivale a t_ed, horas semanales}
#'   \item{t_leisure}{Ocio; suma de t_vsyo_csar + t_vsyo_aa + t_mcm_leer + t_mcm_video +
#'     t_mcm_audio + t_mcm_computador, horas semanales}
#'   \item{t_personal_care}{Cuidados personales fisiologicos (excluye sueno y comidas); equivale
#'     a t_cpaf_cp, horas semanales}
#'   \item{t_meals}{Comer y beber; equivale a t_cpag_comer, horas semanales}
#'   \item{t_sleep}{Dormir; equivale a t_cpag_dormir, ajustado para que la suma sea 168 horas,
#'     horas semanales}
#'   \item{t_commute1}{Traslados asociados a trabajo remunerado, educacion y salud; equivale a
#'     t_tt1, horas semanales}
#'
#'   \item{Tw}{Paid work time (equivalent to t_to / t_paid_work)}
#'   \item{Tf_social}{Social life and recreation time (equivalent to t_vsyo_csar)}
#'   \item{Tf_hobbies}{Hobbies and arts time (equivalent to t_vsyo_aa)}
#'   \item{Tf_read}{Reading time (equivalent to t_mcm_leer)}
#'   \item{Tf_listen}{Audio consumption time (equivalent to t_mcm_audio)}
#'   \item{Tf_watch}{TV and video consumption time (equivalent to t_mcm_video)}
#'   \item{Tf_computer}{Recreational computer/internet use time (equivalent to t_mcm_computador)}
#'   \item{Tc_meals}{Time spent eating and drinking (equivalent to t_cpag_comer)}
#'   \item{Tc_sleep}{Time spent sleeping, adjusted to balance 168 hours (equivalent to t_cpag_dormir)}
#'   \item{Tc_other}{All other time use (job search, domestic work, care, personal care,
#'     volunteering, commuting, education)}
#'
#'   \item{t_total}{Total weekly hours across all activities (should equal 168)}
#'   \item{w}{Hourly wage rate: ing_trab / Tw (thousands CLP per hour)}
#'
#'   \item{Ef_food}{Imputed food expenditure (weekly thousands CLP)}
#'   \item{Ef_recreation}{Imputed recreation and culture expenditure (weekly thousands CLP)}
#'   \item{Ef_restaurants}{Imputed restaurants and food away from home expenditure (weekly thousands CLP)}
#'   \item{Ef_communications}{Imputed communications expenditure (weekly thousands CLP)}
#'   \item{Ef_clothing}{Imputed clothing expenditure (weekly thousands CLP)}
#'   \item{Ec}{Imputed committed and remaining expenditures (weekly thousands CLP, includes
#'     utilities, health, transport, education, household goods, and savings)}
#' }
#'
#' @details
#' The dataset is produced by running \code{data_processing/data_processing.R}. Time
#' variables are normalized to a 168-hour week using a two-step procedure: paid work
#' (\code{t_paid_work}) and sleep (\code{t_sleep}) are treated as reliable anchors;
#' all other activities are scaled proportionally. Weekend time use is imputed via a
#' twin-matching matrix when the respondent's diary day was a weekday. Outliers are
#' removed using Vallejo's method by groups of quintile, employment status, age
#' bracket, and sex. Expenditures are imputed from EPF VIII using a fractional MNL
#' model for budget shares and a linear regression for the savings rate, then
#' allocated to individuals by \code{prop_ing_hogar}. The script writes
#' \code{data/enut-i.dta} plus \code{data/enut-i.csv}. English copies are written as
#' \code{data/enut-i-ENG.dta} and \code{data/enut-i-ENG.csv}.
#'
#' The time-use categories are aggregations of the 25 detailed categories in
#' \code{enut_i_raw}. See \code{agregar_actividades()} in
#' \code{data_processing/processing_functions.R} for the exact aggregation mapping.
#'
#' Compared to ENUT II, this dataset includes \code{t_job_search} as a separate activity
#' category (absent in ENUT II) and a single commute category (\code{t_commute1}) since
#' ENUT I does not distinguish a second commute type. \code{Tc_other} therefore also
#' absorbs job search time.
#'
#' @source <https://www.ine.gob.cl/enut>
#' @source <https://www.ine.gob.cl/estadisticas/sociales/ingresos-y-gastos/encuesta-de-presupuestos-familiares>
#'
#' @docType data
#' @keywords datasets
#' @name enut_i
#' @usage data(enut_i)
#' @format A data frame with 9,497 rows and 98 variables
NULL

#' enut-i-raw
#'
#' Processed raw activity dataset from the first National Time-Use Survey
#' (ENUT I), applied by the Instituto Nacional de Estadisticas de Chile in
#' 2015. Contains the 25 detailed weekly time-use activity categories and
#' household expenditures imputed from the VIII Encuesta de Presupuestos
#' Familiares (EPF). This dataset is the detailed companion to \code{enut_i},
#' which aggregates the same records into model-ready time categories.
#'
#' All income and expenditure variables are expressed in weekly thousands of
#' Chilean pesos. Time variables are expressed in weekly hours and normalized to
#' sum to 168.
#'
#' \describe{
#'   \item{id_persona}{Individual identifier}
#'   \item{id_hogar}{Household identifier}
#'
#'   \item{es_trabajador}{1 if employed, CAE = "Ocupada(o)", and positive paid work time}
#'   \item{es_familia}{1 if es_trabajador, lives with partner, and has children in household}
#'
#'   \item{dia_semana}{Weekday of diary (1 = Monday to 5 = Friday)}
#'   \item{dia_fin_semana}{Weekend day of diary (6 = Saturday, 7 = Sunday)}
#'
#'   \item{parentesco}{Relationship to household head}
#'   \item{n_menores_0_5}{Number of household members under age 6}
#'   \item{n_menores_6_11}{Number of household members aged 6 to 11}
#'   \item{n_menores_0_14}{Number of household members under age 15}
#'   \item{n_menores_12_17}{Number of household members aged 12 to 17}
#'   \item{n_menores}{Number of household members under 18}
#'   \item{n_mayores}{Number of adult household members}
#'   \item{n_tiempo}{Number of household members who reported time use}
#'   \item{n_trabajadores}{Number of employed workers in household}
#'   \item{n_profesionales}{Number of household members with completed university or higher}
#'   \item{n_tercera_edad}{Number of household members aged 60 or older}
#'   \item{hay_tercera_edad}{1 if household contains at least one older person who is not the respondent}
#'   \item{n_personas}{Total household size}
#'   \item{edad_promedio}{Mean age of household members}
#'   \item{tiene_hijos}{1 if respondent has children living in the household}
#'   \item{en_pareja}{1 if respondent is married or cohabiting}
#'   \item{vive_pareja}{1 if respondent lives with their partner}
#'
#'   \item{servicio_domestico}{1 if household receives paid domestic service}
#'   \item{ayuda_cercanos}{1 if household receives unpaid help from relatives or neighbours}
#'   \item{fuentes_externas}{1 if household receives any external care support}
#'
#'   \item{sexo}{Female indicator, recoded from the raw survey variable}
#'   \item{edad_anios}{Age in years}
#'   \item{tramo_edad}{Age bracket: "12-24", "25-44", "45-65", or "66+"}
#'   \item{nivel_escolaridad}{Highest completed education level}
#'   \item{estudia}{1 if currently enrolled in education}
#'   \item{trabaja}{1 if currently employed}
#'   \item{horas_trabajo_contratadas}{Weekly contracted hours in main job}
#'   \item{horas_trabajo_habituales}{Usual weekly hours in main job}
#'   \item{dias_trabajo_semana}{Work days per week in main job}
#'   \item{quintil}{Household income quintile}
#'   \item{macrozona}{Geographic macro-zone}
#'   \item{region}{Region code}
#'   \item{prop_ing_hogar}{Respondent share of total personal income in the household}
#'
#'   \item{cae}{Economic activity category}
#'   \item{cise}{Employment status per CISE classification}
#'   \item{ciuo_agrupada}{Occupation group per CIUO-08 grouped classification}
#'
#'   \item{t11_1_1:t15_1_1}{Original ENUT I well-being and household distribution variables retained from the source survey}
#'
#'   \item{ing_ocuppal}{Income from main occupation}
#'   \item{ing_trab}{Total labor income}
#'   \item{ing_jub_aps}{Pension and AFP income}
#'   \item{ing_g}{Imputed group income component from household total}
#'   \item{ing_mon}{Monetary transfers and other non-labor income}
#'   \item{ing_mon_pc}{Per capita monetary income}
#'   \item{ing_gpp}{Individual share of group income}
#'   \item{ing_personal}{Personal income}
#'   \item{ingreso_hogar}{Total household disposable income}
#'   \item{income_person_week}{Household income divided by number of members}
#'
#'   \item{t_to}{Paid work}
#'   \item{t_to_js}{Job search}
#'   \item{t_tcnr_ce}{Unpaid care work, essential care}
#'   \item{t_tcnr_re}{Unpaid care work, education-related care}
#'   \item{t_tcnr_oac}{Unpaid care work, other care}
#'   \item{t_tdnr_psc}{Unpaid domestic work, meal preparation and service}
#'   \item{t_tdnr_lv}{Unpaid domestic work, dwelling cleaning}
#'   \item{t_tdnr_lrc}{Unpaid domestic work, laundry and clothing repair}
#'   \item{t_tdnr_mrm}{Unpaid domestic work, household maintenance and minor repairs}
#'   \item{t_tdnr_admnhog}{Unpaid domestic work, household administration}
#'   \item{t_tdnr_comphog}{Unpaid domestic work, household shopping}
#'   \item{t_tdnr_cmp}{Unpaid domestic work, pets and plants}
#'   \item{t_tvaoh_tv}{Voluntary work in the community}
#'   \item{t_tvaoh_oh}{Direct help to other households}
#'   \item{t_cpaf_cp}{Personal care, physiological care excluding sleep and meals}
#'   \item{t_cpag_comer}{Eating and drinking}
#'   \item{t_cpag_dormir}{Sleeping, adjusted to balance the 168 hour total}
#'   \item{t_ed}{Education and training}
#'   \item{t_vsyo_csar}{Social life and recreation}
#'   \item{t_vsyo_aa}{Arts and hobbies}
#'   \item{t_mcm_leer}{Reading}
#'   \item{t_mcm_video}{Video and television}
#'   \item{t_mcm_audio}{Audio media}
#'   \item{t_mcm_computador}{Recreational computer and internet use}
#'   \item{t_tt1}{Commute and travel time}
#'
#'   \item{t_total}{Total weekly hours across all detailed activities}
#'   \item{w}{Hourly wage rate: ing_trab / t_to}
#'
#'   \item{savings}{Imputed household savings allocated by income share}
#'   \item{alimentos}{Imputed food expenditure}
#'   \item{vestimenta}{Imputed clothing expenditure}
#'   \item{cuentas}{Imputed utility and fixed household expenditure}
#'   \item{hogar}{Imputed household goods and services expenditure}
#'   \item{salud}{Imputed health expenditure}
#'   \item{transporte}{Imputed transportation expenditure}
#'   \item{comunicaciones}{Imputed communications expenditure}
#'   \item{recreacion}{Imputed recreation and culture expenditure}
#'   \item{educacion}{Imputed education expenditure}
#'   \item{restaurantes}{Imputed restaurants and food away from home expenditure}
#'   \item{total_expenses}{Total imputed expenditure: ing_personal - savings}
#' }
#'
#' @details
#' The dataset is produced by running \code{data_processing/data_processing.R}.
#' The script creates detailed time-use categories with \code{agregar_actividades()},
#' imputes weekend diaries through the twin matrix, removes outliers using
#' Vallejo's method, imputes EPF VIII expenditures with \code{imputacion_gastos()},
#' and writes \code{data/enut-i-raw.dta} plus \code{data/enut-i-raw.csv}. English
#' copies are written as \code{data/enut-i-raw-ENG.dta} and
#' \code{data/enut-i-raw-ENG.csv}.
#'
#' @source <https://www.ine.gob.cl/enut>
#' @source <https://www.ine.gob.cl/estadisticas/sociales/ingresos-y-gastos/encuesta-de-presupuestos-familiares>
#'
#' @docType data
#' @keywords datasets
#' @name enut_i_raw
#' @usage data(enut_i_raw)
#' @format A data frame with 9,497 rows and 108 variables
NULL

#' enut-ii
#'
#' Processed dataset from the second National Time-Use Survey (ENUT II), applied by the
#' Instituto Nacional de Estadísticas de Chile. Contains both the original 11 aggregated
#' time-use activity categories and the new 10-category time allocation structure, along with
#' aggregated household expenditures imputed from the IX Encuesta de Presupuestos
#' Familiares (EPF). Used as the primary input for structural time-use models via
#' \code{get_data()} and \code{get_data_tc()}.
#'
#' All income and expenditure variables are expressed in weekly thousands of Chilean pesos,
#' deflated by IPC adjustment (factor 0.362). Time variables are expressed in weekly hours,
#' normalized to sum to 168.
#'
#' \describe{
#'   \item{id_persona}{Individual identifier}
#'   \item{id_hog}{Household identifier}
#'
#'   \item{es_trabajador}{1 if employed, CAE = "Ocupada(o)", and positive paid work time}
#'   \item{es_familia}{1 if es_trabajador, lives with partner, and has children in household}
#'
#'   \item{dia_semana}{Weekday of diary (1 = Monday to 5 = Friday)}
#'   \item{dia_fin_semana}{Weekend day of diary (6 = Saturday, 7 = Sunday)}
#'
#'   \item{parentesco}{Relationship to household head (from pco)}
#'   \item{n_menores_0_5}{Number of household members under age 6 (ages 0-5)}
#'   \item{n_menores_6_11}{Number of household members aged 6-11}
#'   \item{n_menores_0_14}{Number of household members under age 15 (ages 0-14)}
#'   \item{n_menores_12_17}{Number of household members aged 12-17}
#'   \item{n_menores}{Number of household members under 18, capped at 4}
#'   \item{n_mayores}{Number of adult household members (18+), capped at 6}
#'   \item{n_tiempo}{Number of household members who reported time use}
#'   \item{n_trabajadores}{Number of employed workers in household}
#'   \item{n_profesionales}{Number of household members with completed university or higher}
#'   \item{n_tercera_edad}{Number of household members aged 60+}
#'   \item{hay_tercera_edad}{1 if household contains at least one elderly member aged 60+ who is not the respondent}
#'   \item{n_personas}{Total household size}
#'   \item{edad_promedio}{Mean age of household members}
#'   \item{tiene_hijos}{1 if respondent has children living in the household}
#'   \item{en_pareja}{1 if respondent is in a couple relationship (married or cohabiting)}
#'   \item{vive_pareja}{1 if respondent lives with their partner}
#'
#'   \item{servicio_domestico}{1 if household receives paid domestic service}
#'   \item{ayuda_cercanos}{1 if household receives unpaid help from relatives or neighbours}
#'   \item{fuentes_externas}{1 if household receives any external care support (servicio_domestico or ayuda_cercanos)}
#'
#'   \item{sexo}{Female indicator (1 = female, 0 = male; recoded from raw variable)}
#'   \item{edad_anios}{Age in years}
#'   \item{tramo_edad}{Age bracket: "12-24", "25-44", "45-65", or "66+"}
#'   \item{NSE}{Socioeconomic level from original survey classification: 1 = Bajo, 2 = Medio, 3 = Alto}
#'   \item{nivel_escolaridad}{Highest completed education level: "ninguna", "primaria", "secundaria", "técnica", or "universitaria"}
#'   \item{estudia}{1 if currently enrolled in education}
#'   \item{trabaja}{1 if currently employed (from o1)}
#'   \item{horas_trabajo_contratadas}{Horas semanales contratadas segun contrato de trabajo (from o22c)}
#'   \item{horas_trabajo_habituales}{Horas habituales de trabajo en la ocupacion principal (from o22a);
#'     compare with \code{t_paid_work} (diary-measured) and \code{horas_trabajo_contratadas} (contracted)}
#'   \item{dias_trabajo_semana}{Dias trabajados por semana en la ocupacion principal (from o22b);
#'     used internally to scale weekday diary hours during weekend imputation}
#'   \item{quintil}{Household income quintile (1-5), computed from cumulative survey weights}
#'   \item{macrozona}{Geographic macro-zone: "norte", "metropolitana", "centro", or "sur"}
#'   \item{region_ord}{Ordinal region code}
#'   \item{glosa_region}{Region name label}
#'   \item{prop_ing_hogar}{Respondent's share of total personal income in the household; used to allocate household-level expenditures to individuals}
#'
#'   \item{cae}{Economic activity category: "Ocupada(o)", "Desocupada(o)", "Inactiva(o)", "Menor de 15 años", or "Sin clasificar"}
#'   \item{teletrabaja}{1 if the respondent teleworks (from to2_v_ds or to2_v_fds)}
#'   \item{jornada_laboral}{Work schedule type (from jor_to)}
#'   \item{ocup_form}{Occupation formality: 1 = Formal, 2 = Informal}
#'   \item{cise}{Employment status per CISE-2018 classification}
#'   \item{ciuo_agrupada}{Occupation group per CIUO-08 grouped classification}
#'
#'   \item{bs1}{Life satisfaction: "Que tan satisfecho se siente con su vida en general?"
#'     Ordinal 1-5: 1 = Totalmente insatisfecho(a), 5 = Totalmente satisfecho(a).
#'     Missing: 96.}
#'   \item{bs2}{Satisfaction with distribution of domestic tasks in the household:
#'     "Que tan satisfecho se siente con la forma en que se reparten las tareas domesticas en su hogar?"
#'     Ordinal 1-5 (same scale as bs1). 85 = No aplica. Missing: 96.}
#'   \item{bs3}{Satisfaction with distribution of care tasks in the household:
#'     "Que tan satisfecho se siente con la forma en que se reparten las tareas de cuidado en su hogar?"
#'     Ordinal 1-5 (same scale as bs1). 85 = No aplica. Missing: 96.}
#'   \item{bs4}{Perceived time adequacy for paid work:
#'     "Habitualmente, considera que en el trabajo en su ocupacion..."
#'     Ordinal 1-3: 1 = Me falta tiempo, 2 = El tiempo me alcanza bien, 3 = Me sobra tiempo.
#'     85 = No aplica. Missing: 96.}
#'   \item{bs5}{Perceived time adequacy for study/attending classes.
#'     Same scale as bs4. 85 = No aplica. Missing: 96.}
#'   \item{bs6}{Perceived time adequacy for domestic work (cooking, cleaning, shopping).
#'     Same scale as bs4. 85 = No aplica. Missing: 96.}
#'   \item{bs7}{Perceived time adequacy for caring for household members.
#'     Same scale as bs4. 85 = No aplica. Missing: 96.}
#'   \item{bs8}{Perceived time adequacy for personal care (hygiene, eating, health).
#'     Ordinal 1-3 (same scale as bs4, no No-aplica category). Missing: 96.}
#'   \item{bs9}{Perceived time adequacy for leisure and social life (conversation, TV, friends).
#'     Ordinal 1-3 (same scale as bs4, no No-aplica category). Missing: 96.}
#'   \item{bs10}{Self-assessed share of domestic work relative to what is fair:
#'     "Con respecto al trabajo domestico que realiza en su hogar habitualmente, usted considera que hace..."
#'     Ordinal 1-3: 1 = Menos de lo que corresponde, 2 = Lo que corresponde, 3 = Mas de lo que corresponde.
#'     85 = No aplica. Missing: 96.}
#'   \item{bs11}{Self-assessed share of care work relative to what is fair (same scale as bs10).
#'     85 = No aplica. Missing: 96.}
#'   \item{bs12}{Time scarcity due to caring for a person with dependency (PSDF):
#'     "Que tan frecuentemente piensa que debido al tiempo que dedica a esta persona no tiene suficiente tiempo para usted?"
#'     Ordinal 1-5: 1 = Nunca, 5 = Siempre. Missing: 88, 96, 99. N valid: 1378.}
#'   \item{bs13}{Caregiver burden from reconciling PSDF care with other family/work responsibilities.
#'     Same scale as bs12. Missing: 88, 96, 99. N valid: 1378.}
#'   \item{bs14}{Negative effect of PSDF care on relationships with others.
#'     Same scale as bs12. Missing: 88, 96, 99. N valid: 1378.}
#'   \item{bs15}{Perceived health deterioration due to caring for a PSDF person.
#'     Same scale as bs12. Missing: 88, 96, 99. N valid: 1378.}
#'   \item{bs16}{General sense of overload from caring for a PSDF person.
#'     Same scale as bs12. Missing: 88, 96, 99. N valid: 1378.}
#'   \item{bs17}{Time scarcity due to caring for a child (NNA):
#'     "Piensa que debido al tiempo que le dedica al cuidado de [NOMBRE NNA] no tiene suficiente tiempo para usted?"
#'     Same scale as bs12. Missing: 88, 96, 99. N valid: 7071.}
#'   \item{bs18}{Caregiver burden from reconciling NNA care with work responsibilities.
#'     Same scale as bs12. Missing: 88, 96, 99. N valid: 5162.}
#'   \item{bs19}{Negative effect of NNA care on physical or mental health.
#'     Same scale as bs12. Missing: 88, 96, 99. N valid: 7071.}
#'
#'   \item{ing_ocuppal}{Income from main occupation (weekly, thousands CLP)}
#'   \item{ing_trab}{Total labor income (weekly, thousands CLP)}
#'   \item{ing_jub_aps}{Pension and AFP income (weekly, thousands CLP)}
#'   \item{ing_g}{Imputed group income component from household total (weekly, thousands CLP)}
#'   \item{ing_t_hogar}{Total household income (weekly, thousands CLP)}
#'   \item{ing_t_pc}{Per capita household income (weekly, thousands CLP)}
#'   \item{ing_gpp}{Individual share of group income (weekly, thousands CLP)}
#'   \item{ing_personal}{Personal income: ing_trab + ing_jub_aps + ing_gpp (weekly, thousands CLP)}
#'   \item{ingreso_hogar}{Total household disposable income (weekly, thousands CLP)}
#'   \item{income_person_week}{Household income divided by number of members (weekly, thousands CLP)}
#'
#'   \item{t_paid_work}{Trabajo remunerado; equivale a t_to, horas semanales}
#'   \item{t_domestic_work}{Trabajo domestico no remunerado; suma de t_tdnr_psc + t_tdnr_lv +
#'     t_tdnr_lrc + t_tdnr_mrm + t_tdnr_admnhog + t_tdnr_comphog + t_tdnr_cmp, horas semanales}
#'   \item{t_care_work}{Trabajo de cuidado no remunerado; suma de t_tcnr_ce + t_tcnr_re +
#'     t_tcnr_oac, horas semanales}
#'   \item{t_unpaid_voluntary}{Trabajo voluntario y ayuda a otros hogares; suma de t_tvaoh_tv +
#'     t_tvaoh_oh, horas semanales}
#'   \item{t_education}{Educacion y formacion; equivale a t_ed, horas semanales}
#'   \item{t_leisure}{Ocio; suma de t_vsyo_csar + t_vsyo_aa + t_mcm_leer + t_mcm_video +
#'     t_mcm_audio + t_mcm_computador, horas semanales}
#'   \item{t_personal_care}{Cuidados personales fisiologicos (excluye sueno y comidas); equivale
#'     a t_cpaf_cp, horas semanales}
#'   \item{t_meals}{Comer y beber; equivale a t_cpag_comer, horas semanales}
#'   \item{t_sleep}{Dormir; equivale a t_cpag_dormir, ajustado para que la suma sea 168 horas,
#'     horas semanales}
#'   \item{t_commute1}{Traslados asociados a trabajo remunerado, educacion y salud; equivale a
#'     t_tt1, horas semanales}
#'   \item{t_commute2}{Traslados asociados a tramites del hogar y cuidados; equivale a t_tt2,
#'     horas semanales}
#'
#'   \item{Tw}{Paid work time (equivalent to t_to / t_paid_work)}
#'   \item{Tf_social}{Social life and recreation time (equivalent to t_vsyo_csar)}
#'   \item{Tf_hobbies}{Hobbies and arts time (equivalent to t_vsyo_aa)}
#'   \item{Tf_read}{Reading time (equivalent to t_mcm_leer)}
#'   \item{Tf_listen}{Audio consumption time (equivalent to t_mcm_audio)}
#'   \item{Tf_watch}{TV and video consumption time (equivalent to t_mcm_video)}
#'   \item{Tf_computer}{Recreational computer/internet use time (equivalent to t_mcm_computador)}
#'   \item{Tc_meals}{Time spent eating and drinking (equivalent to t_cpag_comer)}
#'   \item{Tc_sleep}{Time spent sleeping, adjusted to balance 168 hours (equivalent to t_cpag_dormir)}
#'   \item{Tc_other}{All other time use (domestic work, care, other personal care, volunteering, commuting, education)}
#'
#'   \item{t_total}{Total weekly hours across all activities (should equal 168)}
#'   \item{w}{Hourly wage rate: ing_trab / t_paid_work (thousands CLP per hour)}
#'
#'   \item{Ef_food}{Imputed food expenditure (weekly thousands CLP)}
#'   \item{Ef_recreation}{Imputed recreation and culture expenditure (weekly thousands CLP)}
#'   \item{Ef_restaurants}{Imputed restaurants and food away from home expenditure (weekly thousands CLP)}
#'   \item{Ef_communications}{Imputed communications expenditure (weekly thousands CLP)}
#'   \item{Ef_clothing}{Imputed clothing expenditure (weekly thousands CLP)}
#'   \item{Ec}{Imputed committed and remaining expenditures (weekly thousands CLP, includes accounts, health, transport, education, household goods, and savings)}
#' }
#'
#' @details
#' The dataset is produced by running \code{data_processing/data_processing.R}. Time
#' variables are normalized to a 168-hour week using a two-step procedure: paid work
#' (\code{t_paid_work}) and sleep (\code{t_sleep}) are treated as reliable anchors;
#' all other activities are scaled proportionally. Weekend time use is imputed via a
#' twin-matching matrix when the respondent's diary day was a weekday. Outliers are
#' removed using Vallejo's method by groups of quintile, employment status, age
#' bracket, and sex. Expenditures are imputed from EPF IX using a fractional MNL
#' model for budget shares and a linear regression for the savings rate, then
#' allocated to individuals by \code{prop_ing_hogar}. The script writes
#' \code{data/enut-ii.dta} plus \code{data/enut-ii.csv}. English copies are written
#' as \code{data/enut-ii-ENG.dta} and \code{data/enut-ii-ENG.csv}.
#'
#' The time-use categories are aggregations of the 25 detailed categories in
#' \code{enut_ii_raw}. See \code{agregar_actividades()} in
#' \code{data_processing/processing_functions.R} for the exact aggregation mapping.
#'
#' @source <https://www.ine.gob.cl/enut>
#' @source <https://www.ine.gob.cl/estadisticas/sociales/ingresos-y-gastos/encuesta-de-presupuestos-familiares>
#'
#' @docType data
#' @keywords datasets
#' @name enut_ii
#' @usage data(enut_ii)
#' @format A data frame with approximately 4,000-5,000 rows
NULL

#' enut-ii-raw
#'
#' Processed dataset from the second National Time-Use Survey (ENUT II), applied by the
#' Instituto Nacional de Estadísticas de Chile. Contains 25 detailed time-use activity
#' categories and household expenditures imputed from the IX Encuesta de Presupuestos
#' Familiares (EPF). Used as the primary input for structural time-use models via
#' \code{get_data()}.
#'
#' All income and expenditure variables are expressed in weekly thousands of Chilean pesos,
#' deflated by IPC adjustment (factor 0.362). Time variables are expressed in weekly hours,
#' normalized to sum to 168.
#'
#' \describe{
#'   \item{id_persona}{Individual identifier}
#'   \item{id_hog}{Household identifier}
#'
#'   \item{es_trabajador}{1 if employed, CAE = "Ocupada(o)", and positive paid work time}
#'   \item{es_familia}{1 if es_trabajador, lives with partner, and has children in household}
#'
#'   \item{dia_semana}{Weekday of diary (1 = Monday to 5 = Friday)}
#'   \item{dia_fin_semana}{Weekend day of diary (6 = Saturday, 7 = Sunday)}
#'
#'   \item{parentesco}{Relationship to household head (from pco)}
#'   \item{n_menores_0_5}{Number of household members under age 6 (ages 0-5)}
#'   \item{n_menores_6_11}{Number of household members aged 6-11}
#'   \item{n_menores_0_14}{Number of household members under age 15 (ages 0-14)}
#'   \item{n_menores_12_17}{Number of household members aged 12-17}
#'   \item{n_menores}{Number of household members under 18, capped at 4}
#'   \item{n_mayores}{Number of adult household members (18+), capped at 6}
#'   \item{n_tiempo}{Number of household members who reported time use}
#'   \item{n_trabajadores}{Number of employed workers in household}
#'   \item{n_profesionales}{Number of household members with completed university or higher}
#'   \item{n_tercera_edad}{Number of household members aged 60+}
#'   \item{hay_tercera_edad}{1 if household contains at least one elderly member aged 60+ who is not the respondent}
#'   \item{n_personas}{Total household size}
#'   \item{edad_promedio}{Mean age of household members}
#'   \item{tiene_hijos}{1 if respondent has children living in the household}
#'   \item{en_pareja}{1 if respondent is in a couple relationship (married or cohabiting)}
#'   \item{vive_pareja}{1 if respondent lives with their partner}
#'
#'   \item{servicio_domestico}{1 if household receives paid domestic service}
#'   \item{ayuda_cercanos}{1 if household receives unpaid help from relatives or neighbours}
#'   \item{fuentes_externas}{1 if household receives any external care support (servicio_domestico or ayuda_cercanos)}
#'
#'   \item{sexo}{Female indicator (1 = female, 0 = male; recoded from raw variable)}
#'   \item{edad_anios}{Age in years}
#'   \item{tramo_edad}{Age bracket: "12-24", "25-44", "45-65", or "66+"}
#'   \item{NSE}{Socioeconomic level from original survey classification: 1 = Bajo, 2 = Medio, 3 = Alto}
#'   \item{nivel_escolaridad}{Highest completed education level: "ninguna", "primaria", "secundaria", "técnica", or "universitaria"}
#'   \item{estudia}{1 if currently enrolled in education}
#'   \item{trabaja}{1 if currently employed (from o1)}
#'   \item{horas_trabajo_contratadas}{Horas semanales contratadas segun contrato de trabajo (from o22c)}
#'   \item{horas_trabajo_habituales}{Horas habituales de trabajo en la ocupacion principal (from o22a);
#'     compare with \code{t_to} (diary-measured) and \code{horas_trabajo_contratadas} (contracted)}
#'   \item{dias_trabajo_semana}{Dias trabajados por semana en la ocupacion principal (from o22b);
#'     used internally to scale weekday diary hours during weekend imputation}
#'   \item{quintil}{Household income quintile (1-5), computed from cumulative survey weights}
#'   \item{macrozona}{Geographic macro-zone: "norte", "metropolitana", "centro", or "sur"}
#'   \item{region_ord}{Ordinal region code}
#'   \item{glosa_region}{Region name label}
#'   \item{prop_ing_hogar}{Respondent's share of total personal income in the household; used to allocate household-level expenditures to individuals}
#'
#'   \item{cae}{Economic activity category: "Ocupada(o)", "Desocupada(o)", "Inactiva(o)", "Menor de 15 años", or "Sin clasificar"}
#'   \item{teletrabaja}{1 if the respondent teleworks (from to2_v_ds or to2_v_fds)}
#'   \item{jornada_laboral}{Work schedule type (from jor_to)}
#'   \item{ocup_form}{Occupation formality: 1 = Formal, 2 = Informal}
#'   \item{cise}{Employment status per CISE-2018 classification}
#'   \item{ciuo_agrupada}{Occupation group per CIUO-08 grouped classification}
#'
#'   \item{bs1}{Life satisfaction: "Que tan satisfecho se siente con su vida en general?"
#'     Ordinal 1-5: 1 = Totalmente insatisfecho(a), 5 = Totalmente satisfecho(a).
#'     Missing: 96.}
#'   \item{bs2}{Satisfaction with distribution of domestic tasks in the household:
#'     "Que tan satisfecho se siente con la forma en que se reparten las tareas domesticas en su hogar?"
#'     Ordinal 1-5 (same scale as bs1). 85 = No aplica. Missing: 96.}
#'   \item{bs3}{Satisfaction with distribution of care tasks in the household:
#'     "Que tan satisfecho se siente con la forma en que se reparten las tareas de cuidado en su hogar?"
#'     Ordinal 1-5 (same scale as bs1). 85 = No aplica. Missing: 96.}
#'   \item{bs4}{Perceived time adequacy for paid work:
#'     "Habitualmente, considera que en el trabajo en su ocupacion..."
#'     Ordinal 1-3: 1 = Me falta tiempo, 2 = El tiempo me alcanza bien, 3 = Me sobra tiempo.
#'     85 = No aplica. Missing: 96.}
#'   \item{bs5}{Perceived time adequacy for study/attending classes.
#'     Same scale as bs4. 85 = No aplica. Missing: 96.}
#'   \item{bs6}{Perceived time adequacy for domestic work (cooking, cleaning, shopping).
#'     Same scale as bs4. 85 = No aplica. Missing: 96.}
#'   \item{bs7}{Perceived time adequacy for caring for household members.
#'     Same scale as bs4. 85 = No aplica. Missing: 96.}
#'   \item{bs8}{Perceived time adequacy for personal care (hygiene, eating, health).
#'     Ordinal 1-3 (same scale as bs4, no No-aplica category). Missing: 96.}
#'   \item{bs9}{Perceived time adequacy for leisure and social life (conversation, TV, friends).
#'     Ordinal 1-3 (same scale as bs4, no No-aplica category). Missing: 96.}
#'   \item{bs10}{Self-assessed share of domestic work relative to what is fair:
#'     "Con respecto al trabajo domestico que realiza en su hogar habitualmente, usted considera que hace..."
#'     Ordinal 1-3: 1 = Menos de lo que corresponde, 2 = Lo que corresponde, 3 = Mas de lo que corresponde.
#'     85 = No aplica. Missing: 96.}
#'   \item{bs11}{Self-assessed share of care work relative to what is fair (same scale as bs10).
#'     85 = No aplica. Missing: 96.}
#'   \item{bs12}{Time scarcity due to caring for a person with dependency (PSDF):
#'     "Que tan frecuentemente piensa que debido al tiempo que dedica a esta persona no tiene suficiente tiempo para usted?"
#'     Ordinal 1-5: 1 = Nunca, 5 = Siempre. Missing: 88, 96, 99. N valid: 1378.}
#'   \item{bs13}{Caregiver burden from reconciling PSDF care with other family/work responsibilities.
#'     Same scale as bs12. Missing: 88, 96, 99. N valid: 1378.}
#'   \item{bs14}{Negative effect of PSDF care on relationships with others.
#'     Same scale as bs12. Missing: 88, 96, 99. N valid: 1378.}
#'   \item{bs15}{Perceived health deterioration due to caring for a PSDF person.
#'     Same scale as bs12. Missing: 88, 96, 99. N valid: 1378.}
#'   \item{bs16}{General sense of overload from caring for a PSDF person.
#'     Same scale as bs12. Missing: 88, 96, 99. N valid: 1378.}
#'   \item{bs17}{Time scarcity due to caring for a child (NNA):
#'     "Piensa que debido al tiempo que le dedica al cuidado de [NOMBRE NNA] no tiene suficiente tiempo para usted?"
#'     Same scale as bs12. Missing: 88, 96, 99. N valid: 7071.}
#'   \item{bs18}{Caregiver burden from reconciling NNA care with work responsibilities.
#'     Same scale as bs12. Missing: 88, 96, 99. N valid: 5162.}
#'   \item{bs19}{Negative effect of NNA care on physical or mental health.
#'     Same scale as bs12. Missing: 88, 96, 99. N valid: 7071.}
#'
#'   \item{ing_ocuppal}{Income from main occupation (weekly, thousands CLP)}
#'   \item{ing_trab}{Total labor income (weekly, thousands CLP)}
#'   \item{ing_jub_aps}{Pension and AFP income (weekly, thousands CLP)}
#'   \item{ing_g}{Imputed group income component from household total (weekly, thousands CLP)}
#'   \item{ing_t_hogar}{Total household income (weekly, thousands CLP)}
#'   \item{ing_t_pc}{Per capita household income (weekly, thousands CLP)}
#'   \item{ing_gpp}{Individual share of group income (weekly, thousands CLP)}
#'   \item{ing_personal}{Personal income: ing_trab + ing_jub_aps + ing_gpp (weekly, thousands CLP)}
#'   \item{ingreso_hogar}{Total household disposable income (weekly, thousands CLP)}
#'   \item{income_person_week}{Household income divided by number of members (weekly, thousands CLP)}
#'
#'   \item{t_to}{Trabajo remunerado (paid work), horas semanales}
#'   \item{t_tcnr_ce}{Trabajo de cuidado no remunerado - cuidados esenciales: cuidado fisico personal de
#'     miembros del hogar (aseo, alimentacion, vestido), horas semanales}
#'   \item{t_tcnr_re}{Trabajo de cuidado no remunerado - cuidados relativos a la ensenanza: acompanamiento
#'     al establecimiento educacional, apoyo en tareas y lectura, horas semanales}
#'   \item{t_tcnr_oac}{Trabajo de cuidado no remunerado - otros cuidados: acompanamiento medico, traslado
#'     al trabajo u otras actividades, otros cuidados no clasificados, horas semanales}
#'   \item{t_tdnr_psc}{Trabajo domestico no remunerado - preparacion y servicio de comida, horas semanales}
#'   \item{t_tdnr_lv}{Trabajo domestico no remunerado - limpieza de vivienda, horas semanales}
#'   \item{t_tdnr_lrc}{Trabajo domestico no remunerado - limpieza y reparacion de ropa, horas semanales}
#'   \item{t_tdnr_mrm}{Trabajo domestico no remunerado - mantenimiento y reparaciones menores del hogar,
#'     horas semanales}
#'   \item{t_tdnr_admnhog}{Trabajo domestico no remunerado - administracion del hogar, horas semanales}
#'   \item{t_tdnr_comphog}{Trabajo domestico no remunerado - compras del hogar, horas semanales}
#'   \item{t_tdnr_cmp}{Trabajo domestico no remunerado - cuidados de mascotas y plantas, horas semanales}
#'   \item{t_tvaoh_tv}{Trabajo voluntario y ayuda a otros hogares - voluntariado en la comunidad,
#'     horas semanales}
#'   \item{t_tvaoh_oh}{Trabajo voluntario y ayuda a otros hogares - ayuda directa a otros hogares,
#'     horas semanales}
#'   \item{t_cpaf_cp}{Cuidados personales - actividades fisiologicas (higiene, visitas medicas, ejercicio),
#'     excluye sueno y comidas, horas semanales}
#'   \item{t_cpag_comer}{Cuidados personales - comer y beber, horas semanales}
#'   \item{t_cpag_dormir}{Cuidados personales - dormir; ajustado para que la suma de actividades sea
#'     168 horas, horas semanales}
#'   \item{t_ed}{Educacion y formacion (asistencia a clases, estudio, capacitacion), horas semanales}
#'   \item{t_vsyo_csar}{Vida social y ocio - convivencia social y actividades recreativas (reuniones,
#'     fiestas, deportes como espectador), horas semanales}
#'   \item{t_vsyo_aa}{Vida social y ocio - artes y aficiones (artes plasticas, artesania, juegos),
#'     horas semanales}
#'   \item{t_mcm_leer}{Medios de comunicacion y masivos - lectura (libros, prensa, digital), horas semanales}
#'   \item{t_mcm_audio}{Medios de comunicacion y masivos - consumo de audio (radio, musica, podcasts),
#'     horas semanales}
#'   \item{t_mcm_video}{Medios de comunicacion y masivos - consumo de video y television, horas semanales}
#'   \item{t_mcm_computador}{Medios de comunicacion y masivos - uso recreativo de computador e internet,
#'     horas semanales}
#'   \item{t_tt1}{Traslados 1 - traslados asociados a trabajo remunerado, educacion y salud (equivalente
#'     ENUT 2015), horas semanales}
#'   \item{t_tt2}{Traslados 2 - traslados asociados a tramites del hogar y cuidados (adicionales ENUT 2023),
#'     horas semanales}
#'
#'   \item{t_total}{Total weekly hours across all activities (should equal 168)}
#'   \item{w}{Hourly wage rate: ing_trab / t_to (thousands CLP per hour)}
#'
#'   \item{alimentos}{Imputed food expenditure (individual share, weekly thousands CLP)}
#'   \item{vestimenta}{Imputed clothing expenditure (individual share, weekly thousands CLP)}
#'   \item{cuentas}{Imputed utility and fixed household expenditure (individual share, weekly thousands CLP)}
#'   \item{hogar}{Imputed household goods and services expenditure (individual share, weekly thousands CLP)}
#'   \item{salud}{Imputed health expenditure (individual share, weekly thousands CLP)}
#'   \item{transporte}{Imputed transportation expenditure (individual share, weekly thousands CLP)}
#'   \item{comunicaciones}{Imputed communications expenditure (individual share, weekly thousands CLP)}
#'   \item{recreacion}{Imputed recreation and culture expenditure (individual share, weekly thousands CLP)}
#'   \item{educacion}{Imputed education expenditure (individual share, weekly thousands CLP)}
#'   \item{restaurantes}{Imputed restaurants and food away from home expenditure (individual share, weekly thousands CLP)}
#'   \item{savings}{Imputed household savings, allocated by income share (weekly thousands CLP)}
#'   \item{total_expenses}{Total imputed expenditure: ing_personal - savings (weekly thousands CLP)}
#' }
#'
#' @details
#' The dataset is produced by running \code{data_processing/data_processing.R}. Time
#' variables are normalized to a 168-hour week using a two-step procedure: paid work
#' (\code{t_to}) and sleep (\code{t_cpag_dormir}) are treated as reliable anchors;
#' all other activities are scaled proportionally. Weekend time use is imputed via a
#' twin-matching matrix when the respondent's diary day was a weekday. Outliers are
#' removed using Vallejo's method by groups of quintile, employment status, age
#' bracket, and sex. Expenditures are imputed from EPF IX using a fractional MNL
#' model for budget shares and a linear regression for the savings rate, then
#' allocated to individuals by \code{prop_ing_hogar}. The script writes
#' \code{data/enut-ii-raw.dta} plus \code{data/enut-ii-raw.csv}. English copies are
#' written as \code{data/enut-ii-raw-ENG.dta} and
#' \code{data/enut-ii-raw-ENG.csv}.
#'
#' @source <https://www.ine.gob.cl/enut>
#' @source <https://www.ine.gob.cl/estadisticas/sociales/ingresos-y-gastos/encuesta-de-presupuestos-familiares>
#'
#' @docType data
#' @keywords datasets
#' @name enut_ii_raw
#' @usage data(enut_ii_raw)
#' @format A data frame with approximately 4,000-5,000 rows and 112 variables
NULL
