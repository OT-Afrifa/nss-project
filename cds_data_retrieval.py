# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 11:59:52 2021

@author: MET_SAT_RECEIVER
"""


import cdsapi

#'ec_earth3_aerchem','ec_earth3_cc',
#'''
models = ['access_cm2','bcc_csm2_mr','cesm2_fv2','cesm2_waccm_fv2','cmcc_cm2_hr4','cmcc_esm2','cnrm_cm6_1_hr',
'fgoals_f3_l','hadgem3_gc31_ll','iitm_esm','inm_cm5_0','ipsl_cm6a_lr','kiost_esm','miroc6',
'miroc_es2l','mpi_esm1_2_hr','mri_esm2_0','norcpm1','noresm2_mm','taiesm1',
'access_esm1_5','awi_esm_1_1_lr','bcc_esm1','canesm5','cesm2','cesm2_waccm','cmcc_cm2_sr5',
'cnrm_cm6_1','cnrm_esm2_1','ec_earth3_veg_lr','fgoals_g3','gfdl_esm4','hadgem3_gc31_mm','inm_cm4_8',
'ipsl_cm5a2_inca','kace_1_0_g','mpi_esm1_2_lr','nesm3','sam0_unicon','ukesm1_0_ll']
#'''


#models = ['cesm2_waccm_fv2','ec_earth3_aerchem','ec_earth3_cc','mpi_esm1_2_hr','cesm2_fv2']
#'historical', 'ssp5_3_4', 
experiments = ['historical', 'ssp5_3_4', 'ssp3_7_0', 'ssp1_1_9', 'ssp1_2_6', 'ssp4_3_4', 'ssp2_4_5', 'ssp4_6_0', 'ssp5_8_5']
#'ssp4_6_0', 'ssp5_3_4', 'ssp5_8_5', 
#experiments = ['historical', 'ssp5_8_5']

#exp = 'ssp3_7_0'


var = ['precipitation', 'near_surface_air_temperature']

for exp in experiments[0:1]:
    for model in sorted(models)[::-1]:
        c = cdsapi.Client()
    
        try:
            c.retrieve(
                'projections-cmip6',
                {
                    'format': 'zip',
                    'temporal_resolution': 'daily',
                    'experiment': exp,
                    'level': 'single_levels',
                    'variable': var[1],
                    'model': [model],
                    'date': '1850-01-01/2100-12-31',
                    'area': [
                        15, -5, 0, 5,
                    ],
                },
                exp+'/temperature_cmip6_download_'+model+'.zip')
        except:
            pass