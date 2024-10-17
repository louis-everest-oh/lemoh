# Louis Everest's Miscellaneous Functions for Ontario Health Occupational Research

Contains functions for occupational health research.

### Developlemt Notes
 This is a beta version of the lemoh package 

## Getting Started
Install from Github 

    library(devtools) 
    install_github("louis-everest-oh/lemoh")
    

## Functions

### Get Table Rapid

A more efficient version of LTASR::get_table()

    get_table_rapid(persondf, strata = c(), break_yr=5)

### Expected Risk of an Event of Interest

Calculates a person table based on time-varying annual exposure data 

    get_table_annual_exposure(person_data, exposure_data,
                          outcomes, cat_vars, cont_vars, exp_vars,
                          bin = 500, mid_year = 0.4516)


## Authors

  - **Louis Everest** - [aut,cre]



## License

This project is licensed under MIT 

## Acknowledgments

