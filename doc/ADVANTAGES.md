# "ClimDown downscaling R package"

## Features (i.e. why is this better than running our scripts?)

* Easy configuration

    * sensible defaults (see config.R)
    * configurable by user, by run, by environment variable
    * uses R's built in `options` method

* One downscaling method <----> one script

    * all scripts just take 3-4 command line parameters

~~~bash
    $ Rscript downscaling.R [input_gcm] [input_obs] [output_grid] [variable_name]
~~~

* Decomposable

    * all shell scripts are "dumb" wrappers that process command line arguments and delegate
    * most computation is wrapped by a netcdf wrapper that handles chunking and delegates
    * I/O wrapper can *programmatically* determine optimum read size (for some value of optimum)
    * testing is *MUUUUUUUCH* easier (just test with a smaller grid for fast testing)

* It's an actual R package

~~~R
    $ R
    > library(ClimDown)
    > bcca.netcdf.wrapper("my_input.nc", "my_output.nc", "tasmax")
~~~

* Variables are parametrized

    * no redundancy
    * follows DRY principles

* Uses udunits2 and PCICt

    * handles units and time conversion in a systematic way