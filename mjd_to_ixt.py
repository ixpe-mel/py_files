#Takes in a time in MJD and spits out the IXT conversion, ie modified julian date to number of seconds after ixpe epoch

def mjd_to_ixt(mjd_mjdreff,mjdrefi):
    return 86400*(mjd-mjdrefi+mjdreff)
