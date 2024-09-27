#cleaning fits files by cutting unneded energy channels and QUAL!=1
from astropy.io import fits
def load_and_clean(file, Pmin, Pmax):
    with fits.open(str(file)) as hdu:
                data=hdu[1].data
                header=hdu[1].header

    data=data[(Pmin<=data['PI']) & (data['PI']<=Pmax)]
    data=data[data['QUAL']==1]
    return data,header