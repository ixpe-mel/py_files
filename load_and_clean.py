#cleaning fits files by cutting unneded energy channels and QUAL!=1
from astropy.io import fits
def load_and_clean(file, Pmin, Pmax):
    with fits.open(str(file)) as hdu:
                data=hdu[1].data
                print('num of events',len(data['PI']))
                GTI=hdu[2].data
                GTI_header=hdu[2].header
                header=hdu[1].header
                prihdr = hdu[0].header
                prihdu = fits.PrimaryHDU(header=prihdr)

    data=data[(Pmin<=data['PI']) & (data['PI']<=Pmax)]
    data=data[data['QUAL']==1]
    return data,header,GTI,GTI_header,prihdu