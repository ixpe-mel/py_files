
#rxj0440 pd and pa polynomials
import numpy as np
def rxj_PD_poly_mymachine(x): #x is phase
        #y=13.599084715966825 - 77.14279098237614*(-1.0+0.3183098861837907*x) - 58.68027525042033*(-1.0+0.3183098861837907*x)**2+2459.257155713502*(-1.0+0.3183098861837907*x)**3 - 585.7971512178035 * (-1.0+0.3183098861837907*x)**4 - 42389.345465598555*(-1.0+0.3183098861837907*x)**5+31377.701921467713*(-1.0+0.3183098861837907*x)**6+389184.2681757424*(-1.0+0.3183098861837907*x)**7 - 295253.60901547165*(-1.0+0.3183098861837907*x)**8 - 2047779.3589603745*(-1.0+0.3183098861837907*x)**9+1270324.491309807*(-1.0+0.3183098861837907*x)**10+6643443.7902267715*(-1.0+0.3183098861837907*x)**11 - 2966775.4033494988*(-1.0+0.3183098861837907*x)**12 - 13851400.5918384*(-1.0+0.3183098861837907*x)**13 + 3872312.125065608*(-1.0+0.3183098861837907*x)**14 + 18812624.00901786*(-1.0+0.3183098861837907*x)**15 - 2556606.7984732436*(-1.0+0.3183098861837907*x)**16 - 16428501.620406935*(-1.0+0.3183098861837907*x)**17 + 390641.771729522*(-1.0+0.3183098861837907*x)**18+8797093.817109447*(-1.0+0.3183098861837907*x)**19 + 432821.28297364846*(-1.0+0.3183098861837907*x)**20 - 2584215.490575322*(-1.0+0.3183098861837907*x)**21 - 178202.6787844045*(-1.0+0.3183098861837907*x)**22+309558.65888897196*(-1.0+0.3183098861837907*x)**23
        #y=y/100

                                                         
        y_new=2.723e-26*x**47  - 4.368e-25*x**46  + 4.079e-25*x**45  + 1.32e-23*x**44 

        + 4.738e-23*x**43  - 1.649e-22*x**42  - 2.777e-21*x**41  - 1.445e-20*x**40 

        - 9.707e-21*x**39  + 4.632e-19*x**38  + 4.215e-18*x**37  + 1.783e-17*x**36 

        - 1.054e-17*x**35  - 7.876e-16*x**34  - 6.416e-15*x**33  - 2.447e-14*x**32 

        + 4.341e-14*x**31  + 1.351e-12*x**30  + 9.761e-12*x**29  + 2.764e-11*x**28 

        - 1.66e-10*x**27  - 2.466e-09*x**26  - 1.274e-08*x**25  + 2.937e-09*x**24 

        + 5.43e-07*x**23  + 3.717e-06*x**22  + 2.529e-06*x**21  - 0.0001375*x**20 

        - 0.000892*x**19  + 0.001124*x**18  + 0.04198*x**17  + 0.1024*x**16  - 1.489*x**15 

        - 6.236*x**14  + 63.46*x**13  + 134.4*x**12  - 3252*x**11  + 1.771e+04*x**10 

        - 5.507e+04*x**9 + 1.132e+05*x**8 - 1.613e+05*x**7 + 1.61e+05*x**6 - 1.11e+05*x**5

        + 5.075e+04*x**4 - 1.415e+04*x**3 + 2033*x**2 - 102.7*x + 9.976

        

        #y_new
        #print(y_new)
        return y_new

def rxj_PD_poly(x): #x is phase
        #defined over a phase of 0-2pi in percentage
        y=13.599084715966825 - 77.14279098237614*(-1.0+0.3183098861837907*x) - 58.68027525042033*(-1.0+0.3183098861837907*x)**2+2459.257155713502*(-1.0+0.3183098861837907*x)**3 - 585.7971512178035 * (-1.0+0.3183098861837907*x)**4 - 42389.345465598555*(-1.0+0.3183098861837907*x)**5+31377.701921467713*(-1.0+0.3183098861837907*x)**6+389184.2681757424*(-1.0+0.3183098861837907*x)**7 - 295253.60901547165*(-1.0+0.3183098861837907*x)**8 - 2047779.3589603745*(-1.0+0.3183098861837907*x)**9+1270324.491309807*(-1.0+0.3183098861837907*x)**10+6643443.7902267715*(-1.0+0.3183098861837907*x)**11 - 2966775.4033494988*(-1.0+0.3183098861837907*x)**12 - 13851400.5918384*(-1.0+0.3183098861837907*x)**13 + 3872312.125065608*(-1.0+0.3183098861837907*x)**14 + 18812624.00901786*(-1.0+0.3183098861837907*x)**15 - 2556606.7984732436*(-1.0+0.3183098861837907*x)**16 - 16428501.620406935*(-1.0+0.3183098861837907*x)**17 + 390641.771729522*(-1.0+0.3183098861837907*x)**18+8797093.817109447*(-1.0+0.3183098861837907*x)**19 + 432821.28297364846*(-1.0+0.3183098861837907*x)**20 - 2584215.490575322*(-1.0+0.3183098861837907*x)**21 - 178202.6787844045*(-1.0+0.3183098861837907*x)**22+309558.65888897196*(-1.0+0.3183098861837907*x)**23
        y=y/100
        return y


def rxj_PA_poly(x): #x is phase -- this one is plotted in the original polynmials jupter notebook on your machine
        #defined over a phase of 0-2pi in degrees
        y=-47.51371495007874+148.21418784359244*(-1.0+0.3183098861837907*x)+801.6988958213333*(-1.0+0.3183098861837907*x)**2+4903.389784714883*(-1.0+0.3183098861837907*x)**3+4734.255447873519*(-1.0+0.3183098861837907*x)**4 - 97926.81107908147*(-1.0+0.3183098861837907*x)**5 - 216514.69617385636*(-1.0+0.3183098861837907*x)**6+795099.611389075*(-1.0+0.3183098861837907*x)**7+2118171.3736018618*(-1.0+0.3183098861837907*x)**8 - 3579331.684025455*(-1.0+0.3183098861837907*x)**9 - 10530547.628104404*(-1.0+0.3183098861837907*x)**10+9866246.255630272*(-1.0+0.3183098861837907*x)**11+31113675.137438625*(-1.0+0.3183098861837907*x)**12 - 17315717.76681997*(-1.0+0.3183098861837907*x)**13 - 57728035.585747115*(-1.0+0.3183098861837907*x)**14+19410405.351077296*(-1.0+0.3183098861837907*x)**15+67940681.9315697*(-1.0+0.3183098861837907*x)**16 - 13444291.496889476*(-1.0+0.3183098861837907*x)**17 - 49255143.79648739*(-1.0+0.3183098861837907*x)**18 + 5237070.340115794*(-1.0+0.3183098861837907*x)**19 + 20070327.746498633*(-1.0+0.3183098861837907*x)**20 - 876607.7536815401*(-1.0+0.3183098861837907*x)**21 - 3518073.2554108216*(-1.0+0.3183098861837907*x)**22
        y=np.radians(y)
        return y

def rxj_PA_poly_doro(x):
        y=1.0172227098360436e-20 * x**40 + -2.9331357617626305e-19 * x**39 + 2.387391490653428e-18 * x**38 + 3.0135556080558276e-18 * x**37 + -7.538094648074236e-17 * x**36 + -3.182027489654221e-16 * x**35 + 1.4582194843481271e-15 * x**34 + 1.9033506237032367e-14 * x**33 + 4.8068116927573155e-14 * x**32 + -4.589985510569483e-13 * x**31 + -4.721443737632912e-12 * x**30 + -1.1641593816482196e-11 * x**29 + 1.151844611618483e-10 * x**28 + 1.2204070962730063e-09 * x**27 + 3.193014008433967e-09 * x**26 + -2.980734763928548e-08 * x**25 + -3.1921502902234144e-07 * x**24 + -6.993994054645151e-07 * x**23 + 9.278513383639047e-06 * x**22 + 8.060625225411939e-05 * x**21 + 2.1471288502546246e-05 * x**20 + -0.003318292261948712 * x**19 + -0.01448381594620449 * x**18 + 0.09282480341078345 * x**17 + 0.8607043444859448 * x**16 + -2.6656650845296044 * x**15 + -38.22820183915451 * x**14 + 161.57917114937905 * x**13 + 1178.2736055103835 * x**12 + -13512.247156379306 * x**11 + 62344.31295752182 * x**10 + -175141.00302049995 * x**9 + 328804.9215399549 * x**8 + -422294.72996428923 * x**7 + 366823.34662177286 * x**6 + -206230.0926248056 * x**5 + 68120.61813376007 * x**4 + -10299.397709098299 * x**3 + -29.287843271925937 * x**2 + 51.310082959796304 * x + 10.506234490076459
        y=np.radians(y)
        return y

def rxj_PD_poly_doro(x):
        y=-1.352715282850456e-18 * x**36 + 3.763105615573653e-17 * x**35 + -2.983518149253466e-16 * x**34 + -3.2603072529434926e-16 * x**33 + 8.98543603424602e-15 * x**32 + 3.293851928634345e-14 * x**31 + -1.8872024298490652e-13 * x**30 + -1.9692753538038082e-12 * x**29 + -2.7477715243267664e-12 * x**28 + 5.701576441400153e-11 * x**27 + 4.0168642038163623e-10 * x**26 + 7.051898757497217e-11 * x**25 + -1.4365756371296923e-08 * x**24 + -7.765540103691898e-08 * x**23 + 1.3687052987202457e-07 * x**22 + 3.4928110634855975e-06 * x**21 + 1.0219431118432176e-05 * x**20 + -8.755280128352715e-05 * x**19 + -0.0006546946309217261 * x**18 + 0.0015526323378739538 * x**17 + 0.026288058010144093 * x**16 + -0.04818488301330954 * x**15 + -0.973800769129739 * x**14 + 4.643880349458399 * x**13 + 23.420621165398924 * x**12 + -398.4305998217201 * x**11 + 2530.3048827139414 * x**10 + -9940.983734265532 * x**9 + 26335.746060909765 * x**8 + -47881.733821788555 * x**7 + 59144.033815391136 * x**6 + -48060.22396823785 * x**5 + 24172.55232016002 * x**4 + -6720.63733265986 * x**3 + 812.352771373605 * x**2 + -22.17612683713533 * x + 9.477447306186985
        y=y/100
        return y




#Her X-1
def herx1_PD_poly(x):
    y= 0.00504306249 -0.124489318*x + 1.17939217*x**2 -5.33302258*x**3+ 11.6192006*x**4 -11.0294455*x**5 + 3.66090677*x**6-  8.26120543*x**7
    y=y/100
    return y
 

def herx1_PA_poly(x):
        y=-0.0911197026+  2.21771582*x**1 -21.6415664*x**2+  108.172764*x**3 -293.723974*x**4+ 423.281227*x**5 -287.838318*x**6 + 66.4635420*x**7  +32.9456739*x**8
        y=np.radians(y)
        return y
 
