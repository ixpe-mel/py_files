#rxj0440 pd and pa polynomials
def rxj_PD_poly(x): #x is phase
        y=13.599084715966825 - 77.14279098237614*(-1.0+0.3183098861837907*x) - 58.68027525042033*(-1.0+0.3183098861837907*x)**2+2459.257155713502*(-1.0+0.3183098861837907*x)**3 - 585.7971512178035 * (-1.0+0.3183098861837907*x)**4 - 42389.345465598555*(-1.0+0.3183098861837907*x)**5+31377.701921467713*(-1.0+0.3183098861837907*x)**6+389184.2681757424*(-1.0+0.3183098861837907*x)**7 - 295253.60901547165*(-1.0+0.3183098861837907*x)**8 - 2047779.3589603745*(-1.0+0.3183098861837907*x)**9+1270324.491309807*(-1.0+0.3183098861837907*x)**10+6643443.7902267715*(-1.0+0.3183098861837907*x)**11 - 2966775.4033494988*(-1.0+0.3183098861837907*x)**12 - 13851400.5918384*(-1.0+0.3183098861837907*x)**13 + 3872312.125065608*(-1.0+0.3183098861837907*x)**14 + 18812624.00901786*(-1.0+0.3183098861837907*x)**15 - 2556606.7984732436*(-1.0+0.3183098861837907*x)**16 - 16428501.620406935*(-1.0+0.3183098861837907*x)**17 + 390641.771729522*(-1.0+0.3183098861837907*x)**18+8797093.817109447*(-1.0+0.3183098861837907*x)**19 + 432821.28297364846*(-1.0+0.3183098861837907*x)**20 - 2584215.490575322*(-1.0+0.3183098861837907*x)**21 - 178202.6787844045*(-1.0+0.3183098861837907*x)**22+309558.65888897196*(-1.0+0.3183098861837907*x)**23
        
        y=y/100
        return y

def rxj_PA_poly(x): #x is phase
        y=-47.51371495007874+148.21418784359244*(-1.0+0.3183098861837907*x)+801.6988958213333*(-1.0+0.3183098861837907*x)**2+4903.389784714883*(-1.0+0.3183098861837907*x)**3+4734.255447873519*(-1.0+0.3183098861837907*x)**4 - 97926.81107908147*(-1.0+0.3183098861837907*x)**5 - 216514.69617385636*(-1.0+0.3183098861837907*x)**6+795099.611389075*(-1.0+0.3183098861837907*x)**7+2118171.3736018618*(-1.0+0.3183098861837907*x)**8 - 3579331.684025455*(-1.0+0.3183098861837907*x)**9 - 10530547.628104404*(-1.0+0.3183098861837907*x)**10+9866246.255630272*(-1.0+0.3183098861837907*x)**11+31113675.137438625*(-1.0+0.3183098861837907*x)**12 - 17315717.76681997*(-1.0+0.3183098861837907*x)**13 - 57728035.585747115*(-1.0+0.3183098861837907*x)**14+19410405.351077296*(-1.0+0.3183098861837907*x)**15+67940681.9315697*(-1.0+0.3183098861837907*x)**16 - 13444291.496889476*(-1.0+0.3183098861837907*x)**17 - 49255143.79648739*(-1.0+0.3183098861837907*x)**18 + 5237070.340115794*(-1.0+0.3183098861837907*x)**19 + 20070327.746498633*(-1.0+0.3183098861837907*x)**20 - 876607.7536815401*(-1.0+0.3183098861837907*x)**21 - 3518073.2554108216*(-1.0+0.3183098861837907*x)**22
        y=y/100
        return y

#Her X-1
def herx1_PD_poly(x):
    y= 0.00504306249 -0.124489318*x + 1.17939217*x**2 -5.33302258*x**3+ 11.6192006*x**4 -11.0294455*x**5 + 3.66090677*x**6-  8.26120543*x**7
    y=np.radians(y)
    return y
 

def herx1_PA_poly(x):
        y=-0.0911197026+  2.21771582*x**1 -21.6415664*x**2+  108.172764*x**3 -293.723974*x**4+ 423.281227*x**5 -287.838318*x**6 + 66.4635420*x**7  +32.9456739*x**8
        y=np.radians(y)
        return y
 
