import Signals as Sig
import matplotlib.pyplot as plt
import numpy as np

S = Sig.Record()
S.ImportSac('data/CH.OTER2..EH2.D.2013.180.091135.SAC')
S.ImportSac('data/CH.OTER2..EH2.D.2013.180.091135.SAC')
S.ImportSac('data/CH.OTER2..EH2.D.2013.180.091135.SAC')

a = S.CHN[0]
S.Fourier(0)
S.Fourier(1)
b = S.CHN[0]

S.Taper(0.2)
c = S.CHN[0]

S.Filter(0.0,50.,4)
d = S.CHN[0]

plt.close('all')

plt.figure()
plt.plot(a,'b')
plt.plot(b,'r')
plt.plot(c,'g')
plt.plot(d,'y')
plt.show(block=False)

