#!/usr/bin/env python

import os, sys
import numpy as np

class database():
  def __init__(self): 
    self.molecules = ['CO2','CO','H2','H2O','NH3','O2','NO','NO2','CH4','H2O2','N2']
    self.kb        = 0.000086173354613780

  def error_param(self):
    print('\n\nERROR: ',self.var,' not yet parametrized.')
    sys.exit(9)

  def call_var(self,name):
    self.var           = name
    self.interpolation = {}

    if self.var == self.molecules[0]:
      self.ZPVE                             =  0.306
      self.interpolation['m_1.0-50.0']      = -0.0014285714285714288
      self.interpolation['b_1.0-50.0']      =  0.0014285714285714318
      self.interpolation['m_50.0-100.0']    = -0.00182
      self.interpolation['b_50.0-100.0']    =  0.02099999999999999
      self.interpolation['m_100.0-150.0']   = -0.00198
      self.interpolation['b_100.0-150.0']   =  0.03699999999999998
      self.interpolation['m_150.0-200.0']   = -0.0020999999999999994
      self.interpolation['b_150.0-200.0']   =  0.05499999999999988
      self.interpolation['m_200.0-250.0']   = -0.0021799999999999996
      self.interpolation['b_200.0-250.0']   =  0.07099999999999995
      self.interpolation['m_250.0-300.0']   = -0.002240000000000002
      self.interpolation['b_250.0-300.0']   =  0.08600000000000052
      self.interpolation['m_300.0-350.0']   = -0.0023199999999999974
      self.interpolation['b_300.0-350.0']   =  0.1099999999999991
      self.interpolation['m_350.0-400.0']   = -0.0023799999999999997
      self.interpolation['b_350.0-400.0']   =  0.1309999999999999
      self.interpolation['m_400.0-450.0']   = -0.00242
      self.interpolation['b_400.0-450.0']   =  0.14700000000000002
      self.interpolation['m_450.0-500.0']   = -0.002480000000000002
      self.interpolation['b_450.0-500.0']   =  0.17400000000000104
      self.interpolation['m_500.0-550.0']   = -0.00254
      self.interpolation['b_500.0-550.0']   =  0.20399999999999996
      self.interpolation['m_550.0-600.0']   = -0.002559999999999998
      self.interpolation['b_550.0-600.0']   =  0.21499999999999897
      self.interpolation['m_600.0-650.0']   = -0.00262
      self.interpolation['b_600.0-650.0']   =  0.2509999999999999
      self.interpolation['m_650.0-700.0']   = -0.002640000000000002
      self.interpolation['b_650.0-700.0']   =  0.26400000000000134
      self.interpolation['m_700.0-750.0']   = -0.0027
      self.interpolation['b_700.0-750.0']   =  0.3059999999999998
      self.interpolation['m_750.0-800.0']   = -0.002719999999999998
      self.interpolation['b_750.0-800.0']   =  0.3209999999999984
      self.interpolation['m_800.0-850.0']   = -0.0027400000000000002
      self.interpolation['b_800.0-850.0']   =  0.3370000000000002
      self.interpolation['m_850.0-900.0']   = -0.0028000000000000026
      self.interpolation['b_850.0-900.0']   =  0.3880000000000021
      self.interpolation['m_900.0-950.0']   = -0.0028200000000000005
      self.interpolation['b_900.0-950.0']   =  0.40600000000000014
      self.interpolation['m_950.0-1000.0']  = -0.002839999999999998
      self.interpolation['b_950.0-1000.0']  =  0.42499999999999805
    elif self.var == self.molecules[1]:
      self.ZPVE                             =  0.132
      self.interpolation['m_1.0-50.0']      = -0.0012448979591836737
      self.interpolation['b_1.0-50.0']      =  0.0012448979591836717
      self.interpolation['m_50.0-100.0']    = -0.0016200000000000001
      self.interpolation['b_50.0-100.0']    =  0.01999999999999999
      self.interpolation['m_100.0-150.0']   = -0.0018
      self.interpolation['b_100.0-150.0']   =  0.038000000000000006
      self.interpolation['m_150.0-200.0']   = -0.00188
      self.interpolation['b_150.0-200.0']   =  0.04999999999999999
      self.interpolation['m_200.0-250.0']   = -0.0019599999999999995
      self.interpolation['b_200.0-250.0']   =  0.06599999999999989
      self.interpolation['m_250.0-300.0']   = -0.0020400000000000006
      self.interpolation['b_250.0-300.0']   =  0.08600000000000019
      self.interpolation['m_300.0-350.0']   = -0.0020599999999999998
      self.interpolation['b_300.0-350.0']   =  0.09199999999999986
      self.interpolation['m_350.0-400.0']   = -0.0021199999999999995
      self.interpolation['b_350.0-400.0']   =  0.11299999999999977
      self.interpolation['m_400.0-450.0']   = -0.0021599999999999996
      self.interpolation['b_400.0-450.0']   =  0.1289999999999999
      self.interpolation['m_450.0-500.0']   = -0.0021999999999999997
      self.interpolation['b_450.0-500.0']   =  0.1469999999999999
      self.interpolation['m_500.0-550.0']   = -0.002220000000000002
      self.interpolation['b_500.0-550.0']   =  0.15700000000000092
      self.interpolation['m_550.0-600.0']   = -0.00226
      self.interpolation['b_550.0-600.0']   =  0.17899999999999983
      self.interpolation['m_600.0-650.0']   = -0.00226
      self.interpolation['b_600.0-650.0']   =  0.17899999999999983
      self.interpolation['m_650.0-700.0']   = -0.002320000000000002
      self.interpolation['b_650.0-700.0']   =  0.2180000000000013
      self.interpolation['m_700.0-750.0']   = -0.002319999999999993
      self.interpolation['b_700.0-750.0']   =  0.2179999999999951
      self.interpolation['m_750.0-800.0']   = -0.00234
      self.interpolation['b_750.0-800.0']   =  0.23300000000000032
      self.interpolation['m_800.0-850.0']   = -0.0023800000000000045
      self.interpolation['b_800.0-850.0']   =  0.2650000000000037
      self.interpolation['m_850.0-900.0']   = -0.0023800000000000045
      self.interpolation['b_850.0-900.0']   =  0.2650000000000037
      self.interpolation['m_900.0-950.0']   = -0.0023999999999999933
      self.interpolation['b_900.0-950.0']   =  0.2829999999999937
      self.interpolation['m_950.0-1000.0']  = -0.00242
      self.interpolation['b_950.0-1000.0']  =  0.30200000000000005
    elif self.var == self.molecules[2]:
      self.ZPVE                             =  0.265
      self.interpolation['m_1.0-50.0']      = -0.0005918367346938781
      self.interpolation['b_1.0-50.0']      =  0.0015918367346938814
      self.interpolation['m_50.0-100.0']    = -0.0009999999999999998
      self.interpolation['b_50.0-100.0']    =  0.021999999999999964
      self.interpolation['m_100.0-150.0']   = -0.00114
      self.interpolation['b_100.0-150.0']   =  0.035999999999999976
      self.interpolation['m_150.0-200.0']   = -0.00126
      self.interpolation['b_150.0-200.0']   =  0.05399999999999999
      self.interpolation['m_200.0-250.0']   = -0.00132
      self.interpolation['b_200.0-250.0']   =  0.066
      self.interpolation['m_250.0-300.0']   = -0.0013800000000000002
      self.interpolation['b_250.0-300.0']   =  0.08100000000000002
      self.interpolation['m_300.0-350.0']   = -0.00144
      self.interpolation['b_300.0-350.0']   =  0.09899999999999998
      self.interpolation['m_350.0-400.0']   = -0.0014799999999999991
      self.interpolation['b_350.0-400.0']   =  0.11299999999999966
      self.interpolation['m_400.0-450.0']   = -0.0015199999999999992
      self.interpolation['b_400.0-450.0']   =  0.12899999999999967
      self.interpolation['m_450.0-500.0']   = -0.0015600000000000015
      self.interpolation['b_450.0-500.0']   =  0.14700000000000069
      self.interpolation['m_500.0-550.0']   = -0.0015799999999999992
      self.interpolation['b_500.0-550.0']   =  0.15699999999999958
      self.interpolation['m_550.0-600.0']   = -0.0016200000000000014
      self.interpolation['b_550.0-600.0']   =  0.17900000000000083
      self.interpolation['m_600.0-650.0']   = -0.0016399999999999993
      self.interpolation['b_600.0-650.0']   =  0.19099999999999961
      self.interpolation['m_650.0-700.0']   = -0.0016599999999999991
      self.interpolation['b_650.0-700.0']   =  0.20399999999999952
      self.interpolation['m_700.0-750.0']   = -0.0016800000000000016
      self.interpolation['b_700.0-750.0']   =  0.21800000000000108
      self.interpolation['m_750.0-800.0']   = -0.0016999999999999993
      self.interpolation['b_750.0-800.0']   =  0.23299999999999943
      self.interpolation['m_800.0-850.0']   = -0.0017200000000000015
      self.interpolation['b_800.0-850.0']   =  0.24900000000000122
      self.interpolation['m_850.0-900.0']   = -0.0017399999999999948
      self.interpolation['b_850.0-900.0']   =  0.2659999999999956
      self.interpolation['m_900.0-950.0']   = -0.0017400000000000037
      self.interpolation['b_900.0-950.0']   =  0.26600000000000357
      self.interpolation['m_950.0-1000.0']  = -0.0017799999999999995
      self.interpolation['b_950.0-1000.0']  =  0.3039999999999994
    elif self.var == self.molecules[3]:
      self.ZPVE                             =  0.563
      self.interpolation['m_1.0-50.0']      = -0.0010816326530612233
      self.interpolation['b_1.0-50.0']      =  0.0010816326530612291
      self.interpolation['m_50.0-100.0']    = -0.0015400000000000004
      self.interpolation['b_50.0-100.0']    =  0.024000000000000077
      self.interpolation['m_100.0-150.0']   = -0.0017000000000000003
      self.interpolation['b_100.0-150.0']   =  0.04000000000000009
      self.interpolation['m_150.0-200.0']   = -0.0018399999999999994
      self.interpolation['b_150.0-200.0']   =  0.06099999999999994
      self.interpolation['m_200.0-250.0']   = -0.0019199999999999994
      self.interpolation['b_200.0-250.0']   =  0.07699999999999996
      self.interpolation['m_250.0-300.0']   = -0.0019999999999999996
      self.interpolation['b_250.0-300.0']   =  0.09699999999999998
      self.interpolation['m_300.0-350.0']   = -0.002040000000000002
      self.interpolation['b_300.0-350.0']   =  0.10900000000000065
      self.interpolation['m_350.0-400.0']   = -0.0020999999999999994
      self.interpolation['b_350.0-400.0']   =  0.12999999999999978
      self.interpolation['m_400.0-450.0']   = -0.0021399999999999995
      self.interpolation['b_400.0-450.0']   =  0.1459999999999998
      self.interpolation['m_450.0-500.0']   = -0.0021799999999999996
      self.interpolation['b_450.0-500.0']   =  0.16399999999999992
      self.interpolation['m_500.0-550.0']   = -0.0022199999999999998
      self.interpolation['b_500.0-550.0']   =  0.18399999999999994
      self.interpolation['m_550.0-600.0']   = -0.00226
      self.interpolation['b_550.0-600.0']   =  0.20599999999999996
      self.interpolation['m_600.0-650.0']   = -0.0022799999999999977
      self.interpolation['b_600.0-650.0']   =  0.21799999999999864
      self.interpolation['m_650.0-700.0']   = -0.002320000000000002
      self.interpolation['b_650.0-700.0']   =  0.24400000000000155
      self.interpolation['m_700.0-750.0']   = -0.00234
      self.interpolation['b_700.0-750.0']   =  0.25800000000000023
      self.interpolation['m_750.0-800.0']   = -0.0023799999999999997
      self.interpolation['b_750.0-800.0']   =  0.2879999999999998
      self.interpolation['m_800.0-850.0']   = -0.0023799999999999997
      self.interpolation['b_800.0-850.0']   =  0.2879999999999998
      self.interpolation['m_850.0-900.0']   = -0.00242
      self.interpolation['b_850.0-900.0']   =  0.32200000000000006
      self.interpolation['m_900.0-950.0']   = -0.002440000000000002
      self.interpolation['b_900.0-950.0']   =  0.34000000000000186
      self.interpolation['m_950.0-1000.0']  = -0.00246
      self.interpolation['b_950.0-1000.0']  =  0.359
    elif self.var == self.molecules[4]:
      self.ZPVE                             =  0.908
      self.interpolation['m_1.0-50.0']      = -0.0011428571428571438
      self.interpolation['b_1.0-50.0']      =  0.00114285714285714
      self.interpolation['m_50.0-100.0']    = -0.0016199999999999993
      self.interpolation['b_50.0-100.0']    =  0.02499999999999991
      self.interpolation['m_100.0-150.0']   = -0.0017799999999999995
      self.interpolation['b_100.0-150.0']   =  0.040999999999999925
      self.interpolation['m_150.0-200.0']   = -0.0019000000000000017
      self.interpolation['b_150.0-200.0']   =  0.059000000000000274
      self.interpolation['m_200.0-250.0']   = -0.0019999999999999996
      self.interpolation['b_200.0-250.0']   =  0.07899999999999985
      self.interpolation['m_250.0-300.0']   = -0.0020599999999999998
      self.interpolation['b_250.0-300.0']   =  0.09399999999999986
      self.interpolation['m_300.0-350.0']   = -0.0021199999999999995
      self.interpolation['b_300.0-350.0']   =  0.11199999999999977
      self.interpolation['m_350.0-400.0']   = -0.0021799999999999996
      self.interpolation['b_350.0-400.0']   =  0.1329999999999999
      self.interpolation['m_400.0-450.0']   = -0.0022199999999999998
      self.interpolation['b_400.0-450.0']   =  0.1489999999999999
      self.interpolation['m_450.0-500.0']   = -0.002280000000000002
      self.interpolation['b_450.0-500.0']   =  0.17600000000000093
      self.interpolation['m_500.0-550.0']   = -0.00232
      self.interpolation['b_500.0-550.0']   =  0.19599999999999995
      self.interpolation['m_550.0-600.0']   = -0.00234
      self.interpolation['b_550.0-600.0']   =  0.20700000000000007
      self.interpolation['m_600.0-650.0']   = -0.0023999999999999976
      self.interpolation['b_600.0-650.0']   =  0.24299999999999855
      self.interpolation['m_650.0-700.0']   = -0.002440000000000002
      self.interpolation['b_650.0-700.0']   =  0.26900000000000146
      self.interpolation['m_700.0-750.0']   = -0.00246
      self.interpolation['b_700.0-750.0']   =  0.2829999999999999
      self.interpolation['m_750.0-800.0']   = -0.0025
      self.interpolation['b_750.0-800.0']   =  0.31299999999999994
      self.interpolation['m_800.0-850.0']   = -0.00254
      self.interpolation['b_800.0-850.0']   =  0.3450000000000002
      self.interpolation['m_850.0-900.0']   = -0.0025600000000000024
      self.interpolation['b_850.0-900.0']   =  0.3620000000000019
      self.interpolation['m_900.0-950.0']   = -0.0025800000000000003
      self.interpolation['b_900.0-950.0']   =  0.3799999999999999
      self.interpolation['m_950.0-1000.0']  = -0.0026399999999999935
      self.interpolation['b_950.0-1000.0']  =  0.4369999999999936      
    elif self.var == self.molecules[5]:
      self.ZPVE                             =  0.097
      self.interpolation['m_1.0-50.0']      = -0.0013673469387755102
      self.interpolation['b_1.0-50.0']      =  0.0013673469387755072
      self.interpolation['m_50.0-100.0']    = -0.00178
      self.interpolation['b_50.0-100.0']    =  0.021999999999999992
      self.interpolation['m_100.0-150.0']   = -0.00192
      self.interpolation['b_100.0-150.0']   =  0.03600000000000003
      self.interpolation['m_150.0-200.0']   = -0.0020199999999999997
      self.interpolation['b_150.0-200.0']   =  0.050999999999999934
      self.interpolation['m_200.0-250.0']   = -0.0020999999999999994
      self.interpolation['b_200.0-250.0']   =  0.06699999999999995
      self.interpolation['m_250.0-300.0']   = -0.0021799999999999996
      self.interpolation['b_250.0-300.0']   =  0.08699999999999997
      self.interpolation['m_300.0-350.0']   = -0.0021999999999999997
      self.interpolation['b_300.0-350.0']   =  0.09299999999999997
      self.interpolation['m_350.0-400.0']   = -0.00226
      self.interpolation['b_350.0-400.0']   =  0.11399999999999999
      self.interpolation['m_400.0-450.0']   = -0.002300000000000002
      self.interpolation['b_400.0-450.0']   =  0.130000000000001
      self.interpolation['m_450.0-500.0']   = -0.00234
      self.interpolation['b_450.0-500.0']   =  0.1479999999999999
      self.interpolation['m_500.0-550.0']   = -0.0023599999999999975
      self.interpolation['b_500.0-550.0']   =  0.1579999999999988
      self.interpolation['m_550.0-600.0']   = -0.002400000000000002
      self.interpolation['b_550.0-600.0']   =  0.18000000000000127
      self.interpolation['m_600.0-650.0']   = -0.00242
      self.interpolation['b_600.0-650.0']   =  0.19199999999999995
      self.interpolation['m_650.0-700.0']   = -0.00246
      self.interpolation['b_650.0-700.0']   =  0.21799999999999997
      self.interpolation['m_700.0-750.0']   = -0.00246
      self.interpolation['b_700.0-750.0']   =  0.21799999999999997
      self.interpolation['m_750.0-800.0']   = -0.0025
      self.interpolation['b_750.0-800.0']   =  0.248
      self.interpolation['m_800.0-850.0']   = -0.002519999999999998
      self.interpolation['b_800.0-850.0']   =  0.26399999999999824
      self.interpolation['m_850.0-900.0']   = -0.00254
      self.interpolation['b_850.0-900.0']   =  0.28100000000000014
      self.interpolation['m_900.0-950.0']   = -0.0025600000000000024
      self.interpolation['b_900.0-950.0']   =  0.29900000000000215
      self.interpolation['m_950.0-1000.0']  = -0.0025800000000000003
      self.interpolation['b_950.0-1000.0']  =  0.31800000000000006
    elif self.var == self.molecules[6]:
      self.ZPVE                             =  0.119
      self.interpolation['m_1.0-50.0']      = -0.001326530612244898
      self.interpolation['b_1.0-50.0']      =  0.0013265306122449
      self.interpolation['m_50.0-100.0']    = -0.0017
      self.interpolation['b_50.0-100.0']    =  0.01999999999999999
      self.interpolation['m_100.0-150.0']   = -0.0018599999999999999
      self.interpolation['b_100.0-150.0']   =  0.035999999999999976
      self.interpolation['m_150.0-200.0']   = -0.0019799999999999996
      self.interpolation['b_150.0-200.0']   =  0.05399999999999994
      self.interpolation['m_200.0-250.0']   = -0.0020400000000000006
      self.interpolation['b_200.0-250.0']   =  0.06600000000000011
      self.interpolation['m_250.0-300.0']   = -0.0020999999999999986
      self.interpolation['b_250.0-300.0']   =  0.08099999999999963
      self.interpolation['m_300.0-350.0']   = -0.0021600000000000018
      self.interpolation['b_300.0-350.0']   =  0.09900000000000064
      self.interpolation['m_350.0-400.0']   = -0.0021999999999999997
      self.interpolation['b_350.0-400.0']   =  0.11299999999999988
      self.interpolation['m_400.0-450.0']   = -0.00224
      self.interpolation['b_400.0-450.0']   =  0.129
      self.interpolation['m_450.0-500.0']   = -0.00226
      self.interpolation['b_450.0-500.0']   =  0.1379999999999999
      self.interpolation['m_500.0-550.0']   = -0.002320000000000002
      self.interpolation['b_500.0-550.0']   =  0.16800000000000104
      self.interpolation['m_550.0-600.0']   = -0.0023199999999999974
      self.interpolation['b_550.0-600.0']   =  0.16799999999999837
      self.interpolation['m_600.0-650.0']   = -0.0023600000000000023
      self.interpolation['b_600.0-650.0']   =  0.1920000000000015
      self.interpolation['m_650.0-700.0']   = -0.0023799999999999997
      self.interpolation['b_650.0-700.0']   =  0.20499999999999985
      self.interpolation['m_700.0-750.0']   = -0.00242
      self.interpolation['b_700.0-750.0']   =  0.23299999999999987
      self.interpolation['m_750.0-800.0']   = -0.00242
      self.interpolation['b_750.0-800.0']   =  0.23299999999999987
      self.interpolation['m_800.0-850.0']   = -0.00246
      self.interpolation['b_800.0-850.0']   =  0.2649999999999997
      self.interpolation['m_850.0-900.0']   = -0.00246
      self.interpolation['b_850.0-900.0']   =  0.2649999999999999
      self.interpolation['m_900.0-950.0']   = -0.0024999999999999957
      self.interpolation['b_900.0-950.0']   =  0.30099999999999616
      self.interpolation['m_950.0-1000.0']  = -0.0025
      self.interpolation['b_950.0-1000.0']  =  0.30100000000000016      
    elif self.var == self.molecules[7]:
      self.ZPVE                             =  0.233
      self.interpolation['m_1.0-50.0']      = -0.001591836734693878
      self.interpolation['b_1.0-50.0']      =  0.0005918367346938735
      self.interpolation['m_50.0-100.0']    = -0.00206
      self.interpolation['b_50.0-100.0']    =  0.023999999999999994
      self.interpolation['m_100.0-150.0']   = -0.0022400000000000002
      self.interpolation['b_100.0-150.0']   =  0.04199999999999998
      self.interpolation['m_150.0-200.0']   = -0.0023599999999999997
      self.interpolation['b_150.0-200.0']   =  0.05999999999999989
      self.interpolation['m_200.0-250.0']   = -0.00244
      self.interpolation['b_200.0-250.0']   =  0.07599999999999996
      self.interpolation['m_250.0-300.0']   = -0.00252
      self.interpolation['b_250.0-300.0']   =  0.09599999999999997
      self.interpolation['m_300.0-350.0']   = -0.0025800000000000003
      self.interpolation['b_300.0-350.0']   =  0.1140000000000001
      self.interpolation['m_350.0-400.0']   = -0.002659999999999998
      self.interpolation['b_350.0-400.0']   =  0.14199999999999924
      self.interpolation['m_400.0-450.0']   = -0.0026800000000000023
      self.interpolation['b_400.0-450.0']   =  0.15000000000000102
      self.interpolation['m_450.0-500.0']   = -0.002759999999999998
      self.interpolation['b_450.0-500.0']   =  0.18599999999999905
      self.interpolation['m_500.0-550.0']   = -0.0027800000000000047
      self.interpolation['b_500.0-550.0']   =  0.1960000000000024
      self.interpolation['m_550.0-600.0']   = -0.002839999999999998
      self.interpolation['b_550.0-600.0']   =  0.22899999999999876
      self.interpolation['m_600.0-650.0']   = -0.00286
      self.interpolation['b_600.0-650.0']   =  0.24099999999999988
      self.interpolation['m_650.0-700.0']   = -0.002919999999999998
      self.interpolation['b_650.0-700.0']   =  0.2799999999999987
      self.interpolation['m_700.0-750.0']   = -0.0029400000000000003
      self.interpolation['b_700.0-750.0']   =  0.29400000000000004
      self.interpolation['m_750.0-800.0']   = -0.0029800000000000004
      self.interpolation['b_750.0-800.0']   =  0.3240000000000003
      self.interpolation['m_800.0-850.0']   = -0.0029999999999999983
      self.interpolation['b_800.0-850.0']   =  0.3399999999999985
      self.interpolation['m_850.0-900.0']   = -0.0030400000000000028
      self.interpolation['b_850.0-900.0']   =  0.37400000000000233
      self.interpolation['m_900.0-950.0']   = -0.003079999999999998
      self.interpolation['b_900.0-950.0']   =  0.40999999999999837
      self.interpolation['m_950.0-1000.0']  = -0.003100000000000005
      self.interpolation['b_950.0-1000.0']  =  0.4290000000000047    
    elif self.var == self.molecules[8]:
      self.ZPVE                             =  1.182
      self.interpolation['m_1.0-50.0']      = -0.00120408163265306
      self.interpolation['b_1.0-50.0']      =  0.0012040816326530576
      self.interpolation['m_50.0-100.0']    = -0.0016599999999999991
      self.interpolation['b_50.0-100.0']    =  0.02400000000000002
      self.interpolation['m_100.0-150.0']   = -0.0018400000000000016
      self.interpolation['b_100.0-150.0']   =  0.04200000000000026
      self.interpolation['m_150.0-200.0']   = -0.0019599999999999995
      self.interpolation['b_150.0-200.0']   =  0.05999999999999994
      self.interpolation['m_200.0-250.0']   = -0.0020399999999999997
      self.interpolation['b_200.0-250.0']   =  0.07599999999999996
      self.interpolation['m_250.0-300.0']   = -0.0021199999999999995
      self.interpolation['b_250.0-300.0']   =  0.09599999999999997
      self.interpolation['m_300.0-350.0']   = -0.0021799999999999996
      self.interpolation['b_300.0-350.0']   =  0.11399999999999999
      self.interpolation['m_350.0-400.0']   = -0.00224
      self.interpolation['b_350.0-400.0']   =  0.135
      self.interpolation['m_400.0-450.0']   = -0.002280000000000002
      self.interpolation['b_400.0-450.0']   =  0.1510000000000009
      self.interpolation['m_450.0-500.0']   = -0.00234
      self.interpolation['b_450.0-500.0']   =  0.17799999999999994
      self.interpolation['m_500.0-550.0']   = -0.0023999999999999976
      self.interpolation['b_500.0-550.0']   =  0.20799999999999885
      self.interpolation['m_550.0-600.0']   = -0.00242
      self.interpolation['b_550.0-600.0']   =  0.21900000000000008
      self.interpolation['m_600.0-650.0']   = -0.002480000000000002
      self.interpolation['b_600.0-650.0']   =  0.25500000000000145
      self.interpolation['m_650.0-700.0']   = -0.00254
      self.interpolation['b_650.0-700.0']   =  0.29400000000000004
      self.interpolation['m_700.0-750.0']   = -0.002559999999999998
      self.interpolation['b_700.0-750.0']   =  0.3079999999999987
      self.interpolation['m_750.0-800.0']   = -0.0026000000000000025
      self.interpolation['b_750.0-800.0']   =  0.33800000000000185
      self.interpolation['m_800.0-850.0']   = -0.00266
      self.interpolation['b_800.0-850.0']   =  0.3860000000000001
      self.interpolation['m_850.0-900.0']   = -0.002679999999999998
      self.interpolation['b_850.0-900.0']   =  0.40299999999999825
      self.interpolation['m_900.0-950.0']   = -0.0027400000000000002
      self.interpolation['b_900.0-950.0']   =  0.4570000000000003
      self.interpolation['m_950.0-1000.0']  = -0.002759999999999998
      self.interpolation['b_950.0-1000.0']  =  0.4759999999999982
    elif self.var == self.molecules[9]:
      self.ZPVE                             =  0.694
      self.interpolation['m_1.0-50.0']      = -0.0014285714285714275
      self.interpolation['b_1.0-50.0']      =  0.00042857142857141706
      self.interpolation['m_50.0-100.0']    = -0.0018799999999999995
      self.interpolation['b_50.0-100.0']    =  0.02300000000000002
      self.interpolation['m_100.0-150.0']   = -0.0020800000000000007
      self.interpolation['b_100.0-150.0']   =  0.04300000000000015
      self.interpolation['m_150.0-200.0']   = -0.0021999999999999997
      self.interpolation['b_150.0-200.0']   =  0.061
      self.interpolation['m_200.0-250.0']   = -0.0023
      self.interpolation['b_200.0-250.0']   =  0.08100000000000002
      self.interpolation['m_250.0-300.0']   = -0.002380000000000001
      self.interpolation['b_250.0-300.0']   =  0.10100000000000031
      self.interpolation['m_300.0-350.0']   = -0.00246
      self.interpolation['b_300.0-350.0']   =  0.125
      self.interpolation['m_350.0-400.0']   = -0.00254
      self.interpolation['b_350.0-400.0']   =  0.15300000000000002
      self.interpolation['m_400.0-450.0']   = -0.0025800000000000003
      self.interpolation['b_400.0-450.0']   =  0.16900000000000004
      self.interpolation['m_450.0-500.0']   = -0.00266
      self.interpolation['b_450.0-500.0']   =  0.20500000000000007
      self.interpolation['m_500.0-550.0']   = -0.002699999999999996
      self.interpolation['b_500.0-550.0']   =  0.22499999999999787
      self.interpolation['m_550.0-600.0']   = -0.0027400000000000002
      self.interpolation['b_550.0-600.0']   =  0.24700000000000033
      self.interpolation['m_600.0-650.0']   = -0.0028000000000000026
      self.interpolation['b_600.0-650.0']   =  0.2830000000000017
      self.interpolation['m_650.0-700.0']   = -0.002839999999999998
      self.interpolation['b_650.0-700.0']   =  0.3089999999999986
      self.interpolation['m_700.0-750.0']   = -0.0028800000000000023
      self.interpolation['b_700.0-750.0']   =  0.33700000000000196
      self.interpolation['m_750.0-800.0']   = -0.002919999999999998
      self.interpolation['b_750.0-800.0']   =  0.36699999999999866
      self.interpolation['m_800.0-850.0']   = -0.0029800000000000004
      self.interpolation['b_800.0-850.0']   =  0.4150000000000005
      self.interpolation['m_850.0-900.0']   = -0.0029999999999999983
      self.interpolation['b_850.0-900.0']   =  0.4319999999999986
      self.interpolation['m_900.0-950.0']   = -0.003020000000000005
      self.interpolation['b_900.0-950.0']   =  0.4500000000000046
      self.interpolation['m_950.0-1000.0']  = -0.003079999999999998
      self.interpolation['b_950.0-1000.0']  =  0.5069999999999983            
    elif self.var == self.molecules[10]:
      self.ZPVE                             =  0.150
      self.interpolation['m_1.0-50.0']      = -0.0012244897959183673
      self.interpolation['b_1.0-50.0']      =  0.0012244897959183682
      self.interpolation['m_50.0-100.0']    = -0.0016399999999999997
      self.interpolation['b_50.0-100.0']    =  0.021999999999999992
      self.interpolation['m_100.0-150.0']   = -0.00178
      self.interpolation['b_100.0-150.0']   =  0.035999999999999976
      self.interpolation['m_150.0-200.0']   = -0.0018799999999999995
      self.interpolation['b_150.0-200.0']   =  0.050999999999999934
      self.interpolation['m_200.0-250.0']   = -0.0019600000000000017
      self.interpolation['b_200.0-250.0']   =  0.06700000000000039
      self.interpolation['m_250.0-300.0']   = -0.0020199999999999997
      self.interpolation['b_250.0-300.0']   =  0.08199999999999985
      self.interpolation['m_300.0-350.0']   = -0.00208
      self.interpolation['b_300.0-350.0']   =  0.09999999999999998
      self.interpolation['m_350.0-400.0']   = -0.0021199999999999995
      self.interpolation['b_350.0-400.0']   =  0.11399999999999977
      self.interpolation['m_400.0-450.0']   = -0.0021399999999999995
      self.interpolation['b_400.0-450.0']   =  0.12199999999999978
      self.interpolation['m_450.0-500.0']   = -0.002200000000000002
      self.interpolation['b_450.0-500.0']   =  0.1490000000000009
      self.interpolation['m_500.0-550.0']   = -0.0022199999999999998
      self.interpolation['b_500.0-550.0']   =  0.1589999999999998
      self.interpolation['m_550.0-600.0']   = -0.0022399999999999976
      self.interpolation['b_550.0-600.0']   =  0.1699999999999986
      self.interpolation['m_600.0-650.0']   = -0.0022799999999999977
      self.interpolation['b_600.0-650.0']   =  0.19399999999999862
      self.interpolation['m_650.0-700.0']   = -0.0023
      self.interpolation['b_650.0-700.0']   =  0.20700000000000007
      self.interpolation['m_700.0-750.0']   = -0.002320000000000002
      self.interpolation['b_700.0-750.0']   =  0.22100000000000164
      self.interpolation['m_750.0-800.0']   = -0.00234
      self.interpolation['b_750.0-800.0']   =  0.2360000000000002
      self.interpolation['m_800.0-850.0']   = -0.0023600000000000023
      self.interpolation['b_800.0-850.0']   =  0.252000000000002
      self.interpolation['m_850.0-900.0']   = -0.0023799999999999997
      self.interpolation['b_850.0-900.0']   =  0.2689999999999999
      self.interpolation['m_900.0-950.0']   = -0.0023999999999999976
      self.interpolation['b_900.0-950.0']   =  0.2869999999999977
      self.interpolation['m_950.0-1000.0']  = -0.00242
      self.interpolation['b_950.0-1000.0']  =  0.30600000000000005
    else:
      print('\n\nERROR: ',self.var,' not found. It either needs to be parametrized')
      print('                  ... or the element/molecule does not exist.')
      sys.exit(10)

def error_not_recognized(keyword,var,example_list):
  print('\n\nERROR:',keyword,' \"',var,'\" not recognized. Try e.g. ', end='')
  for i in example_list:
    print(str(i)+',', end='')
  print('...')
  sys.exit(2)

def eval_var(string,var,example_list):
  try:
    variable = [float(i) for i in var.split(',')]
    return variable
  except ValueError:
    error_not_recognized(string,var,example_list)

print("\n\n")
print("######################################################################")
print("#                                                                    #")
print("#              Chemical potential \u0394\u03BC(0 \u2192 T) calculator               #")
print("#                                                                    #")
print("#              Parametrized with PBE-D3 using Rigid                  #")
print("#              Rotator/Translator and Harmonic                       #")
print("#              Oscillator approximation                              #")
print("#                                                                    #")
print("#              Reliable for temperatures below 1100K                 #")
print("#                                                                    #")
print("#              Temperatures below 1 K are not supported              #")
print("#                                                                    #")
print("######################################################################")
print("\n")
print("  ___________________________________________________________________")
print(" | INSTRUCTIONS:                                                     |")
print(" |   - State your temperatures and pressures separated by commas     |")
print(" |     This script will create all combinations of your parameters   |")
print(" |     For single temperature/pressure, simply write the number      |")
print(" |                                                                   |")
print(" |                                                                   |")
print(" | Example:                                                          |")
print(" |   - Molecule        : CO2                                         |")
print(" |     Temperature (K) : 150,160,170                                 |")
print(" |     Pressure (bar)  : 1,2,3                                       |")
print(" |                                                                   |")
print(" |   Calculates CO2 chemical potential for 150 K - 1,2,3 bar         |")
print(" |                                         160 K - 1,2,3 bar         |")
print(" |                                         170 K - 1,2,3 bar         |")
print(" |                                                                   |")
print(" |___________________________________________________________________|")
print("  Molecules available:", *database().molecules)
print("\n\n")

molecule  = input("Molecule        : ")
data      = database()

data.call_var(molecule)


temp        = input("Temperature (K) : ")
temperature = eval_var('Temperature',temp,[200,300,400])


press    = input("Pressure (bar)  : ")
pressure = eval_var('Pressure',press,[1,2,3])
ZPVE     = data.ZPVE

print("\nG    = ZPVE + \u0394\u03BC(0 \u2192 T)\n\nZPVE =","{:.3f}".format(round(ZPVE,3))," eV \n")

for i in temperature:
  for j in pressure:
    for val in data.interpolation:
      temprange1=float(val.split('_')[1].split('-')[0])
      temprange2=float(val.split('_')[1].split('-')[1])
      if i >= temprange1 and i <= temprange2:
        if val.split('_')[0] == 'b':
          b = data.interpolation[val]
        elif val.split('_')[0] == 'm':
          m = data.interpolation[val]
      elif i > 1000:
        b = data.interpolation['b_950.0-1000.0']
        m = data.interpolation['m_950.0-1000.0']
 
    chempot_analytical= m*i+b +data.kb*i*np.log(j)
    G                 = ZPVE+chempot_analytical
    print("\u0394\u03BC(0 \u2192 T) at",i,"K,",j,"bar :","{:.3f}".format(round(chempot_analytical,3))," eV,  G :","{:.3f}".format(round(G,3))," eV")
print("")
if max(temperature) > 1000:
    print("\n\nWARNING:   Temperatures above 1000K detected!!!")
    print("           Parametrization done between 0 and 1000K")
    print("           Values between 1000K and 1100K are still acceptable.")
    print("           Values above 1100K are unreliable.\n")


