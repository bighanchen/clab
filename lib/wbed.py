#!/usr/bin/python
#0

""" Module/Scripts Description

Copyright (c) 2008 Yunfei Wang <tszn1984@gmail.com>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  Yunfei Wang
@contact: tszn1984@gmail.com
@modified from: zbed by nimezhu@163.com
"""

# ------------------------------------
# python modules
# ------------------------------------

import os,sys,string,random,cPickle,gzip,copy,wRNA,numpy
from math import log,sqrt
from bisect import bisect
from zSeqIO import *
from wWigIO import *

# ------------------------------------
# constants
# ------------------------------------

_numBins   = 37450;
_binLevels = 7;
_RefSeqAccnToChr={"NC_000001.10":"1","NC_000002.11":"2","NC_000003.11":"3","NC_000004.11":"4","NC_000005.9":"5","NC_000006.11":"6","NC_000007.13":"7","NC_000008.10":"8","NC_000009.11":"9","NC_000010.10":"10","NC_000011.9":"11","NC_000012.11":"12","NC_000013.10":"13","NC_000014.8":"14","NC_000015.9":"15","NC_000016.9":"16","NC_000017.10":"17","NC_000018.9":"18","NC_000019.9":"19","NC_000020.10":"20","NC_000021.8":"21","NC_000022.10":"22","NC_000023.10":"X","NC_000024.9":"Y","NT_113878.1":"HSCHR1_RANDOM_CTG5","NT_167207.1":"HSCHR1_RANDOM_CTG12","NT_113885.1":"HSCHR4_RANDOM_CTG2","NT_113888.1":"HSCHR4_RANDOM_CTG3","NT_113901.1":"HSCHR7_RANDOM_CTG1","NT_113909.1":"HSCHR8_RANDOM_CTG1","NT_113907.1":"HSCHR8_RANDOM_CTG4","NT_113914.1":"HSCHR9_RANDOM_CTG1","NT_113916.2":"HSCHR9_RANDOM_CTG2","NT_113915.1":"HSCHR9_RANDOM_CTG4","NT_113911.1":"HSCHR9_RANDOM_CTG5","NT_113921.2":"HSCHR11_RANDOM_CTG2","NT_113941.1":"HSCHR17_RANDOM_CTG1","NT_113943.1":"HSCHR17_RANDOM_CTG2","NT_113930.1":"HSCHR17_RANDOM_CTG3","NT_113945.1":"HSCHR17_RANDOM_CTG4","NT_113947.1":"HSCHR18_RANDOM_CTG1","NT_113948.1":"HSCHR19_RANDOM_CTG1","NT_113949.1":"HSCHR19_RANDOM_CTG2","NT_113950.2":"HSCHR21_RANDOM_CTG9","NT_113961.1":"HSCHRUN_RANDOM_CTG1","NT_113923.1":"HSCHRUN_RANDOM_CTG2","NT_167208.1":"HSCHRUN_RANDOM_CTG3","NT_167209.1":"HSCHRUN_RANDOM_CTG4","NT_167210.1":"HSCHRUN_RANDOM_CTG5","NT_167211.1":"HSCHRUN_RANDOM_CTG6","NT_167212.1":"HSCHRUN_RANDOM_CTG7","NT_113889.1":"HSCHRUN_RANDOM_CTG9","NT_167213.1":"HSCHRUN_RANDOM_CTG10","NT_167214.1":"HSCHRUN_RANDOM_CTG11","NT_167215.1":"HSCHRUN_RANDOM_CTG13","NT_167216.1":"HSCHRUN_RANDOM_CTG14","NT_167217.1":"HSCHRUN_RANDOM_CTG15","NT_167218.1":"HSCHRUN_RANDOM_CTG16","NT_167219.1":"HSCHRUN_RANDOM_CTG17","NT_167220.1":"HSCHRUN_RANDOM_CTG19","NT_167221.1":"HSCHRUN_RANDOM_CTG20","NT_167222.1":"HSCHRUN_RANDOM_CTG21","NT_167223.1":"HSCHRUN_RANDOM_CTG22","NT_167224.1":"HSCHRUN_RANDOM_CTG23","NT_167225.1":"HSCHRUN_RANDOM_CTG24","NT_167226.1":"HSCHRUN_RANDOM_CTG25","NT_167227.1":"HSCHRUN_RANDOM_CTG26","NT_167228.1":"HSCHRUN_RANDOM_CTG27","NT_167229.1":"HSCHRUN_RANDOM_CTG28","NT_167230.1":"HSCHRUN_RANDOM_CTG29","NT_167231.1":"HSCHRUN_RANDOM_CTG30","NT_167232.1":"HSCHRUN_RANDOM_CTG31","NT_167233.1":"HSCHRUN_RANDOM_CTG32","NT_167234.1":"HSCHRUN_RANDOM_CTG33","NT_167235.1":"HSCHRUN_RANDOM_CTG34","NT_167236.1":"HSCHRUN_RANDOM_CTG35","NT_167237.1":"HSCHRUN_RANDOM_CTG36","NT_167238.1":"HSCHRUN_RANDOM_CTG37","NT_167239.1":"HSCHRUN_RANDOM_CTG38","NT_167240.1":"HSCHRUN_RANDOM_CTG39","NT_167241.1":"HSCHRUN_RANDOM_CTG40","NT_167242.1":"HSCHRUN_RANDOM_CTG41","NT_167243.1":"HSCHRUN_RANDOM_CTG42","NT_167244.1":"HSCHR6_MHC_APD_CTG1","NT_113891.2":"HSCHR6_MHC_COX_CTG1","NT_167245.1":"HSCHR6_MHC_DBB_CTG1","NT_167246.1":"HSCHR6_MHC_MANN_CTG1","NT_167247.1":"HSCHR6_MHC_MCF_CTG1","NT_167248.1":"HSCHR6_MHC_QBL_CTG1","NT_167249.1":"HSCHR6_MHC_SSTO_CTG1","NT_167250.1":"HSCHR4_1_CTG9","NT_167251.1":"HSCHR17_1_CTG5","NC_012920.1":"MT","NW_004070864.1":"HG1472_PATCH","NW_003571030.1":"HG989_PATCH","NW_003871056.1":"HG1292_PATCH","NW_003871055.2":"HG1287_PATCH","NW_003315905.1":"HSCHR1_1_CTG31","NW_003315906.1":"HSCHR1_2_CTG31","NW_003315907.1":"HSCHR1_3_CTG31","NW_004070863.1":"HG1471_PATCH","NW_003871057.1":"HG1293_PATCH","NW_004070865.1":"HG1473_PATCH","NW_003315903.1":"HG999_1_PATCH","NW_003315904.1":"HG999_2_PATCH","NW_003315908.1":"HSCHR2_1_CTG1","NW_003571032.1":"HG686_PATCH","NW_003571033.2":"HSCHR2_2_CTG12","NW_003315909.1":"HSCHR2_1_CTG12","NW_003571031.1":"HG1007_PATCH","NW_003871060.1":"HSCHR3_1_CTG1","NW_003871059.1":"HG325_PATCH","NW_003315910.1":"HG186_PATCH","NW_003315911.1":"HG280_PATCH","NW_003871058.1":"HG1091_PATCH","NW_003315912.1":"HG991_PATCH","NW_003315913.1":"HSCHR3_1_CTG2_1","NW_003315915.1":"HSCHR4_1_CTG6","NW_003315916.1":"HSCHR4_2_CTG9","NW_003571035.1":"HG706_PATCH","NW_003315914.1":"HSCHR4_1_CTG12","NW_003571034.1":"HG1032_PATCH","NW_003315920.1":"HSCHR5_2_CTG1","NW_003571036.1":"HSCHR5_3_CTG1","NW_003315917.2":"HSCHR5_1_CTG1","NW_003315918.1":"HSCHR5_1_CTG2","NW_003871061.1":"HG1063_PATCH","NW_003315919.1":"HSCHR5_1_CTG5","NW_004070866.1":"HG27_PATCH","NW_003871063.1":"HG1322_PATCH","NW_003315921.1":"HSCHR6_1_CTG5","NW_003871062.1":"HG1304_PATCH","NW_003571039.1":"HG736_PATCH","NW_003571038.1":"HG14_PATCH","NW_003871064.1":"HG1257_PATCH","NW_003571041.1":"HG946_PATCH","NW_003571037.1":"HG115_PATCH","NW_003871065.1":"HG1308_PATCH","NW_003315922.2":"HSCHR7_1_CTG6","NW_003571040.1":"HG7_PATCH","NW_003571042.1":"HG19_PATCH","NW_003871066.1":"HG418_PATCH","NW_003315923.1":"HG104_HG975_PATCH","NW_003315924.1":"HG243_PATCH","NW_003315928.1":"HSCHR9_1_CTG1","NW_003871067.1":"HG962_PATCH","NW_003315929.1":"HSCHR9_1_CTG35","NW_003315930.1":"HSCHR9_2_CTG35","NW_003315931.1":"HSCHR9_3_CTG35","NW_004070869.1":"HG1502_PATCH","NW_003315925.1":"HG79_PATCH","NW_004070867.1":"HG1500_PATCH","NW_004070868.1":"HG1501_PATCH","NW_003315926.1":"HG998_1_PATCH","NW_003315927.1":"HG998_2_PATCH","NW_003571043.1":"HG905_PATCH","NW_003871071.1":"HG871_PATCH","NW_003315932.1":"HG544_PATCH","NW_003315934.1":"HSCHR10_1_CTG2","NW_003315935.1":"HSCHR10_1_CTG5","NW_003871068.1":"HG1211_PATCH","NW_003871070.1":"HG339_PATCH","NW_003871069.1":"HG311_PATCH","NW_003315933.1":"HG995_PATCH","NW_004070870.1":"HG1479_PATCH","NW_003871075.1":"HG256_PATCH","NW_003871082.1":"HG873_PATCH","NW_003315936.1":"HSCHR11_1_CTG1_1","NW_003571045.1":"HG281_PATCH","NW_003871073.1":"HG142_HG150_NOVEL_TEST","NW_003871074.1":"HG151_NOVEL_TEST","NW_003571046.1":"HG536_PATCH","NW_004070871.1":"HG865_PATCH","NW_003871081.1":"HG414_PATCH","NW_003871079.1":"HG348_PATCH","NW_003871077.1":"HG305_PATCH","NW_003871080.1":"HG388_HG400_PATCH","NW_003871078.1":"HG306_PATCH","NW_003871072.1":"HG122_PATCH","NW_003871076.1":"HG299_PATCH","NW_003571048.1":"HG858_PATCH","NW_003571049.1":"HSCHR12_1_CTG1","NW_003871083.1":"HG344_PATCH","NW_003571047.1":"HG1133_PATCH","NW_003571050.1":"HSCHR12_2_CTG2","NW_003315938.1":"HSCHR12_1_CTG2","NW_003315939.1":"HSCHR12_1_CTG2_1","NW_003315941.1":"HSCHR12_2_CTG2_1","NW_003315942.2":"HSCHR12_3_CTG2_1","NW_003315940.1":"HSCHR12_1_CTG5","NW_003315937.1":"HG996_PATCH","NW_003571051.1":"HG531_PATCH","NW_003315943.1":"HSCHR15_1_CTG4","NW_003315944.1":"HSCHR15_1_CTG8","NW_003871084.1":"HG971_PATCH","NW_003315945.1":"HSCHR16_1_CTG3_1","NW_003871085.1":"HG1208_PATCH","NW_003315946.1":"HSCHR16_2_CTG3_1","NW_004070872.1":"HG417_PATCH","NW_003315952.1":"HSCHR17_1_CTG1","NW_003315951.1":"HG990_PATCH","NW_003315950.1":"HG987_PATCH","NW_003871090.1":"HG883_PATCH","NW_003315949.1":"HG75_PATCH","NW_003315948.1":"HG745_PATCH","NW_003871091.1":"HSCHR17_4_CTG4","NW_003871093.1":"HSCHR17_6_CTG4","NW_003871092.1":"HSCHR17_5_CTG4","NW_003315953.1":"HSCHR17_1_CTG4","NW_003571052.1":"HG185_PATCH","NW_003871086.1":"HG1146_PATCH","NW_003315947.1":"HG183_PATCH","NW_003871088.1":"HG747_PATCH","NW_003315954.1":"HSCHR17_2_CTG4","NW_003315955.1":"HSCHR17_3_CTG4","NW_003871089.1":"HG748_PATCH","NW_003871087.1":"HG271_PATCH","NW_003315956.1":"HSCHR18_1_CTG1_1","NW_003315959.1":"HSCHR18_2_CTG1_1","NW_003315960.1":"HSCHR18_2_CTG2","NW_003315957.1":"HSCHR18_1_CTG2","NW_003315958.1":"HSCHR18_1_CTG2_1","NW_003315961.1":"HSCHR18_2_CTG2_1","NW_003871094.1":"HG729_PATCH","NW_003571053.1":"HG730_PATCH","NW_003315962.1":"HSCHR19_1_CTG3","NW_003315964.2":"HSCHR19_2_CTG3","NW_003315965.1":"HSCHR19_3_CTG3","NW_003315963.1":"HSCHR19_1_CTG3_1","NW_003571054.1":"HSCHR19LRC_COX1_CTG1","NW_003571055.1":"HSCHR19LRC_COX2_CTG1","NW_003571056.1":"HSCHR19LRC_LRC_I_CTG1","NW_003571057.1":"HSCHR19LRC_LRC_J_CTG1","NW_003571058.1":"HSCHR19LRC_LRC_S_CTG1","NW_003571059.1":"HSCHR19LRC_LRC_T_CTG1","NW_003571060.1":"HSCHR19LRC_PGF1_CTG1","NW_003571061.1":"HSCHR19LRC_PGF2_CTG1","NW_003315966.1":"HSCHR20_1_CTG1","NW_003871095.1":"HG144_PATCH","NW_003571063.2":"HG506_HG507_HG1000_PATCH","NW_003315967.1":"HSCHR21_1_CTG1_1","NW_003315968.1":"HSCHR21_2_CTG1_1","NW_003315969.1":"HSCHR21_3_CTG1_1","NW_003315970.1":"HSCHR21_4_CTG1_1","NW_004070874.1":"HG1487_PATCH","NW_004070873.1":"HG1486_PATCH","NW_004070875.1":"HG1488_PATCH","NW_003871096.1":"HG329_PATCH","NW_003315972.1":"HSCHR22_1_CTG2","NW_003315971.1":"HSCHR22_1_CTG1","NW_004070876.1":"HG497_PATCH","NW_003571064.1":"HG480_HG481_PATCH","NW_003871098.1":"HG1423_PATCH","NW_003871099.1":"HG1424_PATCH","NW_004070879.1":"HG1435_PATCH","NW_004070880.1":"HG1436_HG1432_PATCH","NW_004070877.1":"HG1433_PATCH","NW_004070881.1":"HG1437_PATCH","NW_004070882.1":"HG1438_PATCH","NW_003871100.1":"HG1425_PATCH","NW_003871101.2":"HG1426_PATCH","NW_004070883.1":"HG1439_PATCH","NW_004070884.1":"HG1440_PATCH","NW_004070885.1":"HG1441_PATCH","NW_003871102.1":"HG375_PATCH","NW_004070878.1":"HG1434_PATCH","NW_004070891.1":"HG1462_PATCH","NW_004070892.1":"HG1463_PATCH","NW_004070893.1":"HG1490_PATCH","NW_004070886.1":"HG1442_PATCH","NW_004070887.1":"HG1443_HG1444_PATCH","NW_004070888.1":"HG1453_PATCH","NW_004070889.1":"HG1458_PATCH","NW_004070890.1":"HG1459_PATCH","NW_003871103.2":"HG1497_PATCH"}
_GenBankAccnToChr={"CM000663.1":"1","CM000664.1":"2","CM000665.1":"3","CM000666.1":"4","CM000667.1":"5","CM000668.1":"6","CM000669.1":"7","CM000670.1":"8","CM000671.1":"9","CM000672.1":"10","CM000673.1":"11","CM000674.1":"12","CM000675.1":"13","CM000676.1":"14","CM000677.1":"15","CM000678.1":"16","CM000679.1":"17","CM000680.1":"18","CM000681.1":"19","CM000682.1":"20","CM000683.1":"21","CM000684.1":"22","CM000685.1":"X","CM000686.1":"Y","GL000191.1":"HSCHR1_RANDOM_CTG5","GL000192.1":"HSCHR1_RANDOM_CTG12","GL000193.1":"HSCHR4_RANDOM_CTG2","GL000194.1":"HSCHR4_RANDOM_CTG3","GL000195.1":"HSCHR7_RANDOM_CTG1","GL000196.1":"HSCHR8_RANDOM_CTG1","GL000197.1":"HSCHR8_RANDOM_CTG4","GL000198.1":"HSCHR9_RANDOM_CTG1","GL000199.1":"HSCHR9_RANDOM_CTG2","GL000200.1":"HSCHR9_RANDOM_CTG4","GL000201.1":"HSCHR9_RANDOM_CTG5","GL000202.1":"HSCHR11_RANDOM_CTG2","GL000203.1":"HSCHR17_RANDOM_CTG1","GL000204.1":"HSCHR17_RANDOM_CTG2","GL000205.1":"HSCHR17_RANDOM_CTG3","GL000206.1":"HSCHR17_RANDOM_CTG4","GL000207.1":"HSCHR18_RANDOM_CTG1","GL000208.1":"HSCHR19_RANDOM_CTG1","GL000209.1":"HSCHR19_RANDOM_CTG2","GL000210.1":"HSCHR21_RANDOM_CTG9","GL000211.1":"HSCHRUN_RANDOM_CTG1","GL000212.1":"HSCHRUN_RANDOM_CTG2","GL000213.1":"HSCHRUN_RANDOM_CTG3","GL000214.1":"HSCHRUN_RANDOM_CTG4","GL000215.1":"HSCHRUN_RANDOM_CTG5","GL000216.1":"HSCHRUN_RANDOM_CTG6","GL000217.1":"HSCHRUN_RANDOM_CTG7","GL000218.1":"HSCHRUN_RANDOM_CTG9","GL000219.1":"HSCHRUN_RANDOM_CTG10","GL000220.1":"HSCHRUN_RANDOM_CTG11","GL000221.1":"HSCHRUN_RANDOM_CTG13","GL000222.1":"HSCHRUN_RANDOM_CTG14","GL000223.1":"HSCHRUN_RANDOM_CTG15","GL000224.1":"HSCHRUN_RANDOM_CTG16","GL000225.1":"HSCHRUN_RANDOM_CTG17","GL000226.1":"HSCHRUN_RANDOM_CTG19","GL000227.1":"HSCHRUN_RANDOM_CTG20","GL000228.1":"HSCHRUN_RANDOM_CTG21","GL000229.1":"HSCHRUN_RANDOM_CTG22","GL000230.1":"HSCHRUN_RANDOM_CTG23","GL000231.1":"HSCHRUN_RANDOM_CTG24","GL000232.1":"HSCHRUN_RANDOM_CTG25","GL000233.1":"HSCHRUN_RANDOM_CTG26","GL000234.1":"HSCHRUN_RANDOM_CTG27","GL000235.1":"HSCHRUN_RANDOM_CTG28","GL000236.1":"HSCHRUN_RANDOM_CTG29","GL000237.1":"HSCHRUN_RANDOM_CTG30","GL000238.1":"HSCHRUN_RANDOM_CTG31","GL000239.1":"HSCHRUN_RANDOM_CTG32","GL000240.1":"HSCHRUN_RANDOM_CTG33","GL000241.1":"HSCHRUN_RANDOM_CTG34","GL000242.1":"HSCHRUN_RANDOM_CTG35","GL000243.1":"HSCHRUN_RANDOM_CTG36","GL000244.1":"HSCHRUN_RANDOM_CTG37","GL000245.1":"HSCHRUN_RANDOM_CTG38","GL000246.1":"HSCHRUN_RANDOM_CTG39","GL000247.1":"HSCHRUN_RANDOM_CTG40","GL000248.1":"HSCHRUN_RANDOM_CTG41","GL000249.1":"HSCHRUN_RANDOM_CTG42","GL000250.1":"HSCHR6_MHC_APD_CTG1","GL000251.1":"HSCHR6_MHC_COX_CTG1","GL000252.1":"HSCHR6_MHC_DBB_CTG1","GL000253.1":"HSCHR6_MHC_MANN_CTG1","GL000254.1":"HSCHR6_MHC_MCF_CTG1","GL000255.1":"HSCHR6_MHC_QBL_CTG1","GL000256.1":"HSCHR6_MHC_SSTO_CTG1","GL000257.1":"HSCHR4_1_CTG9","GL000258.1":"HSCHR17_1_CTG5","J01415.2":"MT","JH806574.1":"HG1472_PATCH","GL949741.1":"HG989_PATCH","JH636053.1":"HG1292_PATCH","JH636052.3":"HG1287_PATCH","GL383518.1":"HSCHR1_1_CTG31","GL383519.1":"HSCHR1_2_CTG31","GL383520.1":"HSCHR1_3_CTG31","JH806573.1":"HG1471_PATCH","JH636054.1":"HG1293_PATCH","JH806575.1":"HG1473_PATCH","GL383516.1":"HG999_1_PATCH","GL383517.1":"HG999_2_PATCH","GL383521.1":"HSCHR2_1_CTG1","GL877871.1":"HG686_PATCH","GL582966.2":"HSCHR2_2_CTG12","GL383522.1":"HSCHR2_1_CTG12","GL877870.2":"HG1007_PATCH","JH636055.1":"HSCHR3_1_CTG1","JH159132.1":"HG325_PATCH","GL383523.1":"HG186_PATCH","GL383524.1":"HG280_PATCH","JH159131.1":"HG1091_PATCH","GL383525.1":"HG991_PATCH","GL383526.1":"HSCHR3_1_CTG2_1","GL383528.1":"HSCHR4_1_CTG6","GL383529.1":"HSCHR4_2_CTG9","GL582967.1":"HG706_PATCH","GL383527.1":"HSCHR4_1_CTG12","GL877872.1":"HG1032_PATCH","GL383532.1":"HSCHR5_2_CTG1","GL949742.1":"HSCHR5_3_CTG1","GL339449.2":"HSCHR5_1_CTG1","GL383530.1":"HSCHR5_1_CTG2","JH159133.1":"HG1063_PATCH","GL383531.1":"HSCHR5_1_CTG5","JH806576.1":"HG27_PATCH","JH636057.1":"HG1322_PATCH","GL383533.1":"HSCHR6_1_CTG5","JH636056.1":"HG1304_PATCH","GL582970.1":"HG736_PATCH","GL582969.1":"HG14_PATCH","JH159134.2":"HG1257_PATCH","GL582972.1":"HG946_PATCH","GL582968.1":"HG115_PATCH","JH636058.1":"HG1308_PATCH","GL383534.2":"HSCHR7_1_CTG6","GL582971.1":"HG7_PATCH","GL949743.1":"HG19_PATCH","JH159135.1":"HG418_PATCH","GL383535.1":"HG104_HG975_PATCH","GL383536.1":"HG243_PATCH","GL383539.1":"HSCHR9_1_CTG1","JH636059.1":"HG962_PATCH","GL383540.1":"HSCHR9_1_CTG35","GL383541.1":"HSCHR9_2_CTG35","GL383542.1":"HSCHR9_3_CTG35","JH806579.1":"HG1502_PATCH","GL339450.1":"HG79_PATCH","JH806577.1":"HG1500_PATCH","JH806578.1":"HG1501_PATCH","GL383537.1":"HG998_1_PATCH","GL383538.1":"HG998_2_PATCH","GL877873.1":"HG905_PATCH","JH636060.1":"HG871_PATCH","GL383543.1":"HG544_PATCH","GL383545.1":"HSCHR10_1_CTG2","GL383546.1":"HSCHR10_1_CTG5","JH591181.2":"HG1211_PATCH","JH591183.1":"HG339_PATCH","JH591182.1":"HG311_PATCH","GL383544.1":"HG995_PATCH","JH806580.1":"HG1479_PATCH","JH591184.1":"HG256_PATCH","JH591185.1":"HG873_PATCH","GL383547.1":"HSCHR11_1_CTG1_1","GL582973.1":"HG281_PATCH","JH159136.1":"HG142_HG150_NOVEL_TEST","JH159137.1":"HG151_NOVEL_TEST","GL949744.1":"HG536_PATCH","JH806581.1":"HG865_PATCH","JH159143.1":"HG414_PATCH","JH159141.2":"HG348_PATCH","JH159139.1":"HG305_PATCH","JH159142.2":"HG388_HG400_PATCH","JH159140.1":"HG306_PATCH","JH720443.1":"HG122_PATCH","JH159138.1":"HG299_PATCH","GL582974.1":"HG858_PATCH","GL877875.1":"HSCHR12_1_CTG1","JH720444.1":"HG344_PATCH","GL949745.1":"HG1133_PATCH","GL877876.1":"HSCHR12_2_CTG2","GL383549.1":"HSCHR12_1_CTG2","GL383550.1":"HSCHR12_1_CTG2_1","GL383552.1":"HSCHR12_2_CTG2_1","GL383553.2":"HSCHR12_3_CTG2_1","GL383551.1":"HSCHR12_1_CTG5","GL383548.1":"HG996_PATCH","GL582975.1":"HG531_PATCH","GL383554.1":"HSCHR15_1_CTG4","GL383555.1":"HSCHR15_1_CTG8","JH720445.1":"HG971_PATCH","GL383556.1":"HSCHR16_1_CTG3_1","JH720446.1":"HG1208_PATCH","GL383557.1":"HSCHR16_2_CTG3_1","JH806582.1":"HG417_PATCH","GL383563.1":"HSCHR17_1_CTG1","GL383562.1":"HG990_PATCH","GL383561.1":"HG987_PATCH","JH159145.1":"HG883_PATCH","GL383560.1":"HG75_PATCH","GL383559.1":"HG745_PATCH","JH159146.1":"HSCHR17_4_CTG4","JH159148.1":"HSCHR17_6_CTG4","JH159147.1":"HSCHR17_5_CTG4","GL383564.1":"HSCHR17_1_CTG4","GL582976.1":"HG185_PATCH","JH720447.1":"HG1146_PATCH","GL383558.1":"HG183_PATCH","JH159144.1":"HG747_PATCH","GL383565.1":"HSCHR17_2_CTG4","GL383566.1":"HSCHR17_3_CTG4","JH591186.1":"HG748_PATCH","JH636061.1":"HG271_PATCH","GL383567.1":"HSCHR18_1_CTG1_1","GL383570.1":"HSCHR18_2_CTG1_1","GL383571.1":"HSCHR18_2_CTG2","GL383568.1":"HSCHR18_1_CTG2","GL383569.1":"HSCHR18_1_CTG2_1","GL383572.1":"HSCHR18_2_CTG2_1","JH159149.1":"HG729_PATCH","GL582977.1":"HG730_PATCH","GL383573.1":"HSCHR19_1_CTG3","GL383575.2":"HSCHR19_2_CTG3","GL383576.1":"HSCHR19_3_CTG3","GL383574.1":"HSCHR19_1_CTG3_1","GL949746.1":"HSCHR19LRC_COX1_CTG1","GL949747.1":"HSCHR19LRC_COX2_CTG1","GL949748.1":"HSCHR19LRC_LRC_I_CTG1","GL949749.1":"HSCHR19LRC_LRC_J_CTG1","GL949750.1":"HSCHR19LRC_LRC_S_CTG1","GL949751.1":"HSCHR19LRC_LRC_T_CTG1","GL949752.1":"HSCHR19LRC_PGF1_CTG1","GL949753.1":"HSCHR19LRC_PGF2_CTG1","GL383577.1":"HSCHR20_1_CTG1","JH720448.1":"HG144_PATCH","GL582979.2":"HG506_HG507_HG1000_PATCH","GL383578.1":"HSCHR21_1_CTG1_1","GL383579.1":"HSCHR21_2_CTG1_1","GL383580.1":"HSCHR21_3_CTG1_1","GL383581.1":"HSCHR21_4_CTG1_1","JH806584.1":"HG1487_PATCH","JH806583.1":"HG1486_PATCH","JH806585.1":"HG1488_PATCH","JH720449.1":"HG329_PATCH","GL383583.1":"HSCHR22_1_CTG2","GL383582.1":"HSCHR22_1_CTG1","JH806586.1":"HG497_PATCH","GL877877.1":"HG480_HG481_PATCH","JH720451.1":"HG1423_PATCH","JH720452.1":"HG1424_PATCH","JH806589.1":"HG1435_PATCH","JH806590.1":"HG1436_HG1432_PATCH","JH806587.1":"HG1433_PATCH","JH806591.1":"HG1437_PATCH","JH806592.1":"HG1438_PATCH","JH720453.1":"HG1425_PATCH","JH720454.2":"HG1426_PATCH","JH806593.1":"HG1439_PATCH","JH806594.1":"HG1440_PATCH","JH806595.1":"HG1441_PATCH","JH720455.1":"HG375_PATCH","JH806588.1":"HG1434_PATCH","JH806601.1":"HG1462_PATCH","JH806602.1":"HG1463_PATCH","JH806603.1":"HG1490_PATCH","JH806596.1":"HG1442_PATCH","JH806597.1":"HG1443_HG1444_PATCH","JH806598.1":"HG1453_PATCH","JH806599.1":"HG1458_PATCH","JH806600.1":"HG1459_PATCH","JH159150.2":"HG1497_PATCH"}
_UCSCAccnToChr={"chr1":"1","chr2":"2","chr3":"3","chr4":"4","chr5":"5","chr6":"6","chr7":"7","chr8":"8","chr9":"9","chr10":"10","chr11":"11","chr12":"12","chr13":"13","chr14":"14","chr15":"15","chr16":"16","chr17":"17","chr18":"18","chr19":"19","chr20":"20","chr21":"21","chr22":"22","chrX":"X","chrY":"Y","chrM":"M","chr17_ctg5_hap1":"HSCHR17_1_CTG5","chr1_gl000191_random":"HSCHR1_RANDOM_CTG5","chr1_191_random":"HSCHR1_RANDOM_CTG5","chr1_gl000192_random":"HSCHR1_RANDOM_CTG12","chr1_192_random":"HSCHR1_RANDOM_CTG12","chr4_ctg9_hap1":"HSCHR4_1_CTG9","chr4_gl000193_random":"HSCHR4_RANDOM_CTG2","chr4_193_random":"HSCHR4_RANDOM_CTG2","chr6_apd_hap1":"HSCHR6_MHC_APD_CTG1","chr6_cox_hap2":"HSCHR6_MHC_COX_CTG1","chr6_dbb_hap3":"HSCHR6_MHC_DBB_CTG1","chr6_mann_hap4":"HSCHR6_MHC_MANN_CTG1","chr6_mcf_hap5":"HSCHR6_MHC_MCF_CTG1","chr6_qbl_hap6":"HSCHR6_MHC_QBL_CTG1","chr6_ssto_hap7":"HSCHR6_MHC_SSTO_CTG1","chr7_gl000195_random":"HSCHR7_RANDOM_CTG1","chr7_195_random":"HSCHR7_RANDOM_CTG1","chrUn_gl000212":"HSCHRUN_RANDOM_CTG2","chrUn_212":"HSCHRUN_RANDOM_CTG2","chrUn_gl000214":"HSCHRUN_RANDOM_CTG4","chrUn_214":"HSCHRUN_RANDOM_CTG4","chrUn_gl000219":"HSCHRUN_RANDOM_CTG10","chrUn_219":"HSCHRUN_RANDOM_CTG10","chrUn_gl000220":"HSCHRUN_RANDOM_CTG11","chrUn_220":"HSCHRUN_RANDOM_CTG11","chrUn_gl000221":"HSCHRUN_RANDOM_CTG13","chrUn_221":"HSCHRUN_RANDOM_CTG13","chrUn_gl000222":"HSCHRUN_RANDOM_CTG14","chrUn_222":"HSCHRUN_RANDOM_CTG14","chrUn_gl000223":"HSCHRUN_RANDOM_CTG15","chrUn_223":"HSCHRUN_RANDOM_CTG15","chrUn_gl000228":"HSCHRUN_RANDOM_CTG21","chrUn_228":"HSCHRUN_RANDOM_CTG21","chrUn_gl000241":"HSCHRUN_RANDOM_CTG34","chrUn_241":"HSCHRUN_RANDOM_CTG34","chr17_gl000204_random":"HSCHR17_RANDOM_CTG2","chr17_204_random":"HSCHR17_RANDOM_CTG2","chr17_gl000205_random":"HSCHR17_RANDOM_CTG3","chr17_205_random":"HSCHR17_RANDOM_CTG3","chr19_gl000209_random":"HSCHR19_RANDOM_CTG2","chr19_209_random":"HSCHR19_RANDOM_CTG2","chr4_gl000194_random":"HSCHR4_RANDOM_CTG3","chr4_194_random":"HSCHR4_RANDOM_CTG3","chr9_gl000201_random":"HSCHR9_RANDOM_CTG5","chr9_201_random":"HSCHR9_RANDOM_CTG5","chrUn_gl000211":"HSCHRUN_RANDOM_CTG1","chrUn_211":"HSCHRUN_RANDOM_CTG1","chrUn_gl000213":"HSCHRUN_RANDOM_CTG3","chrUn_213":"HSCHRUN_RANDOM_CTG3","chrUn_gl000218":"HSCHRUN_RANDOM_CTG9","chrUn_218":"HSCHRUN_RANDOM_CTG9","chrUn_gl000227":"HSCHRUN_RANDOM_CTG20","chrUn_227":"HSCHRUN_RANDOM_CTG20","chrUn_gl000229":"HSCHRUN_RANDOM_CTG22","chrUn_229":"HSCHRUN_RANDOM_CTG22","chrUn_gl000237":"HSCHRUN_RANDOM_CTG30","chrUn_237":"HSCHRUN_RANDOM_CTG30","chrUn_gl000243":"HSCHRUN_RANDOM_CTG36","chrUn_243":"HSCHRUN_RANDOM_CTG36","chrUn_gl000247":"HSCHRUN_RANDOM_CTG40","chrUn_247":"HSCHRUN_RANDOM_CTG40","chrUn_215":"HSCHRUN_RANDOM_CTG5","chrUn_224":"HSCHRUN_RANDOM_CTG16","chrUn_238":"HSCHRUN_RANDOM_CTG31","chrUn_242":"HSCHRUN_RANDOM_CTG35","chrUn_244":"HSCHRUN_RANDOM_CTG37","chrUn_249":"HSCHRUN_RANDOM_CTG42"}
_ChrToUCSCAccn={"10":"chr10","11":"chr11","12":"chr12","13":"chr13","14":"chr14","15":"chr15","16":"chr16","17":"chr17","18":"chr18","19":"chr19","1":"chr1","20":"chr20","21":"chr21","22":"chr22","2":"chr2","3":"chr3","4":"chr4","5":"chr5","6":"chr6","7":"chr7","8":"chr8","9":"chr9","HSCHRUN_RANDOM_CTG17":"chrUn_gl000225","HSCHR17_1_CTG5":"chr17_ctg5_hap1","HSCHR17_RANDOM_CTG2":"chr17_gl000204_random","HSCHR17_RANDOM_CTG3":"chr17_gl000205_random","HSCHR19_RANDOM_CTG2":"chr19_gl000209_random","HSCHR1_RANDOM_CTG12":"chr1_gl000192_random","HSCHR1_RANDOM_CTG5":"chr1_gl000191_random","HSCHR4_1_CTG9":"chr4_ctg9_hap1","HSCHR4_RANDOM_CTG2":"chr4_gl000193_random","HSCHR4_RANDOM_CTG3":"chr4_gl000194_random","HSCHR6_MHC_APD_CTG1":"chr6_apd_hap1","HSCHR6_MHC_COX_CTG1":"chr6_cox_hap2","HSCHR6_MHC_DBB_CTG1":"chr6_dbb_hap3","HSCHR6_MHC_MANN_CTG1":"chr6_mann_hap4","HSCHR6_MHC_MCF_CTG1":"chr6_mcf_hap5","HSCHR6_MHC_QBL_CTG1":"chr6_qbl_hap6","HSCHR6_MHC_SSTO_CTG1":"chr6_ssto_hap7","HSCHR7_RANDOM_CTG1":"chr7_gl000195_random","HSCHR9_RANDOM_CTG5":"chr9_gl000201_random","HSCHRUN_RANDOM_CTG10":"chrUn_gl000219","HSCHRUN_RANDOM_CTG11":"chrUn_gl000220","HSCHRUN_RANDOM_CTG13":"chrUn_gl000221","HSCHRUN_RANDOM_CTG14":"chrUn_gl000222","HSCHRUN_RANDOM_CTG15":"chrUn_gl000223","HSCHRUN_RANDOM_CTG16":"chrUn_gl000224","HSCHRUN_RANDOM_CTG1":"chrUn_gl000211","HSCHRUN_RANDOM_CTG20":"chrUn_gl000227","HSCHRUN_RANDOM_CTG21":"chrUn_gl000228","HSCHRUN_RANDOM_CTG22":"chrUn_gl000229","HSCHRUN_RANDOM_CTG2":"chrUn_gl000212","HSCHRUN_RANDOM_CTG30":"chrUn_gl000237","HSCHRUN_RANDOM_CTG31":"chrUn_gl000238","HSCHRUN_RANDOM_CTG34":"chrUn_gl000241","HSCHRUN_RANDOM_CTG35":"chrUn_gl000242","HSCHRUN_RANDOM_CTG36":"chrUn_gl000243","HSCHRUN_RANDOM_CTG37":"chrUn_gl000244","HSCHRUN_RANDOM_CTG3":"chrUn_gl000213","HSCHRUN_RANDOM_CTG40":"chrUn_gl000247","HSCHRUN_RANDOM_CTG42":"chrUn_gl000249","HSCHRUN_RANDOM_CTG4":"chrUn_gl000214","HSCHRUN_RANDOM_CTG5":"chrUn_gl000215","HSCHRUN_RANDOM_CTG9":"chrUn_gl000218","HSCHRUN_RANDOM_CTG6":"chrUn_gl000216","HSCHR9_RANDOM_CTG4":"chr9_gl000200_random","HSCHRUN_RANDOM_CTG7":"chrUn_gl000217","HSCHR9_RANDOM_CTG2":"chr9_gl000199_random","HSCHR19_RANDOM_CTG1":"chr19_gl000208_random","HSCHR9_RANDOM_CTG1":"chr9_gl000198_random","HSCHRUN_RANDOM_CTG26":"chrUn_gl000233","HSCHRUN_RANDOM_CTG23":"chrUn_gl000230","HSCHRUN_RANDOM_CTG29":"chrUn_gl000236","HSCHRUN_RANDOM_CTG33":"chrUn_gl000240","HSCHR17_RANDOM_CTG4":"chr17_gl000206_random","HSCHRUN_RANDOM_CTG25":"chrUn_gl000232","HSCHRUN_RANDOM_CTG27":"chrUn_gl000234","HSCHR11_RANDOM_CTG2":"chr11_gl000202_random","HSCHRUN_RANDOM_CTG41":"chrUn_gl000248","HSCHR8_RANDOM_CTG1":"chr8_gl000196_random","HSCHRUN_RANDOM_CTG39":"chrUn_gl000246","HSCHR17_RANDOM_CTG1":"chr17_gl000203_random","HSCHR8_RANDOM_CTG4":"chr8_gl000197_random","HSCHRUN_RANDOM_CTG38":"chrUn_gl000245","HSCHRUN_RANDOM_CTG28":"chrUn_gl000235","HSCHRUN_RANDOM_CTG32":"chrUn_gl000239","HSCHR21_RANDOM_CTG9":"chr21_gl000210_random","HSCHRUN_RANDOM_CTG24":"chrUn_gl000231","HSCHRUN_RANDOM_CTG19":"chrUn_gl000226","HSCHR18_RANDOM_CTG1":"chr18_gl000207_random","HG1436_HG1432_PATCH":"chrUn_jh806590","HG1438_PATCH":"chrUn_jh806592","HG1439_PATCH":"chrUn_jh806593","HG1441_PATCH":"chrUn_jh806595","HG1459_PATCH":"chrUn_jh806600","":"","M":"chrM","X":"chrX","Y":"chrY"} # bins range in size from 16kb to 512Mb 
# Bin  0          spans 512Mbp,   # Level 1
# Bins 1-8        span 64Mbp,     # Level 2
# Bins 9-72       span 8Mbp,      # Level 3
# Bins 73-584     span 1Mbp       # Level 4
# Bins 585-4680   span 128Kbp     # Level 5
# Bins 4681-37449 span 16Kbp      # Level 6
_binOffsetsExtended = [32678+4096+512+64+8+1, 4096+512+64+8+1, 512+64+8+1, 64+8+1, 8+1, 1, 0]
_binFirstShift = 14;      # How much to shift to get to finest bin. 
_binNextShift  = 3;       # How much to shift to get to next larger bin.

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------


class Bed:
    '''This class is used to process BED format lines. The default BED format is at least three fields.'''
    def __init__(self,x,description=None):
        '''Initiate the bed from either line or list.'''
        try:
            x=x.rstrip("\n\r").split("\t")
        except:
            if isinstance(x[-1],basestring):
                x[-1].rstrip("\n\r")
        self.chr=x[0].strip()
        self.start=int(x[1])
        if self.start<0:
            self.start=0
        self.stop=int(x[2])
        try:
            self.id=x[3]
        except:
            self.id="NONAME"
        try:
            self.score=float(x[4])
        except:
            self.score=1
        try:
            self.strand=x[5]
        except:
            self.strand="."
        try:
            self.otherfields=x[6:]
        except:
            self.otherfields=[]
        self.description=description
    def __str__(self):
        '''Return the bed in basestring format.'''
        string=self.chr+"\t"+str(self.start)+"\t"+str(self.stop)+"\t"+str(self.id)+"\t"+("%-5.2f\t"% self.score)+self.strand
        return string
    def __add__(A,B):
        '''Add Bed A and B together. Please test isOverlap for bed merge!'''
        if not A and not B: return None
        if not A: return B   #B+0=B
        if not B: return A   #A+0=A
        if A.chr==B.chr: #A+B
            start=min(A.start,B.start)
            stop=max(A.stop,B.stop)
            return Bed([A.chr,start,stop,A.id+","+B.id,(A.score*A.length()+B.score*B.length())/(stop-start),A.strand if A.strand==B.strand else "."])
        return A # A and B are not in the same chromosome.
    def __cmp__(self,other):
        return cmp(self.chr,other.chr) or cmp(self.start,other.start) or cmp(self.stop,other.stop) or cmp(other.strand,self.strand) or cmp(self.score,other.score)#'.'<'-'<'+'
    def length(self):
        '''Return the bed range.'''
        return self.stop-self.start
    def isOverlap(self,B):
        '''Return a bool value of the status of if two bed are overlapped.'''
        if not B: return False
        if(self.chr != B.chr) : return False
        if (self.stop <= B.start) : return False
        if (B.stop <= self.start) : return False
        return True
    def overlapLength(self,B):
        '''Return the overlapped length. >0: overlap; =0: neighbour; <0: -distance; -1000000000: not in the same chromosome.'''
        if not B: return -1000000000
        if self.chr==B.chr:
            return self.length()+B.length()-max(self.stop,B.stop)+min(self.start,B.start)
        return -1000000000
    def testBoundary(self):
        '''Swap the start and stop position and reverse the strand if start is bigger than stop.'''
        if self.start>self.stop:
            print >>sys.stderr, "start is bigger than stop, swap them and reverse the strand."
            self.stop,self.start=self.start,self.stop
            if self.strand=="+": 
                self.strand="-"
            elif self.strand=="-":
                self.strand="+"        
    def distance(self,B):
        '''Get the distance between two Beds.'''
        return -self.overlapLenght(B)
        #if(self.chr != bed.chr): return 1000000000 #In different chromosome
        #if(self.isOverlap(bed)): return -1 #Overlapped
        #return min(abs(self.start-bed.stop),abs(self.stop-bed.start))
    def getSeq(self,fn="/home/cxw/Genome/hg19/hg19.2bit"):
        '''Get fasta sequence from 2bit file'''
        if(fn is None): 
            print >>sys.stderr,"2Bit file not specified."
            return ""
        seq=getSeq(fn,self.chr,self.start,self.stop) #Default is lowercase ???
        seq=seq.upper()
        if "-" in self.strand:
            return Utils.rc(seq)
        return seq
    def getWig(self,fn):
        '''get base value from bigWig file.'''
        return getDepthList(fn,self.chr,self.start,self.stop)
    def updownExtend(self,up=0,down=0):
        '''Return the a new Bed with upstream up and downstream down'''
        if self.strand=="-":
            start=self.start-down
            stop=self.stop+up
        else:
            start=self.start-up
            stop=self.stop+down
        tbed=Bed([self.chr,start,stop,self.id+("_up"+str(up) if up!=0 else "")+("_down"+str(down) if down!=0 else ""),self.score,self.strand])
        tbed.testBoundary()
        return tbed
    def stringGraph(self,scale=1):
        '''Illustrate the bed in graph mode.'''
        n=int(self.length()*scale)
        if n==0: n=1
        if self.strand == "+": return ">"*n
        if self.strand == "-": return "<"*n
        return "|"*n
    def setDepth(self,cover):
        '''Add depth for each base.'''
        self.depth=cover
    def strandCmp(self,bed):
        '''Test if the same strand.'''
        if bed.strand == "." or self.strand==".": return "."
        return "+" if self.strand == bed.strand else "-"
    def getBIN(self):
        '''Get the genome BIN.'''
        start=self.start>>_binFirstShift
        stop=(self.stop-1)>>_binFirstShift
        for i in range(_binLevels):
            if start==stop:
                return _binOffsetsExtended[i] + start
            else:
                start >>= _binNextShift
                stop  >>= _binNextShift
        assert "Error! Bed range is out of 512M."
    def toTSS(self):
        '''Change to TSS.'''
        return TSS([self.chr,self.stop-1 if self.strand=='-' else self.start,self.id,self.score,self.strand])
    def toSite(self,stype='TSS'):
        '''Change to Site.'''
        if stype=='TSS':
            return Site([self.chr,self.stop-1 if self.strand=='-' else self.start,self.id,self.score,self.strand],stype)
        if stype=='EnzymeDigest':
            return Site([self.chr,self.stop if self.strand=='-' else self.start,self.id,self.score],stype)
        if stype=='TIS':
            pass
        if stype=='TTS':
            pass
        return None
        
class BedCoverage(Bed):
    '''Bed6 format with coverage infomation.'''
    def __init__(self,x,description=None):
        '''Initiation'''
        Bed.__init__(self,x,description)
        if isinstance(self.otherfields[0],str):
            self.depth=[float(i) for i in self.otherfields[0].split(',')]
            self.otherfields=self.otherfields[1:]
        elif isinstance(self.otherfields[0],list):
            self.depth=self.otherfields[0]
            self.otherfields=self.otherfields[1:]
        else:
            self.depth=[self.score for i in xrange(self.length())]

class BedList(list):
    '''List class for hold line objects.'''
    def __init__(self,data=[],description=None):
        '''Initiate data by list class.'''
        list.__init__(self,data)
        self.sorted=0
        self.description=description
    def readfile(self,infile,format='bed'):
        '''Read data fromfile by  ColumnReader.'''
        for item in ColumnReader(infile,format):
            self.append(item)
    def sort(self):
        '''sort BedList.'''
        list.sort(self)
        self.sorted=1
    def bisect(self,item):
        '''Find the nearest item for comparation.'''
        if not self.sorted:
            self.sort()
        return bisect(self,item)
    def clear(self):
        '''Clear List.'''
        del self[:]
        self.sorted=0
    def mergeSort(bedfiles): #generator
        '''Merge Beds from multiple files. The Bed in each file should be sorted.'''
        print >>sys.stderr, "Make sure each file input is sorted!"
        if isinstance(bedfiles,str):
            bedfiles=[bedfiles]
        if len(bedfiles)==0:
            print >>sys.stderr, "No bed file names provided...."
        elif len(bedfiles)==1: # for single file
            for tbed in ColumnReader(bedfiles[0],'bed'):
                yield tbed
        else: #for multiple file
            cr=[] #List for iteraters
            bedfiledict={} #record bedfile index for iteration
            beds=BedList()
            for index,bedfile in enumerate(bedfiles):
                bedfiledict[bedfile]=index
                cr.append(ColumnReader(bedfile,'bed'))
                try:
                    beds.append(cr[index].next())
                    beds[index].description=bedfile #record the source file number for iteration
                except Exception,e:
                    print >>sys.stderr,"Read Bed error!",e
                    raise
            beds.sort()
            while True:
                if len(beds)>0:
                    yield beds[0] # yield the minimum bed
                    try:
                        tbed=cr[bedfiledict[beds[0].description]].next()
                        tbed.description=beds[0].description
                        beds.insert(beds.bisect(tbed),tbed) # insert(pos,item) insert tbed in  the right position
                    except StopIteration:
                        print >>sys.stderr,bedfiles[bedfiledict[beds[0].description]]+" is finished..."
                    del beds[0]                
                else:
                    break
    mergeSort=staticmethod(mergeSort)
    def mergeBeds(bedfiles,forcestrand=False): #generator
        '''Merge Beds from multiple files. The overlapped beds are combined. The Bed in each file should be sorted.'''
        beds=[None,None]
        bedcount=[0,0]
        bedindex={".":0,"+":0,"-":(1 if forcestrand else 0)}
        for tbed in BedList.mergeSort(bedfiles):
            if tbed.isOverlap(beds[bedindex[tbed.strand]]):
                beds[bedindex[tbed.strand]]+=tbed
            else:
                if beds[bedindex[tbed.strand]]:
                    bedcount[bedindex[tbed.strand]]+=1
                    beds[bedindex[tbed.strand]].description=beds[bedindex[tbed.strand]].id
                    beds[bedindex[tbed.strand]].id="Region_"+str(bedcount[bedindex[tbed.strand]])
                    yield beds[bedindex[tbed.strand]]
                beds[bedindex[tbed.strand]]=tbed
        for tbed in beds:
            if tbed:
                bedcount[bedindex[tbed.strand]]+=1
                tbed.description="Region_"+str(bedcount[bedindex[tbed.strand]])
                tbed.id,tbed.description=tbed.description,tbed.id
                yield tbed
        assert "Reach this line."
    mergeBeds=staticmethod(mergeBeds)
            
class GeneBed(Bed):
    '''UCSC GenePred format.'''
    def __init__(self,x):
        '''Initiate from GeneBed lines. GeneBed names column are not allowed to be numbers.'''
        try:
            self.bin=int(x[0])
            if self.bin<10000:
                x=x[1:]
        except:
            pass
        self.id=x[0]
        if x[1] in _ChrToUCSCAccn.keys():
            self.chr=_ChrToUCSCAccn[x[1]]
        else:
            self.chr=x[1]
        self.strand=x[2]
        self.start=int(x[3])
        self.stop=int(x[4])
        self.txstart=int(x[5])
        self.txstop=int(x[6])
        self.exoncount=int(x[7])
        if isinstance(x[8],basestring):
            self.exonstarts=[int(p) for p in x[8].split(",")[0:-1]]
            self.exonstops=[int(p) for p in x[9].split(",")[0:-1]]
        else:
            self.exonstarts=[int(p) for p in x[8]]
            self.exonstops=[int(p) for p in x[9]]                    
        self.score=int(x[10])
        self.protein_id=x[11]
        self.cdsstartstat=x[12]
        self.cdsendstat=x[13]
        self.exonframes=x[14]
        try:
            self.description=x[15]
        except:
            self.description=""
        
    def __str__(self):
        '''Return GeneBed line.'''
        return "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s,\t%s,\t%d\t%s\t%s\t%s\t%s\t%s" % (self.id,self.chr,self.strand,self.start,self.stop,self.txstart,self.txstop,self.exoncount,",".join([str(p) for p in self.exonstarts]),",".join([str(p) for p in self.exonstops]),self.score,self.protein_id,self.cdsstartstat,self.cdsendstat,self.exonframes,self.description)
    def toBed(self):
        '''Transform to Bed format.'''
        return Bed([self.chr,self.start,self.stop,self.id,self.score,self.strand])
    def getExon(self,i):
        '''Get the Bed format of ith exon.'''
        if i>self.exoncount or i<1:
            return None
        p= (self.exoncount-i if self.strand=='-' else i-1)
        return Bed([self.chr,self.exonstarts[p],self.exonstops[p],self.id+":exon_"+str(i),0,self.strand])
    def getIntron(self,i):
        '''Get the Bed format of ith intron.'''
        if i>0 and i<self.exoncount:
            p=(self.exoncount-i if self.strand=='-' else i)
            #Notice: GeneBed is 1 based and Bed is 0 based.
            return Bed([self.chr,self.exonstops[p-1],self.exonstarts[p],self.id+":intron_"+str(i),0,self.strand])
        else:
            return None
    def _getUTRs(self,end='left'):
        '''Get the UTR.'''
        if end=='left': #Get the UTR located in the left end of chromosome
            if self.start==self.txstart: #No UTR
                return None
            utr=copy.deepcopy(self)
            for i in range(self.exoncount):
                if self.txstart<=self.exonstops[i]: #the txStart site locates in (i+1)th exon.
                    break
            utr.exonstarts=utr.exonstarts[0:i+1]
            utr.exonstops=utr.exonstops[0:i+1]
            utr.stop=utr.txstop=utr.exonstops[i]=self.txstart
            utr.txstart=self.start
            utr.exoncount=i+1
            return utr
        else: #get the UTR located in the right end of chromosome
            if self.stop==self.txstop: #No UTR
                return None
            utr=copy.deepcopy(self)
            for i in range(self.exoncount-1,-1,-1):
                if self.txstop>=self.exonstarts[i]:
                    break
            utr.exonstarts=utr.exonstarts[i:self.exoncount]
            utr.exonstops=utr.exonstops[i:self.exoncount]
            utr.start=utr.txstart=utr.exonstarts[0]=self.txstop
            utr.txstop=self.stop
            utr.exoncount=self.exoncount-i
            return utr
    def getUTR5(self):
        '''Get the 5UTR.'''
        if self.strand=='-':
            utr5=self._getUTRs('right')
        else:
            utr5=self._getUTRs('left')
        if utr5:
            utr5.id+=':UTR5'
        return utr5
    def getUTR3(self):
        '''Get the 3'UTR.'''
        if self.strand=='-':
            utr3=self._getUTRs('left')
        else:
            utr3=self._getUTRs('right')
        if utr3:
            utr3.id+=':UTR3'
        return utr3
    def overlapLength(self,B):
        '''Return the overlap length between genebed and bed. Only exons are considered.'''
        l=0
        for exon in self.exons():
            l+=exon.overlapLength(B)
        return l
    def exons(self):
        '''Iterate all exons.'''
        for i in range(1,self.exoncount+1):
            yield self.getExon(i)
    def introns(self):
        '''Iterate all introns.'''
        for i in range(1,self.exoncount):
            yield self.getIntron(i)
    def getcDNALength(self):
        '''Return cdna length.'''
        l=0
        cdsbed=Bed([self.chr,self.txstart,self.txstop])
        for exon in self.exons():
            l+=exon.length()
        return l
    def getcDNASeq(self,fn="/disk/Genome/hg19/hg19.2bit"):
        '''Get cDNA Sequence.'''
        seq=""
        for i in range(self.exoncount):
            seq+=self.getExon(i+1).getSeq(fn)
        return seq
    def getCDSLength(self):
        '''Return CDS length.'''
        l=0
        for exon in self.exons():
            l+=cdsbed.OverlapLength(exon)
        return l
    def getCDSSeq(self,fn="/disk/Genome/hg19/hg19.2bit"):
        '''get CDS sequence.'''
        seq=""
        cdsbed=Bed([self.chr,self.txstart,self.txstop])
        for exon in self.exons():
            if cdsbed.isOverlap(exon):
                seq+=Bed([self.chr,max(self.txstart,exon.start),min(self.txstop,exon.stop),'',0,self.strand]).getSeq(fn)
        return seq
    def getCDSGenoPos(self):
        '''Reture Genome Position for each base in CDS'''
        pos=[]
        cdsbed=Bed([self.chr,self.txstart,self.txstop])
        for exon in self.exons():
            if cdsbed.isOverlap(exon):
                poslist=range(max(self.txstart,exon.start)+1,min(self.txstop,exon.stop)+1)
                if "-" in self.strand:
                    poslist.reverse()
                    pos+=poslist
                else:
                    pos+=poslist
        return pos
    def getExonGenoPos(self):
        '''Reture Genome Position for each base in Exon'''
        pos=[]
        for exon in self.exons():
            poslist=range(exon.start+1,exon.stop+1)
            if "-" in self.strand:
                poslist.reverse()
                pos+=poslist
            else:
                pos+=poslist
        return pos
    def getGeneLength(self):
        '''Return gene length.'''
        return self.stop-self.start+1

class GeneBedList(BedList): #useless
    '''A list for Genes.'''
    def __init__(self,data=[],description=None):
        BedList.__init__(self,data,description)
        
class Site:
    '''site.'''
    def __init__(self,x,stype=None,description=None):
        '''Initiation'''
        self.chr=x[0]
        self.pos=int(x[1])
        try:
            self.id=x[2]
        except:
            self.id="NONAME"
        try:
            self.score=float(x[3])
        except:
            self.score=1.0
        try:
            self.strand=x[4]
        except:
            self.strand="."
        try:
            self.otherfields=x[5:]
        except:
            self.otherfileds=[]
        self.type=stype
        self.description=description
    def __str__(self):
        '''Return string of Site.'''
        return "%s\t%d\t%s\t%-5.3f\t%s\t%s" %(self.chr,self.pos,self.id,self.score,self.strand,self.type if self.type else '')
    def __cmp__(self,site):
        '''Compare self and other Site.'''
        return cmp(self.chr,site.chr) or cmp(self.pos,site.pos) or cmp(site.strand,self.strand) or cmp(self.score,site.score)

class SiteList(list):
    '''A list of Sites.'''
    def __init__(self,data=[],description=None):
        '''Initiate data by list class.'''
        list.__init__(self,data)
        self.sorted=0
        self.description=description
    def readfile(self,infile):
        '''Read data fromfile by  ColumnReader.'''
        for line in ColumnReader(infile):
            self.append(Site(item))
    def mergeSites(self):
        '''Merge sites in the same position.'''
        if not self.sorted:
            self.sort()
        tsite=None
        for site in self:
            if not tsite:
                tsite=site
            elif site==tsite:
                tsite.score+=site.score
            else:
                yield tsite
                tsite=site
        if tsite:
            yield tsite

class Fasta:
    '''Fasta format.'''
    def __init__(self,name=None,seq=None,description=None):
        '''Initiate the fasta record.'''
        self.id=name
        self.seq=seq
        self.description=description
        self.evalue=None
        self.structure=None
    def length(self):
        '''get the length of sequence.'''
        return len(self.seq)
    def __str__(self):
        '''String for output of Fasta.'''
        return ">"+self.id+"\n"+self.seq
    def Fold(self,constraints="",temperature=20):
        '''Fold the struture for the sequence with stem-loop constraints at fixed temperature.'''
        (self.evalue,self.structure)=wRNA.fold(self.seq,constraints,temperature)

class ColumnReader:
    '''Read column files. Support format:Bed,Gene/GeneBed/Tab/GenePred,Bowtie and SOAP'''
    def __init__(self,infile,format=None):
        '''Initiate file handel. If the line name is combined with reads number, use countsep to split it.'''
        if isinstance(infile,basestring):
            try:
                self.infile=open(infile)
            except Exception,e:
                print >>sys.stderr,"Can't open file.",e
                raise
        else:
            self.infile=infile
        try:
            self.line=self.infile.readline()
        except Exception,e:
            print >>sys.stderr,"Can't read file:",e
            raise
        #Move to the first valid line
        while self.line and self.line[0]=='#':
            self.line=self.infile.readline()
        self.format=format.lower() if format else None
    def next(self):
        '''Read one valid line from file.'''
        while self.line and self.line[0]=="#":
            self.line=self.infile.readline()
        if self.line:
            x=self.line.rstrip('\n\r').split('\t')
            self.line=self.infile.readline()
            if not self.format: return x
            if self.format=='bed':
                return Bed(x)
            if self.format=='bowtie':
                return Utils.BowtieToBed(x)
            if self.format=='soap':
                return Utils.SOAPToBed(x)
            if self.format in ['gene','genebed','genepred','tab']:
                return GeneBed(x)
        else:
            self.infile.close()
            raise StopIteration
    def __iter__(self):
        '''Iterator'''
        return self

class SeqReader:
    '''Read DNA/RNA sequence in Fasta/Fastaq/CSF format.'''
    def __init__(self,infile,format='fasta'):
        '''Initiate file handle.'''
        try:
            self.infile=open(infile)
        except:
            self.infile=infile
        #self.line=self.infile.readline()
        #while self.line and "#" in self.line:
        #    self.line=self.infile.readline()
        self.format=format.lower()
        #self.record=None
    def next(self):
        '''Next Record.'''
        line=self.infile.readline()
        if line:
            if self.format=='fasta' or self.format=='csf':
                description=line.lstrip('>').rstrip("\n\r")
                name=description.split()[0] #split by blanks
                os.linesep='>'
                line=self.infile.readline()
                os.linesep='\n'
                seq=line.rstrip('>').replace('\n','').replace('\r','').upper()
                return Fasta(name,seq,description)
            if self.format=='fastaq':
                name=line.lstrip('@').rstrip("\n\r")
                seq=self.infile.readline().rstrip("\n\r")
                self.infile.readline()
                description=self.infile.readline().rstrip("\n\r")
                return Fasta(name,seq,description)
                return 
            assert False, "File format error!"
        else:
            raise StopIteration
    def __iter__(self):
        return self
                          
class BedMap(dict):
    '''Map Beds with chromsomes and BIN values.'''
    def __init__(self,organism='ce6'):
        '''Default genome is ce6.'''
        dict.__init__(self)
        for chro in Utils.genomeSize(organism).keys():
            self[chro]={}
    def findOverlap(self,bed,fraction=0.5,forcestrand=False):
        '''Find overlaps with bed and put into bedlist.'''
        bedlist=BedList()
        minover=bed.length()*fraction
        maxover=0
        startBin,stopBin= bed.start >> _binFirstShift, (bed.stop-1) >> _binFirstShift
        for i in range(_binLevels):
            offset = _binOffsetsExtended[i]
            for j in xrange(startBin+offset,stopBin+offset+1):
                if not self[bed.chr].has_key(j):
                    continue
                for item in self[bed.chr][j]:
                    overlen=bed.overlapLength(item)
                    if overlen>=minover:
                        if not forcestrand or (forcestrand and bed.strand==item.strand):
                            if maxover<overlen:
                                maxover=overlen
                                bed.description=item
                            bedlist.append(item)
            startBin >>= _binNextShift
            stopBin >>= _binNextShift
        return bedlist
    def intersectBed(self,bed,fraction=0.5,outputoption=1,forcestrand=False):
        '''Intersect bed with items in BedMap.'''
        #outputoption=0 for No overlap,1 for best overlap and 2 for all overlaps}
        bedlist=self.findOverlap(bed,fraction,forcestrand)
        if outputoption==0: #No overlap
            return (True if len(bedlist)==0 else False)
        if outputoption==1: #best overlap
            return bed.description # if description==None, no overlap found.
        #return all overlaps
        return bedlist
    def loadBedToMap(self,bedlist=None,bedtype='bed'):
        '''Load Bed to BedMap from either Bedlist or ColumnReader handle or bedfile. Load one bedlist once.'''
        if bedlist:
            if isinstance(bedlist,str) or isinstance(bedlist,file):
                bedlist=ColumnReader(bedlist,bedtype)
            for bed in bedlist:
                _bin=bed.getBIN()
                self[bed.chr].setdefault(_bin,BedList())
                self[bed.chr][_bin].append(bed)
    def __iter__(self): #useless !!!
        '''Output the bedMap by iteration'''
        for chro in sorted(self.keys()):
            for _bin in sorted(self[chro].keys()):
                for bed in self[chro][_bin]:
                    yield bed

class Utils:
    def rc(seq):
        comps = {'A':"T", 'C':"G", 'G':"C", 'T':"A",
                'B':"V", 'D':"H", 'H':"D", 'K':"M",
                'M':"K", 'R':"Y", 'V':"B", 'Y':"R",
                'W':'W', 'N':'N', 'S':'S'}
        return ''.join([comps[x] for x in seq.upper()[::-1]])
    rc=staticmethod(rc)
    def MW(seq):
        mws={'A':313.21,'C':289.19,'G':329.21,'T':304.2,'I':314.2,'N':308.95,'R':321.21,'Y':296.69,'M':301.2,'K':316.7,'S':309.2,'W':308.71,'H':302.2,'B':307.53,'D':315.54,'V':310.53,'p':79.98,'X':0,'U':290.17}
        return sum([mws[x] for x in seq.upper()])
    MW=staticmethod(MW)
    def TM(seq,Na=100):
        seq=seq.upper()
        if len(seq)<25:
            tm={'A':2,'T':2,'C':4,'G':4}
            return sum([tm[x] for x in seq])
        else:
            N=float(len(seq))
            gc=(seq.count('C')+seq.count('G'))/N
            return 81.5+16.6*log(Na/1000.0,10)+0.41*gc+600.0/N
    TM=staticmethod(TM)
    def toRNA(seq):
        return seq.upper().replace('T','U')
    toRNA=staticmethod(toRNA)
    def toDNA(seq):
        return seq.upper().replace('U','T')
    toDNA=staticmethod(toDNA)
    def toProtein(seq,table):
        '''Translate DNA or RNA to protein according to standard translation table.'''
        seq=seq.upper().rstrip()
        if "U" in seq:
            seq=toDNA(seq)
        if len(seq)%3!=0:
            print >>sys.stderr, "Sequcence length should be 3*N."
            return None
        p=""
        for i in xrange(len(seq)/3):
            p+=table[seq[i*3:(i+1)*3]]
        return p
    toProtein=staticmethod(toProtein)
    def fastaToCSFasta(seq,starter='T'):
        trans= [ ['0','1','2','3'], ['1','0','3','2'], ['2','3','0','1'], ['3','2','1','0']]
        bases= {'A':0,'C':1,'G':2,'T':3}
        csf=starter
        seq=seq.upper()
        tseq=starter+seq
        for i in range(len(seq)):
            csf+=trans[bases[tseq[i]]][bases[tseq[i+1]]]
        return csf
    fastaToCSFasta=staticmethod(fastaToCSFasta)
    def csFastaToFasta(seq):
        trans=[ [0,1,2,3], [1,0,3,2], [2,3,0,1], [3,2,1,0]]
        bases={'A':0,'C':1,'G':2,'T':3}
        csbases={0:'A',1:'C',2:'G',3:'T'}
        f=seq[0]
        for i in range(len(seq)-1):
            f+=csbases[trans[bases[f[i]]][int(seq[i+1])]]
        return f[1:]
    csFastaToFasta=staticmethod(csFastaToFasta)
    def BowtieToBed(x):
        '''Bowtie map result to Bed.Bowtie is 0 based.'''
        x[3]=int(x[3])
        return Bed([x[2],x[3],x[3]+len(x[4]),x[0],1,x[1]])
    BowtieToBed=staticmethod(BowtieToBed)
    def SOAPToBed(x):
        '''SOAP map result to Bed. SOAP is 1 based.'''
        x[8]=int(x[8])-1
        return Bed([x[7],x[8],x[8]+int(x[5]),x[0],1,x[6]])
    SOAPToBed=staticmethod(SOAPToBed)
    def genomeSize(gversion='ce6'):
        '''Genome size dictionary.'''
        genome={}
        genome['ce6']={'chrI':15072421,'chrII':15279323,'chrIII':13783681,'chrIV':17493785,'chrM':13794,'chrV':20919568,'chrX':17718854}
        return genome[gversion]
    genomeSize=staticmethod(genomeSize)
    def translateTables(tabletype="standard"):
        '''Translation tables. @ for stop condons.'''
        tables={}
        tables["standard"]={
            'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
            'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
            'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
            'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
            'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
            'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T',
            'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
            'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A',
            'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D',
            'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
            'GGG': 'G', 'TAA': '@', 'TAG': '@', 'TGA': '@'}
        return tables[tabletype]
    translateTables=staticmethod(translateTables)
    def dump(handle,filename=None):
        '''Dump object into hard disk.'''
        if not filename:
            filename="tmp"+str(random.randint(0,99))+".dmp.gz"
        fh=gzip.open(filename,'wb')
        cPickle.dump(handle,fh)
        fh.close()
    dump=staticmethod(dump)
    def load(filename):
        '''Load dumpped object into memory.'''
        fh=gzip.open(filename,'rb')
        th=cPickle.load(fh)
        fh.close()
        return th
    load=staticmethod(load)
    def GFFReader(filename):
        '''Read GFF format file and transform to GeneBed format.'''
        chroTab={'I':'chrI','II':'chrII','III':'chrIII','IV':'chrIV','MtDNA':'chrM','V':'chrV','X':'chrX'}
        gchr=None
        for line in ColumnReader(open(filename)):
            order=line[8][0:2]
            name=line[8].split(';')[0].split(':')[1]
            if order=='ID':
                if gchr:
                    gexons.sort()
                    if gtxstart-1 in gexons:
                        cur=gexons.index(gtxstart)
                        del gexons[cur]
                        del gexons[cur-1]
                    if gtxstop+1 in gexons:
                        cur=gexons.index(gtxstop+1)
                        del gexons[cur]
                        del gexons[cur-1]
                    exonstarts=[]
                    exonstops=[]
                    exoncount=len(gexons)/2
                    for i in range(exoncount):
                        exonstarts.append(gexons[2*i]-1)
                        exonstops.append(gexons[i*2+1])
                    yield GeneBed([gname,gchr,gstrand,gstart-1,gstop,gtxstart-1,gtxstop,exoncount,exonstarts,exonstops])
                #Initiate new record
                gchr=chroTab[line[0]]
                gname=name
                gstrand=line[6]
                gstart=int(line[3])
                gstop=int(line[4])
                gtxstart=gstop
                gtxstop=gstart
                gexons=[]
            elif order=='Pa':
                start=int(line[3])
                stop=int(line[4])
                gexons.append(start)
                gexons.append(stop)
                if line[2]=='coding_exon':
                    gtxstart=min(gtxstart,start)
                    gtxstop=max(gtxstop,stop)
        gexons.sort()
        if gtxstart-1 in gexons:
            cur=gexons.index(gtxstart)
            del gexons[cur]
            del gexons[cur-1]
        if gtxstop+1 in gexons:
            cur=gexons.index(gtxstop+1)
            del gexons[cur]
            del gexons[cur-1]
        exonstarts=[]
        exonstops=[]
        exoncount=len(gexons)/2
        for i in range(exoncount):
            exonstarts.append(gexons[2*i]-1)
            exonstops.append(gexons[i*2+1])
        yield GeneBed([gname,gchr,gstrand,gstart-1,gstop,gtxstart-1,gtxstop,exoncount,exonstarts,exonstops])
    GFFReader=staticmethod(GFFReader)
    def GFF3Reader(filename):
        '''Read GFF3 format file and transform to GenePred format.'''
        #chroTab={'I':'chrI','II':'chrII','III':'chrIII','IV':'chrIV','MtDNA':'chrM','V':'chrV','X':'chrX'}
        tschr=None
        for line in ColumnReader(open(filename)):
            order=line[8][0:6]
            if line[2] in ["C_gene_segment","D_gene_segment","D_loop","J_gene_segment","match","region","rRNA","tRNA","V_gene_segment"]:
                pass
            elif not cmp(order,"ID=rna"):
                if tschr:
                    exonstarts=[]
                    exonends=[]
                    exoncount=len(tsexons)/2
                    for i in range(exoncount):
                        if not cmp(tsstrand,"+"):
                            exonstarts.append(tsexons[2*i]-1)
                            exonends.append(tsexons[i*2+1])
                        else:
                            exonstarts.append(tsexons[len(tsexons)-2*i-2]-1)
                            exonends.append(tsexons[len(tsexons)-i*2-1])
                    if not cmp(tsstrand, "+"):
                        tsstart=tsexons[0]-1
                        tsend=tsexons[-1]
                        if tscds:
                            cdsstart=tscds[0]-1
                            cdsend=tscds[-1]
                        else:
                            cdsstart=tsend
                            cdsend=tsend
                    else:
                        tsstart=tsexons[-2]-1
                        tsend=tsexons[1]
                        if tscds:
                            cdsstart=tscds[-2]-1
                            cdsend=tscds[1]
                        else:
                            cdsstart=tsend
                            cdsend=tsend
                    description=tsdesc+"|"+genedesc
                    yield GeneBed([tsname,_RefSeqAccnToChr[tschr],tsstrand,tsstart,tsend,cdsstart,cdsend,exoncount,exonstarts,exonends,0,genename,None,None,None,description])
                #Initiate new record
                tschr=line[0]
                tsstrand=line[6]
                tsexons=[]
                tscds=[]
                tsname=line[8].split(';')[1].split('=')[1]
                desc_list=line[8].split(';')
                Name=""
                product=""
                for item in desc_list:
                    item=item.split('=')
                    if not cmp(item[0],"Name"):
                        Name=item[1]
                        break
                    elif not cmp(item[0],"Dbxref"):
                        Dbxref=item[1]
                    elif not cmp(item[0],"product"):
                        product=item[1]
                if Name:
                    tsname=Name
                elif product:
                    tsname=product
                else:
                    tsname=Dbxref
                tsdesc=line[8]
                genename=lastgenename
                genedesc=lastgenedesc
                mark=1
            elif not cmp(line[8][0:7],"ID=gene"):
                lastgenename=line[8].split(';')[1].split('=')[1]
                lastgenedesc=line[8]
                mark=0
            else:
                if (not cmp(line[2],"exon")) and mark:
                    tsexons.append(int(line[3]))
                    tsexons.append(int(line[4]))
                elif (not cmp(line[2],"CDS")) and mark:
                    tscds.append(int(line[3]))
                    tscds.append(int(line[4]))
        exonstarts=[]
        exonends=[]
        exoncount=len(tsexons)/2
        for i in range(exoncount):
            if not cmp(tsstrand,"+"):
                exonstarts.append(tsexons[2*i]-1)
                exonends.append(tsexons[i*2+1])
            else:
                exonstarts.append(tsexons[len(tsexons)-2*i-2]-1)
                exonends.append(tsexons[len(tsexons)-i*2-1])
        if not cmp(tsstrand,"+"):
            tsstart=tsexons[0]-1
            tsend=tsexons[-1]
            if tscds:
                cdsstart=tscds[0]-1
                cdsend=tscds[-1]
            else:
                cdsstart=tsend
                cdsend=tsend
        else:
            tsstart=tsexons[-2]-1
            tsend=tsexons[1]
            if tscds:
                cdsstart=tscds[-2]-1
                cdsend=tscds[1]
            else:
                cdsstart=tsend
                cdsend=tsend
        description=tsdesc+"|"+genedesc
        yield GeneBed([tsname,_RefSeqAccnToChr[tschr],tsstrand,tsstart,tsend,cdsstart,cdsend,exoncount,exonstarts,exonends,0,genename,None,None,None,description])
    GFF3Reader=staticmethod(GFF3Reader)
    def ENSEMBLGTFReader(filename):
        '''Read GTF format file and transform to GeneBed format.'''
        #chroTab={'I':'chrI','II':'chrII','III':'chrIII','IV':'chrIV','MtDNA':'chrM','V':'chrV','X':'chrX'}
        tschr=None
        for line in ColumnReader(open(filename)):
            for i in line[8].split(';'):
                if not cmp(i.split('"')[0]," exon_number "):
                    exon_number=i.split('"')[1]
            if (not cmp(exon_number,"1")) and (not cmp(line[2],"exon")):
                if tschr:
                    exonstarts=[]
                    exonends=[]
                    exoncount=len(tsexons)/2
                    for i in range(exoncount):
                        if not cmp(tsstrand,"+"):
                            exonstarts.append(tsexons[2*i]-1)
                            exonends.append(tsexons[i*2+1])
                        else:
                            exonstarts.append(tsexons[len(tsexons)-2*i-2]-1)
                            exonends.append(tsexons[len(tsexons)-i*2-1])
                    if not cmp(tsstrand, "+"):
                        tsstart=tsexons[0]-1
                        tsend=tsexons[-1]
                        if tscds:
                            cdsstart=tscds[0]-1
                            cdsend=tscds[-1]
                        else:
                            cdsstart=tsend
                            cdsend=tsend
                        if startcodon:
                            cdsstart=startcodon[0]-1
                        if stopcodon:
                            cdsend=stopcodon[1]
                    else:
                        tsstart=tsexons[-2]-1
                        tsend=tsexons[1]
                        if tscds:
                            cdsstart=tscds[-2]-1
                            cdsend=tscds[1]
                        else:
                            cdsstart=tsend
                            cdsend=tsend
                        if startcodon:
                            cdsend=startcodon[1]
                        if stopcodon:
                            cdsstart=stopcodon[0]-1
                    #description="Gene_id="+geneid+";"+"Gene_name="+genename+";"+"Gene_biotype="+genebiotype+";"+"Transcript_id="+tsid+";"+"Transcript_name="+tsname+";"+"Transcript_biotype="+tsbiotype+";"
                    description="Gene_id="+geneid+";"+"Gene_name="+genename+";"+"Transcript_id="+tsid+";"+"Transcript_name="+tsname+";"
                    yield GeneBed([tsid,tschr,tsstrand,tsstart,tsend,cdsstart,cdsend,exoncount,exonstarts,exonends,0,geneid,None,None,None,description])
                #Initiate new record
                tschr=line[0]
                tsstrand=line[6]
                tsexons=[]
                tscds=[]
                startcodon=[]
                stopcodon=[]
                for i in line[8].split(';'):
                    if not cmp(i.split('"')[0]," transcript_id "):
                        tsid=i.split('"')[1]
                    #if not cmp(i.split('"')[0]," linc_name "):
                    if not cmp(i.split('"')[0]," transcript_name "):
                        tsname=i.split('"')[1]
                    if not cmp(i.split('"')[0]," gene_id "):
                        geneid=i.split('"')[1]
                #tsname=line[8].split(';')[4].split('"')[1]
                #tsbiotype=line[1]
                #geneid=line[8].split(';')[0].split('"')[1]
                genename=tsname
                #genename=line[8].split(';')[7].split('"')[1]
                #genebiotype=line[8].split(';')[4].split('"')[1]
                if not cmp(line[2],"exon"):
                    tsexons.append(int(line[3]))
                    tsexons.append(int(line[4]))
                elif not cmp(line[2],"CDS"):
                    tscds.append(int(line[3]))
                    tscds.append(int(line[4]))
                elif not cmp(line[2],"start_codon"):
                    startcodon.append(int(line[3]))
                    startcodon.append(int(line[4]))
                elif not cmp(line[2],"stop_codon"):
                    stopcodon.append(int(line[3]))
                    stopcodon.append(int(line[4]))
            else:
                if not cmp(line[2],"exon"):
                    tsexons.append(int(line[3]))
                    tsexons.append(int(line[4]))
                elif not cmp(line[2],"CDS"):
                    tscds.append(int(line[3]))
                    tscds.append(int(line[4]))
                elif not cmp(line[2],"start_codon"):
                    startcodon.append(int(line[3]))
                    startcodon.append(int(line[4]))
                elif not cmp(line[2],"stop_codon"):
                    stopcodon.append(int(line[3]))
                    stopcodon.append(int(line[4]))
        exonstarts=[]
        exonends=[]
        exoncount=len(tsexons)/2
        for i in range(exoncount):
            if not cmp(tsstrand,"+"):
                exonstarts.append(tsexons[2*i]-1)
                exonends.append(tsexons[i*2+1])
            else:
                exonstarts.append(tsexons[len(tsexons)-2*i-2]-1)
                exonends.append(tsexons[len(tsexons)-i*2-1])
        if not cmp(tsstrand,"+"):
            tsstart=tsexons[0]-1
            tsend=tsexons[-1]
            if tscds:
                cdsstart=tscds[0]-1
                cdsend=tscds[-1]
            else:
                cdsstart=tsend
                cdsend=tsend
            if startcodon:
                cdsstart=startcodon[0]-1
            if stopcodon:
                cdsend=stopcodon[1]
        else:
            tsstart=tsexons[-2]-1
            tsend=tsexons[1]
            if tscds:
                cdsstart=tscds[-2]-1
                cdsend=tscds[1]
            else:
                cdsstart=tsend
                cdsend=tsend
            if startcodon:
                cdsend=startcodon[1]
            if stopcodon:
                cdsstart=stopcodon[0]-1
        #description="Gene_id="+geneid+";"+"Gene_name="+genename+";"+"Gene_biotype="+genebiotype+";"+"Transcript_id="+tsid+";"+"Transcript_name="+tsname+";"+"Transcript_biotype="+tsbiotype+";"
        description="Gene_id="+geneid+";"+"Gene_name="+genename+";"+"Transcript_id="+tsid+";"+"Transcript_name="+tsname+";"
        yield GeneBed([tsid,tschr,tsstrand,tsstart,tsend,cdsstart,cdsend,exoncount,exonstarts,exonends,0,geneid,None,None,None,description])
    ENSEMBLGTFReader=staticmethod(ENSEMBLGTFReader)
    def HumanLincRNACatalogGTFReader(filename):
        '''Read GTF format file and transform to GeneBed format.'''
        #chroTab={'I':'chrI','II':'chrII','III':'chrIII','IV':'chrIV','MtDNA':'chrM','V':'chrV','X':'chrX'}
        tschr=None
        for line in ColumnReader(open(filename)):
            for i in line[8].split(';'):
                if not cmp(i.split('"')[0]," exon_number "):
                    exon_number=i.split('"')[1]
            if (not cmp(exon_number,"1")) and (not cmp(line[2],"exon")):
                if tschr:
                    exonstarts=[]
                    exonends=[]
                    exoncount=len(tsexons)/2
                    for i in range(exoncount):
                        exonstarts.append(tsexons[2*i]-1)
                        exonends.append(tsexons[i*2+1])
                    tsstart=tsexons[0]-1
                    tsend=tsexons[-1]
                    if tscds:
                        cdsstart=tscds[0]-1
                        cdsend=tscds[-1]
                    else:
                        cdsstart=tsend
                        cdsend=tsend
                    if startcodon:
                        cdsstart=startcodon[0]-1
                    if stopcodon:
                        cdsend=stopcodon[1]
                    #description="Gene_id="+geneid+";"+"Gene_name="+genename+";"+"Gene_biotype="+genebiotype+";"+"Transcript_id="+tsid+";"+"Transcript_name="+tsname+";"+"Transcript_biotype="+tsbiotype+";"
                    description="Gene_id="+geneid+";"+"Gene_name="+genename+";"+"Transcript_id="+tsid+";"+"Transcript_name="+tsname+";"
                    yield GeneBed([tsid,tschr,tsstrand,tsstart,tsend,cdsstart,cdsend,exoncount,exonstarts,exonends,0,geneid,None,None,None,description])
                #Initiate new record
                tschr=line[0]
                tsstrand=line[6]
                tsexons=[]
                tscds=[]
                startcodon=[]
                stopcodon=[]
                for i in line[8].split(';'):
                    if not cmp(i.split('"')[0]," transcript_id "):
                        tsid=i.split('"')[1]
                    if not cmp(i.split('"')[0]," linc_name "):
                        tsname=i.split('"')[1]
                    if not cmp(i.split('"')[0],"gene_id "):
                        geneid=i.split('"')[1]
                #tsname=line[8].split(';')[4].split('"')[1]
                #tsbiotype=line[1]
                #geneid=line[8].split(';')[0].split('"')[1]
                genename=tsname
                #genename=line[8].split(';')[7].split('"')[1]
                #genebiotype=line[8].split(';')[4].split('"')[1]
                if not cmp(line[2],"exon"):
                    tsexons.append(int(line[3]))
                    tsexons.append(int(line[4]))
                elif not cmp(line[2],"CDS"):
                    tscds.append(int(line[3]))
                    tscds.append(int(line[4]))
                elif not cmp(line[2],"start_codon"):
                    startcodon.append(int(line[3]))
                    startcodon.append(int(line[4]))
                elif not cmp(line[2],"stop_codon"):
                    stopcodon.append(int(line[3]))
                    stopcodon.append(int(line[4]))
            else:
                if not cmp(line[2],"exon"):
                    tsexons.append(int(line[3]))
                    tsexons.append(int(line[4]))
                elif not cmp(line[2],"CDS"):
                    tscds.append(int(line[3]))
                    tscds.append(int(line[4]))
                elif not cmp(line[2],"start_codon"):
                    startcodon.append(int(line[3]))
                    startcodon.append(int(line[4]))
                elif not cmp(line[2],"stop_codon"):
                    stopcodon.append(int(line[3]))
                    stopcodon.append(int(line[4]))
        exonstarts=[]
        exonends=[]
        exoncount=len(tsexons)/2
        for i in range(exoncount):
            exonstarts.append(tsexons[2*i]-1)
            exonends.append(tsexons[i*2+1])
        tsstart=tsexons[0]-1
        tsend=tsexons[-1]
        if tscds:
            cdsstart=tscds[0]-1
            cdsend=tscds[-1]
        else:
            cdsstart=tsend
            cdsend=tsend
        if startcodon:
            cdsstart=startcodon[0]-1
        if stopcodon:
            cdsend=stopcodon[1]
        #description="Gene_id="+geneid+";"+"Gene_name="+genename+";"+"Gene_biotype="+genebiotype+";"+"Transcript_id="+tsid+";"+"Transcript_name="+tsname+";"+"Transcript_biotype="+tsbiotype+";"
        description="Gene_id="+geneid+";"+"Gene_name="+genename+";"+"Transcript_id="+tsid+";"+"Transcript_name="+tsname+";"
        yield GeneBed([tsid,tschr,tsstrand,tsstart,tsend,cdsstart,cdsend,exoncount,exonstarts,exonends,0,geneid,None,None,None,description])
    HumanLincRNACatalogGTFReader=staticmethod(HumanLincRNACatalogGTFReader)
    def CompareInterval(I1start,I1end,I2start,I2end):
        '''Compute the overlapping of two intervals,0-based'''
        I1start=I1start+1
        I2start=I2start+1
        if I1end<I2start:
            return 0
        elif I1start>I2end:
            return 0
        elif I1start<=I2start and I1end>=I2start and I1end<=I2end:
            return I1end-I2start+1
        elif I1start<=I2start and I1end>=I2end:
            return I2end-I2start+1
        elif I1start>=I2start and I1end<=I2end:
            return I1end-I1start+1
        elif I1start>=I2start and I1start<=I2end and I1end>=I2end:
            return I2end-I1start+1
    CompareInterval=staticmethod(CompareInterval)
    def CompareGenePred(filename1,filename2):
        '''Compare two genePred file, find overlapping genes'''
        current=[]
        for line in ColumnReader(open(filename1)):
            if current:
                if current[1]<line[1]:
                    pass
                elif current[1]>line[1]:
                    continue
                elif (not cmp(current[1],line[1])) and (not cmp(current[2],line[2])):
                    if int(current[4])<int(line[3]):
                        pass
                    elif int(current[3])>int(line[4]):
                        continue
                    else:
                        lexoncount=int(line[7])
                        lstart=line[8].split(",")
                        lend=line[9].split(",")
                        rexoncount=int(current[7])
                        rstart=current[8].split(",")
                        rend=current[9].split(",")
                        overlapsize=0
                        lsize=0
                        rsize=0
                        for i in range(lexoncount):
                            lsize+=int(lend[i])-int(lstart[i])
                        for i in range(rexoncount):
                            rsize+=int(rend[i])-int(rstart[i])
                        for i in range(lexoncount):
                            for j in range(rexoncount):
                                overlapsize+=Utils.CompareInterval(int(lstart[i]),int(lend[i]),int(rstart[j]),int(rend[j]))
                        if float(overlapsize)/min(lsize,rsize)>0.9:
                            return "\t".join(line[0:10])+"\t"+str(float(overlapsize)/lsize)+"\t"+"\t".join(current[0:10])+"\t"+str(float(overlapsize)/rsize)
            for row in ColumnReader(open(filename2)):
                if row[1]<line[1]:
                    continue
                elif row[1]>line[1]:
                    current=row
                    break
                elif (not cmp(row[1],line[1])) and (not cmp(row[2],line[2])):
                    if int(row[4])<int(line[3]):
                        continue
                    elif int(row[3])>int(line[4]):
                        current=row
                        break
                    else:
                        lexoncount=int(line[7])
                        lstart=line[8].split(",")
                        lend=line[9].split(",")
                        rexoncount=int(row[7])
                        rstart=row[8].split(",")
                        rend=row[9].split(",")
                        overlapsize=0
                        lsize=0
                        rsize=0
                        for i in range(lexoncount):
                            lsize+=int(lend[i])-int(lstart[i])
                        for i in range(rexoncount):
                            rsize+=int(rend[i])-int(rstart[i])
                        for i in range(lexoncount):
                            for j in range(rexoncount):
                                overlapsize+=Utils.CompareInterval(int(lstart[i]),int(lend[i]),int(rstart[j]),int(rend[j]))
                        if float(overlapsize)/min(lsize,rsize)>0.9:
                            return "\t".join(line[0:10])+"\t"+str(float(overlapsize)/lsize)+"\t"+"\t".join(row[0:10])+"\t"+str(float(overlapsize)/rsize)
    CompareGenePred=staticmethod(CompareGenePred)
    def GENCODEGTFReader(filename):
        '''Read GTF format file and transform to GeneBed format.'''
        #chroTab={'I':'chrI','II':'chrII','III':'chrIII','IV':'chrIV','MtDNA':'chrM','V':'chrV','X':'chrX'}
        tschr=None
        for line in ColumnReader(open(filename)):
            #exon_number=line[8].split(';')[2].split('"')[1]
            if not cmp(line[2],"transcript"):
                if tschr:
                    exonstarts=[]
                    exonends=[]
                    exoncount=len(tsexons)/2
                    for i in range(exoncount):
                        if not cmp(tsstrand,"+"):
                            exonstarts.append(tsexons[2*i]-1)
                            exonends.append(tsexons[i*2+1])
                        else:
                            exonstarts.append(tsexons[len(tsexons)-2*i-2]-1)
                            exonends.append(tsexons[len(tsexons)-i*2-1])
                    if not cmp(tsstrand, "+"):
                        tsstart=tsexons[0]-1
                        tsend=tsexons[-1]
                        if tscds:
                            cdsstart=tscds[0]-1
                            cdsend=tscds[-1]
                        else:
                            cdsstart=tsend
                            cdsend=tsend
                        if startcodon:
                            cdsstart=startcodon[0]-1
                        if stopcodon:
                            cdsend=stopcodon[1]
                    else:
                        tsstart=tsexons[-2]-1
                        tsend=tsexons[1]
                        if tscds:
                            cdsstart=tscds[-2]-1
                            cdsend=tscds[1]
                        else:
                            cdsstart=tsend
                            cdsend=tsend
                        if startcodon:
                            cdsend=startcodon[1]
                        if stopcodon:
                            cdsstart=stopcodon[0]-1
                    #description="Gene_id:"+geneid+";"+"Gene_name:"+genename+";"+"Gene_type:"+genetype+";"+"Transcript_id:"+tsid+";"+"Transcript_name:"+tsname+";"+"Transcript_type:"+tstype+";"+"Havana_gene:"+havana_gene+";"+"Havana_transcript:"+havana_ts+";"
                    yield GeneBed([tsid,tschr,tsstrand,tsstart,tsend,cdsstart,cdsend,exoncount,exonstarts,exonends,0,geneid,None,None,None,description])
                #Initiate new record
                tschr=line[0]
                tsstrand=line[6]
                tsexons=[]
                tscds=[]
                startcodon=[]
                stopcodon=[]
                '''
                line8=line[8].split(';')
                tsid=line8[1].split('"')[1]
                tstype=line8[5].split('"')[1]
                tsname=line8[7].split('"')[1]
                geneid=line8[0].split('"')[1]
                genename=line8[4].split('"')[1]
                genetype=line8[2].split('"')[1]
                havana_gene=line8[10].split('"')[1]
                havana_ts=line8[11].split('"')[1]
                '''
                line8=line[8].split(';')
                tsid=line8[1].split('"')[1]
                geneid=line8[0].split('"')[1]
                description=line[8]
            else:
                if not cmp(line[2],"exon"):
                    tsexons.append(int(line[3]))
                    tsexons.append(int(line[4]))
                elif not cmp(line[2],"CDS"):
                    tscds.append(int(line[3]))
                    tscds.append(int(line[4]))
                elif not cmp(line[2],"start_codon"):
                    startcodon.append(int(line[3]))
                    startcodon.append(int(line[4]))
                elif not cmp(line[2],"stop_codon"):
                    stopcodon.append(int(line[3]))
                    stopcodon.append(int(line[4]))
        exonstarts=[]
        exonends=[]
        exoncount=len(tsexons)/2
        for i in range(exoncount):
            if not cmp(tsstrand,"+"):
                exonstarts.append(tsexons[2*i]-1)
                exonends.append(tsexons[i*2+1])
            else:
                exonstarts.append(tsexons[len(tsexons)-2*i-2]-1)
                exonends.append(tsexons[len(tsexons)-i*2-1])
        if not cmp(tsstrand,"+"):
            tsstart=tsexons[0]-1
            tsend=tsexons[-1]
            if tscds:
                cdsstart=tscds[0]-1
                cdsend=tscds[-1]
            else:
                cdsstart=tsend
                cdsend=tsend
            if startcodon:
                cdsstart=startcodon[0]-1
            if stopcodon:
                cdsend=stopcodon[1]
        else:
            tsstart=tsexons[-2]-1
            tsend=tsexons[1]
            if tscds:
                cdsstart=tscds[-2]-1
                cdsend=tscds[1]
            else:
                cdsstart=tsend
                cdsend=tsend
            if startcodon:
                cdsend=startcodon[1]
            if stopcodon:
                cdsstart=stopcodon[0]-1
        #description="Gene_id:"+geneid+";"+"Gene_name:"+genename+";"+"Gene_type:"+genetype+";"+"Transcript_id:"+tsid+";"+"Transcript_name:"+tsname+";"+"Transcript_type:"+tstype+";"+"Havana_gene:"+havana_gene+";"+"Havana_transcript:"+havana_ts+";"
        yield GeneBed([tsid,tschr,tsstrand,tsstart,tsend,cdsstart,cdsend,exoncount,exonstarts,exonends,0,geneid,None,None,None,description])
    GENCODEGTFReader=staticmethod(GENCODEGTFReader)
    def interregion(n1,n2,n3,n4):
        if n4<n1 or n3>n2:
            return 0
        elif n3<=n1 and n4>=n2:
            return n2-n1+1
        elif n3<n1 and n4>=n1 and n4<=n2:
            return n4-n1+1
        elif n3>=n1 and n4<=n2:
            return n4-n3+1
        elif n3>=n1 and n3<=n2 and n4>n2:
            return n2-n3+1
    #interregion=staticmethod(interregion)
    def ConvertChr(chro):
        if chro in _RefSeqAccnToChr.keys():
            return _RefSeqAccnToChr[chro]
        elif chro in _GenBankAccnToChr.keys():
            return _GenBankAccnToChr[chro]
        elif chro in _UCSCAccnToChr.keys():
            return _UCSCAccnToChr[chro]
        else:
            return chro
    ConvertChr=staticmethod(ConvertChr)
