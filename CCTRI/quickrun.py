import plotly.express as px
import numpy as np
import math
from pathlib import Path


import matplotlib.pyplot as plt






data = [[64822, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 4368, 1025373, 9039724, 9898648, 0, 0, 8380968, 9898277, 9052604, 0, 4682, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [120758, 9987750, 9987750, 9987750, 9987750, 9987750, 9987750, 9987750, 9987750, 9987750, 9987750, 9987750, 9987750, 9987750, 9987750, 9987750, 9987754, 9987753, 9987752, 9987763, 9987758, 9988048, 9988701, 9902872, 0, 0, 8286378, 9891192, 0, 0, 9988395, 9988357, 9988344, 9988340, 9988339, 9988336, 9988336, 9988336, 9988336, 9988336, 9988336, 0, 0, 9988336, 9988336, 9988336, 9988336, 9988336, 9988336, 9988336], [632280, 9953944, 9953944, 9953944, 9953944, 9953944, 9953944, 9953944, 9953942, 9953942, 9953940, 9953940, 9953937, 9953942, 9953943, 9953944, 9953943, 9953943, 9953910, 9953937, 9954048, 9954710, 9957761, 9897605, 0, 0, 7986742, 9850864, 0, 9958223, 9956976, 9956841, 9956796, 9956791, 9956781, 9956791, 9956793, 9956793, 9956791, 9956793, 9956793, 0, 0, 9956795, 0, 9956795, 9956795, 9956795, 9956795, 9956795], [1336160, 9901322, 9901320, 9901314, 9901314, 9901312, 9901324, 9901322, 9901312, 9901320, 9901298, 9901302, 9901306, 9901306, 9901289, 9901288, 9901346, 0, 9901396, 9901347, 9901614, 9902854, 9907338, 9864132, 0, 0, 7443784, 9759582, 0, 0, 9908904, 9908143, 9908047, 9907994, 9907949, 9907982, 9907974, 9907966, 9907976, 9907963, 9907959, 0, 0, 9907954, 9907960, 9907964, 9907966, 9907964, 9907960, 9907960], [2040510, 9830775, 9830771, 9830765, 9830776, 9830814, 9830726, 9830746, 9830693, 9830752, 9830714, 9830710, 9830651, 9830672, 9830674, 9830778, 9830810, 0, 9831206, 9831405, 9831680, 9833719, 9838285, 9797137, 0, 0, 6632712, 9590698, 0, 0, 9845244, 9843522, 9843144, 9843064, 9842958, 9842881, 9842854, 9842822, 9842860, 9842844, 9842840, 0, 0, 9842816, 0, 9842820, 9842818, 9842810, 9842818, 9842816], [2742470, 9744232, 9744228, 9744210, 9744078, 9744104, 9744137, 9744202, 9744038, 9744106, 9744068, 9743944, 9744192, 9744198, 9744336, 9744314, 9744528, 9744308, 9744490, 9745018, 9746104, 9748448, 9749452, 9690653, 0, 0, 5580074, 9293526, 0, 0, 9766306, 9764073, 9763184, 9762772, 9762312, 9762286, 9762180, 9762128, 9762330, 9762234, 9762172, 0, 0, 9762188, 0, 9762208, 9762216, 9762192, 9762180, 9762178], [2838240, 9639350, 9639166, 9639046, 9639154, 9639086, 9638830, 9638928, 9639106, 9639058, 9639114, 9639656, 9639440, 9639492, 9639470, 9639496, 9639532, 9639902, 9640118, 9640444, 9640920, 9641318, 9630742, 9524874, 0, 4327404, 4408810, 8805386, 0, 0, 9666759, 9665412, 9664074, 9663012, 9662741, 9662420, 9662500, 9662354, 9662226, 9662274, 9662201, 0, 0, 9662226, 0, 9662144, 9662244, 9662166, 9662122, 9662062], [2667868, 9509028, 9509146, 9509324, 9508898, 9509074, 9509198, 9508672, 9509078, 9509210, 9508838, 9509572, 9509258, 9509570, 9509510, 9509496, 9509204, 0, 9508834, 9507750, 9503414, 9494010, 9452312, 9256710, 0, 4347116, 3292880, 8058844, 0, 0, 9527160, 9532908, 9534152, 9534144, 9534546, 9534248, 9534246, 9534516, 9534280, 9534478, 9534634, 9534464, 0, 9534594, 0, 9534838, 9534716, 9534796, 9534720, 9534608], [2303778, 9336922, 9337226, 9337124, 9337044, 9336572, 9336362, 9336376, 9335274, 9335658, 9336194, 9335332, 9335086, 9333638, 9333824, 9332220, 9331018, 9328288, 9324178, 9315998, 9299410, 9262892, 9154662, 0, 0, 3980194, 2375448, 7051576, 8736370, 0, 9297252, 9334046, 9347716, 9352928, 9356766, 9357792, 9359014, 9359836, 9360134, 9360294, 9360946, 9360412, 0, 9360724, 0, 9360988, 9360778, 9360492, 9360376, 9360792], [1817810, 9089372, 9089236, 9089062, 9088876, 9088390, 9088462, 9087924, 9087400, 9085964, 9084366, 9083592, 9083438, 9080108, 9076378, 9071970, 9066304, 9053222, 9039000, 9012514, 8961286, 8860850, 8638070, 0, 0, 3328380, 1676958, 5809544, 0, 0, 8877812, 8990112, 9040298, 9064446, 9078672, 9087170, 9092526, 9095774, 9098770, 9101282, 9102380, 9102806, 0, 9104052, 0, 9103400, 9103336, 9103986, 9103706, 9104370], [1619500, 8707832, 8706830, 8706926, 8706800, 8706988, 8705556, 8702574, 8700966, 8698794, 8695216, 8691846, 8686140, 8678274, 8668228, 8654478, 8634574, 8606468, 8564430, 8494998, 8380152, 8171308, 7772460, 0, 0, 2582128, 1113988, 4434592, 6518454, 0, 8104444, 8357340, 8487528, 8559792, 8605148, 8631558, 8647466, 8660118, 8669770, 8676224, 8680396, 8682760, 0, 8686256, 0, 8688684, 8688700, 8689846, 8690328, 8689330], [1299888, 8108510, 8108476, 8107512, 8106264, 8106032, 8100638, 8096740, 8094628, 8085874, 8078296, 8066376, 8051610, 8034178, 8009136, 7974850, 7930438, 7867312, 7775244, 7631832, 7418612, 7069856, 6503834, 0, 0, 1924196, 617648, 0, 0, 6054160, 6788584, 7225010, 7488038, 7650926, 7755048, 7823180, 7868988, 7903832, 7927134, 7944394, 7954366, 7965566, 0, 7979490, 0, 7987436, 7988866, 7991728, 7993492, 7993264], [1132200, 7192744, 7191352, 7190050, 7185196, 7181184, 7177042, 7170332, 7158412, 7145248, 7123876, 7099484, 7072210, 7037928, 6989780, 6927186, 6845004, 6732298, 6575360, 6358630, 6049316, 5611450, 4985738, 0, 0, 1489692, 124210, 1696806, 0, 4139378, 4923726, 5477660, 5862292, 6127682, 6317806, 6450010, 6542012, 6616658, 6668074, 6707270, 6741086, 0, 0, 6794294, 0, 6814760, 6823118, 6828976, 6834564, 6836618], [1091020, 5914284, 5913760, 5911124, 5908530, 5903492, 5892526, 5881936, 5864988, 5842848, 5817304, 5780822, 5738696, 5686832, 5615422, 5528678, 5418278, 5271684, 5082148, 4831676, 4510658, 4102932, 3577610, 0, 0, 1277934, 338546, 592176, 0, 2217702, 2845848, 3357378, 3752670, 4061930, 4293038, 4473676, 4609146, 4712054, 4796508, 4863510, 4914902, 0, 0, 5010636, 0, 5051700, 5066158, 5080464, 5090734, 5095830], [1122904, 4438332, 4436966, 4436538, 4432532, 4426738, 4417658, 4404646, 4388356, 4367218, 4333962, 4298828, 4250236, 4201156, 4125964, 4033940, 0, 3790734, 3622968, 3418236, 3175366, 2884226, 2539390, 0, 0, 1228230, 734322, 236006, 0, 685910, 1079630, 1414182, 1702996, 1950216, 2147834, 2311584, 2438684, 2543646, 2635660, 2705264, 2752870, 2804842, 0, 2874612, 0, 2929304, 2948030, 2962564, 2977722, 2982744], [1197210, 3101526, 3099416, 3100590, 3100120, 3094300, 3089672, 3077452, 3067176, 3049788, 3025852, 2998808, 2966584, 2922092, 2871596, 2809658, 2733576, 2639764, 2526888, 2397102, 2251996, 2089628, 1908570, 0, 0, 1250350, 1014360, 777106, 0, 324520, 126640, 52044, 216264, 353726, 473574, 573786, 661910, 736428, 799778, 848180, 891238, 0, 0, 983496, 1006532, 1028854, 1046800, 1060372, 1073826, 1080420], [1259422, 2173960, 2174802, 2174276, 2175080, 2173018, 2171982, 2168078, 2161482, 2148386, 2135720, 2121892, 2100708, 2075746, 2050272, 2017170, 1974196, 1930510, 1875732, 1807060, 1736998, 1661262, 1579068, 0, 0, 1292410, 1191740, 0, 0, 898780, 810746, 728004, 659414, 592924, 539516, 491412, 449574, 406862, 375730, 348806, 323124, 0, 0, 269474, 257864, 240004, 226504, 220680, 210620, 209344], [1305918, 1666414, 1668482, 1667484, 1667960, 1670242, 1667620, 1668116, 1664580, 1662480, 1656590, 1651766, 1643176, 1633742, 1620274, 1604964, 1586760, 1571016, 1547924, 1522566, 1491442, 1459472, 1432080, 0, 0, 1325468, 1287522, 1249974, 0, 1177318, 1145824, 1114582, 1088084, 1060894, 1038664, 1018884, 1001768, 989730, 976584, 962218, 953042, 0, 0, 923678, 914238, 907904, 902996, 898558, 894036, 893130], [1323508, 1441642, 1445090, 1446386, 1447532, 1447548, 1445940, 1449222, 1448150, 1447200, 1443298, 1442178, 1440550, 1437768, 1431808, 1427586, 1424196, 1415210, 1406982, 1400074, 1391910, 1378710, 1367380, 0, 0, 1336184, 1321894, 1310194, 0, 1288630, 1277160, 1265460, 1257394, 1251656, 1244140, 1237804, 1232822, 1228764, 1225990, 1222244, 1217544, 0, 0, 1203764, 1202282, 1199998, 1200616, 1196130, 1194730, 1194854], [1329074, 1360060, 1360224, 1359516, 1360032, 1362120, 1362600, 1362998, 1363006, 1362036, 1361998, 1363318, 1359858, 1357896, 1355708, 1354098, 0, 1351298, 1349930, 1348024, 1345102, 1340372, 1338544, 0, 0, 1331638, 1327704, 0, 0, 1318502, 1314142, 1312582, 1310838, 1306480, 1306356, 1305020, 1303780, 1305406, 1302802, 1301240, 1301638, 0, 0, 1298300, 1294372, 1293830, 1293242, 1294106, 1292128, 1293578], [1316476, 1324036, 1324254, 1324504, 1324254, 1323636, 1324624, 1323982, 1323674, 1324132, 1323166, 1322662, 1323046, 1321942, 1322732, 1323622, 1322684, 1321302, 1321654, 1322320, 1318280, 1316486, 1315912, 1317436, 1317492, 1317830, 1317122, 1316756, 1315334, 1315482, 1314638, 1315186, 1315248, 1314370, 1311638, 1311290, 1311176, 1309898, 1310110, 1312802, 1311022, 0, 0, 1308504, 1307912, 1307236, 1308486, 1308864, 1306312, 1308136], [1298062, 1295544, 1295602, 1297738, 1299506, 1298266, 1297328, 1298420, 1298580, 1297154, 1297660, 1296920, 1297226, 1295988, 1296288, 1295610, 1294978, 1295976, 1294000, 1293682, 1294682, 1296516, 0, 1293742, 1295664, 1295234, 1293228, 1294328, 1294040, 1295268, 1296266, 1295226, 1295074, 1293732, 1294002, 1294426, 1294546, 1292280, 1295762, 1295996, 1296124, 0, 0, 1296544, 1294238, 1295342, 1295918, 1297404, 1295906, 1296454], [1264574, 1264692, 1265006, 1266516, 1266226, 1264730, 1266214, 1268114, 1269252, 1267386, 1266124, 1265238, 1264186, 1264710, 1265288, 1264130, 1264108, 1266806, 1263934, 1261984, 1264388, 1264288, 0, 1266728, 1262792, 1264064, 1265796, 1264136, 1266642, 1267358, 1266696, 1265898, 1264808, 1266924, 1265510, 1266590, 1265210, 1265440, 1266368, 1267114, 1266768, 0, 0, 1267750, 1266386, 1266046, 1265262, 1266104, 1266004, 1265182], [1227362, 1230244, 1230584, 1231844, 1229116, 1228028, 1227358, 1230026, 1230076, 1231138, 1230364, 1229602, 1230542, 1229706, 1229436, 1228918, 1228928, 1228256, 1229512, 1227834, 1228316, 1229174, 0, 1228264, 1227644, 1229524, 1227922, 1227730, 1227038, 1228394, 1228090, 1226872, 1228056, 1227668, 1228150, 1229588, 1230276, 1229324, 1229294, 1230274, 1229480, 0, 0, 1231566, 1230846, 1233470, 1231916, 1231790, 1232492, 1232456], [1183568, 1183058, 1184086, 1184434, 1182808, 1184852, 1183884, 1183820, 1183370, 1183952, 1183912, 1185264, 1184948, 1185344, 1184026, 1184792, 1183858, 1183980, 1184684, 1185584, 1185144, 0, 0, 1183750, 1182226, 1182122, 1183616, 1184058, 1185308, 1184034, 1184856, 1185536, 1183712, 1184364, 1185114, 1186066, 1184166, 1183590, 1185672, 1186132, 0, 0, 0, 1185996, 1184288, 1183396, 1184908, 1185308, 1184680, 1186216], [1128756, 1132932, 1131800, 1130962, 1131682, 1131928, 1132932, 1132224, 1132236, 1132444, 1133452, 1132836, 1132828, 1131632, 1131738, 1132308, 1132510, 1131422, 1131100, 1133632, 1133178, 0, 0, 1131132, 1129998, 1129304, 1130354, 1131036, 1131446, 1132104, 1133424, 1132988, 1133378, 1132344, 1131984, 1131162, 1129344, 1131048, 1130726, 1133250, 0, 0, 0, 1133394, 1133648, 1132884, 1132366, 1132734, 1134382, 1133724], [1071828, 1074566, 1074846, 1073800, 1074612, 1074718, 1074376, 1073356, 1074498, 1074446, 1073846, 1072626, 1073090, 1073754, 1073196, 1073380, 1074228, 1074968, 1075664, 1074244, 1073996, 0, 0, 1073906, 1073574, 1074852, 1074128, 1076152, 1075612, 1074844, 1075492, 1075660, 1073692, 1073946, 1074538, 1073148, 1071514, 1070594, 1070822, 1070876, 0, 0, 0, 1073862, 0, 1073826, 1074780, 1073202, 1071954, 1073922], [1012102, 1013958, 1013912, 1014154, 1012930, 1012454, 1012602, 1013240, 1013684, 1013614, 1013924, 1012488, 1012488, 1013340, 1014098, 1014894, 1015466, 1014564, 1013066, 1012778, 1014400, 0, 0, 1015630, 1014324, 1013980, 1012838, 1013594, 1013188, 1014722, 1014168, 1014410, 1013072, 1012100, 1012358, 1011704, 1010494, 1012212, 1011434, 1011806, 0, 0, 0, 1010662, 0, 1010076, 1010626, 1010908, 1012420, 1012866], [948758, 953214, 952596, 951962, 951982, 952018, 952256, 951878, 952204, 951086, 951112, 951490, 951988, 952212, 951438, 951028, 951378, 952546, 951376, 951166, 952680, 0, 0, 951028, 952042, 953064, 953482, 953206, 952970, 954348, 954108, 952246, 951526, 951494, 952220, 951320, 951062, 951654, 951302, 950730, 0, 0, 0, 949882, 0, 949802, 949680, 949736, 950750, 951400], [885210, 893672, 893948, 892308, 892514, 892174, 891840, 892460, 892804, 892620, 892664, 892554, 891528, 891762, 891740, 892136, 892158, 891588, 891312, 890314, 891650, 0, 0, 892774, 892050, 891840, 893656, 894256, 894404, 893608, 893364, 892636, 893200, 894116, 893334, 892978, 890760, 892790, 892076, 892146, 0, 0, 0, 892130, 0, 892022, 893100, 893460, 893444, 892188], [830782, 837160, 837206, 837996, 836958, 837240, 837750, 837560, 837284, 837652, 837374, 837274, 837232, 836726, 835954, 834532, 834118, 835198, 836214, 836602, 836946, 0, 0, 837502, 837632, 837232, 837112, 836834, 836580, 838178, 838620, 836982, 836866, 836182, 835700, 834838, 835934, 835286, 834988, 835878, 0, 0, 0, 837092, 0, 835286, 834884, 835784, 835778, 836834], [778830, 785636, 785204, 784746, 785492, 785560, 785296, 784726, 785586, 786592, 786854, 787038, 786550, 787150, 787252, 786978, 786064, 785934, 785678, 786658, 786340, 0, 0, 786518, 785416, 786428, 786468, 785470, 786140, 786814, 785126, 784918, 786646, 787566, 787198, 786410, 786436, 786510, 786212, 784442, 0, 0, 0, 0, 0, 784556, 783746, 784968, 784724, 784570], [735340, 742726, 741662, 741832, 741808, 741920, 742180, 742338, 742198, 740570, 741186, 741510, 741210, 741362, 741646, 741890, 742256, 741534, 740380, 739798, 738484, 0, 0, 739512, 738896, 740202, 741070, 739924, 740788, 739942, 740254, 740776, 741832, 741734, 741280, 741080, 741546, 740556, 740994, 739958, 0, 0, 0, 0, 0, 738872, 738932, 738980, 738220, 737130], [698110, 702618, 702424, 702350, 702850, 702560, 701962, 701246, 700992, 701894, 701856, 701520, 702148, 701772, 701292, 701174, 701108, 701954, 702456, 701696, 700490, 0, 0, 700378, 700680, 701114, 701442, 699812, 700634, 701082, 701790, 700976, 701692, 701154, 699790, 699598, 699886, 700010, 699700, 700334, 0, 0, 0, 698682, 0, 697164, 697502, 697948, 697878, 697804], [664122, 668034, 667970, 667888, 666982, 667182, 667316, 667390, 667428, 668150, 668612, 668162, 667848, 667166, 667936, 668380, 668320, 667626, 667298, 666288, 664694, 0, 0, 663076, 663772, 663742, 663870, 662884, 662308, 662674, 663294, 663666, 663576, 662774, 662996, 662356, 663300, 662338, 662408, 662180, 0, 0, 0, 662320, 0, 659338, 659546, 659350, 660154, 660124], [627476, 641758, 640988, 640186, 639616, 638840, 638796, 639140, 639242, 638074, 638516, 636782, 636618, 636476, 636328, 635452, 635354, 635396, 634972, 634630, 634312, 0, 0, 630822, 630816, 630124, 630412, 629282, 628322, 627912, 627442, 627668, 627090, 626032, 626478, 625490, 623828, 624138, 622898, 622086, 0, 0, 0, 622512, 0, 622426, 622530, 620788, 620992, 620972], [592872, 619780, 619800, 618812, 618660, 618922, 619596, 618824, 618242, 617136, 615668, 614556, 613674, 611594, 610388, 609092, 607516, 607050, 607484, 605056, 602958, 0, 599978, 598986, 597818, 597298, 596382, 594836, 593518, 591130, 590310, 588030, 586816, 586182, 583512, 582062, 579780, 578692, 578692, 578218, 0, 0, 0, 577618, 578698, 578228, 578266, 578066, 577780, 577334], [559368, 613358, 613708, 613252, 611562, 611194, 609622, 607752, 606196, 603824, 602218, 600370, 598108, 594528, 592240, 589206, 586224, 583770, 581674, 579834, 576552, 0, 571148, 568646, 567292, 564714, 560630, 557898, 554584, 551348, 548514, 544138, 542620, 538438, 535020, 533306, 530110, 527574, 526472, 524460, 0, 0, 521234, 521020, 520522, 520798, 520422, 521234, 520774, 521058], [522498, 626654, 625576, 625020, 623638, 621586, 618696, 615958, 611914, 609360, 605150, 601010, 596576, 591516, 586876, 581384, 575404, 569922, 564970, 558086, 552592, 0, 542942, 537058, 531600, 525064, 518790, 513724, 508500, 501122, 494560, 488438, 483460, 476954, 469984, 464584, 460002, 456402, 453328, 449360, 0, 0, 444240, 444042, 443668, 442704, 441930, 441570, 441844, 442160], [483952, 671162, 669766, 668144, 664856, 661018, 656066, 650214, 644404, 637174, 630464, 622212, 613142, 605960, 596908, 588628, 578524, 568204, 559310, 549374, 538772, 0, 519552, 509206, 497716, 487606, 475954, 463476, 452786, 441028, 429210, 417064, 406124, 396086, 385182, 375550, 367266, 359308, 351662, 345590, 0, 0, 334986, 334288, 333210, 330864, 331302, 330372, 330310, 329914], [443328, 761130, 759282, 756460, 750920, 744000, 735454, 725342, 714374, 702024, 691020, 677970, 664168, 649064, 634206, 619372, 603218, 587606, 569348, 551972, 535644, 518124, 500904, 484150, 465192, 445944, 426602, 408394, 388010, 367276, 346740, 325950, 305596, 286678, 269834, 252752, 236926, 222072, 209762, 200414, 191868, 184196, 179824, 177302, 174076, 171186, 169992, 167972, 167350, 166498], [404930, 923270, 918462, 912748, 905494, 895564, 882400, 867110, 851056, 833234, 813542, 792508, 770692, 745910, 720136, 695054, 668526, 640730, 612320, 584504, 554910, 525520, 496574, 464900, 0, 0, 0, 0, 0, 269900, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16296, 26552, 0, 44112, 0, 0, 60500, 0, 0, 0], [361180, 1177054, 1172760, 1164156, 1152744, 1137770, 1118954, 1098426, 1074586, 1047646, 1017244, 983982, 949506, 912170, 871964, 831218, 786840, 742702, 696794, 650140, 603604, 554172, 506180, 457128, 0, 0, 0, 0, 198144, 143706, 0, 0, 11988, 0, 0, 0, 0, 0, 0, 0, 311052, 330466, 347408, 0, 0, 0, 0, 0, 0, 0], [313032, 1561056, 1552824, 1541876, 1526072, 1505108, 1478738, 1447954, 1412510, 1372504, 1328768, 1281456, 1227336, 1169214, 1108758, 1045144, 978512, 909868, 840904, 768276, 693964, 617780, 542148, 465062, 0, 0, 0, 0, 0, 26194, 0, 0, 263216, 0, 0, 0, 0, 0, 0, 0, 722120, 752864, 778708, 0, 0, 0, 0, 0, 0, 0], [264420, 2117936, 2107950, 2092340, 2069224, 2040818, 2003512, 1960388, 1908880, 1852458, 1788258, 1716738, 1641788, 1560002, 1471592, 1376930, 1276752, 1174326, 1067036, 956908, 844676, 730118, 613296, 493766, 0, 0, 0, 0, 0, 247774, 0, 0, 605860, 0, 0, 0, 0, 0, 0, 0, 1281204, 1327852, 1367692, 0, 0, 0, 0, 0, 0, 0], [213242, 2883080, 2870298, 2850730, 2820148, 2782066, 2734728, 2678652, 2611822, 2535222, 2449214, 2352706, 2247488, 2131688, 2007712, 1873490, 1731026, 1580844, 1424612, 1262564, 1094386, 922330, 746208, 564978, 0, 0, 0, 0, 0, 549178, 728234, 0, 1071976, 0, 0, 0, 0, 0, 1874058, 0, 2042066, 0, 2167870, 0, 0, 0, 0, 0, 0, 0], [160952, 3886910, 3873100, 3849474, 3814392, 3770254, 3713738, 3647574, 3567052, 3472990, 3365372, 3242954, 3105430, 2951726, 2783878, 2600564, 2403758, 2191918, 1967154, 1731684, 1485280, 1230246, 968276, 697902, 0, 145218, 134382, 0, 0, 963948, 0, 0, 1727356, 0, 0, 0, 0, 0, 2839480, 0, 3064370, 0, 3225936, 3289664, 0, 0, 3417458, 0, 0, 0], [108230, 5134820, 5121718, 5097722, 5063860, 5018188, 4962362, 4893758, 4809388, 4709504, 4592092, 4453678, 4295144, 4113988, 3906276, 3672898, 3413752, 3129530, 2820272, 2484348, 2127834, 1750984, 1355882, 945612, 524718, 96294, 0, 0, 0, 1589432, 1979330, 0, 2687568, 0, 0, 0, 0, 0, 4131024, 4273242, 4393114, 4492766, 4575998, 0, 4701354, 0, 0, 0, 0, 0], [58534, 6597568, 6587942, 6571650, 6546344, 6510996, 6466034, 6408672, 6339034, 6254960, 6151462, 6030340, 5883972, 5708692, 5499348, 5252048, 4963052, 4625342, 4235058, 3789900, 3288424, 2730606, 2117328, 1459834, 768460, 52664, 664592, 0, 0, 0, 0, 0, 4198240, 0, 0, 0, 0, 0, 5806970, 5939558, 6049032, 0, 6211624, 0, 0, 0, 0, 0, 0, 0], [15030, 8234560, 8230258, 8221278, 8208778, 8190690, 8166516, 8136074, 8098058, 8051132, 7992310, 7921104, 7835258, 7729378, 7596536, 7430206, 7225846, 6967596, 6639188, 6219032, 5675618, 4967970, 4053584, 2903310, 1534108, 22796, 0, 0, 0, 0, 0, 0, 6653398, 0, 0, 0, 0, 0, 0, 7893944, 7959434, 8013074, 8055680, 0, 0, 8140122, 0, 0, 0, 0]]





plt.figure("heatmap")
plt.imshow(data, cmap='plasma')

plt.title("$\int_0^\infty P(z)dz - \int_{-\infty}^0 P(z)dz$ after 10 iterations")
plt.xlabel("$ \\theta$")
plt.ylabel("$\phi$")
plt.xlim([1,49])
plt.xticks([10,20,30,40,49],["$\\frac{\pi}{10}$","$\\frac{\pi}{5}$","$\\frac{3\pi}{10}$","$\\frac{2\pi}{5}$","$\\frac{\pi}{2}$"])
plt.yticks([10,20,30,40,49],["$\\frac{\pi}{10}$","$\\frac{\pi}{5}$","$\\frac{3\pi}{10}$","$\\frac{2\pi}{5}$","$\\frac{\pi}{2}$"])
plt.ylim([1,49.5])

cbar = plt.colorbar(ticks=[0,2000000,4000000,6000000,8000000,10000000])
cbar .ax.set_yticklabels(["0","0.2","0.4","0.6","0.8","1"])

plt.savefig("heatmap.png", dpi=300,bbox_inches="tight")
plt.show()
#fig = px.imshow(data)
#fig.show()