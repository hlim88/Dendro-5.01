//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
//Automatic code generated by sympy based swsh module in Dendro-GR 
//python module to generate code for far-field energy ratiation extraction (Gravitational Waves).
//(c) 2016 University of Utah, All rights reserved.
//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ifndef DENDRO_SWSH_H
#define DENDRO_SWSH_H
namespace swsh { 
// spin weighted sperical harmonic basis evaluated at the 025 lebedev quadrature points 
static double m2Y2_0_REAL [] = { 
0.38627420202318940579, 0.38627420202318940579, 0.38627420202318940579, 0.38627420202318940579, 0.00000000000000000000, 0.00000000000000000000, 0.25751613468212647406, 0.25751613468212641855, 0.25751613468212647406, 0.25751613468212641855, 0.25751613468212647406, 0.25751613468212641855, 0.25751613468212647406, 0.25751613468212641855, 0.15588841297716954371, 0.15588841297716957146, 0.15588841297716954371, 0.15588841297716957146, 0.15588841297716954371, 0.15588841297716957146, 0.15588841297716954371, 0.15588841297716957146, 0.30832999553460477271, 0.30832999553460477271, 0.30832999553460477271, 0.30832999553460477271, 0.30832999553460477271, 0.30832999553460477271, 0.30832999553460477271, 0.30832999553460477271, 0.30832999553460477271, 0.30832999553460477271, 0.30832999553460477271, 0.30832999553460477271, 0.30832999553460477271, 0.30832999553460477271, 0.30832999553460477271, 0.30832999553460477271, 0.04907624866524754742, 0.04907624866524769314, 0.04907624866524754742, 0.04907624866524769314, 0.04907624866524754742, 0.04907624866524769314, 0.04907624866524754742, 0.04907624866524769314, 0.36173607769056576045, 0.36173607769056576045, 0.36173607769056581596, 0.36173607769056581596, 0.36173607769056576045, 0.36173607769056576045, 0.36173607769056581596, 0.36173607769056581596, 0.36173607769056576045, 0.36173607769056576045, 0.36173607769056581596, 0.36173607769056581596, 0.36173607769056576045, 0.36173607769056576045, 0.36173607769056581596, 0.36173607769056581596, 0.37659432956711380580, 0.37659432956711375029, 0.37659432956711380580, 0.37659432956711375029, 0.37659432956711380580, 0.37659432956711375029, 0.37659432956711380580, 0.37659432956711375029, 0.19797703723963275269, 0.19797703723963275269, 0.19797703723963266942, 0.19797703723963266942, 0.19797703723963275269, 0.19797703723963275269, 0.19797703723963266942, 0.19797703723963266942, 0.19797703723963275269, 0.19797703723963275269, 0.19797703723963266942, 0.19797703723963266942, 0.19797703723963275269, 0.19797703723963275269, 0.19797703723963266942, 0.19797703723963266942, 0.33523894248345925684, 0.33523894248345936786, 0.33523894248345925684, 0.33523894248345936786, 0.33523894248345925684, 0.33523894248345936786, 0.33523894248345925684, 0.33523894248345936786, 0.21865473078145991614, 0.21865473078145991614, 0.21865473078145988839, 0.21865473078145988839, 0.21865473078145991614, 0.21865473078145991614, 0.21865473078145988839, 0.21865473078145988839, 0.21865473078145991614, 0.21865473078145991614, 0.21865473078145988839, 0.21865473078145988839, 0.21865473078145991614, 0.21865473078145991614, 0.21865473078145988839, 0.21865473078145988839, 0.00126001393515946824, 0.00126001393515947106, 0.00126001393515946824, 0.00126001393515947106, 0.00126001393515946824, 0.00126001393515947106, 0.00126001393515946824, 0.00126001393515947106, 0.38564419505560981749, 0.38564419505560981749, 0.38564419505560976198, 0.38564419505560976198, 0.38564419505560981749, 0.38564419505560981749, 0.38564419505560976198, 0.38564419505560976198, 0.38564419505560981749, 0.38564419505560981749, 0.38564419505560976198, 0.38564419505560976198, 0.38564419505560981749, 0.38564419505560981749, 0.38564419505560976198, 0.38564419505560976198, 0.38627420202318940579, 0.38627420202318940579, 0.38627420202318940579, 0.38627420202318940579, 0.38627420202318940579, 0.38627420202318940579, 0.38627420202318940579, 0.38627420202318940579, 0.13101315898097995460, 0.13101315898097998236, 0.13101315898097995460, 0.13101315898097998236, 0.25526104304220958996, 0.25526104304220958996, 0.25526104304220958996, 0.25526104304220958996, 0.13101315898097995460, 0.13101315898097998236, 0.13101315898097995460, 0.13101315898097998236, 0.25526104304220958996, 0.25526104304220958996, 0.25526104304220958996, 0.25526104304220958996, 0.38627420202318940579, 0.38627420202318940579, 0.38627420202318940579, 0.38627420202318940579, 0.38627420202318940579, 0.38627420202318940579, 0.38627420202318940579, 0.38627420202318940579, 0.04856720746496115404, 0.04856720746496121649, 0.04856720746496115404, 0.04856720746496121649, 0.33770699455822833501, 0.33770699455822833501, 0.33770699455822833501, 0.33770699455822833501, 0.04856720746496115404, 0.04856720746496121649, 0.04856720746496115404, 0.04856720746496121649, 0.33770699455822833501, 0.33770699455822833501, 0.33770699455822833501, 0.33770699455822833501, 0.11135413049493768367, 0.11135413049493786408, 0.11135413049493768367, 0.11135413049493786408, 0.11135413049493768367, 0.11135413049493786408, 0.11135413049493768367, 0.11135413049493786408, 0.29486267532652921108, 0.29486267532652932211, 0.29486267532652921108, 0.29486267532652932211, 0.29486267532652921108, 0.29486267532652932211, 0.29486267532652921108, 0.29486267532652932211, 0.11135413049493768367, 0.11135413049493786408, 0.11135413049493768367, 0.11135413049493786408, 0.11135413049493768367, 0.11135413049493786408, 0.11135413049493768367, 0.11135413049493786408, 0.36633159822491201396, 0.36633159822491206947, 0.36633159822491201396, 0.36633159822491206947, 0.36633159822491201396, 0.36633159822491206947, 0.36633159822491201396, 0.36633159822491206947, 0.29486267532652921108, 0.29486267532652932211, 0.29486267532652921108, 0.29486267532652932211, 0.29486267532652921108, 0.29486267532652932211, 0.29486267532652921108, 0.29486267532652932211, 0.36633159822491201396, 0.36633159822491206947, 0.36633159822491201396, 0.36633159822491206947, 0.36633159822491201396, 0.36633159822491206947, 0.36633159822491201396, 0.36633159822491206947
};
static double m2Y2_0_IMAG [] = { 
0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000
};
// spin weighted sperical harmonic basis evaluated at the 025 lebedev quadrature points 
static double m2Y2_1_REAL [] = { 
-0.31539156525251987873, 0.31539156525251987873, -0.00000000000000005794, 0.00000000000000001931, 0.00000000000000000000, 0.00000000000000000000, -0.28722192684951997066, -0.07696088334783987572, -0.28722192684951991515, -0.07696088334783986185, 0.28722192684951985964, 0.07696088334783984797, 0.28722192684951991515, 0.07696088334783987572, -0.25108961253798778301, -0.03226098849053212619, -0.25108961253798772750, -0.03226098849053211232, 0.25108961253798772750, 0.03226098849053211232, 0.25108961253798778301, 0.03226098849053211925, -0.20531647861531393540, -0.20531647861531401866, -0.07803412241320574483, -0.07803412241320577258, 0.20531647861531412969, 0.20531647861531404642, 0.07803412241320581422, 0.07803412241320578646, -0.35298782840188996923, 0.35298782840188991372, -0.13415920435443431935, 0.13415920435443429160, -0.35298782840189002474, 0.35298782840189008025, -0.13415920435443434711, 0.13415920435443434711, -0.15376260475030459340, -0.00522120487176613510, -0.15376260475030453789, -0.00522120487176613423, 0.15376260475030453789, 0.00522120487176613336, 0.15376260475030456565, 0.00522120487176613423, -0.09952719943100231392, -0.09952719943100221678, -0.05945661019106873030, -0.05945661019106867479, 0.09952719943100217515, 0.09952719943100227229, 0.05945661019106865397, 0.05945661019106870254, -0.36894664665270771042, 0.36894664665270771042, -0.22040524677416917232, 0.22040524677416917232, -0.36894664665270759940, 0.36894664665270771042, -0.22040524677416911681, 0.22040524677416917232, -0.25506213824759410569, -0.18534475564333935393, -0.25506213824759405018, -0.18534475564333932618, 0.25506213824759399467, 0.18534475564333929842, 0.25506213824759405018, 0.18534475564333932618, -0.37394743818962561388, -0.37394743818962561388, -0.06645945570130787350, -0.06645945570130788738, 0.37394743818962561388, 0.37394743818962561388, 0.06645945570130787350, 0.06645945570130788738, -0.08478587166617967963, 0.08478587166617954085, -0.01506848906192498512, 0.01506848906192496083, -0.08478587166617983228, 0.08478587166617988780, -0.01506848906192501288, 0.01506848906192502155, -0.28327942900064190246, -0.13224298113692009582, -0.28327942900064179144, -0.13224298113692006806, 0.28327942900064179144, 0.13224298113692006806, 0.28327942900064184695, 0.13224298113692009582, -0.34462193023456644259, -0.34462193023456649810, -0.07090047990299537528, -0.07090047990299538916, 0.34462193023456644259, 0.34462193023456655361, 0.07090047990299537528, 0.07090047990299538916, -0.19015854305724744222, 0.19015854305724730344, -0.03912209519352562170, 0.03912209519352559395, -0.19015854305724735895, 0.19015854305724738671, -0.03912209519352560089, 0.03912209519352561476, -0.02545366340571406441, -0.00002079118031708876, -0.02545366340571405747, -0.00002079118031708876, 0.02545366340571405400, 0.00002079118031708876, 0.02545366340571406094, 0.00002079118031708876, -0.01325162582800474588, -0.01325162582800503037, -0.01222282875802578606, -0.01222282875802604626, 0.01325162582800491762, 0.01325162582800514313, 0.01222282875802594218, 0.01222282875802615035, -0.32759318265681308668, 0.32759318265681308668, -0.30216031043141639012, 0.30216031043141639012, -0.32759318265681308668, 0.32759318265681308668, -0.30216031043141633461, 0.30216031043141639012, -0.18367907418450504786, -0.18367907418450488133, 0.18367907418450474255, 0.18367907418450493684, -0.25638610948949180912, -0.25638610948949175361, 0.25638610948949180912, 0.25638610948949175361, -0.33299430139133234796, -0.03436384697767747021, 0.33299430139133234796, 0.03436384697767747021, -0.40570133669631924800, -0.10707088228266434249, 0.40570133669631924800, 0.10707088228266434249, -0.00000000000000006117, -0.00000000000000000631, 0.00000000000000002039, 0.00000000000000000210, -0.00000000000000007453, -0.00000000000000001967, 0.00000000000000002484, 0.00000000000000000656, -0.11183398203892121192, -0.11183398203892094824, 0.11183398203892083722, 0.11183398203892104539, -0.29489828736998935366, -0.29489828736998929815, 0.29489828736998935366, 0.29489828736998935366, -0.21640129900771382099, -0.00726666507012819780, 0.21640129900771382099, 0.00726666507012819780, -0.39946560433878225416, -0.19033097040119648091, 0.39946560433878225416, 0.19033097040119648091, -0.00000000000000003975, -0.00000000000000000133, 0.00000000000000001325, 0.00000000000000000044, -0.00000000000000007338, -0.00000000000000003496, 0.00000000000000002446, 0.00000000000000001165, -0.13211996672933737362, -0.01120542872631084819, -0.13211996672933731811, -0.01120542872631084125, 0.13211996672933726260, 0.01120542872631083778, 0.13211996672933734587, 0.01120542872631084472, -0.10652417464726557372, -0.03680122080838230808, -0.10652417464726583740, -0.03680122080838239135, 0.10652417464726578189, 0.03680122080838237747, 0.10652417464726579577, 0.03680122080838237747, -0.28286421347339624210, -0.02399042977200736837, -0.28286421347339624210, -0.02399042977200736490, 0.28286421347339624210, 0.02399042977200736490, 0.28286421347339624210, 0.02399042977200736837, -0.18828879854214325418, -0.11856584470325989833, -0.18828879854214330969, -0.11856584470325993996, 0.18828879854214342071, 0.11856584470326000935, 0.18828879854214342071, 0.11856584470326000935, -0.39551273473480352827, -0.13663895103341475168, -0.39551273473480341725, -0.13663895103341472392, 0.39551273473480341725, 0.13663895103341472392, 0.39551273473480347276, 0.13663895103341475168, -0.32653311188562222922, -0.20561857388259577317, -0.32653311188562217371, -0.20561857388259574542, 0.32653311188562217371, 0.20561857388259574542, 0.32653311188562234024, 0.20561857388259582868
};
static double m2Y2_1_IMAG [] = { 
0.00000000000000003862, -0.00000000000000007725, -0.31539156525251987873, 0.31539156525251987873, 0.00000000000000000000, 0.00000000000000000000, -0.28722192684951985964, -0.07696088334783986185, 0.28722192684951991515, 0.07696088334783987572, -0.28722192684951997066, -0.07696088334783988960, 0.28722192684951991515, 0.07696088334783986185, -0.25108961253798772750, -0.03226098849053211232, 0.25108961253798778301, 0.03226098849053211925, -0.25108961253798783853, -0.03226098849053212619, 0.25108961253798772750, 0.03226098849053211925, -0.35298782840189013577, 0.35298782840189008025, -0.13415920435443437486, 0.13415920435443434711, -0.35298782840189002474, 0.35298782840189002474, -0.13415920435443431935, 0.13415920435443434711, -0.20531647861531415744, -0.20531647861531426846, -0.07803412241320582809, -0.07803412241320586973, 0.20531647861531407417, 0.20531647861531404642, 0.07803412241320580034, 0.07803412241320578646, -0.15376260475030453789, -0.00522120487176613336, 0.15376260475030456565, 0.00522120487176613510, -0.15376260475030459340, -0.00522120487176613510, 0.15376260475030456565, 0.00522120487176613423, -0.36894664665270759940, 0.36894664665270765491, -0.22040524677416911681, 0.22040524677416914456, -0.36894664665270765491, 0.36894664665270759940, -0.22040524677416914456, 0.22040524677416914456, -0.09952719943100203637, -0.09952719943100200861, -0.05945661019106856376, -0.05945661019106854989, 0.09952719943100228617, 0.09952719943100202249, 0.05945661019106871642, 0.05945661019106855683, -0.25506213824759399467, -0.18534475564333929842, 0.25506213824759410569, 0.18534475564333935393, -0.25506213824759410569, -0.18534475564333938169, 0.25506213824759405018, 0.18534475564333932618, -0.08478587166617990167, 0.08478587166617981841, -0.01506848906192502328, 0.01506848906192501114, -0.08478587166618004045, 0.08478587166617983228, -0.01506848906192504930, 0.01506848906192501114, -0.37394743818962566939, -0.37394743818962572490, -0.06645945570130788738, -0.06645945570130788738, 0.37394743818962561388, 0.37394743818962561388, 0.06645945570130788738, 0.06645945570130787350, -0.28327942900064179144, -0.13224298113692006806, 0.28327942900064184695, 0.13224298113692009582, -0.28327942900064190246, -0.13224298113692012357, 0.28327942900064184695, 0.13224298113692009582, -0.19015854305724749773, 0.19015854305724741447, -0.03912209519352562864, 0.03912209519352561476, -0.19015854305724746998, 0.19015854305724733120, -0.03912209519352562170, 0.03912209519352560089, -0.34462193023456649810, -0.34462193023456655361, -0.07090047990299538916, -0.07090047990299540304, 0.34462193023456649810, 0.34462193023456649810, 0.07090047990299538916, 0.07090047990299538916, -0.02545366340571405400, -0.00002079118031708876, 0.02545366340571406094, 0.00002079118031708876, -0.02545366340571406441, -0.00002079118031708876, 0.02545366340571405747, 0.00002079118031708876, -0.32759318265681308668, 0.32759318265681308668, -0.30216031043141633461, 0.30216031043141633461, -0.32759318265681308668, 0.32759318265681308668, -0.30216031043141633461, 0.30216031043141633461, -0.01325162582800464700, -0.01325162582800447526, -0.01222282875802569238, -0.01222282875802553626, 0.01325162582800472680, 0.01325162582800449608, 0.01222282875802576697, 0.01222282875802555360, -0.25638610948949164259, 0.25638610948949180912, -0.25638610948949192014, 0.25638610948949175361, -0.18367907418450485357, 0.18367907418450493684, -0.18367907418450485357, 0.18367907418450493684, 0.00000000000000004078, 0.00000000000000000421, -0.00000000000000008156, -0.00000000000000000842, 0.00000000000000004968, 0.00000000000000001311, -0.00000000000000009937, -0.00000000000000002622, -0.33299430139133234796, -0.03436384697767747021, 0.33299430139133234796, 0.03436384697767747021, -0.40570133669631924800, -0.10707088228266434249, 0.40570133669631924800, 0.10707088228266434249, -0.29489828736998924263, 0.29489828736998935366, -0.29489828736998940917, 0.29489828736998929815, -0.11183398203892098988, 0.11183398203892105927, -0.11183398203892097600, 0.11183398203892103151, 0.00000000000000002650, 0.00000000000000000089, -0.00000000000000005300, -0.00000000000000000178, 0.00000000000000004892, 0.00000000000000002331, -0.00000000000000009784, -0.00000000000000004662, -0.21640129900771382099, -0.00726666507012819780, 0.21640129900771382099, 0.00726666507012819780, -0.39946560433878225416, -0.19033097040119648091, 0.39946560433878225416, 0.19033097040119648091, -0.28286421347339624210, -0.02399042977200736490, 0.28286421347339624210, 0.02399042977200736837, -0.28286421347339629762, -0.02399042977200736837, 0.28286421347339624210, 0.02399042977200736490, -0.39551273473480352827, -0.13663895103341475168, 0.39551273473480347276, 0.13663895103341472392, -0.39551273473480347276, -0.13663895103341475168, 0.39551273473480347276, 0.13663895103341472392, -0.13211996672933729036, -0.01120542872631083951, 0.13211996672933734587, 0.01120542872631084645, -0.13211996672933740138, -0.01120542872631084819, 0.13211996672933731811, 0.01120542872631084298, -0.32653311188562234024, -0.20561857388259585644, 0.32653311188562228473, 0.20561857388259582868, -0.32653311188562222922, -0.20561857388259580093, 0.32653311188562222922, 0.20561857388259580093, -0.10652417464726562923, -0.03680122080838232196, 0.10652417464726590679, 0.03680122080838241910, -0.10652417464726594842, -0.03680122080838243298, 0.10652417464726569862, 0.03680122080838234971, -0.18828879854214344847, -0.11856584470326002323, 0.18828879854214353173, 0.11856584470326006486, -0.18828879854214355949, -0.11856584470326009262, 0.18828879854214330969, 0.11856584470325992609
};
// spin weighted sperical harmonic basis evaluated at the 025 lebedev quadrature points 
static double m2Y2_2_REAL [] = { 
0.15769578262625993936, 0.15769578262625993936, -0.15769578262625993936, -0.15769578262625993936, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000012012, 0.00000000000000000862, -0.00000000000000007207, -0.00000000000000000517, -0.00000000000000016817, -0.00000000000000001207, 0.00000000000000002402, 0.00000000000000000172, 0.00000000000000015165, 0.00000000000000000250, -0.00000000000000009099, -0.00000000000000000150, -0.00000000000000021231, -0.00000000000000000350, 0.00000000000000003033, 0.00000000000000000050, -0.16374463519550841450, -0.16374463519550830348, -0.02365312718476659040, -0.02365312718476657305, -0.16374463519550810919, -0.16374463519550822022, -0.02365312718476654877, -0.02365312718476656265, 0.16374463519550808144, 0.16374463519550785939, 0.02365312718476654183, 0.02365312718476651060, 0.16374463519550819246, 0.16374463519550824797, 0.02365312718476655918, 0.02365312718476656959, 0.00000000000000018065, 0.00000000000000000021, -0.00000000000000010839, -0.00000000000000000012, -0.00000000000000025290, -0.00000000000000000029, 0.00000000000000003613, 0.00000000000000000004, -0.21366731842088110271, -0.21366731842088115823, -0.07625262909384315779, -0.07625262909384317167, -0.21366731842088118598, -0.21366731842088113047, -0.07625262909384318555, -0.07625262909384317167, 0.21366731842088129700, 0.21366731842088129700, 0.07625262909384322718, 0.07625262909384322718, 0.21366731842088113047, 0.21366731842088129700, 0.07625262909384315779, 0.07625262909384322718, 0.00000000000000006478, 0.00000000000000003420, -0.00000000000000003887, -0.00000000000000002052, -0.00000000000000009069, -0.00000000000000004789, 0.00000000000000001296, 0.00000000000000000684, 0.41030013381929480998, 0.41030013381929486549, 0.01295967451168810944, 0.01295967451168811117, 0.41030013381929464344, 0.41030013381929486549, 0.01295967451168810423, 0.01295967451168811117, -0.41030013381929503202, -0.41030013381929519856, -0.01295967451168811638, -0.01295967451168812158, -0.41030013381929486549, -0.41030013381929480998, -0.01295967451168811117, -0.01295967451168810944, 0.00000000000000008976, 0.00000000000000001956, -0.00000000000000005385, -0.00000000000000001174, -0.00000000000000012566, -0.00000000000000002739, 0.00000000000000001795, 0.00000000000000000391, 0.23134379264896862138, 0.23134379264896876016, 0.00979197688322411969, 0.00979197688322412663, 0.23134379264896867689, 0.23134379264896892670, 0.00979197688322412142, 0.00979197688322413357, -0.23134379264896870465, -0.23134379264896898221, -0.00979197688322412316, -0.00979197688322413530, -0.23134379264896889894, -0.23134379264896881567, -0.00979197688322413183, -0.00979197688322412836, 0.00000000000000019281, 0.00000000000000000000, -0.00000000000000011568, -0.00000000000000000000, -0.00000000000000026993, -0.00000000000000000000, 0.00000000000000003856, 0.00000000000000000000, -0.17013251366087328575, -0.17013251366087325800, -0.14474129167157406828, -0.14474129167157406828, -0.17013251366087328575, -0.17013251366087325800, -0.14474129167157406828, -0.14474129167157404052, 0.17013251366087331351, 0.17013251366087331351, 0.14474129167157406828, 0.14474129167157409603, 0.17013251366087328575, 0.17013251366087331351, 0.14474129167157406828, 0.14474129167157409603, -0.05072398626174904473, -0.05072398626174925290, -0.05072398626174941944, -0.05072398626174918351, 0.05072398626174927372, 0.05072398626174919739, 0.05072398626174929454, 0.05072398626174919739, 0.51829177655975622319, 0.00551955758077276974, 0.51829177655975622319, 0.00551955758077276974, 0.39486075499302036862, 0.02750260662401041514, 0.39486075499302036862, 0.02750260662401041514, -0.51829177655975622319, -0.00551955758077276974, -0.51829177655975622319, -0.00551955758077276974, -0.39486075499302036862, -0.02750260662401041514, -0.39486075499302036862, -0.02750260662401041514, -0.11804082378591149172, -0.11804082378591168601, -0.11804082378591175539, -0.11804082378591160274, 0.11804082378591164437, 0.11804082378591158886, 0.11804082378591165825, 0.11804082378591161662, 0.59046237320233529999, 0.00066579846235640880, 0.59046237320233529999, 0.00066579846235640880, 0.28935724408535512531, 0.06568928000751306229, 0.28935724408535512531, 0.06568928000751306229, -0.59046237320233529999, -0.00066579846235640880, -0.59046237320233529999, -0.00066579846235640880, -0.28935724408535512531, -0.06568928000751306229, -0.28935724408535512531, -0.06568928000751306229, -0.34401833519866614752, -0.00247458076970680793, -0.34401833519866631406, -0.00247458076970680924, -0.34401833519866642508, -0.00247458076970681010, -0.34401833519866625855, -0.00247458076970680837, -0.30130901842936153034, -0.03596170104083799834, -0.30130901842936130830, -0.03596170104083797753, -0.30130901842936136381, -0.03596170104083797753, -0.30130901842936136381, -0.03596170104083797753, 0.34401833519866636957, 0.00247458076970680967, 0.34401833519866620303, 0.00247458076970680793, 0.34401833519866609201, 0.00247458076970680750, 0.34401833519866631406, 0.00247458076970680880, -0.11897221595382902193, -0.04717542068712516495, -0.11897221595382893866, -0.04717542068712513026, -0.11897221595382881376, -0.04717542068712508169, -0.11897221595382879988, -0.04717542068712507475, 0.30130901842936147483, 0.03596170104083799834, 0.30130901842936125279, 0.03596170104083796365, 0.30130901842936119728, 0.03596170104083796365, 0.30130901842936141932, 0.03596170104083798447, 0.11897221595382878601, 0.04717542068712506781, 0.11897221595382868886, 0.04717542068712503311, 0.11897221595382863335, 0.04717542068712501230, 0.11897221595382896642, 0.04717542068712514414
};
static double m2Y2_2_IMAG [] = { 
-0.00000000000000003862, -0.00000000000000007725, 0.00000000000000005794, 0.00000000000000001931, 0.00000000000000000000, 0.00000000000000000000, 0.39235244860036005976, 0.02816963840300017521, -0.39235244860036005976, -0.02816963840300017521, -0.39235244860036005976, -0.02816963840300017521, 0.39235244860036005976, 0.02816963840300017521, 0.49532390352962796243, 0.00817687077330366170, -0.49532390352962796243, -0.00817687077330366170, -0.49532390352962796243, -0.00817687077330366170, 0.49532390352962796243, 0.00817687077330366170, 0.28788125707206407844, -0.28788125707206413395, 0.04158482492880352654, -0.04158482492880354042, -0.28788125707206424497, 0.28788125707206413395, -0.04158482492880355430, 0.04158482492880354042, 0.28788125707206424497, -0.28788125707206435600, 0.04158482492880355430, -0.04158482492880357512, -0.28788125707206418946, 0.28788125707206413395, -0.04158482492880354736, 0.04158482492880354042, 0.59003221734599142767, 0.00068032391911488591, -0.59003221734599142767, -0.00068032391911488591, -0.59003221734599142767, -0.00068032391911488591, 0.59003221734599142767, 0.00068032391911488591, 0.12432519682809260730, -0.12432519682809251016, 0.04436861561615932509, -0.04436861561615928345, -0.12432519682809245465, 0.12432519682809256567, -0.04436861561615926264, 0.04436861561615930427, 0.12432519682809228811, -0.12432519682809226036, 0.04436861561615920713, -0.04436861561615919325, -0.12432519682809259343, 0.12432519682809227424, -0.04436861561615931121, 0.04436861561615920019, 0.21157475437241354821, 0.11172039364430877417, -0.21157475437241354821, -0.11172039364430877417, -0.21157475437241354821, -0.11172039364430877417, 0.21157475437241354821, 0.11172039364430877417, 0.19613940820789782515, -0.19613940820789765862, 0.00619522802887749729, -0.00619522802887749208, -0.19613940820789813047, 0.19613940820789768638, -0.00619522802887750683, 0.00619522802887749208, 0.19613940820789735331, -0.19613940820789704800, 0.00619522802887748168, -0.00619522802887747213, -0.19613940820789768638, 0.19613940820789779740, -0.00619522802887749295, 0.00619522802887749642, 0.29317115921212094642, 0.06389052096134795189, -0.29317115921212094642, -0.06389052096134795189, -0.29317115921212094642, -0.06389052096134795189, 0.29317115921212094642, 0.06389052096134795189, 0.36706713600375251438, -0.36706713600375240336, 0.01553667323070945065, -0.01553667323070944718, -0.36706713600375245887, 0.36706713600375229234, -0.01553667323070944892, 0.01553667323070944198, 0.36706713600375245887, -0.36706713600375229234, 0.01553667323070944892, -0.01553667323070944024, -0.36706713600375234785, 0.36706713600375240336, -0.01553667323070944371, 0.01553667323070944545, 0.62975391326164520400, 0.00000042017341572184, -0.62975391326164520400, -0.00000042017341572184, -0.62975391326164520400, -0.00000042017341572184, 0.62975391326164520400, 0.00000042017341572184, 0.01378678017726082429, -0.01378678017726111919, 0.01172918878296517683, -0.01172918878296542664, -0.01378678017726100123, 0.01378678017726123542, -0.01172918878296532776, 0.01172918878296552725, 0.01378678017726072021, -0.01378678017726054326, 0.01172918878296508836, -0.01172918878296493744, -0.01378678017726080347, 0.01378678017726056408, -0.01172918878296515949, 0.01172918878296495479, 0.14931522720682752214, -0.14931522720682743888, -0.14931522720682738337, 0.14931522720682746663, 0.14931522720682743888, -0.14931522720682746663, -0.14931522720682743888, 0.14931522720682746663, -0.00000000000000012694, -0.00000000000000000135, -0.00000000000000025389, -0.00000000000000000270, -0.00000000000000009671, -0.00000000000000000674, -0.00000000000000019343, -0.00000000000000001347, 0.00000000000000019042, 0.00000000000000000203, 0.00000000000000006347, 0.00000000000000000068, 0.00000000000000014507, 0.00000000000000001010, 0.00000000000000004836, 0.00000000000000000337, 0.10456731696879298377, -0.10456731696879277560, -0.10456731696879269233, 0.10456731696879285887, 0.10456731696879281723, -0.10456731696879287274, -0.10456731696879278948, 0.10456731696879284499, -0.00000000000000014462, -0.00000000000000000016, -0.00000000000000028924, -0.00000000000000000033, -0.00000000000000007087, -0.00000000000000001609, -0.00000000000000014174, -0.00000000000000003218, 0.00000000000000021693, 0.00000000000000000024, 0.00000000000000007231, 0.00000000000000000008, 0.00000000000000010631, 0.00000000000000002413, 0.00000000000000003544, 0.00000000000000000804, 0.41104158379209204677, 0.00295669007936552640, -0.41104158379209188023, -0.00295669007936552510, -0.41104158379209182472, -0.00295669007936552467, 0.41104158379209199126, 0.00295669007936552597, 0.17499856092673465868, 0.02088635104063267522, -0.17499856092673504726, -0.02088635104063272380, -0.17499856092673496399, -0.02088635104063271339, 0.17499856092673499175, 0.02088635104063271686, 0.41104158379209182472, 0.00295669007936552510, -0.41104158379209199126, -0.00295669007936552597, -0.41104158379209210228, -0.00295669007936552684, 0.41104158379209193575, 0.00295669007936552554, 0.20555259328343897240, 0.08150667770396927136, -0.20555259328343902792, -0.08150667770396928524, -0.20555259328343911118, -0.08150667770396931300, 0.20555259328343911118, 0.08150667770396931300, 0.17499856092673474195, 0.02088635104063268563, -0.17499856092673515828, -0.02088635104063273767, -0.17499856092673521379, -0.02088635104063274461, 0.17499856092673485297, 0.02088635104063269951, 0.20555259328343911118, 0.08150667770396932688, -0.20555259328343916669, -0.08150667770396934075, -0.20555259328343919445, -0.08150667770396935463, 0.20555259328343902792, 0.08150667770396928524
};
}// end of namespace swsh
#endif
