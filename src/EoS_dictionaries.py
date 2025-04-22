# -*- coding: utf-8 -*-
"""
Created on Fri Apr  4 15:31:43 2025

@author: sinjac
"""


Pt = {
			"name" : "Pt",
			"Uo" : 0.0,
			"n" : 1.0,
			"Vo" : 60.3836,#check
			"Ko" : 277.93567,#check
			"kk" : 5.35,
			"QE1" : 177.0,
			"mE1" : 1.5,
			"QE2" : 143.0,
			"mE2" : 1.5,
			"delta" : 0.167,
			"t" : -0.343,
			"ao" : 0.0,
			"eo" : 80.6,
			"m" : 0.0,
			"g" : 0.06,
			"co" : 3.78,
			"c2" : -0.25,
			"To" : 298.15,
			"Z" : 78.0,
			"Therm_alpha_298" : 2.5718816e-05,
			"Therm_diff_temp" : -0.01312,
			"Therm_diff_alpha" : 5.779e-09,
			"Therm_diff_Kprime" : 0.000303
			}

Ag = {
			"name" : "Ag",
			"Uo" : 0.0,
			"n" : 1.0,
			"Vo" : 68.08211,#check
			"Ko" : 1000,#check
			"kk" : 6.15,
			"QE1" : 115.0,
			"mE1" : 1.5,
			"QE2" : 199.0,
			"mE2" : 1.5,
			"delta" : 0.178,
			"t" : 2.21,
			"ao" : 0.0,
			"eo" : 22.1,
			"m" : 0.0,
			"g" : 0.19,
			"co" : 3.75,
			"c2" : 0.98,
			"To" : 298.15,
			"Z" : 47.0,
			"Therm_alpha_298" : 0,
			"Therm_diff_temp" : 0,
			"Therm_diff_alpha" : 0,
			"Therm_diff_Kprime" : 0
			}

Al = {
			"name" : "Al",
			"Uo" : 0.0,
			"n" : 1.0,
			"Vo" : 66.28873,#check
			"Ko" : 728,#check
			"kk" : 4.136,
			"QE1" : 381.0,
			"mE1" : 1.5,
			"QE2" : 202.0,
			"mE2" : 1.5,
			"delta" : -0.242,
			"t" : -0.958,
			"ao" : 0.0,
			"eo" : 64.1,
			"m" : 0.0,
			"g" : 0.33,
			"co" : 1.97,
			"c2" : 0.3,
			"To" : 298.15,
			"Z" : 13.0,
			"Therm_alpha_298" : 6.527e-05,
			"Therm_diff_temp" : -0.00914,
			"Therm_diff_alpha" : 4.395e-8,
			"Therm_diff_Kprime" : 0.000332
			}

Au = {
			"name" : "Au",
			"Uo" : 0.0,
			"n" : 1.0,
			"Vo" : 67.84963,#check
			"Ko" : 1970,#check
			"kk" : 5.662,
			"QE1" : 179.5,
			"mE1" : 1.5,
			"QE2" : 83.0,
			"mE2" : 1.5,
			"delta" : 0.134,
			"t" : 0.087,
			"ao" : 0.0,
			"eo" : 0.0,
			"m" : 0.0,
			"g" : 0.0,
			"co" : 4.1,
			"c2" : 0.25,
			"To" : 298.15,
			"Z" : 79.0,
			"Therm_alpha_298" : 3.911e-05,
			"Therm_diff_temp" : -0.01017,
			"Therm_diff_alpha" : 1.403e-8,
			"Therm_diff_Kprime" : 0.000330
			}

Cu = {
			"name" : "Cu",
			"Uo" : 0.0,
			"n" : 1.0,
			"Vo" : 47.23902,#check
			"Ko" : 1335,#check
			"kk" : 4.744,
			"QE1" : 296,
			"mE1" : 1.5,
			"QE2" : 169,
			"mE2" : 1.5,
			"delta" : -0.07,
			"t" : 1.401,
			"ao" : 0.0,
			"eo" : 27.7,
			"m" : 0.0,
			"g" : 2.18,
			"co" : 3.26,
			"c2" : 0.22,
			"To" : 298.15,
			"Z" : 29,
			"Therm_alpha_298" : 4.744e-05,
			"Therm_diff_temp" : -0.01078,
			"Therm_diff_alpha" : 1.735e-08,
			"Therm_diff_Kprime" : 0.000412
			}

Diamond = {
			"name" : "Diamond",
			"Uo" : 1290,
			"n" : 1.0,
			"Vo" : 45.35398,#check
			"Ko" : 4415,#check
			"kk" : 3.9,
			"QE1" : 684,
			"mE1" : 0.564,
			"QE2" : 1561,
			"mE2" : 2.436,
			"delta" : -0.506,
			"t" : 1.085,
			"ao" : 0.0,
			"eo" : 0.0,
			"m" : 0.0,
			"g" : 0.0,
			"co" : 0.66,
			"c2" : 0.68,
			"To" : 298.15,
			"Z" : 6.0,
			"Therm_alpha_298" : 0,
			"Therm_diff_temp" : 0,
			"Therm_diff_alpha" : 0,
			"Therm_diff_Kprime" : 0
			}

MgO = {
			"name" : "MgO",
			"Uo" : 0,
			"n" : 2.0,
			"Vo" : 74.71098,#check
			"Ko" : 1603,#check
			"kk" : 4.01,
			"QE1" : 748.0,
			"mE1" : 3.0,
			"QE2" : 401.0,
			"mE2" : 3,
			"delta" : -0.235,
			"t" : 0.301,
			"ao" : -17.4,
			"eo" : 0.0,
			"m" : 4.95,
			"g" : 0.0,
			"co" : 1.75,
			"c2" : 0.68,
			"To" : 298.15,
			"Z" : 10.34,
			"Therm_alpha_298" : 3.502e-05,
			"Therm_diff_temp" : -0.00624,
			"Therm_diff_alpha" : 0.775e-08,
			"Therm_diff_Kprime" : 0.000193
			}

Mo = {
			"name" : "Mo",
			"Uo" : 0,
			"n" : 1.0,
			"Vo" : 31.11518,#check
			"Ko" : 2600,#check
			"kk" : 4.2,
			"QE1" : 353.0,
			"mE1" : 1.5,
			"QE2" : 222.0,
			"mE2" : 1.5,
			"delta" : -0.802,
			"t" : -0.791,
			"ao" : 0,
			"eo" : 143.2,
			"m" : 0.0,
			"g" : 2.66,
			"co" : 2.75,
			"c2" : -0.95,
			"To" : 298.15,
			"Z" : 42,
			"Therm_alpha_298" : 1.359e-05,
			"Therm_diff_temp" : -0.00945,
			"Therm_diff_alpha" : 0.491e-08,
			"Therm_diff_Kprime" : 0.000139
			}

Nb = {
			"name" : "Nb",
			"Uo" : 0,
			"n" : 1.0,
			"Vo" : 35.96064,#check
			"Ko" : 1705,#check
			"kk" : 3.65,
			"QE1" : 134.0,
			"mE1" : 1.5,
			"QE2" : 302.0,
			"mE2" : 1.5,
			"delta" : -0.326,
			"t" : -0.763,
			"ao" : 0,
			"eo" : 115.9,
			"m" : 0.0,
			"g" : 0.9,
			"co" : 2.89,
			"c2" : -1.92,
			"To" : 298.15,
			"Z" : 41,
			"Therm_alpha_298" : 2.040e-05,
			"Therm_diff_temp" : -0.00372,
			"Therm_diff_alpha" : 0.358e-08,
			"Therm_diff_Kprime" : 0.000069
			}

Ta = {
			"name" : "Ta",
			"Uo" : 0,
			"n" : 1.0,
			"Vo" : 36.07023,#check
			"Ko" : 1910,#check
			"kk" : 3.83,
			"QE1" : 254.0,
			"mE1" : 1.5,
			"QE2" : 101.0,
			"mE2" : 1.5,
			"delta" : -0.101,
			"t" : -0.148,
			"ao" : 0,
			"eo" : 82.3,
			"m" : 0.0,
			"g" : 0.12,
			"co" : 3.74,
			"c2" : -2.49,
			"To" : 298.15,
			"Z" : 73,
			"Therm_alpha_298" : 1.942e-05,
			"Therm_diff_temp" : -0.00388,
			"Therm_diff_alpha" : 0.206e-08,
			"Therm_diff_Kprime" : 0.000077
			}

W = {
			"name" : "W",
			"Uo" : 0,
			"n" : 1.0,
			"Vo" : 31.72294,#check
			"Ko" : 3080,#check
			"kk" : 4.12,
			"QE1" : 172.0,
			"mE1" : 1.5,
			"QE2" : 309.0,
			"mE2" : 1.5,
			"delta" : -0.686,
			"t" : -0.591,
			"ao" : 0,
			"eo" : 100.1,
			"m" : 0.0,
			"g" : 2.77,
			"co" : 3.49,
			"c2" : -1.81,
			"To" : 298.15,
			"Z" : 74,
			"Therm_alpha_298" : 0,
			"Therm_diff_temp" : 0,
			"Therm_diff_alpha" : 0,
			"Therm_diff_Kprime" : 0
			}
dictionaries = [Ag, Al, Au, Cu, Diamond, MgO, Mo, Pt, Nb, Ta, W]
