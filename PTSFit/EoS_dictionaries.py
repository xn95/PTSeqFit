# -*- coding: utf-8 -*-
"""
Created on Fri Apr  4 15:31:43 2025

@author: sinjac

Required dictionary format:

dictionary = {
			"name" : string defining name to show in GUI,
			"U_0" : reference energy set to 0 if unsure (Jmol−1),
			"n" : atoms per formula unit,
			"V_0" : unit cell volume (Å^3),
			"K_0" : isothermal bulk modulus (GPa),
			"K_prime" : K' pressure derivative of K_0,
			"theta_1" : Einstein temperature 1 (K),
			"Ein_1" : Einstein number 1,
			"theta_2" : Einstein temperature 2 (K),
			"Ein_2" : Einstein number 1,
			"delta" : normalising constant,
			"t" : Grüneisen parameter,
			"a_0" : intrinsic anharmonicity parameter (10−6 K−1),
			"e_0" : free electrons parameter (10−6 K−1),
			"m" : anharmonic analogue of the Grüneisen parameter,
			"g" : electronic analogue of the Grüneisen parameter,
			"c_0" : Holzapfel EoS parameter 0,
			"c_2" : Holzapfel EoS parameter 2,
			"T_0" : reference temperature (K),
			"Z" : average atomic number  ,
			"Therm_alpha_298" : 𝛼(K−1),
			"Therm_diff_temp" : ∂𝐾0/∂𝑇 (GPa 𝐾−1),
			"Therm_diff_alpha" : ∂𝛼𝑇/∂𝑇 (𝐾−2),
			"Therm_diff_Kprime" : ∂𝐾′/∂𝑇 (𝐾−1)
			}



"""
dictionaries = []#append dictionaries to this list

Pt = {
			"name" : "Pt",
			"U_0" : 0.0,
			"n" : 1.0,
			"V_0" : 60.3836,
			"K_0" : 277.93567,
			"K_prime" : 5.35,
			"theta_1" : 177.0,
			"Ein_1" : 1.5,
			"theta_2" : 143.0,
			"Ein_2" : 1.5,
			"delta" : 0.167,
			"t" : -0.343,
			"a_0" : 0.0,
			"e_0" : 80.6,
			"m" : 0.0,
			"g" : 0.06,
			"c_0" : 3.78,
			"c_2" : -0.25,
			"T_0" : 298.15,
			"Z" : 78.0,
			"Therm_alpha_298" : 2.5718816e-05,
			"Therm_diff_temp" : -0.01312,
			"Therm_diff_alpha" : 5.779e-09,
			"Therm_diff_Kprime" : 0.000303
			}
dictionaries.append(Pt)

Ag = {
			"name" : "Ag",
			"U_0" : 0.0,
			"n" : 1.0,
			"V_0" : 68.08211,
			"K_0" : 1000,
			"K_prime" : 6.15,
			"theta_1" : 115.0,
			"Ein_1" : 1.5,
			"theta_2" : 199.0,
			"Ein_2" : 1.5,
			"delta" : 0.178,
			"t" : 2.21,
			"a_0" : 0.0,
			"e_0" : 22.1,
			"m" : 0.0,
			"g" : 0.19,
			"c_0" : 3.75,
			"c_2" : 0.98,
			"T_0" : 298.15,
			"Z" : 47.0,
			"Therm_alpha_298" : 0,
			"Therm_diff_temp" : 0,
			"Therm_diff_alpha" : 0,
			"Therm_diff_Kprime" : 0
			}
dictionaries.append(Ag)

Al = {
			"name" : "Al",
			"U_0" : 0.0,
			"n" : 1.0,
			"V_0" : 66.28873,
			"K_0" : 728,
			"K_prime" : 4.136,
			"theta_1" : 381.0,
			"Ein_1" : 1.5,
			"theta_2" : 202.0,
			"Ein_2" : 1.5,
			"delta" : -0.242,
			"t" : -0.958,
			"a_0" : 0.0,
			"e_0" : 64.1,
			"m" : 0.0,
			"g" : 0.33,
			"c_0" : 1.97,
			"c_2" : 0.3,
			"T_0" : 298.15,
			"Z" : 13.0,
			"Therm_alpha_298" : 6.527e-05,
			"Therm_diff_temp" : -0.00914,
			"Therm_diff_alpha" : 4.395e-8,
			"Therm_diff_Kprime" : 0.000332
			}
dictionaries.append(Al)

Au = {
			"name" : "Au",
			"U_0" : 0.0,
			"n" : 1.0,
			"V_0" : 67.84963,
			"K_0" : 1970,
			"K_prime" : 5.662,
			"theta_1" : 179.5,
			"Ein_1" : 1.5,
			"theta_2" : 83.0,
			"Ein_2" : 1.5,
			"delta" : 0.134,
			"t" : 0.087,
			"a_0" : 0.0,
			"e_0" : 0.0,
			"m" : 0.0,
			"g" : 0.0,
			"c_0" : 4.1,
			"c_2" : 0.25,
			"T_0" : 298.15,
			"Z" : 79.0,
			"Therm_alpha_298" : 3.911e-05,
			"Therm_diff_temp" : -0.01017,
			"Therm_diff_alpha" : 1.403e-8,
			"Therm_diff_Kprime" : 0.000330
			}
dictionaries.append(Au)

Cu = {
			"name" : "Cu",
			"U_0" : 0.0,
			"n" : 1.0,
			"V_0" : 47.23902,
			"K_0" : 1335,
			"K_prime" : 4.744,
			"theta_1" : 296,
			"Ein_1" : 1.5,
			"theta_2" : 169,
			"Ein_2" : 1.5,
			"delta" : -0.07,
			"t" : 1.401,
			"a_0" : 0.0,
			"e_0" : 27.7,
			"m" : 0.0,
			"g" : 2.18,
			"c_0" : 3.26,
			"c_2" : 0.22,
			"T_0" : 298.15,
			"Z" : 29,
			"Therm_alpha_298" : 4.744e-05,
			"Therm_diff_temp" : -0.01078,
			"Therm_diff_alpha" : 1.735e-08,
			"Therm_diff_Kprime" : 0.000412
			}
dictionaries.append(Cu)

Diamond = {
			"name" : "Diamond",
			"U_0" : 1290,
			"n" : 1.0,
			"V_0" : 45.35398,
			"K_0" : 4415,
			"K_prime" : 3.9,
			"theta_1" : 684,
			"Ein_1" : 0.564,
			"theta_2" : 1561,
			"Ein_2" : 2.436,
			"delta" : -0.506,
			"t" : 1.085,
			"a_0" : 0.0,
			"e_0" : 0.0,
			"m" : 0.0,
			"g" : 0.0,
			"c_0" : 0.66,
			"c_2" : 0.68,
			"T_0" : 298.15,
			"Z" : 6.0,
			"Therm_alpha_298" : 0,
			"Therm_diff_temp" : 0,
			"Therm_diff_alpha" : 0,
			"Therm_diff_Kprime" : 0
			}
dictionaries.append(Diamond)

MgO = {
			"name" : "MgO",
			"U_0" : 0,
			"n" : 2.0,
			"V_0" : 74.71098,
			"K_0" : 1603,
			"K_prime" : 4.01,
			"theta_1" : 748.0,
			"Ein_1" : 3.0,
			"theta_2" : 401.0,
			"Ein_2" : 3,
			"delta" : -0.235,
			"t" : 0.301,
			"a_0" : -17.4,
			"e_0" : 0.0,
			"m" : 4.95,
			"g" : 0.0,
			"c_0" : 1.75,
			"c_2" : 0.68,
			"T_0" : 298.15,
			"Z" : 10.0,
			"Therm_alpha_298" : 3.502e-05,
			"Therm_diff_temp" : -0.00624,
			"Therm_diff_alpha" : 0.775e-08,
			"Therm_diff_Kprime" : 0.000193
			}
dictionaries.append(MgO)

Mo = {
			"name" : "Mo",
			"U_0" : 0,
			"n" : 1.0,
			"V_0" : 31.11518,
			"K_0" : 2600,
			"K_prime" : 4.2,
			"theta_1" : 353.0,
			"Ein_1" : 1.5,
			"theta_2" : 222.0,
			"Ein_2" : 1.5,
			"delta" : -0.802,
			"t" : -0.791,
			"a_0" : 0,
			"e_0" : 143.2,
			"m" : 0.0,
			"g" : 2.66,
			"c_0" : 2.75,
			"c_2" : -0.95,
			"T_0" : 298.15,
			"Z" : 42,
			"Therm_alpha_298" : 1.359e-05,
			"Therm_diff_temp" : -0.00945,
			"Therm_diff_alpha" : 0.491e-08,
			"Therm_diff_Kprime" : 0.000139
			}
dictionaries.append(Mo)

Nb = {
			"name" : "Nb",
			"U_0" : 0,
			"n" : 1.0,
			"V_0" : 35.96064,
			"K_0" : 1705,
			"K_prime" : 3.65,
			"theta_1" : 134.0,
			"Ein_1" : 1.5,
			"theta_2" : 302.0,
			"Ein_2" : 1.5,
			"delta" : -0.326,
			"t" : -0.763,
			"a_0" : 0,
			"e_0" : 115.9,
			"m" : 0.0,
			"g" : 0.9,
			"c_0" : 2.89,
			"c_2" : -1.92,
			"T_0" : 298.15,
			"Z" : 41,
			"Therm_alpha_298" : 2.040e-05,
			"Therm_diff_temp" : -0.00372,
			"Therm_diff_alpha" : 0.358e-08,
			"Therm_diff_Kprime" : 0.000069
			}
dictionaries.append(Nb)

Ta = {
			"name" : "Ta",
			"U_0" : 0,
			"n" : 1.0,
			"V_0" : 36.07023,
			"K_0" : 1910,
			"K_prime" : 3.83,
			"theta_1" : 254.0,
			"Ein_1" : 1.5,
			"theta_2" : 101.0,
			"Ein_2" : 1.5,
			"delta" : -0.101,
			"t" : -0.148,
			"a_0" : 0,
			"e_0" : 82.3,
			"m" : 0.0,
			"g" : 0.12,
			"c_0" : 3.74,
			"c_2" : -2.49,
			"T_0" : 298.15,
			"Z" : 73,
			"Therm_alpha_298" : 1.942e-05,
			"Therm_diff_temp" : -0.00388,
			"Therm_diff_alpha" : 0.206e-08,
			"Therm_diff_Kprime" : 0.000077
			}
dictionaries.append(Ta)

W = {
			"name" : "W",
			"U_0" : 0,
			"n" : 1.0,
			"V_0" : 31.72294,
			"K_0" : 3080,
			"K_prime" : 4.12,
			"theta_1" : 172.0,
			"Ein_1" : 1.5,
			"theta_2" : 309.0,
			"Ein_2" : 1.5,
			"delta" : -0.686,
			"t" : -0.591,
			"a_0" : 0,
			"e_0" : 100.1,
			"m" : 0.0,
			"g" : 2.77,
			"c_0" : 3.49,
			"c_2" : -1.81,
			"T_0" : 298.15,
			"Z" : 74,
			"Therm_alpha_298" : 0,
			"Therm_diff_temp" : 0,
			"Therm_diff_alpha" : 0,
			"Therm_diff_Kprime" : 0
			}
dictionaries.append(W)
