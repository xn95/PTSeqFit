# -*- coding: utf-8 -*-
"""
Created on Thu Mar  6 11:52:30 2025

@author: sinjac
P, V calculator based on visual basic code published by Sokolova et. al: https://doi.org/10.1016/j.cageo.2016.06.002
T calculator based on McHardy et. al: https://doi.org/10.1080/08957959.2023.2187294

General use:

This file contains a class handling calculations of EoS via a Birch–Murnaghan EoS
The class requres a calibrant dictionary as 1st arg and a list of PVT as the 2nd arg

The format of PVT is a list of [P,V,T] in GPa, Å^3 and K
Setting one of these values to None will calculate the unknown, eg:
BM(crystal, [10, 60, None]) will calculate the temperature of a crystal system at 10 GPa with unit cell volume of 60 Å^3 based on EoS parameters defined by the crystal dictionary
See EoS_dictionaries for the format of a dictionary

The BM object will calculate various thermodynamic parameters which can be accessed by the appropriate class attribute
"""

import numpy as np
from scipy.optimize import fsolve
from EoS_dictionaries import dictionaries

class BM: #container for EoS parameters and functions
	def __init__(self, calibrant, PVT):

		self.calibrant = calibrant
		self.PVT = PVT
		self.loaded_params = self.calibrant
        
		def f_gamV(x, n, Z, V_0, K_0, K_prime, gamVo, gb, beta):
			fw = -np.log(3 * K_0 / 10 / (1003.6 * (Z * n / (V_0 * 10)) ** (5 / 3)))
			ff = x ** (1 / 3)
			aa = 1.5 * (K_prime - 3) - fw
			KT = K_0 * 1000 / ff ** 6 * np.exp(fw * (1 - ff)) * ((-5 / ff ** 2 + 4 / ff) * (1 + aa * ff - aa * ff ** 2) + (1 / ff - 1) * (1 + aa * ff - aa * ff ** 2) * (-fw) + (1 / ff - 1) * (aa - 2 * aa * ff)) * (-x)
			Px = 3 * K_0 * 1000 * np.exp(fw * (1 - ff)) * (1 / ff ** 5 - 1 / ff ** 4) * (1 + aa * ff - aa * ff ** 2)
			ex = 3 / ff ** 4 * K_0 * 1000 * np.exp(fw * (1 - ff)) * ((-5 / ff ** 2 + 4 / ff) * (1 + aa * ff - aa * ff ** 2) - (1 / ff - 1) * (1 + aa * ff - aa * ff ** 2) * fw + (1 / ff - 1) * (aa - 2 * aa * ff))
			ex1 = 1 / ff ** 3 * K_0 * 1000 * np.exp(fw * (1 - ff)) * fw * ((-5 / ff ** 2 + 4 / ff) * (1 + aa * ff - aa * ff ** 2) - (1 / ff - 1) * (1 + aa * ff - aa * ff ** 2) * fw + (1 / ff - 1) * (aa - 2 * aa * ff))
			ex2 = K_0 * 1000 * np.exp(fw * (1 - ff)) * ((10 / ff ** 3 - 4 / ff ** 2) * (1 + aa * ff - aa * ff ** 2) + (-5 / ff ** 2 + 4 / ff) * (aa - 2 * aa * ff) + fw / ff ** 2 * (1 + aa * ff - aa * ff ** 2) - (1 / ff - 1) * (aa - 2 * aa * ff) * fw - (aa - 2 * aa * ff) / ff ** 2 - 2 * aa * (1 / ff - 1)) / ff ** 3
			kkx = (ex + ex1 - ex2) / (-KT / ff) / 3
			gt = gb - beta * x ** (1 / 3)
			f_gamV = (gamVo + (-3 * KT + 2 * Px * gt + 9 * KT * kkx - 6 * gt * KT) / 6 / (3 * KT - 2 * Px * gt)) / x

			return f_gamV

		def I_gamV(Z, n, b, V_0, K_0, K_prime, gamVo, gb, beta):
			e = 0.0001
			a = 1

			if b < 1:
				a = b
				b = 1

			s2 = 1
			h = b - a
			S = f_gamV(a, n, Z, V_0, K_0, K_prime, gamVo, gb, beta) + f_gamV(b, n, Z, V_0, K_0, K_prime, gamVo, gb, beta)
			x = a + h
			s3 = s2
			h = h / 2
			s1 = 0
			x = a + h

			while x > e:
				s3 = s2
				h = h / 2
				s1 = 0
				x = a + h

				while  x < b:
					s1 += 2 * f_gamV(x, n, Z, V_0, K_0, K_prime, gamVo, gb, beta)
					x += 2 * h

				S = S + s1
				s2 = (S + s1) * h / 3
				x = abs(s3 - s2) / 15

			return -s2 if b > 1 else s2

		def xAP2(n, Z, Tr, V_0, K_0, K_prime, QBo, d, mb, QBo1, d1, mb1, QEo1, Ein_1, QEo2, Ein_2, gamVo, gb, beta, a_0, m, mm, ae, TK, Pbar):
			R = 8.31451
			e = 0.00000000001
			L = 0.4
			Rx = 1.3
			FL = 5000000
			iter_count = 0

			while abs(Rx - L) > e:
				x = (L + Rx) / 2
				iter_count += 1
				fw = -np.log(3 * K_0 / 10 / (1003.6 * (Z * n / (V_0 * 10)) ** (5 / 3)))
				ff = x ** (1 / 3)
				aa = 1.5 * (K_prime - 3) - fw
				Px = 3 * K_0 * 1000 * np.exp(fw * (1 - ff)) * (1 / ff ** 5 - 1 / ff ** 4) * (1 + aa * ff - aa * ff ** 2)
				KT = K_0 * 1000 / ff ** 6 * np.exp(fw * (1 - ff)) * ((-5 / ff ** 2 + 4 / ff) * (1 + aa * ff - aa * ff ** 2) + (1 / ff - 1) * (1 + aa * ff - aa * ff ** 2) * (-fw) + (1 / ff - 1) * (aa - 2 * aa * ff)) * (-x)
				ex = 3 / ff ** 4 * K_0 * 1000 * np.exp(fw * (1 - ff)) * ((-5 / ff ** 2 + 4 / ff) * (1 + aa * ff - aa * ff ** 2) - (1 / ff - 1) * (1 + aa * ff - aa * ff ** 2) * fw + (1 / ff - 1) * (aa - 2 * aa * ff))
				ex1 = 1 / ff ** 3 * K_0 * 1000 * np.exp(fw * (1 - ff)) * fw * ((-5 / ff ** 2 + 4 / ff) * (1 + aa * ff - aa * ff ** 2) - (1 / ff - 1) * (1 + aa * ff - aa * ff ** 2) * fw + (1 / ff - 1) * (aa - 2 * aa * ff))
				ex2 = K_0 * 1000 * np.exp(fw * (1 - ff)) * ((10 / ff ** 3 - 4 / ff ** 2) * (1 + aa * ff - aa * ff ** 2) + (-5 / ff ** 2 + 4 / ff) * (aa - 2 * aa * ff) + fw / ff ** 2 * (1 + aa * ff - aa * ff ** 2) - (1 / ff - 1) * (aa - 2 * aa * ff) * fw - (aa - 2 * aa * ff) / ff ** 2 - 2 * aa * (1 / ff - 1)) / ff ** 3
				kkx = (ex + ex1 - ex2) / (-KT / ff) / 3
				gt = gb - beta * x ** (1 / 3)
				gamV = (-3 * KT + 2 * Px * gt + 9 * KT * kkx - 6 * gt * KT) / 6 / (3 * KT - 2 * Px * gt) + gamVo
				expp = np.exp(I_gamV(Z, n, x, V_0, K_0, K_prime, gamVo, gb, beta) * (1))
				ff2 = (x + 0.00001) ** (1 / 3)
				ex = 3 / ff2 ** 4 * K_0 * 1000 * np.exp(fw * (1 - ff2)) * ((-5 / ff2 ** 2 + 4 / ff2) * (1 + aa * ff2 - aa * ff2 ** 2) - (1 / ff2 - 1) * (1 + aa * ff2 - aa * ff2 ** 2) * fw + (1 / ff2 - 1) * (aa - 2 * aa * ff2))
				ex1 = 1 / ff2 ** 3 * K_0 * 1000 * np.exp(fw * (1 - ff2)) * fw * ((-5 / ff2 ** 2 + 4 / ff2) * (1 + aa * ff2 - aa * ff2 ** 2) - (1 / ff2 - 1) * (1 + aa * ff2 - aa * ff2 ** 2) * fw + (1 / ff2 - 1) * (aa - 2 * aa * ff2))
				ex2 = K_0 * 1000 * np.exp(fw * (1 - ff2)) * ((10 / ff2 ** 3 - 4 / ff2 ** 2) * (1 + aa * ff2 - aa * ff2 ** 2) + (-5 / ff2 ** 2 + 4 / ff2) * (aa - 2 * aa * ff2) + fw / ff2 ** 2 * (1 + aa * ff2 - aa * ff2 ** 2) - (1 / ff2 - 1) * (aa - 2 * aa * ff2) * fw - (aa - 2 * aa * ff2) / ff2 ** 2 - 2 * aa * (1 / ff2 - 1)) / ff2 ** 3
				ff1 = (x - 0.00001) ** (1 / 3)
				ex = 3 / ff1 ** 4 * K_0 * 1000 * np.exp(fw * (1 - ff1)) * ((-5 / ff1 ** 2 + 4 / ff1) * (1 + aa * ff1 - aa * ff1 ** 2) - (1 / ff1 - 1) * (1 + aa * ff1 - aa * ff1 ** 2) * fw + (1 / ff1 - 1) * (aa - 2 * aa * ff1))
				ex1 = 1 / ff1 ** 3 * K_0 * 1000 * np.exp(fw * (1 - ff1)) * fw * ((-5 / ff1 ** 2 + 4 / ff1) * (1 + aa * ff1 - aa * ff1 ** 2) - (1 / ff1 - 1) * (1 + aa * ff1 - aa * ff1 ** 2) * fw + (1 / ff1 - 1) * (aa - 2 * aa * ff1))
				ex2 = K_0 * 1000 * np.exp(fw * (1 - ff1)) * ((10 / ff1 ** 3 - 4 / ff1 ** 2) * (1 + aa * ff1 - aa * ff1 ** 2) + (-5 / ff1 ** 2 + 4 / ff1) * (aa - 2 * aa * ff1) + fw / ff1 ** 2 * (1 + aa * ff1 - aa * ff1 ** 2) - (1 / ff1 - 1) * (aa - 2 * aa * ff1) * fw - (aa - 2 * aa * ff1) / ff1 ** 2 - 2 * aa * (1 / ff1 - 1)) / ff1 ** 3
				QB = QBo * expp
				QB1 = QBo1 * expp
				theta_1 = QEo1 * expp
				theta_2 = QEo2 * expp
				VV = x * V_0
				ggr = d * np.log(1 + QB / Tr / d)
				br = 1 / (np.exp(ggr) - 1)
				PthBr = mb * R * ((Tr * QB * d * br / (Tr * d + QB)) * gamV / VV)
				gg = d * np.log(1 + QB / TK / d)
				b = 1 / (np.exp(gg) - 1)
				ex1 = np.exp(QB / TK)
				PthB = mb * R * ((TK * QB * d * b / (TK * d + QB)) * gamV / VV)
				ex2 = np.exp(QB1 / TK)
				e1 = np.exp(theta_1 / TK)
				PthE1 = Ein_1 * R * ((theta_1 / (e1 - 1)) * gamV / VV)
				e2 = np.exp(theta_2 / TK)
				PthE2 = Ein_2 * R * ((theta_2 / (e2 - 1)) * gamV / VV)
				Pa = 1.5 * n * R * a_0 / 1000000 * x ** (m) * (m) / x / V_0 * (TK ** 2 - Tr ** 2)
				Pe = 1.5 * n * R * ae / 1000000 * x ** (mm) * (mm) / x / V_0 * (TK ** 2 - Tr ** 2)
				Pth = PthB - PthBr + PthE1 - PthE2 + Px + Pe + Pa - Pbar
				F = Pth

				if FL * F > 0:
					L = x

				else:
					Rx = x

				if iter_count > 100000:
					print("Max iterations reached")
					break
					return iter_count

			return x

		def P(Z, n, x, V_0, K_0, K_prime):
			fw = -np.log(3 * K_0 / 10 / (1003.6 * (Z * n / (V_0 * 10)) ** (5 / 3)))
			ff = x ** (1 / 3)
			aa = 1.5 * (K_prime - 3) - fw
			P = 3 * K_0 * 1000 * np.exp(fw * (1 - ff)) * (1 / ff ** 5 - 1 / ff ** 4) * (1 + aa * ff - aa * ff ** 2)
			return P

		def I(Z, nn, b, V_0, K_0, K_prime):
			e = 0.001
			a = 1
			if b < 1:
				a = b
				b = 1
			s2 = 1
			h = b - a
			S = (P(Z, nn, a, V_0, K_0, K_prime) + P(Z, nn, b, V_0, K_0, K_prime)) * V_0
			s3 = s2
			h = h / 2
			s1 = 0
			x = a + h

			while x > e:
				s3 = s2
				h = h / 2
				s1 = 0
				x = a + h

				while x < b:
					s1 += 2 * P(Z, nn, x, V_0, K_0, K_prime) * V_0
					x += 2 * h

				S = S + s1
				s2 = (S + s1) * h / 3
				x = abs(s3 - s2) / 15

			return -s2 if b > 1 else s2

		def F(n, Z, V_0, K_0, K_prime, Tr, x, expp, QBo, d, mb, QB1o, d1, mb1, QE1o, Ein_1, QE2o, Ein_2, TK, gamVo, beta, gb, a_0, m, mm, ae):
			R = 8.31451
			fw = -np.log(3 * K_0 / 10 / (1003.6 * (Z * n / (V_0 * 10)) ** (5 / 3)))
			ff = x ** (1 / 3)
			aa = 1.5 * (K_prime - 3) - fw
			ex = 3 / ff ** 4 * K_0 * 1000 * np.exp(fw * (1 - ff)) * ((-5 / ff ** 2 + 4 / ff) * (1 + aa * ff - aa * ff ** 2) - (1 / ff - 1) * (1 + aa * ff - aa * ff ** 2) * fw + (1 / ff - 1) * (aa - 2 * aa * ff))
			ff2 = (x + 0.00001) ** (1 / 3)
			ex = 3 / ff2 ** 4 * K_0 * 1000 * np.exp(fw * (1 - ff2)) * ((-5 / ff2 ** 2 + 4 / ff2) * (1 + aa * ff2 - aa * ff2 ** 2) - (1 / ff2 - 1) * (1 + aa * ff2 - aa * ff2 ** 2) * fw + (1 / ff2 - 1) * (aa - 2 * aa * ff2))
			ff1 = (x - 0.00001) ** (1 / 3)
			ex = 3 / ff1 ** 4 * K_0 * 1000 * np.exp(fw * (1 - ff1)) * ((-5 / ff1 ** 2 + 4 / ff1) * (1 + aa * ff1 - aa * ff1 ** 2) - (1 / ff1 - 1) * (1 + aa * ff1 - aa * ff1 ** 2) * fw + (1 / ff1 - 1) * (aa - 2 * aa * ff1))
			QB = QBo * expp
			QB1 = QB1o * expp
			theta_1 = QE1o * expp
			theta_2 = QE2o * expp
			gg = d * np.log(1 + QB / TK / d)
			ex = np.exp(QB / TK)
			FB = mb * R * ((d - 1) * QB / 2 / d - TK * np.log(1 + 1 / (np.exp(gg) - 1)))
			gg1 = d1 * np.log(1 + QB1 / TK / d1)
			ex = np.exp(QB1 / TK)
			FB1 = mb1 * R * ((d1 - 1) * QB1 / 2 / d1 - TK * np.log(1 + 1 / (np.exp(gg1) - 1)))
			ex = np.exp(theta_1 / TK)
			FE1 = Ein_1 * R * (theta_1 / 2 + TK * np.log(1 - 1 / ex))
			ex = np.exp(theta_2 / TK)
			FE2 = Ein_2 * R * (theta_2 / 2 + TK * np.log(1 - 1 / np.exp(theta_2 / TK)))
			ggr = d * np.log(1 + QB / Tr / d)
			ex = np.exp(QB / Tr)
			FBr = mb * R * ((d - 1) * QB / 2 / d - Tr * np.log(1 + 1 / (np.exp(ggr) - 1)))
			ggr1 = d1 * np.log(1 + QB1 / Tr / d1)
			ex = np.exp(QB1 / Tr)
			FBr1 = mb1 * R * ((d1 - 1) * QB1 / 2 / d1 - Tr * np.log(1 + 1 / (np.exp(ggr1) - 1)))
			ex = np.exp(theta_1 / Tr)
			FE1r = Ein_1 * R * (theta_1 / 2 + Tr * np.log(1 - 1 / ex))
			ex = np.exp(theta_2 / Tr)
			FE2r = Ein_2 * R * (theta_2 / 2 + Tr * np.log(1 - 1 / ex))
			Fa = -1.5 * n * R * a_0 / 1000000 * x ** (m) * (TK ** 2 - Tr ** 2)
			Fel = -1.5 * n * R * ae / 1000000 * x ** (mm) * (TK ** 2 - Tr ** 2)
			F = FB + FB1 + FE1 + FE2 - FBr - FBr1 - FE1r - FE2r + Fel + Fa

			return F

		def S(n, Z, V_0, K_0, K_prime, Tr, x, expp, QBo, d, mb, QB1o, d1, mb1, QE1o, Ein_1, QE2o, Ein_2, TK, gamVo, beta, gb, a_0, m, mm, ae):
			R = 8.31451
			QB = QBo * expp
			QB1 = QB1o * expp
			theta_1 = QE1o * expp
			theta_2 = QE2o * expp
			gg = d * np.log(1 + QB / TK / d)
			expS = 1 / (np.exp(gg) - 1)
			SB = mb * R * (np.log(1 + expS) + QB * d * expS / (TK * d + QB))
			gg1 = d1 * np.log(1 + QB1 / TK / d1)
			expS1 = 1 / (np.exp(gg1) - 1)
			SB1 = mb1 * R * (np.log(1 + expS1) + QB1 * d1 * expS1 / (TK * d1 + QB1))
			exp1 = np.exp(theta_1 / TK)
			SE1 = Ein_1 * R * (-np.log(1 - 1 / exp1) + theta_1 / TK / (exp1 - 1))
			exp2 = np.exp(theta_2 / TK)
			SE2 = Ein_2 * R * (-np.log(1 - 1 / exp2) + theta_2 / TK / (exp2 - 1))
			Sa = 3 * n * R * a_0 / 1000000 * x ** (m) * (TK)
			See = 3 * n * R * ae / 1000000 * x ** (mm) * (TK)
			S = SB + SB1 + SE1 + SE2 + See + Sa

			return S

		def CCv(n, Z, V_0, K_0, K_prime, Tr, x, expp, QBo, d, mb, QB1o, d1, mb1, QE1o, Ein_1, QE2o, Ein_2, TK, gamVo, beta, gb, a_0, m, mm, ae):
			R = 8.31451
			fw = -np.log(3 * K_0 / 10 / (1003.6 * (Z * n / (V_0 * 10)) ** (5 / 3)))
			ff = x ** (1 / 3)
			aa = 1.5 * (K_prime - 3) - fw
			ex1 = 1 / ff ** 3 * K_0 * 1000 * np.exp(fw * (1 - ff)) * fw * ((-5 / ff ** 2 + 4 / ff) * (1 + aa * ff - aa * ff ** 2) - (1 / ff - 1) * (1 + aa * ff - aa * ff ** 2) * fw + (1 / ff - 1) * (aa - 2 * aa * ff))
			ex2 = K_0 * 1000 * np.exp(fw * (1 - ff)) * ((10 / ff ** 3 - 4 / ff ** 2) * (1 + aa * ff - aa * ff ** 2) + (-5 / ff ** 2 + 4 / ff) * (aa - 2 * aa * ff) + fw / ff ** 2 * (1 + aa * ff - aa * ff ** 2) - (1 / ff - 1) * (aa - 2 * aa * ff) * fw - (aa - 2 * aa * ff) / ff ** 2 - 2 * aa * (1 / ff - 1)) / ff ** 3
			ff2 = (x + 0.00001) ** (1 / 3)
			ex1 = 1 / ff2 ** 3 * K_0 * 1000 * np.exp(fw * (1 - ff2)) * fw * ((-5 / ff2 ** 2 + 4 / ff2) * (1 + aa * ff2 - aa * ff2 ** 2) - (1 / ff2 - 1) * (1 + aa * ff2 - aa * ff2 ** 2) * fw + (1 / ff2 - 1) * (aa - 2 * aa * ff2))
			ex2 = K_0 * 1000 * np.exp(fw * (1 - ff2)) * ((10 / ff2 ** 3 - 4 / ff2 ** 2) * (1 + aa * ff2 - aa * ff2 ** 2) + (-5 / ff2 ** 2 + 4 / ff2) * (aa - 2 * aa * ff2) + fw / ff2 ** 2 * (1 + aa * ff2 - aa * ff2 ** 2) - (1 / ff2 - 1) * (aa - 2 * aa * ff2) * fw - (aa - 2 * aa * ff2) / ff2 ** 2 - 2 * aa * (1 / ff2 - 1)) / ff2 ** 3
			ff1 = (x - 0.00001) ** (1 / 3)
			ex1 = 1 / ff1 ** 3 * K_0 * 1000 * np.exp(fw * (1 - ff1)) * fw * ((-5 / ff1 ** 2 + 4 / ff1) * (1 + aa * ff1 - aa * ff1 ** 2) - (1 / ff1 - 1) * (1 + aa * ff1 - aa * ff1 ** 2) * fw + (1 / ff1 - 1) * (aa - 2 * aa * ff1))
			ex2 = K_0 * 1000 * np.exp(fw * (1 - ff1)) * ((10 / ff1 ** 3 - 4 / ff1 ** 2) * (1 + aa * ff1 - aa * ff1 ** 2) + (-5 / ff1 ** 2 + 4 / ff1) * (aa - 2 * aa * ff1) + fw / ff1 ** 2 * (1 + aa * ff1 - aa * ff1 ** 2) - (1 / ff1 - 1) * (aa - 2 * aa * ff1) * fw - (aa - 2 * aa * ff1) / ff1 ** 2 - 2 * aa * (1 / ff1 - 1)) / ff1 ** 3
			QB = QBo * expp
			QB1 = QB1o * expp
			theta_1 = QE1o * expp
			theta_2 = QE2o * expp
			gg = d * np.log(1 + QB / TK / d)
			b = 1 / (np.exp(gg) - 1)
			CvB = mb * R * ((QB * d / (TK * d + QB)) ** 2 * b * (1 / d + 1 + b))
			gg1 = d1 * np.log(1 + QB1 / TK / d1)
			b1 = 1 / (np.exp(gg1) - 1)
			CvB1 = mb1 * R * ((QB1 * d1 / (TK * d1 + QB1)) ** 2 * b1 * (1 / d1 + 1 + b1))
			ex1 = np.exp(theta_1 / TK)
			CvE1 = Ein_1 * R * (theta_1 ** 2 / TK ** 2 * ex1 / (ex1 - 1) ** 2)
			ex2 = np.exp(theta_2 / TK)
			CvE2 = Ein_2 * R * (theta_2 ** 2 / TK ** 2 * ex2 / (ex2 - 1) ** 2)
			Cva = 3 * n * R * a_0 / 1000000 * x ** (m) * (TK)
			CvE = 3 * n * R * ae / 1000000 * x ** (mm) * (TK)
			CCv = CvB + CvB1 + CvE1 + CvE2 + CvE + Cva

			return CCv

		def Pth(n, Z, V_0, K_0, K_prime, Tr, x, expp, V, QBo, d, mb, QB1o, d1, mb1, QE1o, Ein_1, QE2o, Ein_2, TK, gamVo, beta, gb, a_0, m, mm, ae):
			R = 8.31451
			fw = -np.log(3 * K_0 / 10 / (1003.6 * (Z * n / (V_0 * 10)) ** (5 / 3)))
			ff = x ** (1 / 3)
			aa = 1.5 * (K_prime - 3) - fw
			Px = 3 * K_0 * 1000 * np.exp(fw * (1 - ff)) * (1 / ff ** 5 - 1 / ff ** 4) * (1 + aa * ff - aa * ff ** 2)
			KT = K_0 * 1000 / ff ** 6 * np.exp(fw * (1 - ff)) * ((-5 / ff ** 2 + 4 / ff) * (1 + aa * ff - aa * ff ** 2) + (1 / ff - 1) * (1 + aa * ff - aa * ff ** 2) * (-fw) + (1 / ff - 1) * (aa - 2 * aa * ff)) * (-x)
			ex = 3 / ff ** 4 * K_0 * 1000 * np.exp(fw * (1 - ff)) * ((-5 / ff ** 2 + 4 / ff) * (1 + aa * ff - aa * ff ** 2) - (1 / ff - 1) * (1 + aa * ff - aa * ff ** 2) * fw + (1 / ff - 1) * (aa - 2 * aa * ff))
			ex1 = 1 / ff ** 3 * K_0 * 1000 * np.exp(fw * (1 - ff)) * fw * ((-5 / ff ** 2 + 4 / ff) * (1 + aa * ff - aa * ff ** 2) - (1 / ff - 1) * (1 + aa * ff - aa * ff ** 2) * fw + (1 / ff - 1) * (aa - 2 * aa * ff))
			ex2 = K_0 * 1000 * np.exp(fw * (1 - ff)) * ((10 / ff ** 3 - 4 / ff ** 2) * (1 + aa * ff - aa * ff ** 2) + (-5 / ff ** 2 + 4 / ff) * (aa - 2 * aa * ff) + fw / ff ** 2 * (1 + aa * ff - aa * ff ** 2) - (1 / ff - 1) * (aa - 2 * aa * ff) * fw - (aa - 2 * aa * ff) / ff ** 2 - 2 * aa * (1 / ff - 1)) / ff ** 3
			kkx = (ex + ex1 - ex2) / (-KT / ff) / 3
			gt = gb - beta * x ** (1 / 3)
			gamV = (-3 * KT + 2 * Px * gt + 9 * KT * kkx - 6 * gt * KT) / 6 / (3 * KT - 2 * Px * gt) + gamVo
			ff2 = (x + 0.00001) ** (1 / 3)
			ex = 3 / ff2 ** 4 * K_0 * 1000 * np.exp(fw * (1 - ff2)) * ((-5 / ff2 ** 2 + 4 / ff2) * (1 + aa * ff2 - aa * ff2 ** 2) - (1 / ff2 - 1) * (1 + aa * ff2 - aa * ff2 ** 2) * fw + (1 / ff2 - 1) * (aa - 2 * aa * ff2))
			ex1 = 1 / ff2 ** 3 * K_0 * 1000 * np.exp(fw * (1 - ff2)) * fw * ((-5 / ff2 ** 2 + 4 / ff2) * (1 + aa * ff2 - aa * ff2 ** 2) - (1 / ff2 - 1) * (1 + aa * ff2 - aa * ff2 ** 2) * fw + (1 / ff2 - 1) * (aa - 2 * aa * ff2))
			ex2 = K_0 * 1000 * np.exp(fw * (1 - ff2)) * ((10 / ff2 ** 3 - 4 / ff2 ** 2) * (1 + aa * ff2 - aa * ff2 ** 2) + (-5 / ff2 ** 2 + 4 / ff2) * (aa - 2 * aa * ff2) + fw / ff2 ** 2 * (1 + aa * ff2 - aa * ff2 ** 2) - (1 / ff2 - 1) * (aa - 2 * aa * ff2) * fw - (aa - 2 * aa * ff2) / ff2 ** 2 - 2 * aa * (1 / ff2 - 1)) / ff2 ** 3
			ff1 = (x - 0.00001) ** (1 / 3)
			ex = 3 / ff1 ** 4 * K_0 * 1000 * np.exp(fw * (1 - ff1)) * ((-5 / ff1 ** 2 + 4 / ff1) * (1 + aa * ff1 - aa * ff1 ** 2) - (1 / ff1 - 1) * (1 + aa * ff1 - aa * ff1 ** 2) * fw + (1 / ff1 - 1) * (aa - 2 * aa * ff1))
			ex1 = 1 / ff1 ** 3 * K_0 * 1000 * np.exp(fw * (1 - ff1)) * fw * ((-5 / ff1 ** 2 + 4 / ff1) * (1 + aa * ff1 - aa * ff1 ** 2) - (1 / ff1 - 1) * (1 + aa * ff1 - aa * ff1 ** 2) * fw + (1 / ff1 - 1) * (aa - 2 * aa * ff1))
			ex2 = K_0 * 1000 * np.exp(fw * (1 - ff1)) * ((10 / ff1 ** 3 - 4 / ff1 ** 2) * (1 + aa * ff1 - aa * ff1 ** 2) + (-5 / ff1 ** 2 + 4 / ff1) * (aa - 2 * aa * ff1) + fw / ff1 ** 2 * (1 + aa * ff1 - aa * ff1 ** 2) - (1 / ff1 - 1) * (aa - 2 * aa * ff1) * fw - (aa - 2 * aa * ff1) / ff1 ** 2 - 2 * aa * (1 / ff1 - 1)) / ff1 ** 3
			QB = QBo * expp
			QB1 = QB1o * expp
			theta_1 = QE1o * expp
			theta_2 = QE2o * expp
			ggr = d * np.log(1 + QB / Tr / d)
			br = 1 / (np.exp(ggr) - 1)
			PBr = mb * R * ((QB * (d - 1) / 2 / d + Tr * QB * d * br / (Tr * d + QB)) * gamV / V)
			gg = d * np.log(1 + QB / TK / d)
			b = 1 / (np.exp(gg) - 1)
			ex1 = np.exp(QB / TK)
			PB = mb * R * ((QB * (d - 1) / 2 / d + TK * QB * d * b / (TK * d + QB)) * gamV / V)
			gg1r = d1 * np.log(1 + QB1 / Tr / d1)
			b1r = 1 / (np.exp(gg1r) - 1)
			PB1r = mb1 * R * ((QB1 * (d1 - 1) / 2 / d1 + Tr * QB1 * d1 * b1r / (Tr * d1 + QB1)) * gamV / V)
			gg1 = d1 * np.log(1 + QB1 / TK / d1)
			b1 = 1 / (np.exp(gg1) - 1)
			ex2 = np.exp(QB1 / TK)
			PB1 = mb1 * R * ((QB1 * (d1 - 1) / 2 / d1 + TK * QB1 * d1 * b1 / (TK * d1 + QB1)) * gamV / V)
			e1 = np.exp(theta_1 / TK)
			PE1 = Ein_1 * R * ((theta_1 / 2 + theta_1 / (e1 - 1)) * gamV / V)
			e2 = np.exp(theta_2 / TK)
			PE2 = Ein_2 * R * ((theta_2 / 2 + theta_2 / (e2 - 1)) * gamV / V)
			e1r = np.exp(theta_1 / Tr)
			PE1r = Ein_1 * R * ((theta_1 / 2 + theta_1 / (e1r - 1)) * gamV / V)
			e2r = np.exp(theta_2 / Tr)
			PE2r = Ein_2 * R * ((theta_2 / 2 + theta_2 / (e2r - 1)) * gamV / V)
			Pea = 3 / 2 * n * R * a_0 / 1000000 * x ** (m) * (m) / V * (TK ** 2 - Tr ** 2)
			Pee = 3 / 2 * n * R * ae / 1000000 * x ** (mm) * (mm) / V * (TK ** 2 - Tr ** 2)
			Pth = PB + PB1 + PE1 + PE2 - PBr - PB1r - PE1r - PE2r + Pee + Pea

			return Pth

		def KTth(n, Z, V_0, K_0, K_prime, Tr, x, expp, V, QBo, d, mb, QB1o, d1, mb1, QE1o, Ein_1, QE2o, Ein_2, TK, gamVo, beta, gb, a_0, m, mm, ae):
			R = 8.31451

			fw = -np.log(3 * K_0 / 10 / (1003.6 * (Z * n / (V_0 * 10)) ** (5 / 3)))
			ff = x ** (1 / 3)
			aa = 1.5 * (K_prime - 3) - fw
			Px = 3 * K_0 * 1000 * np.exp(fw * (1 - ff)) * (1 / ff ** 5 - 1 / ff ** 4) * (1 + aa * ff - aa * ff ** 2)
			KT = K_0 * 1000 / ff ** 6 * np.exp(fw * (1 - ff)) * ((-5 / ff ** 2 + 4 / ff) * (1 + aa * ff - aa * ff ** 2) + (1 / ff - 1) * (1 + aa * ff - aa * ff ** 2) * (-fw) + (1 / ff - 1) * (aa - 2 * aa * ff)) * (-x)
			ex = 3 / ff ** 4 * K_0 * 1000 * np.exp(fw * (1 - ff)) * ((-5 / ff ** 2 + 4 / ff) * (1 + aa * ff - aa * ff ** 2) - (1 / ff - 1) * (1 + aa * ff - aa * ff ** 2) * fw + (1 / ff - 1) * (aa - 2 * aa * ff))
			ex1 = 1 / ff ** 3 * K_0 * 1000 * np.exp(fw * (1 - ff)) * fw * ((-5 / ff ** 2 + 4 / ff) * (1 + aa * ff - aa * ff ** 2) - (1 / ff - 1) * (1 + aa * ff - aa * ff ** 2) * fw + (1 / ff - 1) * (aa - 2 * aa * ff))
			ex2 = K_0 * 1000 * np.exp(fw * (1 - ff)) * ((10 / ff ** 3 - 4 / ff ** 2) * (1 + aa * ff - aa * ff ** 2) + (-5 / ff ** 2 + 4 / ff) * (aa - 2 * aa * ff) + fw / ff ** 2 * (1 + aa * ff - aa * ff ** 2) - (1 / ff - 1) * (aa - 2 * aa * ff) * fw - (aa - 2 * aa * ff) / ff ** 2 - 2 * aa * (1 / ff - 1)) / ff ** 3
			kkx = (ex + ex1 - ex2) / (-KT / ff) / 3
			gt = gb - beta * x ** (1 / 3)
			gtx = -beta / 3 * x ** (-2 / 3)
			gamV = (-3 * KT + 2 * Px * gt + 9 * KT * kkx - 6 * gt * KT) / 6 / (3 * KT - 2 * Px * gt) + gamVo
			ff2 = (x + 0.00001) ** (1 / 3)
			KT2 = K_0 * 1000 / ff2 ** 6 * np.exp(fw * (1 - ff2)) * ((-5 / ff2 ** 2 + 4 / ff2) * (1 + aa * ff2 - aa * ff2 ** 2) + (1 / ff2 - 1) * (1 + aa * ff2 - aa * ff2 ** 2) * (-fw) + (1 / ff2 - 1) * (aa - 2 * aa * ff2)) * (-(x + 0.00001))
			ex = 3 / ff2 ** 4 * K_0 * 1000 * np.exp(fw * (1 - ff2)) * ((-5 / ff2 ** 2 + 4 / ff2) * (1 + aa * ff2 - aa * ff2 ** 2) - (1 / ff2 - 1) * (1 + aa * ff2 - aa * ff2 ** 2) * fw + (1 / ff2 - 1) * (aa - 2 * aa * ff2))
			ex1 = 1 / ff2 ** 3 * K_0 * 1000 * np.exp(fw * (1 - ff2)) * fw * ((-5 / ff2 ** 2 + 4 / ff2) * (1 + aa * ff2 - aa * ff2 ** 2) - (1 / ff2 - 1) * (1 + aa * ff2 - aa * ff2 ** 2) * fw + (1 / ff2 - 1) * (aa - 2 * aa * ff2))
			ex2 = K_0 * 1000 * np.exp(fw * (1 - ff2)) * ((10 / ff2 ** 3 - 4 / ff2 ** 2) * (1 + aa * ff2 - aa * ff2 ** 2) + (-5 / ff2 ** 2 + 4 / ff2) * (aa - 2 * aa * ff2) + fw / ff2 ** 2 * (1 + aa * ff2 - aa * ff2 ** 2) - (1 / ff2 - 1) * (aa - 2 * aa * ff2) * fw - (aa - 2 * aa * ff2) / ff2 ** 2 - 2 * aa * (1 / ff2 - 1)) / ff2 ** 3
			kkx2 = (ex + ex1 - ex2) / (-KT2 / ff2) / 3
			ff1 = (x - 0.00001) ** (1 / 3)
			KT1 = K_0 * 1000 / ff1 ** 6 * np.exp(fw * (1 - ff1)) * ((-5 / ff1 ** 2 + 4 / ff1) * (1 + aa * ff1 - aa * ff1 ** 2) + (1 / ff1 - 1) * (1 + aa * ff1 - aa * ff1 ** 2) * (-fw) + (1 / ff1 - 1) * (aa - 2 * aa * ff1)) * (-(x - 0.00001))
			ex = 3 / ff1 ** 4 * K_0 * 1000 * np.exp(fw * (1 - ff1)) * ((-5 / ff1 ** 2 + 4 / ff1) * (1 + aa * ff1 - aa * ff1 ** 2) - (1 / ff1 - 1) * (1 + aa * ff1 - aa * ff1 ** 2) * fw + (1 / ff1 - 1) * (aa - 2 * aa * ff1))
			ex1 = 1 / ff1 ** 3 * K_0 * 1000 * np.exp(fw * (1 - ff1)) * fw * ((-5 / ff1 ** 2 + 4 / ff1) * (1 + aa * ff1 - aa * ff1 ** 2) - (1 / ff1 - 1) * (1 + aa * ff1 - aa * ff1 ** 2) * fw + (1 / ff1 - 1) * (aa - 2 * aa * ff1))
			ex2 = K_0 * 1000 * np.exp(fw * (1 - ff1)) * ((10 / ff1 ** 3 - 4 / ff1 ** 2) * (1 + aa * ff1 - aa * ff1 ** 2) + (-5 / ff1 ** 2 + 4 / ff1) * (aa - 2 * aa * ff1) + fw / ff1 ** 2 * (1 + aa * ff1 - aa * ff1 ** 2) - (1 / ff1 - 1) * (aa - 2 * aa * ff1) * fw - (aa - 2 * aa * ff1) / ff1 ** 2 - 2 * aa * (1 / ff1 - 1)) / ff1 ** 3
			kkx1 = (ex + ex1 - ex2) / (-KT1 / ff1) / 3
			dkkdx = (kkx2 - kkx1) / (0.00002)
			qV = -1 / 2 * (-6 * KT * kkx ** 2 * Px * gt / x + 6 * KT * dkkdx * Px * gt + 6 * KT * kkx * KT * gt / x - 6 * KT * kkx * Px * gtx + 4 * gt ** 2 * KT * kkx * Px / x - 4 * gt ** 2 * KT ** 2 / x - 9 * KT ** 2 * dkkdx + 6 * gtx * KT ** 2) / (3 * KT - 2 * Px * gt) ** 2 * x / gamV
			QB = QBo * expp
			QB1 = QB1o * expp
			theta_1 = QE1o * expp
			theta_2 = QE2o * expp
			VV = V_0 * x
			gg = d * np.log(1 + QB / Tr / d)
			ex = np.exp(QB / Tr)
			b = 1 / (np.exp(gg) - 1)
			KTB = mb * R * ((Tr * QB * d * b / (Tr * d + QB)) * gamV / VV * (1 + gamV - qV) - gamV ** 2 * Tr / VV * (QB * d / (Tr * d + QB)) ** 2 * b * (1 / d + 1 + b))
			dPthBr = (KTB)
			gg = d * np.log(1 + QB / TK / d)
			ex = np.exp(QB / TK)
			b = 1 / (np.exp(gg) - 1)
			KTB = mb * R * ((TK * QB * d * b / (TK * d + QB)) * gamV / VV * (1 + gamV - qV) - gamV ** 2 * TK / VV * (QB * d / (TK * d + QB)) ** 2 * b * (1 / d + 1 + b))
			dPthB = (KTB)
			gg = d1 * np.log(1 + QB1 / Tr / d1)
			ex = np.exp(QB1 / Tr)
			b = 1 / (np.exp(gg) - 1)
			KTB = mb1 * R * ((Tr * QB1 * d1 * b / (Tr * d1 + QB1)) * gamV / VV * (1 + gamV - qV) - gamV ** 2 * Tr / VV * (QB1 * d1 / (Tr * d1 + QB1)) ** 2 * b * (1 / d1 + 1 + b))
			dPthBBr = (KTB)
			gg = d1 * np.log(1 + QB1 / TK / d1)
			ex = np.exp(QB1 / TK)
			b = 1 / (np.exp(gg) - 1)
			KTB = mb1 * R * ((TK * QB1 * d1 * b / (TK * d1 + QB1)) * gamV / VV * (1 + gamV - qV) - gamV ** 2 * TK / VV * (QB1 * d1 / (TK * d1 + QB1)) ** 2 * b * (1 / d1 + 1 + b))
			dPthBB = (KTB)
			e1 = np.exp(theta_1 / TK)
			KTB = Ein_1 * R * ((theta_1 / (e1 - 1)) * gamV / VV * (1 + gamV - qV) - gamV ** 2 * TK / VV * (theta_1 / TK) ** 2 * e1 / (e1 - 1) ** 2)
			dPthE1 = (KTB)
			e2 = np.exp(theta_2 / TK)
			KTB = Ein_2 * R * ((theta_2 / (e2 - 1)) * gamV / VV * (1 + gamV - qV) - gamV ** 2 * TK / VV * (theta_2 / TK) ** 2 * e2 / (e2 - 1) ** 2)
			dPthE2 = (KTB)
			e1 = np.exp(theta_1 / Tr)
			KTB = Ein_1 * R * ((theta_1 / (e1 - 1)) * gamV / VV * (1 + gamV - qV) - gamV ** 2 * Tr / VV * (theta_1 / Tr) ** 2 * e1 / (e1 - 1) ** 2)
			dPthE1r = (KTB)
			e2 = np.exp(theta_2 / Tr)
			KTB = Ein_2 * R * ((theta_2 / (e2 - 1)) * gamV / VV * (1 + gamV - qV) - gamV ** 2 * Tr / VV * (theta_2 / Tr) ** 2 * e2 / (e2 - 1) ** 2)
			dPthE2r = (KTB)
			KTa = 3 / 2 * n * R * a_0 / 1000000 * x ** (m) * (m) / VV * (TK ** 2 - Tr ** 2) * (1 - m)
			KTe = 3 / 2 * n * R * ae / 1000000 * x ** (mm) * (mm) / VV * (TK ** 2 - Tr ** 2) * (1 - mm)
			dPdV = dPthB - dPthBr + dPthBB - dPthBBr + dPthE1 - dPthE1r + dPthE2 - dPthE2r + KTe + KTa
			KTth = dPdV

			return KTth

		def dPdTth(n, Z, V_0, K_0, K_prime, Tr, x, expp, V, QBo, d, mb, QB1o, d1, mb1, QE1o, Ein_1, QE2o, Ein_2, TK, gamVo, beta, gb, a_0, m, mm, ae):
			R = 8.31451
			fw = -np.log(3 * K_0 / 10 / (1003.6 * (Z * n / (V_0 * 10)) ** (5 / 3)))
			ff = x ** (1 / 3)
			aa = 1.5 * (K_prime - 3) - fw
			Px = 3 * K_0 * 1000 * np.exp(fw * (1 - ff)) * (1 / ff ** 5 - 1 / ff ** 4) * (1 + aa * ff - aa * ff ** 2)
			KT = K_0 * 1000 / ff ** 6 * np.exp(fw * (1 - ff)) * ((-5 / ff ** 2 + 4 / ff) * (1 + aa * ff - aa * ff ** 2) + (1 / ff - 1) * (1 + aa * ff - aa * ff ** 2) * (-fw) + (1 / ff - 1) * (aa - 2 * aa * ff)) * (-x)
			ex = 3 / ff ** 4 * K_0 * 1000 * np.exp(fw * (1 - ff)) * ((-5 / ff ** 2 + 4 / ff) * (1 + aa * ff - aa * ff ** 2) - (1 / ff - 1) * (1 + aa * ff - aa * ff ** 2) * fw + (1 / ff - 1) * (aa - 2 * aa * ff))
			ex1 = 1 / ff ** 3 * K_0 * 1000 * np.exp(fw * (1 - ff)) * fw * ((-5 / ff ** 2 + 4 / ff) * (1 + aa * ff - aa * ff ** 2) - (1 / ff - 1) * (1 + aa * ff - aa * ff ** 2) * fw + (1 / ff - 1) * (aa - 2 * aa * ff))
			ex2 = K_0 * 1000 * np.exp(fw * (1 - ff)) * ((10 / ff ** 3 - 4 / ff ** 2) * (1 + aa * ff - aa * ff ** 2) + (-5 / ff ** 2 + 4 / ff) * (aa - 2 * aa * ff) + fw / ff ** 2 * (1 + aa * ff - aa * ff ** 2) - (1 / ff - 1) * (aa - 2 * aa * ff) * fw - (aa - 2 * aa * ff) / ff ** 2 - 2 * aa * (1 / ff - 1)) / ff ** 3
			kkx = (ex + ex1 - ex2) / (-KT / ff) / 3
			gt = gb - beta * x ** (1 / 3)
			gamV = (-3 * KT + 2 * Px * gt + 9 * KT * kkx - 6 * gt * KT) / 6 / (3 * KT - 2 * Px * gt) + gamVo
			ff2 = (x + 0.00001) ** (1 / 3)
			ex = 3 / ff2 ** 4 * K_0 * 1000 * np.exp(fw * (1 - ff2)) * ((-5 / ff2 ** 2 + 4 / ff2) * (1 + aa * ff2 - aa * ff2 ** 2) - (1 / ff2 - 1) * (1 + aa * ff2 - aa * ff2 ** 2) * fw + (1 / ff2 - 1) * (aa - 2 * aa * ff2))
			ex1 = 1 / ff2 ** 3 * K_0 * 1000 * np.exp(fw * (1 - ff2)) * fw * ((-5 / ff2 ** 2 + 4 / ff2) * (1 + aa * ff2 - aa * ff2 ** 2) - (1 / ff2 - 1) * (1 + aa * ff2 - aa * ff2 ** 2) * fw + (1 / ff2 - 1) * (aa - 2 * aa * ff2))
			ex2 = K_0 * 1000 * np.exp(fw * (1 - ff2)) * ((10 / ff2 ** 3 - 4 / ff2 ** 2) * (1 + aa * ff2 - aa * ff2 ** 2) + (-5 / ff2 ** 2 + 4 / ff2) * (aa - 2 * aa * ff2) + fw / ff2 ** 2 * (1 + aa * ff2 - aa * ff2 ** 2) - (1 / ff2 - 1) * (aa - 2 * aa * ff2) * fw - (aa - 2 * aa * ff2) / ff2 ** 2 - 2 * aa * (1 / ff2 - 1)) / ff2 ** 3
			ff1 = (x - 0.00001) ** (1 / 3)
			ex = 3 / ff1 ** 4 * K_0 * 1000 * np.exp(fw * (1 - ff1)) * ((-5 / ff1 ** 2 + 4 / ff1) * (1 + aa * ff1 - aa * ff1 ** 2) - (1 / ff1 - 1) * (1 + aa * ff1 - aa * ff1 ** 2) * fw + (1 / ff1 - 1) * (aa - 2 * aa * ff1))
			ex1 = 1 / ff1 ** 3 * K_0 * 1000 * np.exp(fw * (1 - ff1)) * fw * ((-5 / ff1 ** 2 + 4 / ff1) * (1 + aa * ff1 - aa * ff1 ** 2) - (1 / ff1 - 1) * (1 + aa * ff1 - aa * ff1 ** 2) * fw + (1 / ff1 - 1) * (aa - 2 * aa * ff1))
			ex2 = K_0 * 1000 * np.exp(fw * (1 - ff1)) * ((10 / ff1 ** 3 - 4 / ff1 ** 2) * (1 + aa * ff1 - aa * ff1 ** 2) + (-5 / ff1 ** 2 + 4 / ff1) * (aa - 2 * aa * ff1) + fw / ff1 ** 2 * (1 + aa * ff1 - aa * ff1 ** 2) - (1 / ff1 - 1) * (aa - 2 * aa * ff1) * fw - (aa - 2 * aa * ff1) / ff1 ** 2 - 2 * aa * (1 / ff1 - 1)) / ff1 ** 3
			QB = QBo * expp
			QB1 = QB1o * expp
			theta_1 = QE1o * expp
			theta_2 = QE2o * expp
			gg = d * np.log(1 + QB / TK / d)
			ex = np.exp(QB / TK)
			b = 1 / (np.exp(gg) - 1)
			dPB = mb * R * (gamV / V * (QB * d / (TK * d + QB)) ** 2 * b * (1 / d + 1 + b))
			gg = d1 * np.log(1 + QB1 / TK / d1)
			ex = np.exp(QB1 / TK)
			b = 1 / (np.exp(gg) - 1)
			dPB1 = mb1 * R * (gamV / V * (QB1 * d1 / (TK * d1 + QB1)) ** 2 * b * (1 / d1 + 1 + b))
			e1 = np.exp(theta_1 / TK)
			dP1 = Ein_1 * R * (gamV / V * (theta_1 / TK) ** 2 * e1 / (e1 - 1) ** 2)
			e2 = np.exp(theta_2 / TK)
			dP2 = Ein_2 * R * (gamV / V * (theta_2 / TK) ** 2 * e2 / (e2 - 1) ** 2)
			dPdTa = 3 * n * R * a_0 / 1000000 * x ** (m) * (TK) * m / V
			dPdTe = 3 * n * R * ae / 1000000 * x ** (mm) * (TK) * mm / V
			dPdTth = dPB + dPB1 + dP1 + dP2 + dPdTe + dPdTa

			return dPdTth

		def derive_params():
			self.Ex = I(self.loaded_params["Z"],
					self.loaded_params["n"],
					self.x,
					self.loaded_params["V_0"]* 6.02214 / 400,
					self.loaded_params["K_0"]*10,
					self.loaded_params["K_prime"]
					)
			
			self.exp = np.exp(I_gamV(self.loaded_params["n"],
					self.loaded_params["Z"],
					self.x,
					self.loaded_params["V_0"]* 6.02214 / 400,
					self.loaded_params["K_0"]*10,
					self.loaded_params["K_prime"],
					self.loaded_params["delta"],
					self.loaded_params["t"],
					0))
			
			self.free_energy = F(self.loaded_params["n"],
					self.loaded_params["Z"],
					self.loaded_params["V_0"]* 6.02214 / 400,
					self.loaded_params["K_0"]*10,
					self.loaded_params["K_prime"],
					self.loaded_params["T_0"],
					self.x,
					self.exp,
					1, 1, 0, 1, 1, 0,
					self.loaded_params["theta_1"],
					self.loaded_params["Ein_1"],
					self.loaded_params["theta_2"],
					self.loaded_params["Ein_2"],
					self.T,
					self.loaded_params["delta"],
					0,
					self.loaded_params["t"],
					self.loaded_params["a_0"],
					self.loaded_params["m"],
					self.loaded_params["g"],
					self.loaded_params["e_0"]
					) + self.Ex + self.loaded_params["U_0"]
			
			self.entropy = S(self.loaded_params["n"],
					self.loaded_params["Z"],
					self.loaded_params["V_0"]* 6.02214 / 400,
					self.loaded_params["K_0"]*10,
					self.loaded_params["K_prime"],
					self.loaded_params["T_0"],
					self.x,
					self.exp,
					1, 1, 0, 1, 1, 0,
					self.loaded_params["theta_1"],
					self.loaded_params["Ein_1"],
					self.loaded_params["theta_2"],
					self.loaded_params["Ein_2"],
					self.T,
					self.loaded_params["delta"],
					0,
					self.loaded_params["t"],
					self.loaded_params["a_0"],
					self.loaded_params["m"],
					self.loaded_params["g"],
					self.loaded_params["e_0"]
					)
			
			self.in_energy = self.free_energy + self.entropy * self.T

			self.p_x = P(self.loaded_params["Z"],
					self.loaded_params["n"],
					self.x,
					self.loaded_params["V_0"]* 6.02214 / 400,
					self.loaded_params["K_0"]*10,
					self.loaded_params["K_prime"])
			
			self.Cv = CCv(self.loaded_params["n"],
				    self.loaded_params["Z"],
					self.loaded_params["V_0"]* 6.02214 / 400,
					self.loaded_params["K_0"]*10,
					self.loaded_params["K_prime"],
					self.loaded_params["T_0"],
					self.x,
					self.exp,
					1, 1, 0, 1, 1, 0,
					self.loaded_params["theta_1"],
					self.loaded_params["Ein_1"],
					self.loaded_params["theta_2"],
					self.loaded_params["Ein_2"],
					self.T,
					self.loaded_params["delta"],
					0,
					self.loaded_params["t"],
					self.loaded_params["a_0"],
					self.loaded_params["m"],
					self.loaded_params["g"],
					self.loaded_params["e_0"]
					)
			
			self.p_thermal = Pth(self.loaded_params["n"],
					self.loaded_params["Z"],
					self.loaded_params["V_0"]* 6.02214 / 400,
					self.loaded_params["K_0"]*10,
					self.loaded_params["K_prime"],
					self.loaded_params["T_0"],
					self.x,
					self.exp,
					(self.Vmol/10),
					1, 1, 0, 1, 1, 0,
					self.loaded_params["theta_1"],
					self.loaded_params["Ein_1"],
					self.loaded_params["theta_2"],
					self.loaded_params["Ein_2"],
					self.T,
					self.loaded_params["delta"],
					0,
					self.loaded_params["t"],
					self.loaded_params["a_0"],
					self.loaded_params["m"],
					self.loaded_params["g"],
					self.loaded_params["e_0"])
			
			self.p_real = self.p_x + self.p_thermal

			self.gibbs_energy = self.free_energy + self.p_real * (self.Vmol/10)

			self.enthalpy = self.in_energy + self.p_real * (self.Vmol / 10)

			self.K_T_thermal = KTth(self.loaded_params["n"],
					self.loaded_params["Z"],
					self.loaded_params["V_0"]* 6.02214 / 400,
					self.loaded_params["K_0"]*10,
					self.loaded_params["K_prime"],
					self.loaded_params["T_0"],
					self.x,
					self.exp,
					(self.Vmol/10),
					1, 1, 0, 1, 1, 0,
					self.loaded_params["theta_1"],
					self.loaded_params["Ein_1"],
					self.loaded_params["theta_2"],
					self.loaded_params["Ein_2"],
					self.T,
					self.loaded_params["delta"],
					0,
					self.loaded_params["t"],
					self.loaded_params["a_0"],
					self.loaded_params["m"],
					self.loaded_params["g"],
					self.loaded_params["e_0"]
					)

			self.ff = self.x ** (1/3)

			self.KTx = self.loaded_params["K_0"]*10 /self.ff ** 6 * np.exp(self.loaded_params["c_0"] * (1 - self.ff)) * ((-5 / self.ff ** 2 + 4 / self.ff)*(1 + self.loaded_params["c_2"] * self.ff - self.loaded_params["c_2"] * self.ff**2)+(1 / self.ff - 1)*(1 + self.loaded_params["c_2"] * self.ff - self.loaded_params["c_2"] * self.ff ** 2) * (- self.loaded_params["c_0"]) + (1 / self.ff - 1) * (self.loaded_params["c_2"] - 2 * self.loaded_params["c_2"] * self.ff)) * (- self.x) * 1000

			self.k_t = self.K_T_thermal + self.KTx

			self.dPdT = dPdTth(self.loaded_params["n"],
					self.loaded_params["Z"],
					self.loaded_params["V_0"]* 6.02214 / 400,
					self.loaded_params["K_0"]*10,
					self.loaded_params["K_prime"],
					self.loaded_params["T_0"],
					self.x,
					self.exp,
					self.Vmol,
					1, 1, 0, 1, 1, 0,
					self.loaded_params["theta_1"],
					self.loaded_params["Ein_1"],
					self.loaded_params["theta_2"],
					self.loaded_params["Ein_2"],
					self.T,
					self.loaded_params["delta"],
					0,
					self.loaded_params["t"],
					self.loaded_params["a_0"],
					self.loaded_params["m"],
					self.loaded_params["g"],
					self.loaded_params["e_0"]
					)

			self.alpha = self.dPdT / (self.K_T_thermal + self.KTx)

			self.Cp = self.Cv + self.alpha ** 2 * self.Vmol * self.T *(self.K_T_thermal + self.KTx)

			self.Ks = self.T + self.T * self.Vmol / self.Cv * self.dPdT ** 2

			self.gamma = self.alpha * self.Vmol * self.T / self.Cv


		self.P = self.PVT[0]
		self.V = self.PVT[1]
		self.T = self.PVT[2]

		if self.P == None and self.V != None and self.T != None:

			self.Vmol = self.V * 6.02214 / 40
			self.Vo_mol = self.loaded_params["V_0"] * 6.02214 / 400
			self.x = (self.Vmol /10) / self.Vo_mol

			derive_params()

			self.P = self.p_real / 10000

		if self.P != None and self.V == None and self.T != None:

			self.Vo_mol = self.loaded_params["V_0"] * 6.02214 / 400
			self.x = xAP2(self.loaded_params["n"],
					self.loaded_params["Z"],
					298.15,
					self.Vo_mol,
					self.loaded_params["K_0"]*10,
					self.loaded_params["K_prime"],
					1, 1, 0, 1, 1, 0,
					self.loaded_params["theta_1"],
					self.loaded_params["Ein_1"],
					self.loaded_params["theta_2"],
					self.loaded_params["Ein_2"],
					self.loaded_params["delta"],
					self.loaded_params["t"],
					0,
					self.loaded_params["a_0"],
					self.loaded_params["m"],
					self.loaded_params["g"],
					self.loaded_params["e_0"],
					self.T,
					self.P * 10000)
			
			V_Jbar = self.x * self.Vo_mol
			self.Vmol = V_Jbar * 10
			self.V = self.Vmol *40/6.02214
			
			derive_params()

		if self.P != None and self.V != None and self.T == None:
			self.Vo_mol = self.loaded_params["V_0"] * 6.02214 / 400

			guess_low = 0
			guess_high = 2000
			V = self.V
			V_0 = self.loaded_params["V_0"]
			self.x = V_0/V
			K_0 = self.loaded_params["K_0"]
			K_prime = self.loaded_params["K_prime"]
			alpha_298 = self.loaded_params["Therm_alpha_298"]
			diff_temp = self.loaded_params["Therm_diff_temp"]
			diff_alpha = self.loaded_params["Therm_diff_alpha"]
			Ko_298 = K_0

			x = self.x
			P_BM = (3*K_0 /2)*(x**(7/3)-x**(5/3))*(1+(3/4)*(K_prime-4)*(x**(2/3)-1))
			P_real = self.P
			y = P_real - P_BM

			a = alpha_298 * Ko_298
			b = alpha_298 * diff_temp
			c = Ko_298 * diff_alpha
			d = diff_alpha * diff_temp

			def equ(x):
				return ((a + b*x + c*x + d*x**2) * x) - y


			def solve_x(y):

				initial_guesses = list(range(guess_low, guess_high, 100)) #not actually necessary here, as there is only one solution
				x_roots = fsolve(equ, initial_guesses)
				real_roots = [root.real for root in x_roots if np.isclose(root.imag, 0)]
				return sorted(set(real_roots))

			y_value = y
			solutions = solve_x(y_value)

			self.T =  np.average(solutions)+298.15

			V_Jbar = self.x * self.Vo_mol
			self.Vmol = V_Jbar * 10
			
			derive_params()

