# Description: Compute minimum m vs n for MRB success on McEliece or HQC instances.

import numpy as np
from math import *
import sys

if __name__ == "__main__":
	scheme = str(sys.argv[1]) if len(sys.argv) > 1 else "mceliece"
	kappa = float(sys.argv[2]) if len(sys.argv) > 2 else 0
	succ_probability = int(sys.argv[3]) if len(sys.argv) > 3 else 90

	# set d string for file name
	if kappa == 0:
		d_str = "d_0"
	elif kappa == 0.5:
		d_str = "d_thalf"
	elif kappa == 1:
		d_str = "d_t"
	
	# set parameters dependent on the scheme
	if scheme == "hqc":
		omega = 132 / sqrt(18669 * 2)
		R = 0.5
		c_list = [1.0 + 0.1 * i for i in range(0, 16)]
		parameters_list = [[
			1000 + (j - 1) * 1000, 
			c_list, 
			int(omega * sqrt(1000 + (j - 1) * 1000)), 
			int((1000 + (j - 1) * 1000) * R), int(omega * sqrt(1000 + (j - 1) * 1000)) * kappa
		] for j in range(1, 200, 2)]
	elif scheme == "mceliece":
		omega = 64 / 3488 * log2(3488)
		R = (3488 - 12 * 64) / 3488
		c_list = [1.9 + 0.1 * i for i in range(0, 26)]
		parameters_list = [[
			1000 + (j - 1) * 1000, 
			c_list, 
			int(omega * (1000 + (j - 1) * 1000) / log2(1000 + (j - 1) * 1000)), 
			int((1000 + (j - 1) * 1000) * R), int(omega * (1000 + (j - 1) * 1000) / log2(1000 + (j - 1) * 1000)) * kappa
		] for j in range(1, 50, 2)]

	rng = np.random.default_rng()
	with open("min_m_for_success_" + scheme + "_" + d_str + "_prob_" + str(succ_probability) + ".txt", "a") as f:
		num_repetitions = 5
		num_trials = 100
		max_errors = int(100 - succ_probability)
		best_m = []
		f.write(f"scheme = {scheme}, num_trials = {num_trials}\n")
		for parameters in parameters_list:
			n, c_list, t, k, d = parameters
			for c in c_list:
				m = int(c * (t + 2 * d) * log2(t))
				for _ in range(num_repetitions):
					count_errors = 0
					for T in range(num_trials):

						# generate instance
						scores_D0 = -d * m + rng.binomial(m*(t + 2*d), 1/2, size = n - t)
						scores_D1 = -(d - 1) * m + rng.binomial(m*(t - 1 + 2*d), 1/2, size = t)
						scores = np.concatenate((scores_D1, scores_D0))

						# sort the probabilities
						sorted_scores = np.flip(np.sort(scores))
						sorting = np.flip(np.argsort(scores))

						# check success
						success = np.isin(range(0, t), sorting[:n - k]).all()
						# DEBUG
						# 	print(f"Trial {T}: success = {success}")
						if not success:
							count_errors += 1
							if count_errors > max_errors:
								break

					print(f"n = {n}, k = {k}, t = {t}, m = {c} * t * log2(t) = {m}, success rate: {num_trials - count_errors}/{num_trials}")
					if count_errors <= max_errors:
						break
				if count_errors <= max_errors:
					print(f"Achieved desired probability for m = {m}")
					break
			if not success:
				print(f"For n = {n}, m <= {c} * t * log2(t) = {m} was not sufficient for the desired successes.")
				best_m.append((n, t, np.nan, d))
				print("----------")
				continue
			print(f"For n = {n}, minimum m with desired successes: m = {c} * t * log2(t) = {m}")
			f.write(f"({n}, {m})\n")
			best_m.append((n, t, m, d))
			print("----------")

		f.write("Corollary bound data:\n")
		for item in best_m:
			if not np.isnan(item[2]):
				f.write(f"({item[0]}, {int(2.0 * (item[1] + 2 * item[3]) * log2(item[1]))})\n")