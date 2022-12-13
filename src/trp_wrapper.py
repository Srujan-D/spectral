import numpy as np
import os
import optuna

from ctypes import c_double, c_int, POINTER, Structure, CDLL, byref

class Metrics(Structure):
	_fields_ = [
		("s_avg_acc", c_double),
		("l_avg_acc", c_double),
		("s_max_acc", c_double),
		("l_max_acc", c_double),
		("a_cost", c_double)
	]

	def __str__(self):
		return "\ns_avg_acc: {0}\nl_avg_acc: {1}".format(self._fields_[0][1].value, self._fields_[1][1].value)

class Params(Structure):
	_fields_ = [
		("s_acc_weight", c_double),
		("s_jerk_weight", c_double),
		("l_acc_weight", c_double),
		("l_jerk_weight", c_double),
		("weight_s_ref", c_double),
		("weight_ds_ref", c_double),
		("weight_l_ref", c_double),
		("weight_dl_ref", c_double),
		("weight_end_s", c_double),
		("weight_end_l", c_double),
		("iteration", c_int)
	]
	# {
	# double s_acc_weight;
	# double s_jerk_weight;
	# double l_acc_weight;
	# double l_jerk_weight;
# 
	# double weight_s_ref;
	# double weight_ds_ref;
	# double weight_l_ref;
	# double weight_dl_ref;
# };

cdll = CDLL("/home/srujan_d/RISS/code/btrapz/src/libtrp.so")

# _c_metrics_ =  POINTER(Metrics)

_run_btrapz = cdll.find_traj
_run_btrapz.argtypes = (
   POINTER(Params),
)

_run_btrapz.restype = c_double

def run_btrapz(trial): #, path = "/home/srujan_d/RISS/code/riss/src/btrapz/src/c1.txt"):
	# a = [trial.suggest_float('param_ub_'+str(i),0,100) for i in range(15)]

	# s_acc_weight = trial.suggest_float('s_acc_weight',a[0],a[1])
	# s_jerk_weight = trial.suggest_float('s_jerk_weight',a[2],a[3])
	# l_acc_weight = trial.suggest_float('l_acc_weight',a[4],a[5])
	# l_jerk_weight = trial.suggest_float('l_jerk_weight',a[6],a[7])

	# weight_s_ref = trial.suggest_float('weight_s_ref',a[8],a[9])
	# weight_ds_ref = trial.suggest_float('weight_ds_ref',a[10],a[11])
	# weight_l_ref = trial.suggest_float('weight_l_ref',a[12],a[13])
	# weight_dl_ref = trial.suggest_float('weight_dl_ref',a[14],a[15])

	s_acc_weight = trial.suggest_float('s_acc_weight',0,50)
	s_jerk_weight = trial.suggest_float('s_jerk_weight',0,50)
	l_acc_weight = trial.suggest_float('l_acc_weight',0,50)
	l_jerk_weight = trial.suggest_float('l_jerk_weight',0,50)

	weight_s_ref = trial.suggest_float('weight_s_ref',0,50)
	weight_ds_ref = trial.suggest_float('weight_ds_ref',0,50)
	weight_l_ref = trial.suggest_float('weight_l_ref',0,50)
	weight_dl_ref = trial.suggest_float('weight_dl_ref',0,50)

	weight_end_s = trial.suggest_float('weight_end_s',0,50)
	weight_end_l = trial.suggest_float('weight_end_l',0,50)



	# weight_file = '/home/srujan_d/RISS/code/btrapz/src/weights.txt'
	# f = open(weight_file)
	# lines = f.readlines()
	# s_acc_weight = float(lines[0])
	# s_jerk_weight = float(lines[1])
	# l_acc_weight = float(lines[2])
	# l_jerk_weight = float(lines[3])
	# weight_s_ref = float(lines[4])
	# weight_ds_ref = float(lines[5])
	# weight_l_ref = float(lines[6])
	# weight_dl_ref = float(lines[7])
	iteration = 1
	m = _run_btrapz(Params(s_acc_weight, s_jerk_weight, l_acc_weight, l_jerk_weight, weight_s_ref, weight_ds_ref, weight_l_ref, weight_dl_ref, weight_end_s, weight_end_l, iteration))
	return m

def find_traj():
	weight_file = '/home/srujan_d/RISS/code/btrapz/src/weights.txt'
	f = open(weight_file)
	lines = f.readlines()
	print(lines)
	lines = lines[0].split('\t')

	s_acc_weight = float(lines[0])
	s_jerk_weight = float(lines[1])
	l_acc_weight = float(lines[2])
	l_jerk_weight = float(lines[3])
	weight_s_ref = float(lines[4])
	weight_ds_ref = float(lines[5])
	weight_l_ref = float(lines[6])
	weight_dl_ref = float(lines[7])
	weight_end_s = float(lines[8])
	weight_end_l = float(lines[9])
	iteration = 3
	a = _run_btrapz(Params(s_acc_weight, s_jerk_weight, l_acc_weight, l_jerk_weight, weight_s_ref, weight_ds_ref, weight_l_ref, weight_dl_ref, weight_end_s, weight_end_l, iteration))
	if a == 100000000000:
		return False
	else:
		return True
	# print(a)





# if __name__ == '__main__':
# 	print("LETS START")

# 	# study = optuna.create_study(direction='minimize')
# 	# study.optimize(run_btrapz, n_trials=100)
# 	# trial = study.best_trial

# 	# # m = run_btrapz()

# 	# weight_file = '/home/srujan_d/RISS/code/btrapz/src/weights.txt'
# 	# f = open(weight_file, 'w')

# 	# all_weights = '/home/srujan_d/RISS/code/btrapz/src/all_weights.txt'
# 	# f2 = open(all_weights, 'a')

# 	# for key in trial.params:
# 	# 	f.write(str(round(trial.params[key], 2)))
# 	# 	f.write('\t')
# 	# 	f2.write(str(round(trial.params[key], 2)))
# 	# 	f2.write('\t')

# 	# f2.write('\n')
# 	# f.close()
# 	# f2.close()
# 	# print(round(trial.value,2), trial.params)


# 	weight_file = '/home/srujan_d/RISS/code/btrapz/src/weights.txt'
# 	f = open(weight_file)
# 	lines = f.readlines()
# 	print(lines)
# 	lines = lines[0].split('\t')

# 	s_acc_weight = float(lines[0])
# 	s_jerk_weight = float(lines[1])
# 	l_acc_weight = float(lines[2])
# 	l_jerk_weight = float(lines[3])
# 	weight_s_ref = float(lines[4])
# 	weight_ds_ref = float(lines[5])
# 	weight_l_ref = 0 #float(lines[6])
# 	weight_dl_ref = 0 #float(lines[7])
# 	weight_end_s = float(lines[8])
# 	weight_end_l = float(lines[9])
# 	print("values of params being used:")
# 	print(s_acc_weight, s_jerk_weight, l_acc_weight, l_jerk_weight, weight_s_ref, weight_ds_ref, weight_l_ref, weight_dl_ref, weight_end_s, weight_end_l)
# 	print("running btrapz")
# 	iteration = 31
# 	a = _run_btrapz(Params(s_acc_weight, s_jerk_weight, l_acc_weight, l_jerk_weight, weight_s_ref, weight_ds_ref, weight_l_ref, weight_dl_ref, weight_end_s, weight_end_l, iteration))
# 	print(a)