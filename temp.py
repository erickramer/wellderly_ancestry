import os
import time
import numpy as np

def temp_file(temp_dir="./temp"):


	if not os.path.exists(temp_dir):
		os.mkdir(temp_dir)

	base = abs(hash(time.time() + np.random.randint(0, 10000)))
	return "%s/%s" % (temp_dir, base)