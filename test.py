import itertools as it

if __name__ == "__main__":

	x = it.imap(lambda x: x +1, [1,2,4])

	for a in x:
		print a