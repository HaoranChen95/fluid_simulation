from math import exp

def C(gh):
	result = 2. * gh - 3. + 4. * exp(-gh) - exp(-2. * gh)
	return result

def G(gh):
	result = exp(gh) - 2. * gh - exp(-gh)
	return result

def E_(gh):
    return -C(gh)*C(-gh)- G(gh) **2

def E(gh):
	result = (
    	16. * (exp(gh) + exp(-gh)) -
		4. * (exp(2. * gh) + exp(-2. * gh)) -
		4. * gh * (exp(gh) - exp(-gh)) +
		2. * gh * (exp(2. * gh) - exp(-2. * gh)) - 24.
	)
	return result

for i in range(-8, 3):
	print(i)
	print(E(10 ** i))
	print(E_(10 ** i))