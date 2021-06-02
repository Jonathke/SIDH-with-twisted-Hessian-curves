p = 8975115180965832301039146585311033125034032407066470151659120543496400979689183960070681253859283762908345682362367

K.<i> = GF(p**2, name='i', modulus=x^2+1)

x = PolynomialRing(K, 'x').gen()

j = 5078601030025936612951094422122457493757846004081090863385038905475160431769027333596329668623319709636198670796335 + 730435898302185493450854427674126914241218266865584295975053215582853664936059337902183419108185317207844682910681*i
a = 1
f = (x^4 + 216*a*x)^3 - j*a*((x^3-27*a)^3)
print(f.roots())
d = f.roots()[-1][0]
print("a = {}".format(a))
print("d = {}".format(d))

def genPoint(P):
	l = []
	while not l:
		px = randint(0, p) + randint(0, p)*i
		f = a*px^3 + x^3 + 1 - d*px*x
		l = f.roots()
	y = l[0][0]
	print("{}.X = {}".format(P, px))
	print("{}.Y = {}".format(P, y))

print()
genPoint("aliceP")
print()
genPoint("aliceQ")
print()
genPoint("bobP")
print()
genPoint("bobQ")
