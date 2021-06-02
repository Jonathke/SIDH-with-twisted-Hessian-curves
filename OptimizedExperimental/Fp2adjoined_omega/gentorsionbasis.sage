p = 6143

K.<w> = GF(p**2, name='w', modulus=x^2+x+1)

x = PolynomialRing(K, 'x').gen()

j = 5736 + 2534*w
a = 1
f = (x^4 + 216*a*x)^3 - j*a*((x^3-27*a)^3)
print(f.roots())
d = f.roots()[-1][0]
print("a = {}".format(a))
print("d = {}".format(d))

def genPoint(P):
	l = []
	while not l:
		px = randint(0, p) + randint(0, p)*w
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
