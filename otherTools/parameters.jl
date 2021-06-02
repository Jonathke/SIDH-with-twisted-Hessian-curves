using Primes

# These primes are good example primes, as they have p = ±4 (mod 9), so taking cube roots is easy in F_p^2
plist = [5, 13, 23, 31, 41, 59, 67]
bigplist = [big(2)^107 - 1, big(2)^521 - 1]

# The prime used in the OG SIDH paper. Does not have above property^
OGp = big(2)^63*big(3)^41*11-1

# If p = 2 (mod 3) we have an easy supersingular curve
# E = Affine_EC_TwiHes(F, inv(4*one(F)), zero(F))
# Or twisted to
# E = Affine_EC_TwiHes(F, sqrt(F(b))*inv(4*one(F)), zero(F))
# This is isomorphic to the Weierstrass curve y^2 = x^3 + 1

# Similarly (slightly worse though) If p = 3 (mod 4) we have an easy supersingular curve
# y^2 = x^3 - x
# Translating this to twisted Hessian form will be future work

# If p = 1 (mod 12), just run away, dont look back.


#Scratch all this, much easier than that...
#If p = 2 (mod 3), then 0 is a supersingular j-invariant
#If p = 3 (mod 4), then 1728 is a supersingular j-invariant


#Example:
OGp = big(2)^63*big(3)^41*11-1
F = Fpsquare(OGp)
E = TwiHesEC(F, inv(4*one(F)), zero(F))
#w = F([1850222081870264162797200520152908562431, 2673953603825912900316077551548884232036])
#MAYBE: all points of order two are of the form (X : 1 : 1)
#Want to find roots of x^3-a = 0, where a = 3700444163740528325594401040305817124855
#One lies in F_p, the other two are the roots of x^2 + 3700444163740528325594401040305817124861*x + 4
#Roots are r1 = 3700444163740528325594401040305817124862 + 1647463043911297475037754062791951339209*i
#and r2 = 3700444163740528325594401040305817124862 + 2052981119829230850556646977513865785654*i
#r1 = F([3700444163740528325594401040305817124862, 1647463043911297475037754062791951339209])
#P = E(-r1, one(F), one(F))
#E0 = Isogeny(E,P).codomain
#c = F([3398654792751657841130893587213062545659])*w
#Q = E0(zero(F), -w, one(F))
#E1 = Isogeny(E0, Q).codomain

#Simple example with isogeny over F_p
F = Zmod(31)
E = TwiHesEC(F, F(1), F(7))
listP = allPointsNaive(E)
phi = Isogeny(E, E(F(13), F(1), F(1)))


#Simple example with isogeny over F_p^2
F = Fpsquare(31)
E = TwiHesEC(F, F([4, 5]), F([1,2]))
listP = allPointsNaive(E)
phi = Isogeny(E, E(F([9,19]), F([1]), F([1])))
