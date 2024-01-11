--Check rank 4 and 5 in C^3*^3*C^3
restart
kk=ZZ/3001
S=kk[a,b,c]
vects=apply(3, i->random(S^1,S^{3: 0}))
A=(transpose vects_1)*vects_2;
A1=sub(matrix{{1,0,0},{0,0,0},{0,0,0}},S);
A2=sub(matrix{{0,0,0},{0,1,0},{0,0,0}},S);
A3=sub(matrix{{0,0,0},{0,0,0},{0,0,1}},S);
u=vects_0;
--pick a general rank 4 tensor: e1*e1*e1+e2*e2*e2+e3*e3*e3+v*w*z
T4=a*A1+b*A2+c*A3+u_(0,0)*a*A+u_(0,1)*b*A+u_(0,2)*c*A--tensor
L=ideal flatten entries T
rel=syz(gens L, DegreeLimit=>1);
R=kk[r_0..r_8]
relR=sub(rel,R)
Eq=saturate ideal ((vars R)*relR);--embedding of the tensor in P^8
M=matrix{{r_0..r_2},{r_3..r_5},{r_6..r_8}}
ElipCurve=(Eq+ideal det M);--determinantal elliptic curve 
mp=flatten entries gens minors(2,M);
Rl=R/Eq
mingens sub(ideal mp,Rl)--6 dimensional--Veronses Embedding of P^2



T=kk[s_0..s_8]
f=map(R,T,mp)--the invesrion map
J=preimage_f(Eq)-- the image surface in P^8
dim J, degree J--(3,4)
betti res J, betti res radical J
time Delta=minors(6,jacobian gens J)+J;--singular locus of the surface
dim Delta, degree Delta
----intersection with locus of rk<2 locus
N=matrix{{s_0..s_2},{s_3..s_5},{s_6..s_8}};
rk2=det N;
J2=J+rk2;
dim J2, degree J2, genus J2
betti res J2, betti res radical J2
--substitution in P^5
linear=ideal (gens J2)_{0,1,2}
TL=T/sub(linear,T)
TLL=kk[flatten entries basis(1,TL)]
J2L=sub(sub(J2,TL),TLL)
betti res J2L, betti res radical J2L
dim J2L, degree J2L, genus J2L
time delta=minors(4,jacobian gens J2L)+J2L;--singular
dim delta, degree delta, genus delta
betti res delta
#decompose delta
---intersection with rk=1 locus
rk1=minors(2,N);
J1=J+rk1;--is image of ElipCurve in P^2
dim J1, degree J1, genus J1
betti res J1, betti res radical J1
J1==sub(radical J2, ring J1)--the two curves obtaied from the intersection are the same, set-theoretically
--substitution in P^5
linear=ideal (gens J)_{0,1,2}
TL=T/sub(linear,T)
TLL=kk[flatten entries basis(1,TL)]
J1L=sub(sub(J,TL),TLL)
betti res J1L, betti res radical J1L
J1L==sub(radical J2L, ring J1L)--true
dim J1L, degree J1L, genus J1L
time delta=minors(4,jacobian gens J1L)+J1L;--singular
dim delta, degree delta

----compute the degree for a general rk 5 tensor
vect1=apply(3, i->random(S^1,S^{3: 0}))
vect2=apply(3, i->random(S^1,S^{3: 0}))
A=(transpose vect1_1)*vect1_2;
B=(transpose vect2_1)*vect2_2;
A1=sub(matrix{{1,0,0},{0,0,0},{0,0,0}},S);
A2=sub(matrix{{0,0,0},{0,1,0},{0,0,0}},S);
A3=sub(matrix{{0,0,0},{0,0,0},{0,0,1}},S);
u=vect1_0;
v=vect2_0;
---pick rk 5 tensor
T5=a*A1+b*A2+c*A3+u_(0,0)*a*A+u_(0,1)*b*A+u_(0,2)*c*A+v_(0,0)*a*B+v_(0,1)*b*B+v_(0,2)*c*B
L=ideal flatten entries T5;
rel=syz(gens L, DegreeLimit=>1);
R=kk[r_0..r_8]
relR=sub(rel,R)
Eq5=ideal ((vars R)*relR);--embedding of the tensor in P^8
M=matrix{{r_0..r_2},{r_3..r_5},{r_6..r_8}}
ElipCurve=(Eq+ideal det M);--determinantal elliptic curve 
mp=flatten entries gens minors(2,M);
Rl=R/Eq5
mingens sub(ideal mp,Rl)--6 dimensional--Veronses Embedding of P^2



M=matrix{{r_0..r_2},{r_3..r_5},{r_6..r_8}}
mp=flatten entries gens minors(2,M);
T=kk[s_0..s_8]
f=map(R,T,mp)--inversion map
F=saturate(preimage_f(Eq5));--the image surface
dim F, degree F
betti res F, betti res radical F
time Delta=minors(6,jacobian gens J)+J;--singular locus
dim Delta, degree Delta

********************************C4*C4*C4
restart
kk=ZZ/101
S=kk[a,b,c,d]
UVW=apply(6, i->apply(3, i->random(S^1,S^{4: 0})));
ABCD=apply(6, i->transpose (UVW_i)_1*(UVW_i)_2);
T6=sum apply(6, i->(matrix{{a,b,c,d}}*(transpose (UVW_i)_0))_(0,0)*ABCD_i);
L=ideal flatten entries T6;
rel=syz(gens L, DegreeLimit=>1);
R=kk[r_0..r_15]
relR=sub(rel,R);
Eq=ideal ((vars R)*relR);--embedding of the tensor in P^8
M=matrix{{r_0..r_3},{r_4..r_7},{r_8..r_11},{r_12..r_15}}
mp=flatten entries gens minors(3,M);
T=kk[s_0..s_15]
f=map(R,T,mp)--inversion map
time J=preimage_f(Eq);--the image surface
-- used 38.2436 seconds
dim J, degree J--(4,27)
time res(J, FastNonminimal => true)
elapsedTime C = minimalBetti(B, DegreeLimit=>2, LengthLimit=>3)
---time Delta=minors(12,jacobian gens J)+J;--singular locus
---dim Delta, degree Delta
     -- 121.506 seconds elapsed

             0  1   2    3     4
o24 = total: 1 52 321 2562 10674
          0: 1  .   .    .     .
          1: . 52 236  210     .
          2: .  .  85 2352 10674
	  


---rk 7
vect1=apply(3, i->random(S^1,S^{4: 0}));
vect2=apply(3, i->random(S^1,S^{4: 0}));
vect3=apply(3, i->random(S^1,S^{4: 0}));
A=(transpose vect1_1)*vect1_2;
B=(transpose vect2_1)*vect2_2;
C=(transpose vect3_1)*vect3_2;
A1=sub(matrix{{1,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},S);
A2=sub(matrix{{0,0,0,0},{0,1,0,0},{0,0,0,0},{0,0,0,0}},S);
A3=sub(matrix{{0,0,0,0},{0,0,0,0},{0,0,1,0},{0,0,0,0}},S);
A4=sub(matrix{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,1}},S);
u=vect1_0;
v=vect2_0;
w=vect3_0;
---pick rk 7 tensor
T7=a*A1+b*A2+c*A3+d*A4+u_(0,0)*a*A+u_(0,1)*b*A+u_(0,2)*c*A+
u_(0,3)*d*A+v_(0,0)*a*B+v_(0,1)*b*B+v_(0,2)*c*B+v_(0,3)*d*B+
w_(0,0)*a*C+w_(0,1)*b*C+w_(0,2)*c*C+w_(0,3)*d*C;
L=ideal flatten entries T7;
rel=syz(gens L, DegreeLimit=>1);
R=kk[r_0..r_15]
relR=sub(rel,R)
Eq7=ideal ((vars R)*relR);--embedding of the tensor in P^15
M=matrix{{r_0..r_3},{r_4..r_7},{r_8..r_11},{r_12..r_15}}
mp=flatten entries gens minors(3,M);
T=kk[s_0..s_15]
f=map(R,T,mp)--inversion map
time F=preimage_f(Eq7);
dim F, degree F
elapsedTime C = minimalBetti(F, DegreeLimit=>2,  LengthLimit=>3)

             0  1   2    3     4
o27 = total: 1 52 256 2497 10674
          0: 1  .   .    .     .
          1: . 52 236  145     .
          2: .  .  20 2352 10674
	  
Nmp=apply(4,i->sub(product flatten entries M^{i},{M_(i,0)=>1}))
Nmpc=apply(4,i->sub(product flatten entries M_{i},{M_(3,i)=>1})) 
d1=toList product apply(4,i-> M_(i,i))
d2=toList product apply(4,i-> M_(i,3-i))	  
T=kk[s_0..s_7]	  
f=map(R,T,Nmp|Nmpc)	  

st=set flatten entries vars T
sb=subsets(st,3);
Li={}

time for i from 0 to 559 do(
    if dim(F+ideal toList sb_i)==1 then
    Li=append(Li, degree radical (F+ideal toList sb_i));
    )
Li

    
Inj=J+ideal(s_0,s_13,s_8);
dim Inj, degree Inj
degree radical Inj


Inf=F+ideal(s_0,s_13,s_8);
dim Inf, degree radical Inf
	  

L=ideal random(R^1,R^{1:-1});
P2=Eq+L;
(dim P2, degree P2)	  
Rl=R/P2
mingens sub(ideal mp,Rl)--is 3-uple embedding
V3P2=mingens preimage_f(P2);
betti V3P2
degree ideal V3P2, dim ideal V3P2

V3P2l=(V3P2)_{0..5};
betti oo
isSubset(F,V3P2)

(syz((mingens F)_{27}|V3P2, DegreeLimit=>2))_{15}

P=mingens F;
apply(rank source P, k->rank matrix apply(16,i->apply(16, j->((coefficients(P_{k}_(0,0),Monomials=>{s_i*s_j}))_1)_(0,0))))


####  intersection with rk 1 locus   ####
N=matrix{{s_0..s_3},{s_4..s_7},{s_8..s_11},{s_12..s_15}};
time rk1=minors(2,N);
F1=F+rk1;--is a surface
dim F1, degree F1
time minimalBetti(F1, DegreeLimit=>1,  LengthLimit=>5)
--pick a random hyperplane to intersect
B1=basis(1,T)
L1=ideal flatten entries (random(T^1,T^{16:0})*transpose B1)
IntCurve=F1+L1;
genus IntCurve
time minimalBetti(IntCurve, DegreeLimit=>1,  LengthLimit=>4)

             0  1   2   3   4
o36 = total: 1 63 402 953 901
          0: 1  1   .   .   .
          1: . 62 402 953 901
	  


--rk7
             0  1   2   3
o48 = total: 1 52 236 145
          0: 1  .   .   .
          1: . 52 236 145

--rk6
elapsedTime C = minimalBetti(J, DegreeLimit=>1, LengthLimit=>5)
 -- 10708.7 seconds elapsed
 
             0  1   2    3     4     5     6     7     8     9    10   11   12  13 14 15
o29 = total: 1 52 321 2562 10674 25960 42823 51700 47256 33188 17941 7384 2250 480 64  4
          0: 1  .   .    .     .     .     .     .     .     .     .    .    .   .  .  .
          1: . 52 236  210     .     .     .     .     .     .     .    .    .   .  .  .
          2: .  .  85 2352 10674 25960 42823 51700 47256 33188 17941 7384 2250 480 64  4

--rk 5
             0  1   2   3    4    5
o27 = total: 1 22 171 736 2046 3930
          0: 1  6  15  20   15    6
          1: . 11  81 255  445  465
          2: .  5  75 461 1586 3459
	  
	  
	  

----------------------C5*C5*C5
restart
kk=ZZ/1009
S=kk[a,b,c,d,e]
vect1=apply(3, i->random(S^1,S^{5: 0}))
vect2=apply(3, i->random(S^1,S^{5: 0}))
vect3=apply(3, i->random(S^1,S^{5: 0}))
A=(transpose vect1_1)*vect1_2;
B=(transpose vect2_1)*vect2_2;
C=(transpose vect3_1)*vect3_2;
A1=sub(matrix{{1,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},S);
A2=sub(matrix{{0,0,0,0,0},{0,1,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},S);
A3=sub(matrix{{0,0,0,0,0},{0,0,0,0,0},{0,0,1,0,0},{0,0,0,0,0},{0,0,0,0,0}},S);
A4=sub(matrix{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,1,0},{0,0,0,0,0}},S)
A5=sub(matrix{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,1}},S)
u=vect1_0;
v=vect2_0;
w=vect3_0;
---pick rk 8 tensor
T8=a*A1+b*A2+c*A3+d*A4+e*A5+u_(0,0)*a*A+u_(0,1)*b*A+u_(0,2)*c*A+
u_(0,3)*d*A+u_(0,4)*e*A+v_(0,0)*a*B+v_(0,1)*b*B+v_(0,2)*c*B+
v_(0,3)*d*B+v_(0,4)*e*B+w_(0,0)*a*C+w_(0,1)*b*C+w_(0,2)*c*C+w_(0,3)*d*C+w_(0,4)*e*C
L=ideal flatten entries T8;
rel=syz(gens L, DegreeLimit=>1);
R=kk[r_0..r_24]
relR=sub(rel,R)
Eq8=ideal ((vars R)*relR);
M=matrix{{r_0..r_4},{r_5..r_9},{r_10..r_14},{r_15..r_19},{r_20..r_24}}
mp=flatten entries gens minors(4,M);
T=kk[s_0..s_24]
f=map(R,T,mp)--inversion map
time J=preimage_f(Eq);
-- used 38.2436 seconds
dim J, degree J


P2=Eq+ideal random(R^1,R^{1:-1});	  
Rl=R/P2
mingens sub(ideal mp,Rl)--is 3-uple embedding
time V4P2=preimage_f(P2);
betti mingens V4P2
---rank 9 tensor
b3=sum flatten entries basis(3,R);
gaus3=flatten entries diff(vars R, b3);
Ru=kk[p_0..p_24]
g3=map(R,Ru, gaus3);
time I8=preimage_g3(Eq8);
betti I8
Il=ideal (mingens I8)_{0..10};
Rl=Ru/Il;
I8l=sub(I8,Rl);
betti oo

----n=5, rank 8
after cutting with 4 linear forms
     -- used 600.901 seconds

             0   1    2
o51 = total: 1 140 5014
          0: 1   .    .
          1: .  80  205
          2: .  60 4809



--cutting with three linear and getting curve--not good--all same

i14 : time C8=minimalBetti(B, DegreeLimit=>2,LengthLimit=>2)
    -- used 7230.1 seconds

             0   1    2     3
o15 = total: 1 140 4605 37345
          0: 1   .    .     .
          1: .  80  205     .
          2: .  60 4400 37345

--8--curve-- used 5696.69 seconds

             0   1    2     3
o10 = total: 1 140 5205 49739
          0: 1   .    .     .
          1: .  80  205     .
          2: .  60 5000 49739




--surface -- used 25566.4 seconds

             0   1    2     3
o10 = total: 1 140 4000 34449
          0: 1   .    .     .
          1: .  80  205     .
          2: .  60 3795 34449



 time C10=minimalBetti(G, DegreeLimit=>2, LengthLimit=>3)--original variety
     -- used 7476.34 seconds

            0   1    2     3
o6 = total: 1 140 5200 49629
         0: 1   .    .     .
         1: .  80  205     .
         2: .  60 4995 49629




F=I92+I93;
Il=ideal random(R^1,R^{3:-1});
Rn=R/Il
Fn=sub(F,Rn);
W=ZZ/10009[flatten entries mingens ideal vars Rn]
B=sub(Fn,W);
time C9=minimalBetti(B, DegreeLimit=>2, LengthLimit=>3)



*****************gaus maps***degree 3
load"I81gaus3.m2"
load"I82gaus3.m2"
F=I81+I82;
time C8gaus3=minimalBetti(F, DegreeLimit=>2, LengthLimit=>4)
    -- used 1888.62 seconds

            0  1   2    3     4
o4 = total: 1 60 825 5935 27708
         0: 1 10  45  120   210
         1: . 50 780 5815 27498


*****under gaus3 8,9,10 are the same with the following Betti tables
     -- used 8.87622 seconds

             0  1   2   3    4    5   6   7   8  9 10
o13 = total: 1 50 280 765 1248 1260 790 335 126 40  5
          0: 1  .   .   .    .    .   .   .   .  .  .
          1: . 50 280 765 1248 1260 720 175   .  .  .
          2: .  .   .   .    .    .  70 160 126 40  5




kk=ZZ/32003
R=kk[a,b,c,d,e]
T8=matrix{{13822*a + 4172*b - 6139*c + 1659*d - 11470*e, -15347*a + 4006*b + 14034*c + 10227*d + 7086*e,  1424*a + 15806*b - 15852*c + 7989*d - 13095*e,   10878*a + 9166*b - 10548*c - 3171*d - 6772*e , -15367*a + 7246*b - 1187*c - 6563*d - 14161*e}
,{(-1)*13970*a + 14989*b - 6495*c + 14928*d - 6034*e ,  5745*a - 2451*b - 15267*c + 14198*d - 1506*e,  6804*a + 15102*b + 13654*c + 9790*d + 14393*e ,    8380*a - 2002*b + 9755*c + 5385*d - 9750*e, 9941*a + 14180*b - 11028*c - 13139*d - 13433*e}
,{  -1865*a + 7398*b - 173*c + 12093*d + 14381*e,    2717*a - 9231*b - 12951*c - 7957*d + 1903*e,   5997*a + 7577*b - 8375*c + 11490*d - 13855*e, 15860*a - 13161*b + 12245*c + 7725*d - 14624*e,    3959*a + 2954*b + 9879*c - 5514*d + 12782*e}
,{  13130*a + 4376*b + 1254*c - 7422*d + 12973*e,   (-1)*13476*a - 9094*b + 2042*c - 440*d - 15168*e ,  -8972*a - 6940*b - 7235*c + 11632*d + 7980*e   , 6229*a - 9383*b + 9981*c + 11252*d - 4452*e   ,3583*a - 11663*b - 8857*c - 4693*d + 14214*e}
,{   4603*a - 12400*b - 5195*c + 5389*d + 8105*e ,(-1)*13763*a + 3083*b - 12510*c - 7546*d + 14424*e,   -2111*a - 14198*b + 10079*c - 203*d + 9444*e    ,  874*a + 51*b - 11700*c + 14369*d + 7448*e   ,11444*a + 14185*b - 2374*c - 5084*d - 9555*e}}

L=ideal flatten entries T8;
rel=syz(gens L, DegreeLimit=>1);
R=kk[r_0..r_24]
relR=sub(rel,R)
Eq8=ideal ((vars R)*relR);
M=matrix{{r_0..r_4},{r_5..r_9},{r_10..r_14},{r_15..r_19},{r_20..r_24}}
mp=flatten entries gens minors(4,M);
Sbs=flatten apply(#mp, i->{r_i=>mp_i})
Tns=sub(Eq8,Sbs);
time minimalBetti(Tns, DegreeLimit=>3, LengthLimit=>)



T9=matrix{{ 8123*a + 10309*b + 6834*c - 10520*d - 10885*e , 6203*a - 3539*b + 13796*c + 12221*d - 12704*e ,  -8560*a + 10941*b - 3795*c + 4860*d + 8449*e ,  -4635*a - 10059*b - 5616*c - 4969*d - 9050*e, -6566*a - 15773*b + 8992*c - 14515*d - 10224*e}
,{    9288*a + 15233*b - 2116*c - 317*d - 6022*e,  -12796*a - 3629*b + 15734*c + 1955*d - 4738*e ,   12157*a + 14695*b - 9507*c + 8119*d - 997*e, -11915*a + 7218*b - 10201*c + 3558*d + 15694*e ,  -10683*a + 471*b - 3082*c - 5977*d + 12297*e}
,{ -5143*a - 9318*b + 10850*c + 3259*d + 13839*e ,   7769*a + 13878*b + 112*c + 4386*d + 12330*e  , 11170*a + 5469*b + 11935*c + 1397*d + 6532*e ,  -8209*a - 11505*b - 8813*c - 4579*d - 8053*e, 15832*a - 10303*b + 15195*c + 14222*d - 6898*e}
,{  2901*a + 13913*b - 10760*c - 2283*d - 8519*e  ,   6100*a - 7822*b + 4601*c + 2699*d + 8144*e   ,  7720*a - 8264*b - 12805*c - 2717*d + 481*e  ,  14591*a - 2142*b + 7441*c + 8614*d - 2958*e,  -5888*a - 15148*b + 5905*c - 8814*d + 12926*e}
,{   -334*a - 10352*b + 5613*c - 334*d - 2833*e, 14750*a + 15020*b - 13373*c + 6774*d + 15616*e , -12297*a - 4593*b + 14376*c + 2845*d + 2116*e   , -3504*a - 9876*b + 6153*c - 1452*d - 4025*e ,  11902*a - 105*b + 9817*c + 15434*d + 10760*e}}




----n=6, rk10
kk=ZZ/1009
S=kk[a,b,c,d,e,f]

T9=matrix{{ -175*a - 155*b + 144*c + 471*d - 235*e - 65*f,  -427*a + 241*b + 135*c - 130*d + 336*e + 96*f ,  393*a - 344*b + 476*c + 224*d + 90*e + 231*f ,     7*a + 169*b - 64*c - 464*d - 82*e - 312*f ,   -303*a - 17*b - 359*c + 15*d + 44*e - 220*f ,   -309*a - 212*b + 329*c - 90*d + 37*e - 37*f},
{-265*a + 144*b - 344*c + 147*d + 375*e - 141*f , -48*a + 333*b - 400*c - 423*d + 389*e + 235*f ,-423*a + 351*b - 444*c - 319*d - 424*e + 115*f,  -162*a + 89*b + 135*c + 114*d - 472*e - 357*f  ,  315*a + 136*b - 233*c - 41*d + 390*e - 72*f , -25*a - 153*b + 273*c + 484*d - 416*e + 221*f},
{  35*a - 279*b + 109*c + 132*d + 433*e + 420*f , -465*a - 396*b - 338*c - 483*d + 292*e - 31*f,  151*a + 494*b + 500*c + 352*d - 295*e + 137*f ,   66*a + 446*b + 152*c + 159*d + 11*e + 370*f , -105*a + 58*b + 232*c + 176*d + 458*e + 279*f   ,-180*a + 93*b - 362*c - 330*d - 119*e + 15*f},
{ 10*a + 414*b - 405*c - 170*d - 51*e + 302*f ,   147*a + 74*b - 178*c - 141*d + 44*e + 198*f,     43*a + 6*b + 390*c - 175*d + 344*e + 149*f,   24*a - 326*b - 149*c - 369*d - 352*e + 297*f ,   336*a - 92*b - 460*c + 376*d - 144*e + 26*f    ,267*a - 62*b - 269*c + 145*d - 25*e + 118*f},
{-278*a + 114*b + 345*c - 240*d + 189*e + 446*f,     25*a + 61*b + 379*c - 488*d + 50*e + 308*f  ,-231*a + 385*b - 481*c + 202*d - 58*e + 246*f, -491*a + 346*b + 200*c + 378*d + 299*e + 245*f  , 271*a - 222*b + 442*c + 79*d + 277*e - 367*f , -366*a - 313*b - 323*c + 77*d - 449*e - 438*f},
{    71*a + 200*b + 14*c - 66*d - 382*e + 62*f  ,-222*a - 101*b + 89*c - 257*d - 371*e + 194*f ,    435*a - 397*b + 16*c + 435*d - 61*e + 37*f,   -252*a - 179*b - 157*c - 244*d + 45*e + 14*f   , 93*a + 61*b - 341*c + 423*d + 150*e + 213*f  , 490*a - 276*b - 70*c + 149*d + 453*e - 202*f}}

L=ideal flatten entries T9;
rel=syz(gens L, DegreeLimit=>1);
R=kk[r_0..r_35]
relR=sub(rel,R);
Eq10=ideal ((vars R)*relR);
M=matrix{{r_0..r_5},{r_6..r_11},{r_12..r_17},{r_18..r_23},{r_24..r_29},{r_30..r_35}}
N2=minors(2,M);
time Int=intersect(N2,Eq10);
betti Int
time K1=syz(gens Int, DegreeLimit=>3);
time K2=syz(K1, DegreeLimit=>4)

**************Koszul module construction for tensors in C3*C3 with V=C3+C3
kk=ZZ/101
Z=kk[x_0..x_2,y_0..y_2, Degrees=>{3:{1,0},3:{0,1}}] 
m=ideal vars Z;
F=res(m)--kosul complex

S=kk[a,b,c]
vects=apply(3, i->random(S^1,S^{3: 0}))
A=(transpose vects_1)*vects_2;
A1=sub(matrix{{1,0,0},{0,0,0},{0,0,0}},S);
A2=sub(matrix{{0,0,0},{0,1,0},{0,0,0}},S);
A3=sub(matrix{{0,0,0},{0,0,0},{0,0,1}},S);
u=vects_0;
--pick a general rank 4 tensor: e1*e1*e1+e2*e2*e2+e3*e3*e3+v*w*z
T4=a*A1+b*A2+c*A3+u_(0,0)*a*A+u_(0,1)*b*A+u_(0,2)*c*A--tensor
L=ideal flatten entries T4
rel=syz(gens L, DegreeLimit=>1);
R=kk[r_0..r_8]
relR=sub(rel,R)
Eq=saturate ideal ((vars R)*relR)--embedding of the tensor in P^8	  
Rr=R/Eq	  
mingens ideal vars Rr-- giving the coordinates of P^2 as tensor
EmdOfTensor=submatrix(F.dd_2,{7,11,12})	--pick the submatrix of syzygies corresponding to P^2 (tensor) 
im=image EmdOfTensor	  
Im=image(F.dd_2)	  
isSubset(im,Im)
K=ker(F.dd_1)
isSubset(im,K)
KozMod4=K/im
betti res KozMod4



              0  1  2 3 4
o98 = total: 15 23 15 6 1
          1:  .  3  . . .
          2: 15 20 15 6 1


---koszul module of rk5 tensor
vect1=apply(3, i->random(S^1,S^{3: 0}))
vect2=apply(3, i->random(S^1,S^{3: 0}))
A=(transpose vect1_1)*vect1_2;
B=(transpose vect2_1)*vect2_2;
A1=sub(matrix{{1,0,0},{0,0,0},{0,0,0}},S);
A2=sub(matrix{{0,0,0},{0,1,0},{0,0,0}},S);
A3=sub(matrix{{0,0,0},{0,0,0},{0,0,1}},S);
u=vect1_0;
v=vect2_0;
T5=a*A1+b*A2+c*A3+u_(0,0)*a*A+u_(0,1)*b*A+u_(0,2)*c*A+v_(0,0)*a*B+v_(0,1)*b*B+v_(0,2)*c*B
L=ideal flatten entries T5;
rel=syz(gens L, DegreeLimit=>1);
R=kk[r_0..r_8]
use R
relR=sub(rel,R)
Eq=ideal ((vars R)*relR);--embedding of the tensor in P^8
Rr=R/Eq
mingens ideal vars Rr-- giving the coordinates of P^2 as tensor
EmdOfTensor=submatrix(F.dd_2,{10,11,12})	  
im=image EmdOfTensor	  
Im=image(F.dd_2)	  
isSubset(im,Im)
K=ker(F.dd_1)
isSubset(im,K)
KozMod5=K/im
betti res KozMod5


              0  1  2 3 4
o98 = total: 15 23 15 6 1
          1:  .  3  . . .
          2: 15 20 15 6 1

----kosul module for rk 7,6--are same
restart
kk=ZZ/101
Z=kk[x_0..x_3,y_0..y_3, Degrees=>{4:{1,0},4:{0,1}}] 
m=ideal vars Z;
F=res(m)
 
 
S=kk[a,b,c,d]
vect1=apply(3, i->random(S^1,S^{4: 0}))
vect2=apply(3, i->random(S^1,S^{4: 0}))
A=(transpose vect1_1)*vect1_2;
B=(transpose vect2_1)*vect2_2;
A1=sub(matrix{{1,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},S);
A2=sub(matrix{{0,0,0,0},{0,1,0,0},{0,0,0,0},{0,0,0,0}},S);
A3=sub(matrix{{0,0,0,0},{0,0,0,0},{0,0,1,0},{0,0,0,0}},S);
A4=sub(matrix{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,1}},S)
u=vect1_0;
v=vect2_0;
---pick rk 6 tensor
T6=a*A1+b*A2+c*A3+d*A4+u_(0,0)*a*A+u_(0,1)*b*A+u_(0,2)*c*A+
u_(0,3)*d*A+v_(0,0)*a*B+v_(0,1)*b*B+v_(0,2)*c*B+u_(0,3)*d*B
L=ideal flatten entries T6;
rel=syz(gens L, DegreeLimit=>1);
R=kk[r_0..r_15]
relR=sub(rel,R)
Eq=ideal((vars R)*relR);--embedding of the tensor in P^15 
Rr=R/Eq
mingens ideal vars Rr-- giving the coordinates of P^2 as tensor
EmdOfTensor=submatrix(F.dd_2,{18,22,23,24})	  
im=image EmdOfTensor	  
Im=image(F.dd_2)	  
isSubset(im,Im)
K=ker(F.dd_1)
isSubset(im,K)
KozMod=K/im
betti res KozMod


              0  1  2  3  4 5 6
o35 = total: 28 60 70 56 28 8 1
          1:  .  4  .  .  . . .
          2: 28 56 70 56 28 8 1


********************************


