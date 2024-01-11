--Gorenstein
---A5/P4---
kk=ZZ/1009
R=kk[x,y,z,w]
m=matrix{{w^2+x*y+z^2}}
I=inverseSystem(m, DividedPowers=>true);
degree I
apply(5,i->hilbertFunction(i,I))
degree Hom(I,R/I)
---A5/P^4
kk=ZZ/1009
R=kk[x,y,z,w,u]
m=matrix{{w^2+x*y+z^4+u^2}}
I=inverseSystem(m, DividedPowers=>true)
degree I
apply(5,i->hilbertFunction(i,I))

degree Hom(I,R/I)--dim of tangent space
----

loadPackage"VersalDeformations"
F0=mingens I;
F1=gens ker F0;
T1=cotangentCohomology1(F0);
(F,R)=firstOrderDeformations(F0,F1,T1);
Idef=sum F;
