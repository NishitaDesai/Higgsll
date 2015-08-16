Dimension 4;
vector p8,p9,p10,p11;
index  mu1,nu1,mu2,nu2,i,j,k,l;
symbols s,t,u,a,b,c,mw;
index si1,la1,mu2,nu2,si2,la2;

CFunctions ep(antisymmetric);

* The diagram:
*                           W^+(k2)  __ l^+(p10)
* l^-(p8) ----|  W-*(k1) |----------<__
*             |----------|              nu(p11)
* Nu (p9) ----|          |--------- h(p5)
*
* Trace follows Peskin-Schroeder notation: 
*            opposite to charge flow

* --------- Definitions -------------

Local k1(i)=p8(i)+p9(i);
Local k2(i)=p10(i)+p11(i);

* ----- Trace for quark fields ---
Local T1A(mu1,mu2)=g_(1,p9)*g_(1,mu1)*g_(1,p8)*g_(1,mu2);
Local T1B(mu1,mu2)=g_(1,p9)*g_(1,mu1)*g_(1,p8)*g_(1,mu2)*g5_(1);

* ------ Trace for Lepton Fields
Local T2A(nu1,nu2)=g_(2,p10)*g_(2,nu1)*g_(2,p11)*g_(2,nu2);
Local T2B(nu1,nu2)=g_(2,p10)*g_(2,nu1)*g_(2,p11)*g_(2,nu2)*g5_(2);

Local T1(mu1,mu2)=T1A(mu1,mu2)-4*i_*ep(p9,mu1,p8,mu2);
Local T2(nu1,nu2)=T2A(nu1,nu2)-4*i_*ep(p10,nu1,p11,nu2);

* ------- eps(mu,nu,k1,k2)  ------------
Local eps1(mu1,nu1) = ep(mu1,nu1,p8,p10) + ep(mu1,nu1,p8,p11)
		    + ep(mu1,nu1,p9,p10) + ep(mu1,nu1,p9,p11);
Local eps2(mu2,nu2) = ep(mu2,nu2,p8,p10) + ep(mu2,nu2,p8,p11)
		    + ep(mu2,nu2,p9,p10) + ep(mu2,nu2,p9,p11);


* ----- Total hWW vertex ------

Local vertex(mu1,nu1)=d_(mu1,nu1)+(b+i_*c)*eps1(mu1,nu1);

Local vertexbar(mu2,nu2)=d_(mu2,nu2)+(b-i_*c)*eps2(mu2,nu2);
Local v(mu1,mu2,nu1,nu2)=vertex(mu1,nu1)*vertexbar(mu2,nu2);

*  ----- Standard model only ----
* Local v(mu1,mu2,nu1,nu2) = d_(mu1,nu1)*d_(mu2,nu2);


* ---- Full |M|^2 (leaving out propagator denominators) ----
Local mesq=T1(mu1,mu2)*T2(nu1,nu2)*v(mu1,mu2,nu1,nu2);


* ---- Commands ----

trace4,1;
trace4,2;
contract;
contract;
.sort

* -------- Substitutions -------
id ep(mu1?,nu1?,la1?,si1?)*ep(mu2?,nu2?,la2?,si2?)=
    -d_(mu1,mu2)*(d_(nu1,nu2)*d_(la1,la2)*d_(si1,si2)
                + d_(si1,nu2)*d_(nu1,la2)*d_(la1,si2)
                + d_(la1,nu2)*d_(si1,la2)*d_(nu1,si2)
                - d_(si1,nu2)*d_(la1,la2)*d_(nu1,si2)
                - d_(la1,nu2)*d_(nu1,la2)*d_(si1,si2)
                - d_(nu1,nu2)*d_(si1,la2)*d_(la1,si2))
    +d_(nu1,mu2)*(d_(mu1,nu2)*d_(la1,la2)*d_(si1,si2)
                + d_(si1,nu2)*d_(mu1,la2)*d_(la1,si2)
                + d_(la1,nu2)*d_(si1,la2)*d_(mu1,si2)
                - d_(si1,nu2)*d_(la1,la2)*d_(mu1,si2)
                - d_(la1,nu2)*d_(mu1,la2)*d_(si1,si2)
                - d_(mu1,nu2)*d_(si1,la2)*d_(la1,si2))
    -d_(la1,mu2)*(d_(mu1,nu2)*d_(nu1,la2)*d_(si1,si2)
                + d_(si1,nu2)*d_(mu1,la2)*d_(nu1,si2)
                + d_(nu1,nu2)*d_(si1,la2)*d_(mu1,si2)
                - d_(si1,nu2)*d_(nu1,la2)*d_(mu1,si2)
                - d_(nu1,nu2)*d_(mu1,la2)*d_(si1,si2)
                - d_(mu1,nu2)*d_(si1,la2)*d_(nu1,si2))
    +d_(si1,mu2)*(d_(mu1,nu2)*d_(nu1,la2)*d_(la1,si2)
                + d_(la1,nu2)*d_(mu1,la2)*d_(nu1,si2)
                + d_(nu1,nu2)*d_(la1,la2)*d_(mu1,si2)
                - d_(la1,nu2)*d_(nu1,la2)*d_(mu1,si2)
                - d_(nu1,nu2)*d_(mu1,la2)*d_(la1,si2)
                - d_(mu1,nu2)*d_(la1,la2)*d_(nu1,si2));

id ep(mu1?,nu1?,la1?,si1?)*e_(mu2?,nu2?,la2?,si2?)=
    -d_(mu1,mu2)*(d_(nu1,nu2)*d_(la1,la2)*d_(si1,si2)
                + d_(si1,nu2)*d_(nu1,la2)*d_(la1,si2)
                + d_(la1,nu2)*d_(si1,la2)*d_(nu1,si2)
                - d_(si1,nu2)*d_(la1,la2)*d_(nu1,si2)
                - d_(la1,nu2)*d_(nu1,la2)*d_(si1,si2)
                - d_(nu1,nu2)*d_(si1,la2)*d_(la1,si2))
    +d_(nu1,mu2)*(d_(mu1,nu2)*d_(la1,la2)*d_(si1,si2)
                + d_(si1,nu2)*d_(mu1,la2)*d_(la1,si2)
                + d_(la1,nu2)*d_(si1,la2)*d_(mu1,si2)
                - d_(si1,nu2)*d_(la1,la2)*d_(mu1,si2)
                - d_(la1,nu2)*d_(mu1,la2)*d_(si1,si2)
                - d_(mu1,nu2)*d_(si1,la2)*d_(la1,si2))
    -d_(la1,mu2)*(d_(mu1,nu2)*d_(nu1,la2)*d_(si1,si2)
                + d_(si1,nu2)*d_(mu1,la2)*d_(nu1,si2)
                + d_(nu1,nu2)*d_(si1,la2)*d_(mu1,si2)
                - d_(si1,nu2)*d_(nu1,la2)*d_(mu1,si2)
                - d_(nu1,nu2)*d_(mu1,la2)*d_(si1,si2)
                - d_(mu1,nu2)*d_(si1,la2)*d_(nu1,si2))
    +d_(si1,mu2)*(d_(mu1,nu2)*d_(nu1,la2)*d_(la1,si2)
                + d_(la1,nu2)*d_(mu1,la2)*d_(nu1,si2)
                + d_(nu1,nu2)*d_(la1,la2)*d_(mu1,si2)
                - d_(la1,nu2)*d_(nu1,la2)*d_(mu1,si2)
                - d_(nu1,nu2)*d_(mu1,la2)*d_(la1,si2)
                - d_(mu1,nu2)*d_(la1,la2)*d_(nu1,si2));


 id p8.p8=0;
 id p9.p9=0;
 id p10.p10=0;
 id p11.p11=0;
* id b=0;
* id c=0;
* id b^2=0;
* id c^2=0;

* ---------- Printing --------
 Format Fortran ;
 print +s;
 .end
