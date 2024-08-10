SquareFactorization := function(n)
    // input: n = p^2*q, where p and q are primes of roughly the same size
    // output: p^2

    q2 := n^(1/3);  // we assume that q is roughly equal to p
    e := Sqrt(Log(q2)/Log(Log(q2)));
    B := 2*Round(q2^(1/(2*e)));  // we took B somewhat bigger,
    // not sure if it is better

    prod := 1;  // prod is the integer k
    for p in PrimesUpTo(B) do
        prod *:= p^(Floor(Log(p,B)));  // smaller exponents is better in practice
    end for;

    B2 := 2*B*Round(Log(B));
    stage2primes := PrimesInInterval(B,B2);

    largeprime := NextPrime(Round(10*n^(1/6)));
    // this is the prime r from Algorithm 2.
    // if p/q <= 10, then this r is large enough

    s := 0;

    while true do
        s +:= 1;
        if IsSquarefree(s) eq false then
            continue;
        end if;

        // STAGE 1
        m := n*s;
        Q := QuadraticForms(-4*m);
        Q2 := QuadraticForms(-4*m*largeprime^2);

        for r in PrimesInInterval(3,1000) do
        // we search for an element in C(-4*m)
        // we expect this to work within a few tries
            if LegendreSymbol(-m,r) ne 1 then
                continue;
            end if;
            b := Integers() ! Sqrt(Integers(4*r) ! (-4*m));
            f := Q ! [r,b,(4*m+b^2) div (4*r)];
            break;
        end for;

        g := f^prod;  // this is the main part of stage 1
        // we hope that g is the identity element in 
        // the underlying group C(-4*q*s)

        h := Q2 ! [g[1]*largeprime^2, g[2]*largeprime, g[3]];
        // we lift g to the larger group C(-4mr^2)
        // this makes it possible to retrieve p^2
        h := Reduction(h);
        k := h^(largeprime - LegendreSymbol(-m,largeprime));

        d := GCD(k[1],n);
        if d gt 1 and d lt n then  // if a factor is found
            return(d);  // then return it
        end if;
        
        // STAGE 2
        k2 := k^2;
        smallpowers := [k2];
        k3 := k2;
        smallbound := Round(10*Log(B2));  // largest prime gap
        // for our tables this bound was large enough
        for i in [1..smallbound] do
            k2 *:= k3;
            Append(~smallpowers, k2);
        end for;

        previousprime := 0;
        // for each prime z in [B,B_2] we will consider k^z
        for p1 in stage2primes do
            currentprime := p1;
            if previousprime eq 0 then
                l := k^currentprime;
            else
                smallstep := (currentprime - previousprime) div 2;
                if smallstep gt smallbound then
                    "help", smallbound, smallstep;
                    break;
                else
                    // we try the next prime
                    l *:= smallpowers[smallstep];
                end if;
           end if;

            d := GCD(l[1],n);
            if d gt 1 and d lt n then
                return(d);
            end if;
            previousprime := currentprime;
        end for;

        // if stage 2 was also unsuccessful,
        // then take the next s and try again

    end while;

end function;


// testing the function
p := NextPrime(10^20);
q := NextPrime(10^21);
n := p^2*q;
SquareFactorization(n);
