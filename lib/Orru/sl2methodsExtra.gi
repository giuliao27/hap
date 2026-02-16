InstallMethod(CosetPosFunction,
    "Returns cosetPos(g) function for the congruence subgroup G",
    [ IsIntegerMatrixGroup and IsHAPCongruenceSubgroupGamma0 ],
    function(G)
        local cosetPos, canonicalRep, n, ProjLine;

        if DimensionOfMatrixGroup(G) > 2 then
            TryNextMethod();
        fi;

        n := LevelOfCongruenceSubgroup(G);

        ProjLine := FiniteProjectiveLine(n);

        canonicalRep := function(g)
            local v, vv, U, d, dd;
            v := [g[1][1], g[2][1]];
            vv := List(v, x -> x mod n);
            U := Units(Integers mod n);
            if vv[1] mod n = 0 then
                return [0,1];
            elif ZmodnZObj(vv[1],n) in U then
                return [1,(Inverse(vv[1]) mod n)*vv[2] mod n];
            else
                d := Gcd(vv[1],n);
                dd := n/d;
                return [vv[1], vv[2] mod dd];
            fi;
        end;

        cosetPos := function(g)
            local w;
            w := canonicalRep(g);
            return Position(ProjLine,w);
        end;

        return cosetPos;
    end);