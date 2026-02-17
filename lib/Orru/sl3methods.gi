##########################################################################
##
## Methods for 3x3 congruence subgroups of SL3

##########################################################################
##
## CosetPosFunction( <G> )
##
## Returns a function cosetPos(g) giving the position of the coset gG in 
## the ambient group. 
InstallMethod(CosetPosFunction,
"Returns cosetPos(g) function for the congruence subgroup G",
[ IsIntegerMatrixGroup and IsHAPCongruenceSubgroupGamma0 ],
    function(G)
        local cosetPos, n, ProjPlane;
        if DimensionOfMatrixGroup(G) <> 3 then
            TryNextMethod();
        fi;

        n := LevelOfCongruenceSubgroup(G);
    
        ProjPlane := FiniteProjectivePlane(n);

        cosetPos := function(g)
            local v, vv, U, u, w;
            v := [g[1][1], g[2][1], g[3][1]];
            vv := List(v, x -> x mod n);
            U := Units(Integers mod n);
            for u in U do
                w := List(vv, x -> (Int(u)*x) mod n);
                if w in ProjPlane.Reps then
                    return Position(ProjPlane.Reps,w);
                fi;
            od;
        end;

        return cosetPos;
    end);
##########################################################################
##
## CosetRepFunction( <G> )
##
## Returns a function cosetPos(g) giving a canonical rpresentative of the 
## coset gG in the ambient group. 
InstallMethod(CosetRepFunction,
"Returns cosetRep(g) function for the congruence subgroup G",
[ IsIntegerMatrixGroup and IsHAPCongruenceSubgroupGamma0 ],
     function(G)
        local MatrixInSL3_Hermite, cosetOfInt, cosetRep, n, ProjPlane, cosetPos;
        if DimensionOfMatrixGroup(G) <> 3 then
            TryNextMethod();
        fi;

        n := LevelOfCongruenceSubgroup(G);

        MatrixInSL3_Hermite := function(v)
            local Herm;
            Herm := HermiteNormalFormIntegerMatTransform([[v[1]],[v[2]],[v[3]]]);
            return Inverse(Herm!.rowtrans);
        end;

        ProjPlane := FiniteProjectivePlane(n);

        cosetOfInt:=function(i)
            local x,y,z;
            x := ProjPlane.Reps[i][1];
            y := ProjPlane.Reps[i][2];
            z := ProjPlane.Reps[i][3];

            return MatrixInSL3_Hermite([x,y,z]);
        end;

        cosetPos := CosetPosFunction(G);

        cosetRep:=function(g);
            return cosetOfInt(cosetPos(g));
        end;

        return cosetRep;
     end);