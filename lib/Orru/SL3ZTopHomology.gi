
###########################################################
###########################################################
InstallGlobalFunction(Gamma0_SL3ZTopRationalHomology,
function(n)
    local FirstCondition, SecondCondition, ThirdCondition_1, ThirdCondition_2, A, B, C, i, j, M, SolutionSpace,FiniteProjectivePlane,P;

#Returns the dimension of the degree 3 rational homology of the 
#congruence subgroup Gamma_0(n) in SL(3,Z). 

##########################################
##########################################
FiniteProjectivePlane := function(n)
    local UnitEls, x, y, z, i, c, d, u, UnitsAction, Representatives, 
          RepOf, r, v, w, m, min;


    UnitEls := Units(Integers mod n);

    UnitsAction := function(c, u)
    local uu;
        uu:=Int(u);
        return List(c, x -> (uu * x) mod n);
    end;

    Representatives:=[];

    RepOf:=[];
    for x in [1..n] do
      RepOf[x]:=[];
      for y in [1..n] do
        RepOf[x][y]:=[];
      od;
    od;

    for x in [0..n-1] do
      for y in [0..n-1] do
        for z in [0..n-1] do
          if Gcd(x,y,z,n) = 1 then
            v := [x,y,z];
            if not IsBound(RepOf[x+1][y+1][z+1]) then
              m:=Orbit(UnitEls,v,UnitsAction);
              min:=Minimum(m);
              AddSet(Representatives,min);
              for w in m+1 do
                RepOf[w[1]][w[2]][w[3]]:=min;
              od;
            fi;
          fi;
        od;
      od;
    od;


    return rec(
        Reps := Set(Representatives),
        RepOf:= RepOf
    );
end;
##########################################
##########################################

P:=FiniteProjectivePlane(n);
    FirstCondition := function(p)
        local v;
        v := [-p[2],p[1],p[3]] mod n;
        v:=v+1;
        return P.RepOf[v[1]][v[2]][v[3]];
    end;

    SecondCondition := function(p)
        local v;
        v := [p[3],p[1],p[2]] mod n;
        v:=v+1;
        return P.RepOf[v[1]][v[2]][v[3]];
    end;

    ThirdCondition_1 := function(p)
        local v;
        v := [-p[2],p[1]-p[2],p[3]] mod n;
        v:=v+1;
        return P.RepOf[v[1]][v[2]][v[3]];
    end;

    ThirdCondition_2 := function(p)
        local v;
        v := [p[2]-p[1],p[1],p[3]] mod n;
        v:=v+1;
        return P.RepOf[v[1]][v[2]][v[3]];
    end;

    A := SparseIdentityMat(Length(P.Reps));
    B := SparseIdentityMat(Length(P.Reps));
    C := SparseIdentityMat(Length(P.Reps));

    for i in [1..Length(P.Reps)] do
        if FirstCondition(P.Reps[i]) in P.Reps then
            j := Position(P.Reps,FirstCondition(P.Reps[i]));
            SparseMatAddToEntry(A,i,j,1);
        fi;
        if SecondCondition(P.Reps[i]) in P.Reps then
            j := Position(P.Reps,SecondCondition(P.Reps[i]));
            SparseMatAddToEntry(B,i,j,-1);
        fi;
        if ThirdCondition_1(P.Reps[i]) in P.Reps then
            j := Position(P.Reps,ThirdCondition_1(P.Reps[i]));
            SparseMatAddToEntry(C,i,j,1);
        fi;
        if ThirdCondition_2(P.Reps[i]) in P.Reps then
            j := Position(P.Reps,ThirdCondition_2(P.Reps[i]));
            SparseMatAddToEntry(C,i,j,1);
        fi;
    od;

    M:=SparseMatConcatenation(A,B);Unbind(A);Unbind(B);
    M:=SparseMatConcatenation(M,C);Unbind(C);
    return Length(P.Reps) - RankMatDestructive(M);
end);
###########################################################
###########################################################

