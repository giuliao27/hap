
########################################################
########################################################
InstallMethod(FiniteProjectiveLine,
"Finite projective line for the ring Z/nZ",
[IsInt],
function(n)
return HAP_FiniteProjectiveLineIntegers(n);
end);
########################################################
########################################################

########################################################
########################################################
InstallMethod(FiniteProjectivePlane,
"Finite projective plane for the ring Z/nZ",
[IsInt],
function(n)
return HAP_FiniteProjectivePlaneIntegers(n);
end);
########################################################
########################################################
InstallMethod(FiniteProjectiveLine_alt,
"Finite projective line for the ring Z/nZ",
[IsInt],
function(n)
return HAP_FiniteProjectiveLineIntegers_alt(n);
end);




InstallGlobalFunction(HAP_FiniteProjectiveLineIntegers_alt,
function(n)
    local UnitEls, x, y, i, c, d, u, UnitsAction, Representatives, 
          RepOf, r, v, w, m, min;

    UnitEls := Units(Integers mod n);

    UnitsAction := function(c, u)
      local uu;
      uu :=Int (u);
      return List(c, x -> (uu * x) mod n);
    end;

    Representatives := [];

    RepOf := [];
    for x in [1..n] do
      RepOf[x] := [];
    od;

    for x in [0..n-1] do
        for y in [0..n-1] do
          if Gcd(x,y,n) = 1 then
            v := [x,y];
            if not IsBound(RepOf[x+1][y+1]) then
              m := Orbit(UnitEls,v,UnitsAction);
              min := Minimum(m);
              AddSet(Representatives,min);
              for w in m + 1 do
                RepOf[w[1]][w[2]] := min;
              od;
            fi;
          fi;
        od;
    od;

    return rec(
        Reps := Set(Representatives),
        RepOf:= RepOf
    );
end);

InstallGlobalFunction(HAP_FiniteProjectivePlaneIntegers,
function(n)
    local UnitEls, x, y, z, i, c, d, u, UnitsAction, Representatives, 
          RepOf, r, v, w, m, min;

    UnitEls := Units(Integers mod n);

    UnitsAction := function(c, u)
      local uu;
      uu :=Int (u);
      return List(c, x -> (uu * x) mod n);
    end;

    Representatives := [];

    RepOf := [];
    for x in [1..n] do
      RepOf[x] := [];
      for y in [1..n] do
        RepOf[x][y] := [];
      od;
    od;

    for x in [0..n-1] do
      for y in [0..n-1] do
        for z in [0..n-1] do
          if Gcd(x,y,z,n) = 1 then
            v := [x,y,z];
            if not IsBound(RepOf[x+1][y+1][z+1]) then
              m := Orbit(UnitEls,v,UnitsAction);
              min := Minimum(m);
              AddSet(Representatives,min);
              for w in m + 1 do
                RepOf[w[1]][w[2]][w[3]] := min;
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
end);

InstallGlobalFunction(HAP_FiniteProjectiveLineIntegers,
function(n)
    local Rep, i, factors, Znd, d, z, p, l, toFill, t;

    Rep := [[0,1]];

    for i in [1..n] do
        Add(Rep,[1,i-1]);
    od;
    
    factors := List(DivisorsInt(n));
    Remove(factors,1);
    Remove(factors);

    for p in factors do
        d := n/p;

        Znd := [0..d-1];
        toFill := [];

        for z in Znd do
            if Gcd(p,z,d) = 1 then
                Add(toFill, z);
            fi;
        od;

        t := 1;
        while not IsEmpty(toFill) do
            if Gcd(p,t) = 1 then
                if (t mod d) in toFill then
                    Add(Rep,[p,t]);
                    Remove(toFill, Position(toFill,t mod d));
                fi;
            fi;
            t := t+1;
        od;
    od;

    return Rep;
end);


