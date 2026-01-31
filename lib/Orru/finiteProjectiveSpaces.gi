InstallGlobalFunction(FiniteProjectiveLine,
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

InstallGlobalFunction(FiniteProjectivePlane,
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