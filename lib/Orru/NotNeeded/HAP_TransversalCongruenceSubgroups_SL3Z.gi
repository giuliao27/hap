InstallGlobalFunction(HAP_TransversalCongruenceSubgroups_SL3Z,
function ( G, H )
    local tree, InH, S, T, U, v, p, g, s, n, q, vv, gens, nodes, nodesinv, leaves, 
    ambientGenerators, InLowDegreeNodesModH, one, poscan, nind;
    S := [ [1,0,1], [0,-1,-1], [0,1,0] ];
    T := [ [0,1,0], [0,0,1], [1,0,0]];
    U := [[0,1,0], [1,0,0], [-1,-1,-1]];
    one := IdentityMat( 3 );
    ambientGenerators := [S, T, U];
    tree := [ 1 ];
    cnt := 1;
    leaves := NewDictionary( S, true, SL( 3, Integers ) );
    nodes := [ one ];
    AddDictionary( leaves, one, 1 );
    InH := H!.membershipLight;
    InLowDegreeNodesModH := function ( g )
          local x, gg, B1, B2;
          gg := g ^ -1;
          for x in nodes do
              if InH( x * gg ) then
                  return x;
              fi;
          od;
          return false;
      end;
    while Size( leaves ) > 0 do
        vv := leaves!.entries[1];
        v := vv[1];
        for s in [ 1 .. Length( ambientGenerators ) ] do
            g := v * ambientGenerators[s];
            q := InLowDegreeNodesModH( g );
            if q = false then
                Add( nodes, g );
                AddDictionary( leaves, g, Length( nodes ) );
                p := LookupDictionary( leaves, v );
                Add( tree, [ p, s ] );
            fi;
        od;
        RemoveDictionary( leaves, v );
    od;
    #nodes := Filtered( nodes, function ( g )
    #        return g in G;
    #    end );
    nodesinv := List( nodes, function ( g )
            return g ^ -1;
        end );
    nind := [ 1 .. Length( nodes ) ];
    poscan := function ( x )
          local i;
          for i in nind do
              if InH( x * nodesinv[i] ) then
                  return i;
              fi;
          od;
          return fail;
      end;
    return 
     Objectify( 
       NewType( FamilyObj( G ), 
         IsHapRightTransversalSL3ZSubgroup and IsList and IsDuplicateFreeList 
          and IsAttributeStoringRep ), rec(
          group := G,
          subgroup := H,
          cosets := nodes,
          poscan := poscan ) );
end);