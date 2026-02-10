#descrive come i generatori di un gruppo agiscono sui coset di un sottogruppo.

InstallGlobalFunction(HAP_SL3ZSubgroupTree_fast,
function(G)
    local ambientGenerators, tree, InGmodU, Ugrp, v, p, g, s, n, q, vv, leaves, genTriples, generators, InLowDegreeNodesModG, csts, cnt, vertex2word, one, triple2word, i, j, u,c, a, b;
    
    ambientGenerators:=G!.ambientGenerators;
    one := IdentityMat(3); #to change? stabilizers to consider
    Ugrp := G!.ugrp;
    Ugrp := Elements( Ugrp );
    tree := [  ];
    genTriples := [  ];
    cnt := 1;
    leaves := NewDictionary( ambientGenerators[1], true, SL( 3, Integers ) );
    csts := [  ];
    csts[G!.cosetPos( one )] := 1;
    AddDictionary( leaves, one, G!.cosetPos( one ) );
    InLowDegreeNodesModG := function ( g )
          local pos;
          pos := G!.cosetPos( g );
          if not IsBound( csts[pos] ) then
              return false;
          else
              ;
              return pos;
          fi;
          return;
      end;
    while Size( leaves ) > 0 do
        vv := 1 * leaves!.entries[1];
        v := vv[1];
        p := G!.cosetPos( v );
        for s in [ 1 .. Length( ambientGenerators ) ] do
            g := ambientGenerators[s] * v;
            q := InLowDegreeNodesModG( g );
            if q = false then
                q := G!.cosetPos( g );
                AddDictionary( leaves, g, q );
                csts[q] := 1;
                tree[q] := [ p, s ];
            else
                Add( genTriples, [ p, s, v, g ] );
            fi;
        od;
        RemoveDictionary( leaves, v );
    od;
    vertex2word := function ( v )
          local word;
          word := one;
          while IsBound( tree[v] ) do
              word := word * ambientGenerators[tree[v][2]];
              v := tree[v][1];
          od;
          return word;
      end;
    triple2word := function ( x )
          local u, uu, g, q, c;
          for u in Ugrp do
              c := x[4] ^ -1 * u * vertex2word( G!.cosetPos( x[4] ) );
              if G!.membership( c ) then
                  return c;
              fi;
          od;
          return fail;
      end;
    genTriples := List( genTriples, function ( x )
            return triple2word( x );
        end );
    genTriples := List( genTriples, function ( x )
            return Minimum( x, x ^ -1 );
        end );
    genTriples := SSortedList( genTriples );
    G!.tree := tree;
    G!.GeneratorsOfMagmaWithInverses := genTriples;
    return tree;
end);
