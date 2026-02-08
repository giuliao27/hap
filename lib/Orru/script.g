# SchreierGens := [];
# n:=3;
# ProjPlane := FiniteProjectivePlane(n);
# G:=HAPCongruenceSubgroupGamma0(3,n);
# S := [
#   [[1,1,0],[0,1,0],[0,0,1]],
#   [[1,0,0],[1,1,0],[0,0,1]],
#   [[0,-1,0],[1,0,0],[0,0,1]],
#   [[1,0,0],[0,0,-1],[0,1,0]]
# ];
# for i in [1..Length(ProjPlane.Reps)] do
#   for s in S do
#     r := G!.cosetRep(i);
#     j := G!.cosetPos(r*s);
#     g := r * s * Inverse(G!.cosetRep(j));
#     if G!.membershipLight(g) then
#       Add(SchreierGens, g);
#     fi;
#   od;
# od;
# 
# SchreierGens := Set(SchreierGens);
# Print(SchreierGens);

#function ( G )
#    local ambientGenerators, tree, InGmodU, Ugrp, S, T, U, v, p, g, s, n, q, #vv, leaves, 
#    genTriples, generators, InLowDegreeNodesModG, csts, cnt, vertex2word, one, #triple2word, i, 
#    j, u, c, a, b;
#    S := [ [ 0, -1 ], [ 1, 0 ] ];
#    T := [ [ 1, 1 ], [ 0, 1 ] ];
#    U := S * T;
#    one := IdentityMat( 2 );
#    ambientGenerators := [ S * U, S * U ^ 2 ];
#    Ugrp := G!.ugrp;
#    Ugrp := Elements( Ugrp );
#    tree := [  ];
#    genTriples := [  ];
#    cnt := 1;
#    leaves := NewDictionary( S, true, SL( 2, Integers ) );
#    csts := [  ];
#    csts[G!.cosetPos( one )] := 1;
#    AddDictionary( leaves, one, G!.cosetPos( one ) );
#    InLowDegreeNodesModG := function ( g )
#          local pos;
#          pos := G!.cosetPos( g );
#          if not IsBound( csts[pos] ) then
#              return false;
#          else
#              ;
#              return pos;
#          fi;
#          return;
#      end;
#    while Size( leaves ) > 0 do
#        vv := 1 * leaves!.entries[1];
#        v := vv[1];
#        p := G!.cosetPos( v );
#        for s in [ 1 .. Length( ambientGenerators ) ] do
#            g := ambientGenerators[s] * v;
#            q := InLowDegreeNodesModG( g );
#            if q = false then
#                q := G!.cosetPos( g );
#                AddDictionary( leaves, g, q );
#                csts[q] := 1;
#                tree[q] := [ p, s ];
#            else
#                Add( genTriples, [ p, s, v, g ] );
#            fi;
#        od;
#        RemoveDictionary( leaves, v );
#    od;
#    vertex2word := function ( v )
#          local word;
#          word := one;
#          while IsBound( tree[v] ) do
#              word := word * ambientGenerators[tree[v][2]];
#              v := tree[v][1];
#          od;
#          return word;
#      end;
#    triple2word := function ( x )
#          local u, uu, g, q, c;
#          for u in Ugrp do
#              c := x[4] ^ -1 * u * vertex2word( G!.cosetPos( x[4] ) );
#              if G!.membership( c ) then
#                  return c;
#              fi;
#          od;
#          return fail;
#      end;
#    genTriples := List( genTriples, function ( x )
#            return triple2word( x );
#        end );
#    genTriples := List( genTriples, function ( x )
#            return Minimum( x, x ^ -1 );
#        end );
#    genTriples := SSortedList( genTriples );
#    G!.tree := tree;
#    G!.GeneratorsOfMagmaWithInverses := genTriples;
#    return;
#end

K := ContractibleGcomplex("SL(3,Z)s");
KK:=BarycentricSubdivision(K);
G := HAPCongruenceSubgroupGamma0(3,2);
HAPCongruenceSubgroupTree(G);
#gens := G!.GeneratorsOfMagmaWithInverses;
#G := Group(gens);
SetGeneratorsOfGroup(G,G!.GeneratorsOfMagmaWithInverses);
Y := GComplexToRegularCWComplex(KK,G);

#prova := function ( R, N )
#    local Y, dim, dim1, Cells, boundaries, StandardForm, Tiles, tile, pos, #NeighbourGens, 
#    TileReps, w, NewLeaves, Leaves, i, e, f, n, k, g, b, x, T;
#    T := RightTransversal( R!.group, N );
#    dim := 0;
#    for n in [ 1 .. Length( R ) ] do
#        if R!.dimension( n ) > 0 then
#            dim := dim + 1;
#        else
#            break;
#        fi;
#    od;
#    dim1 := dim - 1;
#    StandardForm := function ( n, x )
#          local e, g, gc, stab;
#          e := AbsInt( x[1] );
#          g := x[2];
#          stab := R!.stabilizer( n, e );
#          gc := g * Elements( stab );
#          gc := List( gc, function ( a )
#                  return T[T!.poscan( a )];
#              end );
#          gc := Minimum( gc );
#          return [ e, gc ];
#      end;
#    Tiles := [  ];
#    for k in [ 1 .. R!.dimension( dim ) ] do
#        tile := 1 * R!.boundary( dim, k );
#        Apply( tile, function ( x )
#              return [ x[1], R!.elts[x[2]] ];
#          end );
#        Apply( tile, function ( e )
#              return StandardForm( dim1, e );
#          end );
#        Add( Tiles, tile );
#    od;
#    TileReps := List( [ 1 .. Length( Tiles ) ], function ( i )
#            return T;
#        end );
#    Cells := List( [ 1 .. dim + 1 ], function ( i )
#            return [  ];
#        end );
#    for k in [ 1 .. R!.dimension( dim ) ] do
#        for g in TileReps[k] do
#            Add( Cells[dim + 1], StandardForm( dim, [ k, g ] ) );
#        od;
#    od;
#    Cells[dim + 1] := SSortedList( Cells[dim + 1] );
#    for n in Reversed( [ 1 .. dim ] ) do
#        for x in Cells[n + 1] do
#            b := 1 * R!.boundary( n, x[1] );
#            Apply( b, function ( z )
#                  return [ z[1], R!.elts[z[2]] ];
#              end );
#            Apply( b, function ( y )
#                  return [ y[1], x[2] * y[2] ];
#              end );
#            Apply( b, function ( y )
#                  return StandardForm( n - 1, y );
#              end );
#            Append( Cells[n], b );
#        od;
#        Cells[n] := SSortedList( Cells[n] );
#    od;
#    boundaries := List( [ 1 .. dim + 2 ], function ( n )
#            return [  ];
#        end );
#    boundaries[1] := List( Cells[1], function ( x )
#            return [ 1, 0 ];
#        end );
#    for n in [ 1 .. dim ] do
#        for k in [ 1 .. Length( Cells[n + 1] ) ] do
#            x := Cells[n + 1][k];
#            b := 1 * R!.boundary( n, x[1] );
#            Apply( b, function ( z )
#                  return [ z[1], R!.elts[z[2]] ];
#              end );
#            Apply( b, function ( y )
#                  return [ y[1], x[2] * y[2] ];
#              end );
#            Apply( b, function ( y )
#                  return StandardForm( n - 1, y );
#              end );
#            Apply( b, function ( c )
#                  return PositionSorted( Cells[n], c );
#              end );
#            b := Concatenation( [ Length( b ) ], b );
#            boundaries[n + 1][k] := b;
#        od;
#    od;
#    Y := RegularCWComplex( boundaries );
#    return Y;
#end;

N := HAP_CongruenceSubgroupGamma0(2);;
gens:=GeneratorsOfGroup(N);;
x:=gens[1]*gens[2];;
x in N;

N := HAPCongruenceSubgroupGamma0(3,2);;
HAPCongruenceSubgroupTree(N);
SetGeneratorsOfGroup(N,N!.GeneratorsOfMagmaWithInverses);
gens:=GeneratorsOfGroup(N);;
x:=gens[1]*gens[2];;
x in N;
#Error, no method found! For debugging hints type ?Recovery from NoMethodFound
#Error, no 3rd choice method found for `Enumerator' on 1 arguments at /usr/share/#gap-4.14.0/lib/methsel2.g:250 called from
#Enumerator( D ) at /usr/share/gap-4.14.0/lib/domain.gi:285 called from
#<function "in for a domain, and an element">( <arguments> )
# called from read-eval loop at *stdin*:139
#type 'quit;' to quit to outer loop