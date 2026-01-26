
##############################################################
##############################################################
InstallGlobalFunction(ResolutionIsomorphismGroup,
function(arg)
local K,hom,fn,G,HH, STAB, ROT, n, k, ln, x,g;
#Inputs a (possibly non-free) resolution K for G and isomorphism hom:G->H.
#Destructively converts K to a resolution for H.

K:=arg[1]; hom:=arg[2];

if Length(arg)=2 then
   fn:=function(x) return ImageElm(hom,x); end;
   HH:=Range(hom);
else
   fn:=hom;
   HH:=arg[3];
fi;

##################################


if IsHapResolution(K) then
   K!.group:=Range(hom);
   if IsMutable(K!.elts) then
      Apply(K!.elts,fn);
   else
      K!.elts:=List(K!.elts,fn);
   fi;
   return true;
fi;


#################################

#################################
K!.group:=HH;

ln:=0;  #The dimension of the complex
for n in [1..100] do
if K!.dimension(n)=0 then break; fi;
ln:=ln+1;
od;

STAB:=[];
ROT:=[];
   for n in [0..ln] do
   STAB[n+1]:=[];
   ROT[n+1]:=[];
      for k in [1..K!.dimension(n)] do
if n>0 then
      x:=Elements(K!.stabilizer(n,k));
      for g in x do
          if not g in K!.elts then Add(K!.elts,g); fi;
      od;
      x:=Filtered(x,g->K!.action(n,k,Position(K!.elts,g))=1);
      x:=List(x,fn);
      x:=Group(x);
      Add(ROT[n+1],x);
fi;
      x:=GeneratorsOfGroup(K!.stabilizer(n,k));
      x:=List(x,fn);
      x:=Group(x);
      Add(STAB[n+1],x);
   od;
od;

##############################
K!.stabilizer:=function(n,k);
    return STAB[n+1][k];
end;
##############################

   if IsMutable(K!.elts) then
      Apply(K!.elts,fn);
   else
      K!.elts:=List(K!.elts,fn);
   fi;


##############################
K!.action:=function ( n, k, g )
    local id, r, u, ans, abk, H;
    if n = 0 then
        return 1;
    fi;
    abk := AbsInt( k );
    H := K!.stabilizer( n, abk );
    if Order( H ) = infinity then
        return 1;
    fi;
    id := CanonicalRightCosetElement( H, Identity( H ) );
    r := CanonicalRightCosetElement( H, K!.elts[g] ^ -1 );
    r := id ^ -1 * r;
    u := r * K!.elts[g];
    if u in ROT[n+1][abk] then
        ans := 1;
    else
        ans := -1;
    fi;
    return ans;
end;

##############################


return true;
#################################


end);
##############################################################
##############################################################

