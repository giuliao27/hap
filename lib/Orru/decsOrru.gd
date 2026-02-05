IsHapCongruenceSubgroup:=NewFilter("IsHapCongruenceSubgroup");;
IsHapRightTransversalSL3ZSubgroup:=NewFilter("IsHapRightTransversalSL3ZSubgroup");;
DeclareGlobalFunction("HAP_GenericCongruenceSubgroup");
DeclareGlobalFunction("FiniteProjectiveLine");
DeclareGlobalFunction("FiniteProjectivePlane");
DeclareGlobalFunction("HAP_SL3ZSubgroupTree_fast");
DeclareGlobalFunction("HAP_TransversalCongruenceSubgroups_SL3Z");
DeclareOperation("HAPCongruenceSubgroupGamma0",[IsInt,IsInt]);
DeclareOperation("HAPCongruenceSubgroupTree",[IsHapCongruenceSubgroup]);
DeclareOperation("in",[IsMatrix, IsHapCongruenceSubgroup and IsGroup]);
DeclareOperation("RightTransversal",[IsMatrixGroup, IsHapCongruenceSubgroup]);

InstallMethod( ViewObj,
"for HapCongruenceSubgroup",
[IsHapCongruenceSubgroup and IsGroup],
10000000,  #Ensures that this method is chosen
function(G)
Print(G!.name," of ",G!.fam);
end);

InstallMethod( PrintObj,
"for HapCongruenceSubgroup",
[IsHapCongruenceSubgroup and IsGroup],
100000000, #Ensures that this method is chosen
function(G)
Print(G!.name," of ",G!.fam);
end);

