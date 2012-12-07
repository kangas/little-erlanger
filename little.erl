%% ---------------------------------------------------------------
%% The Little Schemer – implemented in Erlang
%% Copyright 2008 Matt Kangas. All rights reserved.
%% ---------------------------------------------------------------

-module(little).
-compile(export_all).

test_all() ->
	chapter2_test(),
	chapter3_test(),
	chapter4_test(),
	chapter5_test(),
	chapter6_test(),
	chapter7_test(),
	chapter8_test(),
	chapter9_test(),
	chapter10_test(),
    test_ok.

%% 
%% FIXME: Suggested by jared@alloycode.com (blog: http://alloycode.com/)
%%
%% "It would definitely be interesting to see the 10 commandments 
%% and the 5 rules cast into Erlang, since asking: 
%% (null? lat) doesn't intuitively translate into writing: 
%% lat([]) -> true; lat([H|T]) -> ..."
%%

%% ---------------------------------------------------------------
%% Chapter 1. Little Scheme Primitives
%% ---------------------------------------------------------------

% car/cdr/cons is performed with pattern-matching syntax in Erlang
% "[H|T]" splits a list into "Head" and "Tail"
% Also: "hd()" returns head of list, "tl()" returns tail

car([Head|_]) -> Head.

cdr([_|Tail]) -> Tail.

cons(S, L) -> [S|L].

% null?() and atom?() are Erlang built-in functions: is_null(), is_atom()
% eq?() is expressed as:
% 	=:=	"is identical to"
%  	==	"is equal to", ONLY useful for comparing floats with integers
is_eq(X, Y) -> X =:= Y.


%% ---------------------------------------------------------------
%% Chapter 2: Recurring on lists of atoms ("lats")
%% ---------------------------------------------------------------

% Unit tests: "badmatch" exception is thrown if the return value is wrong
chapter2_test() ->
	Lat = [jack,sprat,could],
	true 	= is_lat(Lat), 
	false 	= is_lat([Lat | [eat,no]]),
	true 	= is_lat2(Lat), 
	false 	= is_lat2([Lat | [eat,no]]),
	true 	= is_lat3(Lat), 
	false 	= is_lat3([Lat | [eat,no]]),
	true 	= is_member(sprat,Lat), 
	false 	= is_member(foo,Lat),
	true 	= is_member2(sprat,Lat), 
	false 	= is_member2(foo,Lat),
	test_ok.

% -------

% "is_list_of_atoms", version 1:
% We express our base case as a separate function clause 
% (think of it as a form of polymorphism).
is_lat([]) -> true;
is_lat([H|T]) -> 
	case is_atom(H) of
		true -> is_lat(T);
		false -> false
	end.

% "is_list_of_atoms", version 2:
% Closer in spirit to Scheme's "cond()"
is_lat2([]) -> true;
is_lat2([H|T]) -> 
	if 
		is_atom(H) -> is_lat2(T);
		true -> false	% looks wierd, but is an "else" clause
	end.

% "is_list_of_atoms", version 3:
% Better Erlang style: put the test into a "guard clause".
% However, not possible when test is a custom function
is_lat3([]) -> true;
is_lat3([H|T]) when is_atom(H) -> is_lat3(T);
is_lat3(_) -> false.


% "is_member", using pattern-matching for our tests
% Will stop testing list at first matching atom (short-circuit behavior)
% Offers same functionality as "lists:member()" library function
is_member(_, []) -> false;
is_member(A, [A|_]) -> true;			% When (A == H)
is_member(A, [_|T]) -> is_member(A, T). % Only when H != A

% "is_member", version 2, without pattern matching.
% The "or" test evaluates *both* sides of the expresson, so unlike above, 
% we always traverse the whole list (needlessly).
% For correct short-circuit behavior, change "or" to "orelse"
is_member2(_, []) -> false;
is_member2(A, [H|T]) -> 
	(A == H) or is_member2(A,T).

% OR: just use lists:member()! :-)

%% ---------------------------------------------------------------
%% Chapter 3: Producing "lats" (lists of atoms)
%% ---------------------------------------------------------------

chapter3_test() ->
	Lat = [bacon,lettuce,'and',tomato],		% need to quote keyword 'and'
	[bacon,lettuce,tomato] = rember('and', Lat),
	[a,c,e] = firsts([[a,b],[c,d],[e,f]]),
	L2 = [ice,cream,with,fudge,for,dessert],
	[ice,cream,with,fudge,topping,for,dessert] = insertR(topping, fudge, L2),
	[ice,cream,with,topping,fudge,for,dessert] = insertL(topping, fudge, L2),
	[ice,cream,with,topping,for,dessert] = subst(topping, fudge, L2),	
	L3 = [banana,ice,cream,with,chocolate,topping],
	[vanilla,ice,cream,with,chocolate,topping] = subst2(vanilla, chocolate, banana, L3),
	L4 = [coffee,cup,tea,cup,'and',hick,cup],
	[coffee,tea,'and',hick] = multirember(cup, L4),
	L5 = [chips,'and',fish,'or',fish],
	[chips,'and',fish,fried,'or',fish,fried] = multiinsertR(fried,fish,L5),
	[chips,'and',fried,fish,'or',fried,fish] = multiinsertL(fried,fish,L5),
	[chips,'and',fried,'or',fried] = multisubst(fried,fish,L5),
	test_ok.

% -------

rember(_, []) -> [];
rember(A, [A|T]) -> T;          % When (H =:= A)
rember(A, [H|T]) -> [H | rember(A,T)].

firsts([]) -> [];
firsts([H|T]) -> 
	[hd(H) | firsts(T)]. % use hd() to replace another [H|T] pattern-match

insertR(_, _, []) -> [];
insertR(New, Old, [Old|T]) -> [Old | [New | T]];
insertR(New, Old, [H|T])   -> [H | insertR(New, Old, T)].

insertL(_, _, []) -> [];
insertL(New, Old, [Old|T]) -> [New | [Old | T]];
insertL(New, Old, [H|T])   -> [H | insertL(New, Old, T)].

subst(_, _, []) -> [];
subst(New, Old, [Old|T]) -> [New | T];
subst(New, Old, [H|T])   -> [H | subst(New, Old, T)].

% Use a guard clause ("when") to enhance pattern-matching
subst2(_, _, _, []) -> [];
subst2(New, O1, O2, [H|T]) when (H =:= O1) or (H =:= O2) -> [New | T];
subst2(New, O1, O2, [H|T]) -> [H | subst2(New, O1, O2, T)].

multirember(_, []) -> [];
multirember(A, [A|T]) -> multirember(A,T);          % If (H == A)
multirember(A, [H|T]) -> [H | multirember(A,T)].    % If (H != A)

multiinsertR(_, _, []) -> [];
multiinsertR(New, Old, [Old|T]) -> [Old | [New | multiinsertR(New, Old, T)]];
multiinsertR(New, Old, [H|T])   -> [H | multiinsertR(New, Old, T)].

multiinsertL(_, _, []) -> [];
multiinsertL(New, Old, [Old|T]) -> [New | [Old | multiinsertL(New, Old, T)]];
multiinsertL(New, Old, [H|T])   -> [H | multiinsertL(New, Old, T)].

multisubst(_, _, []) -> [];
multisubst(New, Old, [Old|T]) -> [New | multisubst(New, Old, T)];
multisubst(New, Old, [H|T])   -> [H | multisubst(New, Old, T)].



%% ---------------------------------------------------------------
%% Chapter 4: Recurring on numbers
%%
%% "Tuple" has a different meaning in Erlang than a list of numbers,
%% but we follow "The Little Schemer" here
%%
%% Erlang has many built-in shortcuts for these tasks,
%% often in the "lists" module. See "Prefer:" comments below
%% ---------------------------------------------------------------

chapter4_test() ->
	5             = plus(3,2),
	1             = minus(3,2),
	10            = addtup([4,3,2,1]),
	10            = mul(2,5),
	[3,7,11]      = tupplus([1,3,5], [2,4,6]),
	[9,11,4]      = tupplus([2,3,4], [7,8]),
	[true,false,false] = lists:map(fun(X) -> lt(X,5) end, [4,5,6]),
	[false,false,true] = lists:map(fun(X) -> gt(X,5) end, [4,5,6]),
	[false,true,false] = lists:map(fun(X) -> equal(X,5) end, [4,5,6]),
	[8, 27, 64]   = lists:map(fun(X) -> pow(X,3) end, [2,3,4]),
	[1,3,3]       = lists:map(fun(X) -> quotient(X,4) end, [6,12,13]),
	5             = len([a,b,c,d,e]),
	c             = pick(3, [a,b,c,d,e]),
	[a,b,d,e]     = rempick(3, [a,b,c,d,e]),
	[pears,prunes,dates] = no_nums([5,pears,6,prunes,9,dates]),
	[5,6,9]       = all_nums([5,pears,6,prunes,9,dates]),
	[true,false,false] = lists:map(fun(X)->is_eqan(X,4) end, [4,5,foo]),
	3             = occur(a, [a,b,a,c,a,b]),
	3             = occur2(a, [a,b,a,c,a,b]),
	[false,true,false] = lists:map(fun is_one/1, [0,1,2]),
	test_ok.

%% Primitives used for this chapter
add1(N) -> N + 1.
sub1(N) -> N - 1.
is_zero(N) -> N =:= 0.
% Note: is_number is a built-in function (BIF)

% -------

plus(N, 0) -> N;
plus(N, M) -> add1(plus(N, sub1(M))).

minus(N, 0) -> N;
minus(N, M) -> sub1(minus(N, sub1(M))).

% Prefer lists:sum() normally in Erlang. This follows the Scheme version.
addtup([]) -> 0;
addtup([H|T]) -> plus(H, addtup(T)).

% Prefer: "*" operator
mul(_, 0) -> 0;
mul(N, M) -> plus(N, mul(N, sub1(M))).

tupplus([], L) -> L;
tupplus(L, []) -> L;
tupplus([H|T], [H2|T2]) ->  [plus(H,H2) | tupplus(T,T2)].

% Prefer: ">" operator
gt(0, _) -> false;
gt(_, 0) -> true;
gt(N, M) -> gt(sub1(N), sub1(M)).

% Prefer: "<" operator
lt(_, 0) -> false;
lt(0, _) -> true;
lt(N, M) -> lt(sub1(N), sub1(M)).

% Prefer: "=:=" operator
equal(0, 0) -> true;
equal(0, _) -> false;
equal(_, 0) -> false;
equal(N, M) -> equal(sub1(N), sub1(M)).

pow(_, 0) -> 1;
pow(N, M) -> mul(N, pow(N, sub1(M))).

% Prefer: "div" operator for quotient, "/" for floating-point division
quotient(N, M) when N < M -> 0;
quotient(N, M) -> add1(quotient(minus(N,M), M)).

% Prefer: erlang:length() BIF
len([]) -> 0;
len([_|T]) -> add1(len(T)).

% Prefer: lists:nth()
pick(1, [H|_]) -> H;
pick(N, [_|T]) -> pick(sub1(N), T).

rempick(1, [_|T]) -> T;
rempick(N, [H|T]) -> [H | rempick(sub1(N),T)].

% Prefer: lists:filter(fun(X)->not is_number(X) end, ...)
no_nums([]) -> [];
no_nums([H|T]) when is_number(H) -> no_nums(T);
no_nums([H|T]) -> [H | no_nums(T)].

% Prefer: lists:filter(fun(X)->is_number(X) end, ...)
all_nums([]) -> [];
all_nums([H|T]) when is_number(H) -> [H | all_nums(T)];
all_nums([_|T]) -> all_nums(T).

% Prefer: "=:=" operator, of course
is_eqan(A1, A2) when (is_number(A1) and is_number(A2)) -> equal(A1, A2);
is_eqan(A1, A2) when (is_number(A1) or is_number(A2)) -> false;
is_eqan(A1, A2) -> is_eq(A1, A2).

occur(_, []) -> 0;
occur(A, [H|T]) ->		% Can't use is_eqan() in guard clause
	case is_eqan(A,H) of
		true -> add1(occur(A, T));
 		false -> occur(A, T)
	end.

% Alternate version using operator "=:=" to test
occur2(_, []) -> 0;
occur2(A, [H|T]) when A =:= H -> add1(occur2(A, T));
occur2(A, [_|T]) -> occur(A, T).

is_one(0) -> false;
is_one(N) -> is_zero(sub1(N)).


%% ---------------------------------------------------------------
%% Chapter 5: Recurring on S-expressions
%% ---------------------------------------------------------------

chapter5_test() ->
	[[coffee], [[tea]]] = remberstar(cup, [[coffee], cup, [[tea], cup]]),
    [[a,[b],c,d], [[c,d]]] = insertRstar(d, c, [[a,[b],c], [[c]]]),
    5      = occurstar(a, [[a],[b,[[[a,c]]],[d,[a]]],[a],[e],[a,e]]),
    [[x],[b,[[[x,c]]]]] = subststar(x, a, [[a],[b,[[[a,c]]]]]),
    [[a,[b],d,c], [[d,c]]] = insertLstar(d, c, [[a,[b],c], [[c]]]),
    true   = memberstar(chips, [[potato],[chips,[[with],fish]]]),
    potato = leftmost([[potato],[chips,[[with],fish]]]),
    true   = is_eqlist([strawberry,ice,cream], [strawberry,ice,cream]),
    false  = is_eqlist([strawberry,ice,cream], [strawberry,cream,ice]),
    true   = is_equal([[strawberry],[ice,cream]], [[strawberry],[ice,cream]]),
    false  = is_equal([[strawberry],[ice,cream]], [[strawberry],[cream,ice]]),
    true   = is_eqlist2([strawberry,ice,cream], [strawberry,ice,cream]),
    false  = is_eqlist2([strawberry,ice,cream], [strawberry,cream,ice]),
    [[bacon,lettuce],tomato] = rember2([foo], [[bacon,lettuce],[foo],tomato]),
	test_ok.

% -------
% "..star" functions recur on *both* H and T (car and cdr)

remberstar(_, []) -> [];
remberstar(A, [H|T]) when is_atom(H) ->
    case H =:= A of
        true -> remberstar(A, T);
        false -> [H | remberstar(A,T)]
    end;
remberstar(A, [H|T]) -> 
    [remberstar(A, H) | remberstar(A, T)].

insertRstar(_, _, []) -> [];
insertRstar(New, Old, [H|T]) when is_atom(H) -> 
    case H =:= Old of
        true -> [Old | [New | insertRstar(New, Old, T)]];
        false -> [H | insertRstar(New, Old, T)]
    end;
insertRstar(New, Old, [H|T]) -> 
    [insertRstar(New, Old, H) | insertRstar(New, Old, T)].

occurstar(_, []) -> 0;
occurstar(A, [H|T]) when is_atom(H) -> 
    case H =:= A of
        true -> add1(occurstar(A, T));
        false -> occurstar(A, T)
    end;
occurstar(A, [H|T]) -> 
    plus(occurstar(A, H), occurstar(A, T)).     % prefer "+" operator

subststar(_, _, []) -> [];
subststar(New, Old, [H|T]) when is_atom(H) -> 
    case H =:= Old of
        true -> [New | subststar(New, Old, T)];
        false -> [H | subststar(New, Old, T)]
    end;
subststar(New, Old, [H|T]) -> 
    [subststar(New, Old, H) | subststar(New, Old, T)].

insertLstar(_, _, []) -> [];
insertLstar(New, Old, [H|T]) when is_atom(H) -> 
    case H =:= Old of
        true -> [New | [Old | insertLstar(New, Old, T)]];
        false -> [H | insertLstar(New, Old, T)]
    end;
insertLstar(New, Old, [H|T]) -> 
    [insertLstar(New, Old, H) | insertLstar(New, Old, T)].

memberstar(_, []) -> false;
memberstar(A, [H|T]) when is_atom(H) -> 
    (H =:= A) or memberstar(A, T);
memberstar(A, [H|T]) -> 
    memberstar(A, H) or memberstar(A, T).

leftmost([]) -> [];
leftmost([H|_]) when is_atom(H) -> H;
leftmost([H|_]) -> leftmost(H).

% First attempt at comparing lists; input *must* be lists
is_eqlist([], []) -> true;
is_eqlist(L1, L2) when ([]=:=L1) or ([]=:=L2) -> false;
is_eqlist([H1|T1], [H2|T2]) when is_atom(H1) and is_atom(H2) -> 
    (H1 =:= H2) and is_eqlist(T1, T2);
is_eqlist([H1|_], [H2|_]) when is_atom(H1) or is_atom(H2) -> false;
is_eqlist([H1|T1], [H2|T2]) -> 
    is_eqlist(H1, H2) and is_eqlist(T1, T2).

% Input may be S-exps
is_equal(S1, S2) when is_atom(S1) and is_atom(S2) -> S1 =:= S2;
is_equal(S1, S2) when is_atom(S1) or is_atom(S2) -> false;
is_equal(S1, S2) -> is_eqlist(S1, S2).

% Second attempt; mutually recur with is_equal()
is_eqlist2([], []) -> true;
is_eqlist2(L1, L2) when ([]=:=L1) or ([]=:=L2) -> false;
is_eqlist2([H1|T1], [H2|T2]) -> 
    is_equal(H1, H2) and is_eqlist2(T1, T2).

% Can't put a custom function (is_equal) into function guard clause,
% so we must use "case"
rember2(_, []) -> [];
rember2(S, [H|T]) -> 
    case is_equal(S, H) of
        true -> T;
        false -> [H | rember2(S, T)]
    end.


%% ---------------------------------------------------------------
%% Chapter 6: Calculations
%% ---------------------------------------------------------------

chapter6_test() ->
    true    = is_numbered(1),
    true    = is_numbered([2,plus,[4,mul,5]]),
    false   = is_numbered([2,mul,sausage]),
    13      = value(13),
    4       = value([1,'+',3]),
    82      = value([1,'+',[3,'^',4]]),
    13      = valueP(13),
    4       = valueP(['+',1,3]),
    82      = valueP(['+',1,['^',3,4]]),

    % Make sure valueP properly coughs up bad_operator exceptions
    try valueP(['+',1,['/',3,4]]) of
        _ -> erlang:error("failed to detect bad operator~n")
    catch 
        throw:{bad_operator,'/'} -> ok
    end,

    true    = is_sero([]),
    false   = is_sero([[]]),
    [[]]    = edd1([]),
    [[]]    = zub1([[],[]]),
    [[],[],[]] = splus([[]], [[],[]]),
    test_ok.

% -------
% "is_atom()" and "is_number()" are mutually-exclusive in Erlang
% e.g. is_atom(1) -> false, is_number(1) -> true

%% Infix notation (standard arithmetic: "2 + 3")
is_numbered(Aexp) when is_number(Aexp) -> true;
is_numbered(Aexp) when is_atom(Aexp) -> false;
is_numbered([]) -> true;
is_numbered([H|T]) -> 
    is_numbered(H) and is_numbered(hd(tl(T))).

value(Nexp) when is_number(Nexp) -> Nexp;
value(Nexp) when is_atom(Nexp) -> none;
value([Nexp1 | T]) -> 
    [Op,Nexp2] = T,
    case Op of
        '+' -> plus(value(Nexp1), value(Nexp2));
        '-' -> minus(value(Nexp1), value(Nexp2));
        '^' -> pow(value(Nexp1), value(Nexp2));
        _ -> throw({bad_operator, Op})
    end.

%% Now use prefix notation ("+ 2 3")
first_subexp([_|T]) -> hd(T).
second_subexp([_|T]) -> hd(tl(T)).
operator([H|_]) -> H.

valueP(Nexp) when is_number(Nexp) -> Nexp;
valueP(Nexp) when is_atom(Nexp) -> none;
valueP(Nexp) -> 
    case operator(Nexp) of
        '+' -> plus(valueP(first_subexp(Nexp)), valueP(second_subexp(Nexp)));
        '-' -> minus(valueP(first_subexp(Nexp)), valueP(second_subexp(Nexp)));
        '^' -> pow(valueP(first_subexp(Nexp)), valueP(second_subexp(Nexp)));
        _ -> throw({bad_operator, operator(Nexp)})
    end.

%% With this version of "value", we could redefine our helper functions
%% (first_subexp, second_subexp, operator) to handle infix notation!

%% Now use a braces-counting notation
is_sero([]) -> true;
is_sero(_) -> false.

edd1(N) -> [[] | N].

zub1([_|T]) -> T.

splus(N, M) -> 
    case is_sero(M) of
        true -> N;
        false -> edd1(splus(N, zub1(M)))
    end.


%% ---------------------------------------------------------------
%% Chapter 7: Sets, Pairs/Relations/Functions
%%
%% See the "sets" module for built-in helpers for many of these operations
%% ---------------------------------------------------------------

chapter7_test() ->
	false   = is_set([apple,peaches,apple,plum]),
	true    = is_set([apple,peaches,pears,plum]),
	true    = is_set([]),
	Notset  = [apple,peach,pear,peach,plum,apple,lemon,peach],
	[pear,plum,apple,lemon,peach] = makeset(Notset),
	[apple,peach,pear,plum,lemon] = makeset2(Notset),
	true    = is_subset([a,b,c], [a,j,b,k,c,l,m]),
	false   = is_subset([a,b,c], [a,j,b,k,l,m,n,o]),
	true    = is_eqset([6,large,chickens,with,wings], [6,chickens,with,large,wings]),
	false   = is_eqset([6,large,chickens,with,wings], [6,large,wings]),
	true    = is_intersect([stewed,tomatoes,'and',macaroni], [macaroni,'and',cheese]),
	false   = is_intersect([a,b,c], [d,e,f]),
	[a,c]   = intersect([a,b,c], [d,a,f,c]),
	[]      = intersect([a,b,c], [d,e,f]),
	[b,d,a,f,c] = union([a,b,c], [d,a,f,c]),
	[a,b,c,d,f] = union_fwd([a,b,c], [d,a,f,c]),
	[a,b,c,d,e,f] = union([a,b,c], [d,e,f]),
	[b]     = difference([a,b,c], [d,a,f,c]),   % Set1 atoms *not* in Set2
	[a,b,c] = difference([a,b,c], [d,e,f]),
	[a]     = intersectall([[a,b,c],[c,a,d,e],[e,f,g,h,a,b]]),
	true    = is_a_pair([pear,pair]),
	true    = is_a_pair([3,7]),
	true    = is_a_pair([[2],[pair]]),
	false   = is_a_pair([[1],2,[3]]),
	a       = first([a,b]),
	b       = second([a,b]),
	[a,b]   = build(a,b),
	c       = third([a,b,c]),
	false   = is_fun([[4,3], [4,2], [7,6], [6,2], [3,4]]),
	true    = is_fun([[8,3], [4,2], [7,6], [6,2], [3,4]]),
	[a,8]   = revpair([8,a]),
	[[a,8],[pie,pumpkin],[sick,got]] = revrel([[8,a], [pumpkin,pie], [got,sick]]),
	[[a,8],[pie,pumpkin],[sick,got]] = revrel2([[8,a], [pumpkin,pie], [got,sick]]),
	[3,2,6,2,4] = seconds([[8,3], [4,2], [7,6], [6,2], [3,4]]),
	false   = is_fullfun([[8,3], [4,2], [7,6], [6,2], [3,4]]),
	true    = is_fullfun([[8,3], [4,8], [7,6], [6,2], [3,4]]),
	true    = is_1to1([[grape,raisin], [plum,prune], [stewed,grape]]),
	test_ok.

% -------
% "is_member()" was defined in Chapter 2. 
% We *could* use that below, but I'm going to use lists:member() 
% for better Erlang style...

is_set([]) -> true;
is_set([H|T]) ->
    case lists:member(H, T) of
        true -> false;
        false -> is_set(T)
    end.

makeset([]) -> [];
makeset([H|T]) -> 
    case lists:member(H, T) of
        true -> makeset(T);
        false -> [H | makeset(T)]
    end.

makeset2([]) -> [];
makeset2([H|T]) -> [H | makeset2( multirember(H, T) )].

is_subset([], _) -> true;
is_subset([H|T], Set2) ->
    lists:member(H, Set2) and is_subset(T, Set2).

is_eqset(Set1, Set2) ->
    is_subset(Set1, Set2) and is_subset(Set2, Set1).

is_intersect([], _) -> false;
is_intersect([H|T], Set2) ->
    lists:member(H, Set2) or is_intersect(T, Set2).

intersect([], _) -> [];
intersect([H|T], Set2) ->
    case lists:member(H, Set2) of
        true -> [H | intersect(T, Set2)];
        false -> intersect(T, Set2)
    end.

union([], Set2) -> Set2;
union([H|T], Set2) ->
    case lists:member(H, Set2) of
        true -> union(T, Set2);
        false -> [H | union(T, Set2)]
    end.

union_fwd([], Set2) -> Set2;
union_fwd([H|T], Set2) ->
    [H | union_fwd(T, multirember(H, Set2))].

% Returns atoms in Set1 *not* in Set2. It's somewhat misnamed...
difference([], _) -> [];
difference([H|T], Set2) ->
    case lists:member(H, Set2) of
        true -> difference(T, Set2);
        false -> [H | difference(T, Set2)]
    end.

%% Prefer lists:merge(ListOfLists) when input lists are sorted
intersectall([]) -> [];
intersectall([H | []]) -> H;
intersectall([H|T]) ->
    intersect(H, intersectall(T)).  % recursively apply intersect() to subsets

is_a_pair(L) -> 2 =:= length(L).

first([H|_]) -> H.

second([_|T]) -> hd(T).

build(S1, S2) -> [S1 | [S2 | []]].

third([_|T]) -> hd(tl(T)).

is_fun(Rel) ->
    is_set(firsts(Rel)).

revrel([]) -> [];
revrel([H|T]) ->
    [build(second(H),first(H)) | revrel(T)].

revpair(P) -> 
    build(second(P), first(P)).

revrel2([]) -> [];
revrel2([H|T]) -> 
    [revpair(H) | revrel2(T)].

% Let's make this a one-liner for kicks
seconds(L) -> lists:map(fun(X) -> second(X) end, L).

is_fullfun(Fun) ->
    is_set(seconds(Fun)).

is_1to1(Fun) ->
    is_fun(revrel(Fun)).


%% ---------------------------------------------------------------
%% Chapter 8: Higher-order functions
%% ---------------------------------------------------------------

chapter8_test() ->
    % "fun is_eq/2" creates a function pointer for the named function 
    % "is_eq" (arity=2)
    [6,2,3] = remberF(fun is_eq/2, 5, [6,2,5,3]),
    [beans,are,good] = 
        remberF(fun is_eq/2, jelly, [jelly,beans,are,good]),
    [lemonade,'and',[cake]] = 
        remberF(fun is_eq/2, [pop,corn], [lemonade,[pop,corn],'and',[cake]]),
    
    % Call curried function, then call result
    Is_eq_salad = is_eqC(salad),
    true    = Is_eq_salad(salad),
    true    = (is_eqC(salad)) (salad),
    
    Rember_eq = remberFC(fun is_eq/2),
    [salad,is,good] = Rember_eq(tuna, [tuna,salad,is,good]),
    [salad,is,good] = (remberFC(fun is_eq/2)) (tuna, [tuna,salad,is,good]),
    
    L2 = [ice,cream,with,fudge,for,dessert],
    Result_InsertL = [ice,cream,with,topping,fudge,for,dessert],
    Result_InsertR = [ice,cream,with,fudge,topping,for,dessert],
    
    Result_InsertL = (insertL_F(fun is_eq/2)) (topping, fudge, L2),
    Result_InsertR = (insertR_F(fun is_eq/2)) (topping, fudge, L2),
    
    % Test our "generic insert"
    InsertL = insertG(fun sis_eqL/3),
    Result_InsertL = InsertL(topping, fudge, L2),
    % just call the lambda, don't name it
    Result_InsertL = (insertG(fun sis_eqL/3)) (topping, fudge, L2),
    Result_InsertR = (insertG(fun sis_eqR/3)) (topping, fudge, L2),
    % test wrappers
    Result_InsertL = insertL2(topping, fudge, L2),
	Result_InsertR = insertR2(topping, fudge, L2),
    % try without helpers or wrappers!
    Result_InsertL = (insertG(fun (New,Old,L) -> [New | [Old | L]] end)) 
                        (topping, fudge, L2),
    
	[ice,cream,with,topping,for,dessert] = substG(topping, fudge, L2),	
    [salad,is,good] = remberG(tuna, [tuna,salad,is,good]),
    
    82 = valueG(['+',1,['^',3,4]]),
    
    % bad_operator test
    try valueG(['+',1,['/',3,4]]) of
        _ -> erlang:error("failed to detect bad operator~n")
    catch 
        throw:{bad_operator,'/'} -> ok
    end,
    
    L4 = [shrimp,salad,tuna,salad,'and',tuna],
	[shrimp,salad,salad,'and'] = (multiremberF(fun is_eq/2)) (tuna, L4),
    [shrimp,salad,salad,'and'] = multirember_iseq(tuna, L4),
    [shrimp,salad,salad,'and'] = multiremberT(fun is_eq_tuna/1, L4),
    
    true  = multiremberCO(tuna, [], fun a_friend/2),
    false = multiremberCO(tuna, [tuna], fun a_friend/2),
    false = multiremberCO(tuna, ['and',tuna], fun a_friend/2),
    3 = multiremberCO(tuna, [strawberries,tuna,'and',swordfish], fun last_friend/2),
    
    [chips,salty,'and',salty,fish,'or',salty,fish,'and',chips,salty] =
        multiinsertLR(salty, fish, chips, 
            [chips,'and',fish,'or',fish,'and',chips]),
    
    [[], 0, 0] 
        = multiinsertLRCO([cranberries], [fish], [chips], [], fun(X,Y,Z)->[X,Y,Z] end),
    
    [[chips,salty,'and',salty,fish,'or',salty,fish,'and',chips,salty], 2, 2] 
        = multiinsertLRCO(salty, fish, chips, 
            [chips,'and',fish,'or',fish,'and',chips], fun(X,Y,Z)->[X,Y,Z] end),
    
    [[2, 8], 10, [[], 6], 2]
        = evensonly_star([[9,1,2,8],3,10,[[9,9],7,6],2]),
    
    [38, 1920, [2, 8], 10, [[], 6], 2]
        = evensonly_starCO( [[9,1,2,8],3,10,[[9,9],7,6],2],
            fun (NewL,Product,Sum) -> [Sum | [Product | NewL]] end),
    
    test_ok.

% -------

remberF(_, _, []) -> [];
remberF(Test, A, [H|T]) -> 
    case Test(A,H) of
        true -> T;
        false -> [H | remberF(Test, A, T)]
    end.

% A "curried" function: returns another function with one argument pre-set
is_eqC(A) -> fun(X) -> X =:= A end.

% Woo hoo!! We can use pattern-matching in anonymous functions too!
remberFC(Test) -> 
    fun (_,[]) -> [];
        (A,[A|T]) -> T;
        (A,[H|T]) -> [H | (remberFC(Test)) (A,T)]
    end.

% Revisit some old friends - transform ala "remberF()"
insertL_F(Test) ->
    fun (_,_,[]) -> [];
        (New,Old,[Old|T]) -> [New | [Old | T]];
        (New,Old,[H|T])   -> [H | (insertL_F(Test)) (New,Old,T)]
    end.

% Oops, above didn't actually use Test(). Let's do it right (but ugler)
insertR_F(Test) ->
    fun (_,_,[]) -> [];
        (New,Old,[H|T]) -> 
            case Test(Old, H) of
                true -> [Old | [New | T]];
                false -> [H | (insertR_F(Test)) (New,Old,T)]
            end
     end.

% In above functions, only one line differs! Abstract this difference
sis_eqL(New, Old, L) -> [New | [Old | L]].

sis_eqR(New, Old, L) -> [Old | [New | L]].

% Now write a generic "insert" function
insertG(SisEq) -> 
    fun (_, _, []) -> [];
        (New,Old,[Old|T]) -> SisEq(New, Old, T);
        (New,Old,[H|T])   -> [H | (insertG(SisEq)) (New,Old,T)]
    end.

% We *could* create specific "inserter-functions" like this:
%   InsertL = insertG(fun sis_eqL/3).
%   InsertR = insertG(fun sis_eqR/3).
% But variables are not allowed in a module, so we have to type a bit more
insertL2(New,Old,L) -> (insertG(fun sis_eqL/3)) (New,Old,L).

insertR2(New,Old,L) -> (insertG(fun sis_eqR/3)) (New,Old,L).

% Helpers for subst(), rember()
sis_eqS(New, _, L) -> [New | L].
%sis_eqRem(_, _, L) -> L.

% Define subst(), rember() using our generic insert function
substG(New, Old, L) -> (insertG(fun sis_eqS/3)) (New,Old,L).
remberG(A, L)       -> (insertG(fun(_,_,X) -> X end)) (false,A,L).

% Rewrite value() from Chapter 6 using these techniques
atom_to_function('+') -> fun plus/2;
atom_to_function('-') -> fun minus/2;
atom_to_function('^') -> fun pow/2;
atom_to_function(X) -> throw({bad_operator, X}).

valueG(Aexp) when is_number(Aexp) -> Aexp;
valueG(Aexp) when is_atom(Aexp) -> none;
valueG(Aexp) -> 
    F = atom_to_function( operator(Aexp) ),
    F( valueG(first_subexp(Aexp)), valueG(second_subexp(Aexp)) ).

% Rewrite our old friend multirember() similarly
multiremberF(Test) -> 
    fun (_, []) -> [];
        (A, [H|T]) -> 
            case Test(A,H) of
                true  -> (multiremberF(Test)) (A,T);
                false -> [H | (multiremberF(Test)) (A,T)]
            end
    end.

% Use this to define a multirember instance
multirember_iseq(A, L) -> (multiremberF(fun is_eq/2)) (A,L).

% Neither the test nor atom (A) change while multirember is running.
% Let's try combining these into one helper function
is_eq_tuna(X) -> X =:= tuna.

% rewrite to accept a test like is_eq_tuna()
multiremberT(_, []) -> [];
multiremberT(Test, [H|T]) -> 
    case Test(H) of
        true -> multiremberT(Test, T);
        false -> [H | multiremberT(Test, T)]
    end.

% Rewrite to accept a *collector* (sometimes called a "continuation")
% 
% What this function really does:
% - traverses a list of atoms, looking for matches with A
% - calls Col() with two lists: (NewLat, Seen)
%   "NewLat":   the input list sans matches of A
%   "Seen":     list of matched atoms
%
% Instead of "cons"ing lists of return values, we "cons" arguments
% to the collector function, which is invoked when we reach our base case
%
multiremberCO(_, [], Col) -> 
    Col([], []);
multiremberCO(A, [A|T], Col) -> 
    multiremberCO(A, T, fun (NewLat,Seen) -> Col(NewLat, [A | Seen]) end);
multiremberCO(A, [H|T], Col) -> 
    multiremberCO(A, T, fun (NewLat,Seen) -> Col([H | NewLat], Seen) end).

% some collector functions to use with multiremberCO()
a_friend(_X, []) -> true;
a_friend(_, _) -> false.

last_friend(X, _) -> length(X).

% Combine two old friends, "multiinsertL" and "multiinsertR"
multiinsertLR(_, _, _, []) -> [];
multiinsertLR(New, OldL, OldR, [OldL | T]) -> 
    [New | [OldL | multiinsertLR(New, OldL, OldR, T)]];
multiinsertLR(New, OldL, OldR, [OldR | T]) -> 
    [OldR | [New | multiinsertLR(New, OldL, OldR, T)]];
multiinsertLR(New, OldL, OldR, [H|T]) -> 
    [H | multiinsertLR(New, OldL, OldR, T)].

% Now rework above to add a collector
multiinsertLRCO(_, _, _, [], Col) -> 
    Col([], 0, 0);
multiinsertLRCO(New, OldL, OldR, [OldL | T], Col) -> 
    multiinsertLRCO(New, OldL, OldR, T, 
        fun (NewLat,L,R) -> Col([New | [OldL | NewLat]], L+1, R) end);
multiinsertLRCO(New, OldL, OldR, [OldR | T], Col) -> 
    multiinsertLRCO(New, OldL, OldR, T,
        fun (NewLat,L,R) -> Col([OldR | [New | NewLat]], L, R+1) end);
multiinsertLRCO(New, OldL, OldR, [H|T], Col) -> 
    multiinsertLRCO(New, OldL, OldR, T,
        fun (NewLat,L,R) -> Col([H | NewLat], L, R) end).

% Remember "star" functions"? They operate on arbitary S-exps
% (empty, atom consed onto list, or list consed onto list)
% Let's try a collector in a star function, so we can process a "tree"

is_even(N) -> (N rem 2) == 0.

evensonly_star([]) -> [];
evensonly_star([H|T]) when is_number(H) -> 
    case is_even(H) of
        true -> [H | evensonly_star(T)];
        false -> evensonly_star(T)
    end;
evensonly_star([H|T]) -> 
    [evensonly_star(H) | evensonly_star(T)].

% Now using a collector. We doubly-recurse in the "is_list" case
evensonly_starCO([], Col) -> 
    Col([], 1, 0);
evensonly_starCO([H|T], Col) when is_number(H) -> 
    case is_even(H) of
        true -> evensonly_starCO(T, 
            fun(NewL,P,S) -> Col([H | NewL], H*P, S) end);
        false -> evensonly_starCO(T, 
            fun(NewL,P,S) -> Col(NewL, P, H+S) end)
    end;
evensonly_starCO([H|T], Col) when is_list(H) -> 
    evensonly_starCO(H, 
        fun (AL,AP,AS) -> evensonly_starCO(T, 
            fun (DL,DP,DS) -> Col([AL|DL], AP*DP, AS+DS) end) end
        );
evensonly_starCO([H|_], _) -> 
    erlang:error({bad_list_element, H}).


%% ---------------------------------------------------------------
%% Chapter 9: Partial Functions and the Y Combinator.
%% 
%% You can write recursive anonymous functions (lambdas) in Erlang
%% *without* Y combinator, but you must assign the lambda to a variable
%% before invoking.
%% 
%% Y Combinator allows *unnamed* anonymous lambdas, but adds complexity.
%% This chapter is important for theoretical understanding,
%% but it largely violates Joe Armstrong's admonition:
%%
%% "You must constantly strive to write 'mind boggling simple code'"
%%
%% (http://www.erlang.org/ml-archive/erlang-questions/200301/msg00062.html)
%% ---------------------------------------------------------------

chapter9_test() ->

    true = looking(caviar, [6,2,4,caviar,5,7,3]),
    false = looking(caviar, [6,2,grits,caviar,5,7,3]),
    
    % Do not run the following! It will never return
    %% looking(caviar, [7,1,2,caviar,5,6,3]),
    
    [a,[b,c]]     = shift( [[a,b],c] ),
    [a,[b,[c,d]]] = shift( [[a,b],[c,d]] ),
    
    io:format("# chapter8: need test for align()~n"),
    
    [a,[b,c]]     = shuffle([a,[b,c]]),
    [a,b]         = shuffle([a,b]),
    % shuffle() is NOT total... will loop forever on [[a,b],[c,d]]
    
    2 = ackermann(1, 0),
    3 = ackermann(1, 1),
    7 = ackermann(2, 2),
    
    % Computational cost of Ackermann() increases *very quickly*.
    % ackermann(4,3) is effectively uncomputable!
    
    % Example recursive lambda without using Y Combinator
    % We must assign the lambda to a variable before calling.
    Fac  = fun(0, _) -> 1; 
              (N, F) -> N * F(N-1, F) end,
    720 = Fac(6, Fac),
    
    % Using Y Combinator
    5   = len_Y([a,b,c,d,e]),
    true.

% -------

% Let's scan a list of atoms. If the value is a number, use it as a pointer
% to the next item. (Sort of like a CPU's program counter)
looking(A, Lat) ->
    keep_looking(A, pick(1, Lat), Lat).

% Example of "unnatural recursion. ("Sorn" is "symbol or number")
% Why "unnatural? It does not always get closer to its goal when recurring
%
% Functions that do not return for some values (like this one)
% are called "partial functions"
%
keep_looking(A, Sorn, Lat) when is_number(Sorn) ->
    keep_looking(A, pick(Sorn, Lat), Lat);
keep_looking(A, Sorn, _) -> 
    Sorn == A.

% This is the *most unnatural* recursion possible,
% and the *most partial* function
eternity(X) -> eternity(X).

% For below... 
% Use functions from chapter 7
shift(Pair) ->
    build(first(first(Pair)), 
        build(second(first(Pair)),
            second(Pair))).

% This is actually a "total" function: it always approaches its goal
% with each recurson step, even though it's hard to tell that's true.
%
% This violates the "Seventh Commandment": recur on subparts of the same nature
%
% (PorA: "pair or atom")
align(PorA) when is_atom(PorA) or is_number(PorA) -> PorA;
align(PorA) -> 
    case is_a_pair(first(PorA)) of
        true -> align( shift(PorA) );
        false -> build( first(PorA), align(second(PorA)) )
    end.

% This function is NOT total: will never yield a result
% if input is [[a,b],[c,d]]
shuffle(PorA) when is_atom(PorA) or is_number(PorA) -> PorA;
shuffle(PorA) -> 
    case is_a_pair(first(PorA)) of
        true -> shuffle( revpair(PorA) );
        false -> build( first(PorA), shuffle(second(PorA)) )
    end.

% Ackermann's function IS total, but unnatural
% (args do not necessarily decrease for the recursion)
ackermann(0, M) -> M+1;
ackermann(N, 0) -> ackermann(N-1, 1);
ackermann(N, M) -> ackermann(N-1, ackermann(N, M-1)).

%%%%%%%%% FIXME ADD MORE %%%%%%%%%%%

%% Y-combinator: simplest version, ONLY works with arity=1 lambdas
%% http://www.erlang.org/ml-archive/erlang-questions/200301/msg00053.html
y_arity1(X) ->
  F = fun (P) -> X(fun (Arg) -> (P(P))(Arg) end) end,
  F(F).

%% Y-Combinator: multi-arity version by Atilla Babo, for arity=0..5 lambdas
%% http://www.erlang.org/pipermail/erlang-questions/2008-February/032673.html
y(F) when is_function(F) ->
	{arity, Arity} = erlang:fun_info(F(F), arity),
	G = case Arity of
		0 -> fun(H) -> F(fun() -> (H(H))() end) end;
		1 -> fun(H) -> F(fun(A) -> (H(H))(A) end) end;
		2 -> fun(H) -> F(fun(A, B) -> (H(H))(A, B) end) end;
		3 -> fun(H) -> F(fun(A, B, C) -> (H(H))(A, B, C) end) end;
		4 -> fun(H) -> F(fun(A, B, C, D) -> (H(H))(A, B, C, D) end) end;
		5 -> fun(H) -> F(fun(A, B, C, D, E) -> (H(H))(A, B, C, D, E) end) end;
        true -> erlang:error(unsupported_arity, Arity)
	end,
	G(G);
y(_) -> erlang:error(not_function).

%% A dynamically-generated, arbitrary-arity version was created by John Webb.
%% http://www.erlang.org/pipermail/erlang-questions/2008-February/032690.html
%% ... Too complex to include here, but good food for thought!

%% Examples: "dots()" runs forever, don't try it...
% dots() ->
%	 (y(fun(F) -> fun() -> io:fwrite("."), F() end end))().

len_Y(I) when is_list(I) ->
	(y(fun(F) -> fun(X, N) -> case X of [] -> N; [_|R] -> F(R, N+1) end end end))(I, 0).


%% ---------------------------------------------------------------
%% Chapter 10: Build a "Little Scheme" interpreter
%% ---------------------------------------------------------------

chapter10_test() ->
	erlang:error(not_implemented).

% -------

