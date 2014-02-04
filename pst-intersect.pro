/tx@IntersectDict 200 dict def
tx@IntersectDict begin

/VecAdd {
    3 -1 roll add 3 1 roll add exch
} bind def
/VecSub {
    neg 3 -1 roll add 3 1 roll neg add exch
} bind def
/VecScale {
  dup 4 -1 roll mul 3 1 roll mul
} bind def
/ToVec {
    [ 3 1 roll ]
} bind def

/MaxPrecision 1e-6 def
/H1Interval [0 0.5] def
/H2Interval [0.5 MaxPrecision add 1] def
/Epsilon 1e-4 def
/MinClippedSizeThreshold 0.8 def

/CurveToString {
    (CurveToString) DebugBegin
    aload pop ([) 3 -1 roll 20 string cvs strcat (, ) strcat exch 20 string cvs strcat (]) strcat
    DebugEnd
} bind def
%
/IntersectBeziers {
    2 copy [0 1] [0 1] MaxPrecision IterateIntersection
    exch TArray dup /lt exch quicksort
    exch TArray dup /lt exch quicksort
    3 -1 roll exch LoadIntersectionPoints 
    3 1 roll LoadIntersectionPoints
} bind def

%
/TArray {
    [ exch
    { dup type /nulltype eq { pop exit } if
	aload pop add 0.5 mul
    } forall ] % -> array of tvalues
} bind def
%
% Curve t
/LoadIntersectionPoints {
    % prepare Curve for use with tx@Func
    exch [ exch { aload pop } forall ]
    tx@Dict begin tx@FuncDict begin 2 dict begin
	dup length 2 idiv 1 sub /BezierType exch def /Points exch def
	[ exch {
	    GetBezierCoor
	} forall ]
    end end end
} bind def
%
% Iteration procedure to compute all intersections of CurveA and CurveB.
% This contains the
%
% [CurveA] [CurveB] [intervalA] [intervalB] precision -> domsA domsB
/IterateIntersection {
    (IterateIntersection) DebugBegin
    11 dict begin
	/precision exch def
	% in order to limit recursion
        /counter 0 def
	/depth 0 def
	/domsA InitInterval
	/domsB InitInterval
	/domsA /domsB 6 2 roll _IterateIntersection
	domsB domsA
    end
    DebugEnd
} bind def
%
% This is the iteration part which is called recursively.
% Before this one must push domA domB CurveA CurveB iter on the stack, so that
% after the recursive call, the previous state of these five variables can be restored.
%
% /domsA /domsB CurveA CurveB domA domB
/_IterateIntersection {
    (_IterateIntersection) DebugBegin
    CloneVec /domB exch def
    CloneVec /domA exch def
    CloneCurve /CurveB exch def
    CloneCurve /CurveA exch def
    /iter 0 def
    /depth depth 1 add def
    /dom null def
    /counter counter 1 add def

    CheckIT {
	(>> curve subdivision performed: dom(A) = ) domA CurveToString strcat
	(, dom(B) = ) strcat domB CurveToString strcat ( <<) strcat ==
    } if
    CurveA IsConstant CurveB IsConstant and {
	CurveA MiddlePoint ToVec
	CurveB MiddlePoint ToVec AreNear {
	    domA domB 4 -1 roll exch PutInterval PutInterval
	} {
	    pop pop
	} ifelse
    }{
	counter 100 lt {
	    % use a loop to simulate some kind of return to exit at different
	    % positions.
	    {
		/iter iter 1 add def
		iter 100 lt
		domA Extent precision ge
		domB Extent precision ge or and not {
		    iter 100 ge {
			false 
		    } {
			CurveA MiddlePoint ToVec
			CurveB MiddlePoint ToVec AreNear {
			    domA domB true
			}{
			    false
			} ifelse
		    } ifelse
		    exit
		} if
		% iter < 100 && (dompA.extent() >= precision || dompB.extent() >= precision)
		CheckIT {
		    (counter: ) counter 20 string cvs strcat
		    (, iter: ) iter 20 string cvs strcat strcat
		    (, depth: ) depth 20 string cvs strcat strcat ==
		} if
	
		CurveA CurveB ClipCurve /dom exch def
	
		CheckIT {(dom : ) dom CurveToString strcat == } if		
		dom IsEmptyInterval {
		    CheckIT { (empty interval, exit) == } if
		    false exit
		} if
		% dom[0] > dom[1], invalid. How to handle?	
		dom aload pop 2 copy min 3 1 roll max gt {
		    CheckIT {
			(dom[0] > dom[1], invalid!) ==
		    } if
		    false exit
		} if

		domB dom MapTo /domB exch def
		CurveB dom Portion

		CurveB IsConstant CurveA IsConstant and {
		    CheckIT {
          		(both curves are constant: ) ==	
			(C1: [ ) CurveA { CurveToString ( ) strcat strcat } forall (]) strcat ==
			(C2: [ ) CurveB { CurveToString ( ) strcat strcat } forall (]) strcat ==
		    } if
		    CurveA MiddlePoint ToVec
		    CurveB MiddlePoint ToVec AreNear {
			domA domB true
		    } {
			false
		    } ifelse
		    exit
		} if
		% if we have clipped less than 20%, we need to subdivide the
		% curve with the largest domain into two sub-curves.
		dom Extent MinClippedSizeThreshold gt {
		    CheckIT {
			(clipped less than 20% : ) ==
			(angle(A) = ) CurveA dup length 1 sub get aload pop
				      CurveA 0 get aload pop VecSub
   				      exch 2 copy 0 eq exch 0 eq and {
					  pop pop (NaN)
				      } {
					  atan 20 string cvs
				      } ifelse strcat ==
		        (angle(B) = ) CurveB dup length 1 sub get aload pop
		                      CurveB 0 get aload pop VecSub
				      exch 2 copy 0 eq exch 0 eq and {
					  pop pop (NaN)
				      } {
					  atan 20 string cvs
				      } ifelse strcat ==
		        (dom : ) == dom == (domB :) == domB ==
		    } if

		    CurveA CurveB domA domB iter% leave those values on the stack to revert to this state after the recursive calls
     		    7 -2 roll 2 copy 9 2 roll 2 copy % /domsA /domsB CurveA CurveB domA domB iter /domsA /domsB /domsA /domsB
	    
		    domA Extent domB Extent gt {
			CurveA CloneCurve dup H1Interval Portion % pC1
			CurveA CloneCurve dup H2Interval Portion % pC2
			domA H1Interval MapTo                    % dompC1
			domA H2Interval MapTo                    % dompC2
			% need: /domsA /domsB pC2 CurveB dompC2 domB   /domsA /domsB pC1 CurveB dompC1 domB
			3 -1 roll exch % /domsA /domsB /domsA /domsB pC1 dompC1 pC2 dompC2
			CurveB exch domB 8 4 roll % /domsA /domsB pC2 CurveB dompC2 domB /domsA /domsB pC1 dompC1
			CurveB exch domB % /domsA /domsB pC2 CurveB dompC2 domB /domsA /domsB pC1 CurveB dompC1 domB
		    } {
			CurveB CloneCurve dup H1Interval Portion % pC1
			CurveB CloneCurve dup H2Interval Portion % pC2
			domB H1Interval MapTo                    % dompC1
			domB H2Interval MapTo                    % dompC2
			% need: /domsB /domsA pC2 CurveA dompC2 domA   /domsB /domsA pC1 CurveA dompC1 domA
			8 -2 roll exch 8 2 roll 6 -2 roll exch 6 2 roll % /domsB /domsA /domsB /domsA pC1 pC2 dompC1 dompC2
			3 -1 roll exch % /domsB /domsA /domsB /domsA pC1 dompC1 pC2 dompC2
			CurveA exch domA 8 4 roll % /domsB /domsA pC2 CurveA dompC2 domA /domsB /domsA pC1 dompC1
			CurveA exch domA          % /domsB /domsA pC2 CurveA dompC2 domA /domsB /domsA pC1 CurveA dompC1 domA
		    } ifelse

		    _IterateIntersection
		    _IterateIntersection
		    
		    % restore the state before the recursive calls.
		    /iter exch def
		    /domB exch def
		    /domA exch def
		    /CurveB exch def
		    /CurveA exch def
		    false exit
		} if
		CurveA CurveB /CurveA exch def /CurveB exch def
		domA domB /domA exch def /domB exch def
		% exchange /domsA and /domsB on the stack!
		exch
	    } loop	

	    % boolean on stack
	    {
		4 -1 roll exch PutInterval PutInterval
		CheckIT {
		    (found an intersection ============================) ==
		} if
	    } { pop pop } ifelse
	} {
	    pop pop
	} ifelse
    } ifelse
    /depth depth 1 sub def
    DebugEnd
} bind def
%
% Add a new interval [newinterval] to the array stored in /Intervals
% /Intervals [newinterval] PutInterval
/PutInterval {
    % create new array to be pushed into /Intervals
    CloneVec exch
    dup load 3 -1 roll exch
    dup dup length dup 3 1 roll 1 sub get % Arr L Arr[L-1]
    2 copy 2 add eq { % array is full, create a new one which is twice as long
	% Arr L Arr[L-1]
	exch 2 mul dup % Arr Arr[L-1] 2L 2L
	array dup 0 6 -1 roll putinterval % Arr[L-1] 2L Arr2
	exch 3 -1 roll % Arr2 2L Arr[L-1]
	3 copy exch 1 sub exch put % Arr2 2L Arr[L-1]
    } if
    % [newinterval] Arr L Arr[L-1]
    1 add 3 copy % [newinterval] Arr L Arr[L-1]+1 Arr L Arr[L-1]+1
    exch 1 sub exch put % [newinterval] Arr L Arr[L-1]+1
    exch pop exch dup 4 2 roll exch put
    def
} bind def
%
% /IntervalName InitInterval
/InitInterval {
    10 array dup 9 -1 put def
} bind def
%
% CheckIT if an interval is empty, which is represented by a [1 0] interval.
%
% [interval] -> empty?
/IsEmptyInterval {
    aload pop 0 eq exch 1 eq and
} bind def
%
% Does a deep copy of the array [Curve].
%
% [Curve] CloneCurve -> [newCurve]
/CloneCurve {
    [ exch {
	CloneVec
    } forall ]
} bind def
%
% Does a deep copy of the vector [X Y]
%
% [X Y] CloneVec -> [Xnew Ynew]
/CloneVec {
    aload pop ToVec
} bind def
%
% Map the sub-interval I in [0,1] into the interval J. Returns a new array.
%
% [J] [I] MapTo -> [Jnew]
/MapTo {
    (MapTo) DebugBegin
    exch aload 0 get 3 1 roll exch sub 2 copy % [I] J0 Jextent J0 Jextent
    5 -1 roll aload aload pop % J0 Jextent J0 Jextent I0 I1 I0 I1
    min 4 -1 roll mul % J0 Jextent J0 I0 I1 min(I0,I1)*Jextent
    4 -1 roll add [ exch % J0 Jextent I0 I1 [ J0new
    6 2 roll max mul add ]
    DebugEnd
} bind def
%
% Compute the portion of the Bezier curve "B" wrt the interval "I"
%
% [CurveB] [I] Portion
/Portion {
    (Portion) DebugBegin
    dup Min 0 eq { % [CurveB] [I]
	% I.min() == 0
	Max dup 1 eq {% [CurveB] I.max()
	    % I.max() == 1
	    pop pop	    
	} { % [CurveB] I.max()
	    LeftPortion
	} ifelse
    } { % [CurveB] [I]
	2 copy Min % [CurveB] [I] [CurveB] I.min()
	RightPortion
	dup Max 1 eq {
	    % I.max() == 1
	    pop pop
	} {% [CurveB] [I]
	    dup aload pop exch sub 1 3 -1 roll Min sub div % [CurveB] (I1-I0)/(1-I.min())
	    LeftPortion
	} ifelse
    } ifelse
    DebugEnd
} bind def
%
% Compute the portion of the Bezier curve "B" wrt the interval [0,t].
%
% [CurveB] t LeftPortion
/LeftPortion {
    (LeftPortion) DebugBegin
    exch dup length 1 sub dup 4 1 roll % L-1 t [CurveB] L-1
    1 1 3 -1 roll { % L-1 t [CurveB] i
	4 -1 roll dup 5 1 roll % L-1 t [CurveB] i L-1
	-1 3 -1 roll % L-1 t [CurveB] L-1 -1 i
	{ % L-1 t [CurveB] j
	    2 copy 5 copy % L-1 t [CurveB] j [CurveB] j t [CurveB] j [CurveB] j 
	    1 sub get 3 1 roll get % L-1 t [CurveB] j [CurveB] j t B[j-1] B[j]
	    Lerp put pop % L-1 t [CurveB]
	} for
    } for
    pop pop pop
    DebugEnd
} bind def
%
% Compute the portion of the Bezier curve "B" wrt the interval [t,1].
%
% [CurveB] t RightPortion
/RightPortion {
    (RightPortion) DebugBegin
    exch dup length 1 sub dup 4 1 roll % L-1 t [CurveB] L-1
    1 1 3 -1 roll {% L-1 t [CurveB] i
	4 -1 roll dup 5 1 roll % L-1 t [CurveB] i L-1
	exch sub 0 1 3 -1 roll  % L-1 t [CurveB] 0 1 L-i-1
	{% L-1 t [CurveB] j
	    2 copy 5 copy
	    get 3 1 roll 1 add get Lerp put pop
	} for
    } for
    pop pop pop
    DebugEnd
} bind def
%
% Given two points and a parameter t \in [0, 1], return a point
% proportionally from A to B by t. Akin to 1 degree bezier.
%
% t [A] [B] Lerp -> [newpoint]
/Lerp {
    (Lerp) DebugBegin
    3 -1 roll dup 1 exch sub 3 1 roll % [A] (1-t) [B] t
    exch aload pop 3 -1 roll VecScale % [A] (1-t) B.x*t B.y*t
    4 2 roll
    exch aload pop 3 -1 roll VecScale VecAdd ToVec % [A.x*(1-t)+B.x*t A.y*(1-t)+B.y*t]
    DebugEnd
} bind def
%
% Test if all points of a curve are near to each other.
%
% [Curve]
/IsConstant {
    aload length [ exch 1 roll ] true 3 1 roll % true [P0] [[P1] ... [Pn]]
    {
	exch dup 4 1 roll % [P0] near? [P0] [Pi]
	AreNear and exch
    } forall
    pop
} bind def
%
% [P1] [P2] AreNear -> bool
/AreNear {
    (AreNear) DebugBegin
    aload pop 3 -1 roll aload pop
    4 copy abs 3 { exch abs max } repeat Epsilon mul
    dup 6 2 roll VecSub abs 4 -1 roll lt exch abs 3 -1 roll lt and
    DebugEnd
} bind def
%
% [P] Min -> min(P.x, P.y)
/Min {
    aload pop min
} bind def
% [P] Max -> max(P.x, P.y)
/Max {
    aload pop max
} bind def
%
% [P] Extent -> P1 - P0
/Extent {
    aload pop exch sub
} bind def
%
% Compute the middle point of the first and last point of [Curve]
% [Curve] -> X Y
/MiddlePoint {
    dup dup length 1 sub get aload pop % [Curve] XN YN
    3 -1 roll 0 get aload pop
    VecAdd 0.5 VecScale
} bind def
%
% MiddlePointA [CurveB] -> A B C
/OrthogonalOrientationLine {
    (OrthogonalOrientationLine) DebugBegin
    dup dup length 1 sub get aload pop 3 -1 roll 0 get aload pop VecSub
    neg exch % rotate by +90 degrees
    4 2 roll 2 copy 6 2 roll VecAdd % MiddlePointA CalcPoint
    ImplicitLine
    DebugEnd
} bind def
%
% Pick an orientation line for a Bezier curve. This uses the first pair of coordinates which are not near.
%
% [Curve] -> A B C
/PickOrientationLine {
    (PickOrientationLine) DebugBegin
    dup dup length 1 sub exch 0 get% [Curve] L-1 P0
    exch -1 1 {% [Curve] P0 i
	3 -1 roll dup 4 1 roll exch get % [Curve] P0 Pi
	2 copy AreNear {
	    pop
	} {
	    exit
	} ifelse
    } for
    3 -1 roll pop
    exch aload pop 3 -1 roll aload pop ImplicitLine
    DebugEnd
} bind def
%
%
% Compute the coefficients A, B, C of the normalized implicit equation
% of the line which goes through the points (X1, Y1) and (X2, Y2).
% (equivalent to orientation_line).
%
% Xi Yi Xj Yj ImplicitLine -> A B C
/ImplicitLine {
    4 copy % Xi Yi Xj Yj Xi Yi Xj Yj
    3 -1 roll sub 7 1 roll sub 5 1 roll % Yj-Yi Xi-Xj Xi Yi Xj Yj
    % Yi*Xj - Xi*Yj
    4 -1 roll mul neg % Yj-Yi Xi-Xj Yi Xj -Yj*Xi
    3 1 roll mul add % Yj-Yi Xi-Xj Yi*Xj-Yj*Xi | l0 l1 l2
    3 1 roll 2 copy tx@Dict begin Pyth end dup dup % l2 l0 l1 L L L
    5 -1 roll exch % l2 l1 L L l0 L
    div 5 1 roll % l0/L l2 l1 L L
    3 1 roll div % l0/L l2 L l1/L
    3 1 roll div % l0/L l1/L l2/L
} bind def
%
% Compute the distance of point (X, Y) from the implicit line given
% by A*x + B*y + C = 0, (A²+B² = 1)
%
% X Y A B C
/distance {
    5 1 roll 3 -1 roll mul 3 1 roll mul add add
} bind def
%
% convert [A.x A.y ... N.x N.y] to [[A.x A.y] ... [N.x N.y]]
/ArrayToPointArray {
    aload length dup 2 idiv {
	3 1 roll [ 3 1 roll ] exch % A.x A.y ... [N.x N.y] L
	dup 1 sub 3 1 roll 1 roll
    } repeat 1 add [ exch 1 roll ]
} bind def
%
% convert [[A.x A.y] ... [N.x N.y]] to [A.x A.y ... N.x N.y]
/PointArrayToArray {
    aload length dup {
	1 add dup 3 -1 roll aload pop 4 -1 roll 1 add 2 roll
    } repeat 1 add [ exch 1 roll ]
} bind def
%
% Clip the Bezier curve B with respect to the Bezier curve A for
% individuating intersection points. The new parameter interval for the
% clipped curve is pushed on the stack.
%
% [CurveA] [CurveB] ClipCurve -> [newinterval]
/ClipCurve {
    (ClipCurve) DebugBegin
    4 dict begin 
    /CurveB exch def /CurveA exch def
    CurveA IsConstant {
    	CurveA MiddlePoint CurveB OrthogonalOrientationLine
    } {
	CurveA PickOrientationLine
    } ifelse
    CheckIT {
	3 copy exch 3 -1 roll (OrientationLine : )
	3 { exch 20 string cvs ( ) strcat strcat } repeat ==
    } if
    CurveA FatLineBounds
    CheckIT { dup (FatLineBounds : ) exch aload pop exch 20 string cvs (, ) strcat exch 20 string cvs strcat strcat == } if
    CurveB ClipCurveInterval
    end % end local dictionary
    DebugEnd
} bind def
%
% A B C [Curve] -> A B C [dmin dmax]
/FatLineBounds {
    (FatLineBounds) DebugBegin
    /dmin 0 def /dmax 0 def
    { 
	4 copy aload pop 5 2 roll distance
	dup dmin lt { dup /dmin exch def } if
	dup dmax gt { dup /dmax exch def } if
	pop pop
    } forall
    [dmin dmax]
    DebugEnd
} bind def
%
% Clip the Bezier curve wrt the fat line defined by the orientation
% line (given by A, B, C) and the interval range "bound". The new parameter interval
% for the clipped curve is pushed on the stack.
%
% A B C [bound] [curve] -> [newinterval]
/ClipCurveInterval {
    (ClipCurveInterval) DebugBegin
    15 dict begin
    /curve exch def
    aload pop 2 copy min /boundMin exch def max /boundMax exch def
    [ 4 1 roll ] cvx /fatline exch def
    % number of sub-intervals
    /n curve length 1 sub def
    % distance curve control points
    /D n 1 add array def
    0 1 n { % i
	dup curve exch get aload pop % i Pi.x Pi.y
	fatline distance % distance d of Point i from the orientation line, on stack; i d
	exch dup n div % d i i/n
	[ exch 4 -1 roll ] % i [ i/n d ]
	D 3 1 roll put 
    } for
    D ConvexHull /D exch def
    % get the x-coordinate of the i-th point, i getX -> D[i][X]
    /getX { D exch get 0 get } def
    % get the y-coordinate of the i-th point, i getY -> D[i][Y]
    /getY { D exch get 1 get } def
    /tmin 1 def /tmax 0 def
    0 getY dup
    boundMin lt /plower exch def
    boundMax gt /phigher exch def
    plower phigher or not {
	% inside the fat line
	tmin 0 getX gt { /tmin 0 getX def } if
	tmax 0 getX lt { /tmax 0 getX def } if	
    } if
    1 1 D length 1 sub {
	/i exch def
	/clower i getY boundMin lt def
	/chigher i getY boundMax gt def
	clower chigher or not {
	    % inside the fat line
	    tmin i getX gt { /tmin i getX def } if
	    tmax i getX lt { /tmax i getX def } if
	} if
	clower plower eq not {
	    % cross the lower bound
	    boundMin i 1 sub i D Intersect % t on stack
	    dup tmin lt { dup /tmin exch def } if
	    dup tmax gt { dup /tmax exch def } if
	    pop 
	    /plower clower def
	} if
	chigher phigher eq not {
	    % cross the upper bound
	    boundMax i 1 sub i D Intersect
	    dup tmin lt { dup /tmin exch def } if
	    dup tmax gt { dup /tmax exch def } if
	    pop 
	    /phigher chigher def
	} if
    } for
    % we have to test the closing segment for intersection
    /i D length 1 sub def
    /clower 0 getY boundMin lt def
    /chigher 0 getY boundMax gt def
    clower plower eq not {
	% cross the lower bound
	boundMin i 0 D Intersect
	dup tmin lt { dup /tmin exch def } if
	dup tmax gt { dup /tmax exch def } if
	pop
    } if
    chigher phigher eq not {
	% cross the lower bound
	boundMax i 0 D Intersect
	dup tmin lt { dup /tmin exch def } if
	dup tmax gt { dup /tmax exch def } if
	pop
    } if
    [tmin tmax]
    end % end of local dictionary
    DebugEnd
} bind def
%
% Get the x component of the intersection point between the line passing
% through points (Xi, Yi) and (Xj, Yj) and the horizonal line Y = "y"
%
% y i j [Curve] -> Xisect
/Intersect {
    % load the coordinates of the i-th and j-th point
    dup 4 -1 roll get aload pop % y j [Curve] Xi Yi
    4 2 roll exch get aload pop % y Xi Yi Xj Yj
    % (Xj - Xi) * (y - Yi)/(Yj - Yi) + Xi
    % We are sure, that Yi != Yj, because this procedure is called only
    % when the lower or upper bound is crossed.
    4 2 roll 2 copy 6 2 roll VecSub % y Xi Yi (Xj-Xi) (Yj-Yi)
    5 2 roll % (Xj-Xi) (Yj-Yi) y Xi Yi
    neg 3 -1 roll add % (Xj-Xi) (Yj-Yi) Xi (y - Yi)
    3 -1 roll div % (Xj-Xi) Xi (y-Yi)/(Yj-Yi)
    3 -1 roll mul add
} bind def

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPLEX HULL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Graham Scal algorithm to compute the convex hull of a set of
% points. Code written by Bill Casselman,
% http://www.math.ubc.ca/~cass/graphics/text/www/
%
% [[X1 Y1] [X2 Y2] ... [Xn Yn]] hull -> [[...] ... [...]]
%
/hulldict 32 dict def

hulldict begin

% u - v 

/vsub { 2 dict begin
/v exch def
/u exch def
[ 
  u 0 get v 0 get sub
  u 1 get v 1 get sub
]
end } def

% u - v rotated 90 degrees

/vperp { 2 dict begin
/v exch def
/u exch def
[ 
  v 1 get u 1 get sub
  u 0 get v 0 get sub
]
end } def

/dot { 2 dict begin
/v exch def
/u exch def
  v 0 get u 0 get mul
  v 1 get u 1 get mul
  add
end } def 

% P Q
% tests whether P < Q in lexicographic order
% i.e xP < xQ, or yP < yQ if xP = yP

/comp { 2 dict begin
/Q exch def
/P exch def
P 0 get Q 0 get lt 
  P 0 get Q 0 get eq
  P 1 get Q 1 get lt 
  and 
or 
end } def

end

% args: an arrya of points C
% effect: returns the array of points on the boundary of
%     the convex hull of C, in clockwise order 

/ConvexHull {
(ConvexHull) DebugBegin
hulldict begin
/C exch def
/comp C quicksort
/n C length def
% Q might circle around to the start
/Q n 1 add array def
Q 0 C 0 get put
Q 1 C 1 get put
/i 2 def
/k 2 def
% i is next point in C to be looked at
% k is next point in Q to be added
% [ Q[0] Q[1] ... ]
% scan the points to make the top hull
n 2 sub {
  % P is the current point at right
  /P C i get def
  /i i 1 add def
  {
    % if k = 1 then just add P 
    k 2 lt { exit } if
    % now k is 2 or more
    % look at Q[k-2] Q[k-1] P: a left turn (or in a line)?
    % yes if (P - Q[k-1])*(Q[k-1] - Q[k-2])^perp >= 0
    P Q k 1 sub get vsub 
    Q k 1 sub get Q k 2 sub get vperp 
    dot 0 lt {
      % not a left turn
      exit
    } if
    /k k 1 sub def
  } loop
  Q k P put
  /k k 1 add def
} repeat

% done with top half
% K is where the right hand point is
/K k 1 sub def

/i n 2 sub def
Q k C i get put
/i i 1 sub def
/k k 1 add def
n 2 sub {
  % P is the current point at right
  /P C i get def
  /i i 1 sub def
  {
    % in this pass k is always 2 or more
    k K 2 add lt { exit } if
    % look at Q[k-2] Q[k-1] P: a left turn (or in a line)?
    % yes if (P - Q[k-1])*(Q[k-1] - Q[k-2])^perp >= 0
    P Q k 1 sub get vsub 
    Q k 1 sub get Q k 2 sub get vperp 
    dot 0 lt {
      % not a left turn
      exit
    } if
    /k k 1 sub def
  } loop
  Q k P put
  /k k 1 add def
} repeat

% strip Q down to [ Q[0] Q[1] ... Q[k-2] ]
% excluding the doubled initial point
[ 0 1 k 2 sub {
  Q exch get
} for ] 
end
DebugEnd
} def

/qsortdict 8 dict def

qsortdict begin

% args: /comp a L R x
% effect: effects a partition into two pieces [L j] [i R]
%     leaves i j on stack

/partition { 8 dict begin
/x exch def
/j exch def
/i exch def
/a exch def
load /comp exch def
{
  {
    a i get x comp exec not {
      exit
    } if
    /i i 1 add def
  } loop
  {
    x a j get comp exec not {
      exit
    } if
    /j j 1 sub def
  } loop
  
  i j le {
    % swap a[i] a[j]
    a j a i get
    a i a j get 
    put put
    /i i 1 add def
    /j j 1 sub def
  } if
  i j gt {
    exit
  } if
} loop
i j
end } def

% args: /comp a L R
% effect: sorts a[L .. R] according to comp

/subsort {
% /c a L R
[ 3 1 roll ] 3 copy
% /c a [L R] /c a [L R]
aload aload pop 
% /c a [L R] /c a L R L R
add 2 idiv
% /c a [L R] /c a L R (L+R)/2
3 index exch get
% /c a [L R] /c a L R x
partition
% /c a [L R] i j
% if j > L subsort(a, L, j)
dup 
% /c a [L R] i j j
3 index 0 get gt {
  % /c a [L R] i j
  5 copy 
  % /c a [L R] i j /c a [L R] i j
  exch pop
  % /c a [L R] i j /c a [L R] j
  exch 0 get exch
  % ... /c a L j 
  subsort
} if
% /c a [L R] i j
pop dup
% /c a [L R] i i
% if i < R subsort(a, i, R)
2 index 1 get lt {
  % /c a [L R] i
  exch 1 get 
  % /c a i R
  subsort
}{
  4 { pop } repeat
} ifelse
} def

end % qsortdict

% args: /comp a
% effect: sorts the array a 
% comp returns truth of x < y for entries in a

/quicksort { qsortdict begin
dup length 1 gt {
% /comp a
dup 
% /comp a a 
length 1 sub 
% /comp a n-1
0 exch subsort
} {
pop pop
} ifelse
end } def

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEBUGGING STUFF %%%%%%%%%%%%%%%%%
/debug {
    dup 1 add copy {==} repeat pop
} bind def
/DebugIT false def
/CheckIT false def
/DebugDepth 0 def
/DebugBegin {
  DebugIT {
    /DebugProcName exch def
    DebugDepth 2 mul string
    0 1 DebugDepth 2 mul 1 sub {
      dup 2 mod 0 eq { (|) }{( )} ifelse
      3 -1 roll dup 4 2 roll
      putinterval
    } for
    DebugProcName strcat ==
    /DebugDepth DebugDepth 1 add def
  }{
    pop
  } ifelse
} bind def
/DebugEnd {
  DebugIT {
    /DebugDepth DebugDepth 1 sub def
    DebugDepth 2 mul 2 add string
    0 1 DebugDepth 2 mul 1 sub {
      dup 2 mod 0 eq { (|) }{ ( ) } ifelse
      3 -1 roll dup 4 2 roll
      putinterval
    } for
    dup DebugDepth 2 mul (+-) putinterval
    ( done) strcat ==
  } if
} bind def
/strcat {
    exch 2 copy
    length exch length add
    string dup dup 5 2 roll
    copy length exch
    putinterval
} bind def
/ShowCurve {
    { aload pop } forall
    8 -2 roll moveto curveto
} bind def
% A B C
/ShowLine {
    exch neg dup 3 1 roll div dup 0 exch moveto 3 1 roll div 200 mul add 200 exch lineto
} bind def
end % tx@IntersectDict
