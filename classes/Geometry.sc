//information about cross product and normal vector:
//"http://ef.gy/linear-algebra:normal-vectors-in-higher-dimensional-spaces"
/*
Line2D{
	var <>p1,<>p2,<>seg;
	*new { arg p1,p2,seg=false;
		^super.newCopyArgs(p1,p2,seg);
	}
	lineIntersection2D{arg  aLine;
		var y43 = aLine.p2.y - aLine.p1.y;
		var x21 = p2.x - p1.x;
		var x43 = aLine.p2.x - aLine.p1.x;
		var y21 = p2.y - p1.y;
		var denom = (y43 * x21) - (x43 * y21);

        var y13 = p1.y - aLine.p1.y;
		var x13 = p1.x - aLine.p1.x;

		var ua = ((x43 * y13) - (y43 * x13)) / denom;

		var x = p1.x + (ua * x21);
		var y = p1.y + (ua * y21);

		aLine.seg=aLine.seg;

		if (denom == 0, {
			^false;
		});
		if (this.seg && ((ua < 0) || (ua > 1)), {
			^false;
		});
		if (aLine.seg, {
			var ub = ((x21 * y13) - (y21 * x13)) / denom;
			if ((ub < 0) || (ub > 1), {
				^false;
			})
		});

		^(x@y)
	}

	projectPoint2Line{arg aPoint;
		^(p1+(p2-p1).asRealVector.proj((aPoint-p1).asRealVector2D));
	}
	closestPointOnLine{arg aPoint;
		^(p1+(p2-p1).asRealVector.proj((aPoint-p1).asRealVector2D));
	}

	closestPointOnSegment{arg aPoint;
		var v1x=(p2-p1).asRealVector.proj((aPoint-p1).asRealVector2D);
		var aLine=Line2D(aPoint,p1+v1x);
		^this.lineIntersection2D(aLine);
	}


	//deprecated
	*lineIntersection2D{arg  p1,  p2, p3,  p4;
		var seg1=false;//if line 1 is a segment.
		var seg2=false;
		var y43 = p4.y - p3.y;
		var x21 = p2.x - p1.x;
		var x43 = p4.x - p3.x;
		var y21 = p2.y - p1.y;
		var denom = (y43 * x21) - (x43 * y21);
        var y13 = p1.y - p3.y;
		var x13 = p1.x - p3.x;
		var ua = ((x43 * y13) - (y43 * x13)) / denom;
		var x = p1.x + (ua * x21);
		var y = p1.y + (ua * y21);
		if (denom == 0, {
			^false;
		});
		if (seg1 && ((ua < 0) || (ua > 1)), {
				^false;
		});
		if (seg2, {
			var ub = ((x21 * y13) - (y21 * x13)) / denom;
			if ((ub < 0) || (ub > 1), {
				^false;
			})
		});
		^(x@y)
	}
}
*/

LineND{
	var <>positionVector,<>directionVector,<>segment;
	*new { arg positionVector,directionVector,segment=false;
		^super.newCopyArgs(positionVector,directionVector,segment);
	}
	*newUsingPoints{arg positionVector1,positionVector2,segment=false;
		^super.newCopyArgs(positionVector1,(positionVector2-positionVector1),segment);
	}

	projectionScalar{arg aPoint;
		^directionVector.projectionScalar(aPoint.asRealVector-positionVector)
	}
	projectPoint2Line{arg aPoint;
		^positionVector+directionVector.proj(aPoint.asRealVector-positionVector);
	}
	closestsPointOnSegment{arg aPoint;
		if(segment,
			{
				^this.skewLinesNearestPoints(
					LineND(aPoint.asRealVector,this.projectPoint2Line(aPoint.asRealVector)-aPoint.asRealVector,true)
				)[0];
			},{
				^this.projectPoint2Line(aPoint.asRealVector);
		})
	}

	isOrthogonal{arg aLine; ^directionVector.isOrthogonal(aLine.directionVector);}

	isParallel{arg aLine;
		var a=this.directionVector.angle(aLine.directionVector);
		if((a==pi)||(a==0),{^true},{^false});
	}



	skewLinesNearestPoints{arg aLine;
		var p1 = this.positionVector;
		var p2 =  aLine.positionVector;
		var d1 = this.directionVector;
		var d2 =  aLine.directionVector;
		var delta = (sum(d1**2)*sum(d2**2)) - ((d1 <|> d2)**2);
		var delta1 = ((d2<|> (p1-p2)) * (d1<|> d2)) - ((d1<|> (p1-p2)) * sum(d2**2));
		var delta2 = ((d2<|> (p1-p2)) * sum(d1**2)) - ((d1<|> d2) * (d1<|> (p1-p2)));
		var t1 = delta1 / delta;
		var t2 = delta2 / delta;

		if (this.segment, {
			t1=t1.min(1).max(0);
		});
		if (aLine.segment, {
			t2=t2.min(1).max(0);
		});
		^[p1+(t1*d1), p2+(t2*d2)];
	}

	intersection{arg aLine;
		var p1,p2;
		#p1,p2 = this.skewLinesNearestPoints(aLine);
		if(p1.dist(p2)<0.00001,{^p1},{^false});
	}

	prPointAt{arg index,value;
		var r=(value/directionVector[index])-(positionVector[index]/directionVector[index]);
		^this.getPointOnLine(r);
	}

	getPointOnLine{arg s;
		^(this.positionVector+(s*this.directionVector));
	}
	//only for compatibility with plane !!! very stupid
	vectorFromPointOnPlane{|dot|
		^this.getPointOnLine(dot.x);
	}

}

Plane{
	var <>positionVector,<>spanA,<>spanB;
	*new { arg positionVector,spanA,spanB;
		^super.newCopyArgs(positionVector,
			spanA,
			spanB)
	}

	*newUsingPoints{arg positionVector1,positionVector2,positionVector3;
		^super.newCopyArgs(positionVector1,
			(positionVector2-positionVector1),
			(positionVector3-positionVector1));
	}
	pointOnPlane{|scalarA=0,scalarB=0|
		^positionVector+(spanA*scalarA)+(spanB*scalarB);
	}

	projectionScalars{arg aPoint;
		var d1 = this.spanA;
		var	d2 = this.spanB;
		var bp = this.positionVector;
		var  a = sum(d1**2);
		var  b = (d1<|> d2);
		var  d = b;
		var  e = sum(d2**2);
		var  c = ((aPoint.asRealVector-bp) <|> d1);
		var  f =  ((aPoint.asRealVector-bp) <|> d2);
		var delta = (a*e)-(b*d) + 0.0;
		var alpha = ((c*e)-(b*f)) / delta;
		var beta = ((a*f)-(c*d))/delta;
		^alpha@beta;
	}

	projectionPoint2Plane{arg aPoint;
		var ppoint;
		ppoint=this.projectionScalars(aPoint.asRealVector);
		^this.vectorFromPointOnPlane(ppoint);
	}
	vectorFromPointOnPlane{|dot|
		^(this.positionVector + (dot.x * this.spanA) + (dot.y * this.spanB));
	}

	dimensions{
		^positionVector.size;
	}


	aPerpendicularToPlane{
		var point=positionVector*2;
		var pointOnPlane=this.projectionPoint2Plane(point);
		^(point-pointOnPlane);
	}

	//special methods

	//strange help functions to get 2d points from vectors
	//when x is someting y whill create a [0,1] at position :[a,b];
	// fiy0{|i,x|((positionVector[i] * -1)-(x*spanA[i]))/spanB[i]}
	// fiy1{|i,x|^((positionVector[i] * -1)-(x*spanA[i])+1)/spanB[i]}
	// fix0{|i,y|^((positionVector[i] * -1)-(y*spanB[i]))/spanA[i]}
	// fix1{|i,y|^((positionVector[i] * -1)-(y*spanB[i])+1)/spanA[i]}

	// fiy00{|i|^(positionVector[i] * -1)/spanB[i]}
	// fiy01{|i|^((positionVector[i] * -1)-spanA[i])/spanB[i]}
	// fiy10{|i|^((positionVector[i] * -1) + 1)/spanB[i]}
	// fiy11{|i|^((positionVector[i] * -1) - spanA[i] + 1)/spanB[i]}
	// fix00{|i|^(positionVector[i] * -1) / spanA[i]}
	// fix01{|i|^((positionVector[i] * -1) - spanB[i])/spanA[i]}
	// fix10{|i|^((positionVector[i] * -1) + 1)/spanA[i]}
	// fix11{|i|^((positionVector[i] * -1) - spanB[i] + 1)/spanA[i]}


	linePlaneIntersection{|aLine|
		var p0=positionVector;
		var l0=aLine.positionVector;
		var n=this.aPerpendicularToPlane.normalize;
		var l=aLine.directionVector;
		var d=((p0-l0)<|>n)/(l<|>n);
		^l0+(d*l);
	}

	isCrossing{arg aLine;
		var p=this.linePlaneIntersection(aLine);
		var pp=this.projectionPoint2Plane(p);
		^(p-pp).postln.sum.abs.postln<0.000001;
	}

	//TODO
	//crossingLine{arg aPlane;}//returns a LineND
	// planePlaneIntersection{|aPlane|
	// 	var n1=this.aPerpendicularToPlane.normalize;
	// 	var n2=aPlane.aPerpendicularToPlane.normalize;
	// 	var r=1;//?i dont know what r is;
	// 	var h2=n2<|>r;
	// 	var h1=n1<|>r;
	// 	var c1=(h1-(h2*(n1<|>n2)))/(1-(n1<|>n2).pow(2));
	// 	var c2=(h2-(h1*(n1<|>n2)))/(1-(n1<|>n2).pow(2));
	// }



	//isParallel{arg aPlane;}
	//isOrthogonal{arg aLine;}

}


