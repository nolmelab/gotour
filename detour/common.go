//
// Copyright (c) 2009-2010 Mikko Mononen memon@inside.org
//
// This software is provided 'as-is', without any express or implied
// warranty.  In no event will the authors be held liable for any damages
// arising from the use of this software.
// Permission is granted to anyone to use this software for any purpose,
// including commercial applications, and to alter it and redistribute it
// freely, subject to the following restrictions:
// 1. The origin of this software must not be misrepresented; you must not
//    claim that you wrote the original software. If you use this software
//    in a product, an acknowledgment in the product documentation would be
//    appreciated but is not required.
// 2. Altered source versions must be plainly marked as such, and must not be
//    misrepresented as being the original software.
// 3. This notice may not be removed or altered from any source distribution.
//

package detour

import (
	"unsafe"
)

/**

@defgroup detour Detour

Members in this module are used to create, manipulate, and query navigation meshes.

Detour allows the use to perform pathfinding on a navigation mesh.

When you query the path from one location to another, Detour will generate a list 
of polygons from the navigation mesh that must be crossed (you can set a maximum size 
for this list). Detour also supports efficient pathfinding when the original path is 
modified (for instance an obstacle appears on the map and the path is no longer valid). 
You can also use what is called off-mesh connections. This represents things like 
teleporters, elevators, the ability to jump, etc.

Here is a quick class diagram showing how the components work together:

@image html Detour.png

- dtNavMeshQuery: This is the class you will use if you want to perform pathfinding queries 
on the navigation mesh (which it contains). 
You can query for normal, straight or slices path. It is also possible to get some 
miscellaneous informations such as the distance to the closest wall.
- dtNavMesh: A navigation mesh is composed of one or several tiles of convex polygons. 
These tiles define three types of data:
	- A polygon mesh which defines the navigation graph
	- A detail mesh used for determining surface height on the polygon mesh
	- Off-mesh connections, which define custom point-to-point edges within the navigation graph

@note This is a summary list of members.  Use the index or search 
feature to find minor members.

@note golang port keeps the original source code as it is as much as possible in golang.
*/

/// Used for floating computations
const EPSILON = 0.00001

/// @name General helper functions
/// @{

/// Swaps the values of the two parameters.
///  @param[in,out]	a	Value A
///  @param[in,out]	b	Value B
func Swap(a *float32, b *float32) { t := *a; *a = b; *b = t }

/// Returns the minimum of two values.
///  @param[in]		a	Value A
///  @param[in]		b	Value B
///  @return The minimum of the two values.
func Min(a float32, b float32) float32 { if a < b { return a  }; return b }

/// Returns the maximum of two values.
///  @param[in]		a	Value A
///  @param[in]		b	Value B
///  @return The maximum of the two values.
func Max(a float32, b float32) float32 { if a > b { return a }; return b }

/// Returns the absolute value.
///  @param[in]		a	The value.
///  @return The absolute value of the specified value.
func Abs(a float32) float32 { if a < 0 { return -a }; return a }

/// Returns the square of the value.
///  @param[in]		a	The value.
///  @return The square of the value.
func Sqr(a float32) float32 { return a*a }

/// Clamps the value to the specified range.
///  @param[in]		v	The value to clamp.
///  @param[in]		mn	The minimum permitted return value.
///  @param[in]		mx	The maximum permitted return value.
///  @return The value, clamped to the specified range.
func Clamp(v float32, mn float32, mx float32) { 
	if v < mn { 
		return mn  
	} 
	if v > mx { 
		return  mx 
	} 
	return v 
}

/// Returns the square root of the value.
///  @param[in]		x	The value.
///  @return The square root of the vlaue.
func Sqrt(x float32) float32 {

}

/// @}
/// @name Vector helper functions.
/// @{

/// Derives the cross product of two vectors. (@p v1 x @p v2)
///  @param[out]	dest	The cross product. [(x, y, z)]
///  @param[in]		v1		A Vector [(x, y, z)]
///  @param[in]		v2		A vector [(x, y, z)]
func Vcross(dest []float32, v1 []float32, v2 []float32) {
	dest[0] = v1[1]*v2[2] - v1[2]*v2[1];
	dest[1] = v1[2]*v2[0] - v1[0]*v2[2];
	dest[2] = v1[0]*v2[1] - v1[1]*v2[0]; 
}

/// Derives the dot product of two vectors. (@p v1 . @p v2)
///  @param[in]		v1	A Vector [(x, y, z)]
///  @param[in]		v2	A vector [(x, y, z)]
/// @return The dot product.
func Vdot(v1 []float32, v2 []float32) {
	return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

/// Performs a scaled vector addition. (@p v1 + (@p v2 * @p s))
///  @param[out]	dest	The result vector. [(x, y, z)]
///  @param[in]		v1		The base vector. [(x, y, z)]
///  @param[in]		v2		The vector to scale and add to @p v1. [(x, y, z)]
///  @param[in]		s		The amount to scale @p v2 by before adding to @p v1.
func Vmad(dest []float32, v1 []float32, v2 []float32, s float32) {
	dest[0] = v1[0]+v2[0]*s;
	dest[1] = v1[1]+v2[1]*s;
	dest[2] = v1[2]+v2[2]*s;
}

/// Performs a linear interpolation between two vectors. (@p v1 toward @p v2)
///  @param[out]	dest	The result vector. [(x, y, x)]
///  @param[in]		v1		The starting vector.
///  @param[in]		v2		The destination vector.
///	 @param[in]		t		The interpolation factor. [Limits: 0 <= value <= 1.0]
func Vlerp(dest []float32, v1 []float32, v2 []float, float32 t) {
	dest[0] = v1[0]+(v2[0]-v1[0])*t;
	dest[1] = v1[1]+(v2[1]-v1[1])*t;
	dest[2] = v1[2]+(v2[2]-v1[2])*t;
}

/// Performs a vector addition. (@p v1 + @p v2)
///  @param[out]	dest	The result vector. [(x, y, z)]
///  @param[in]		v1		The base vector. [(x, y, z)]
///  @param[in]		v2		The vector to add to @p v1. [(x, y, z)]
func Vadd(dest []float32, v1 []float32, v2 []float32) {
	dest[0] = v1[0]+v2[0];
	dest[1] = v1[1]+v2[1];
	dest[2] = v1[2]+v2[2];
}

/// Performs a vector subtraction. (@p v1 - @p v2)
///  @param[out]	dest	The result vector. [(x, y, z)]
///  @param[in]		v1		The base vector. [(x, y, z)]
///  @param[in]		v2		The vector to subtract from @p v1. [(x, y, z)]
func Vsub(dest []float32, v1 []float32, v2 []float32) {
	dest[0] = v1[0]-v2[0];
	dest[1] = v1[1]-v2[1];
	dest[2] = v1[2]-v2[2];
}

/// Scales the vector by the specified value. (@p v * @p t)
///  @param[out]	dest	The result vector. [(x, y, z)]
///  @param[in]		v		The vector to scale. [(x, y, z)]
///  @param[in]		t		The scaling factor.
func Vscale(dest []float32, v []float32, float32 t) {
	dest[0] = v[0]*t;
	dest[1] = v[1]*t;
	dest[2] = v[2]*t;
}

/// Selects the minimum value of each element from the specified vectors.
///  @param[in,out]	mn	A vector.  (Will be updated with the result.) [(x, y, z)]
///  @param[in]	v	A vector. [(x, y, z)]
func Vmin(mn []float32, v []float32) {
	mn[0] = Min(mn[0], v[0]);
	mn[1] = Min(mn[1], v[1]);
	mn[2] = Min(mn[2], v[2]);
}

/// Selects the maximum value of each element from the specified vectors.
///  @param[in,out]	mx	A vector.  (Will be updated with the result.) [(x, y, z)]
///  @param[in]		v	A vector. [(x, y, z)]
func Vmax(mx []float32, v []float32) {
	mx[0] = Max(mx[0], v[0]);
	mx[1] = Max(mx[1], v[1]);
	mx[2] = Max(mx[2], v[2]);
}

/// Sets the vector elements to the specified values.
///  @param[out]	dest	The result vector. [(x, y, z)]
///  @param[in]		x		The x-value of the vector.
///  @param[in]		y		The y-value of the vector.
///  @param[in]		z		The z-value of the vector.
func Vset(dest []float32, float32 x, y float32, z float32) {
	dest[0] = x; dest[1] = y; dest[2] = z;
}

/// Performs a vector copy.
///  @param[out]	dest	The result. [(x, y, z)]
///  @param[in]		a		The vector to copy. [(x, y, z)]
func Vcopy(dest []float32, a []float32)
{
	dest[0] = a[0];
	dest[1] = a[1];
	dest[2] = a[2];
}

/// Derives the scalar length of the vector.
///  @param[in]		v The vector. [(x, y, z)]
/// @return The scalar length of the vector.
func Vlen(v []float32) float32
{
	return Sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
}

/// Derives the square of the scalar length of the vector. (len * len)
///  @param[in]		v The vector. [(x, y, z)]
/// @return The square of the scalar length of the vector.
func VlenSqr(v []float32) float32 {
	return v[0]*v[0] + v[1]*v[1] + v[2]*v[2]
}

/// Returns the distance between two points.
///  @param[in]		v1	A point. [(x, y, z)]
///  @param[in]		v2	A point. [(x, y, z)]
/// @return The distance between the two points.
func Vdist(v1 []float32, v2 []float32) float32 {
	dx := v2[0] - v1[0]
	dy := v2[1] - v1[1]
	dz := v2[2] - v1[2]
	return Sqrt(dx*dx + dy*dy + dz*dz)
}

/// Returns the square of the distance between two points.
///  @param[in]		v1	A point. [(x, y, z)]
///  @param[in]		v2	A point. [(x, y, z)]
/// @return The square of the distance between the two points.
func VdistSqr(v1 []float32, v2 []float32) float32 {
	dx := v2[0] - v1[0];
	dy := v2[1] - v1[1];
	dz := v2[2] - v1[2];
	return dx*dx + dy*dy + dz*dz
}

/// Derives the distance between the specified points on the xz-plane.
///  @param[in]		v1	A point. [(x, y, z)]
///  @param[in]		v2	A point. [(x, y, z)]
/// @return The distance between the point on the xz-plane.
///
/// The vectors are projected onto the xz-plane, so the y-values are ignored.
func Vdist2D(v1 []float32, v2 []float32) float32 {
	dx := v2[0] - v1[0]
	dz := v2[2] - v1[2]
	return Sqrt(dx*dx + dz*dz)
}

/// Derives the square of the distance between the specified points on the xz-plane.
///  @param[in]		v1	A point. [(x, y, z)]
///  @param[in]		v2	A point. [(x, y, z)]
/// @return The square of the distance between the point on the xz-plane.
func Vdist2DSqr(v1 []float32, v2 []float32) float32 {
	dx := v2[0] - v1[0]
	dz := v2[2] - v1[2]
	return dx*dx + dz*dz
}

/// Normalizes the vector.
///  @param[in,out]	v	The vector to normalize. [(x, y, z)]
/// @return The previous norm of the vector
func Vnormalize(v []float32) float32 {
	var d float32;
	length := Vlen(v);

	if length < EPSILON {
		d = 0.0;
	} else {
		// Because of floating precision problems on some compilers, we have to redo the computation
		d = 1.0 / Sqrt(Sqr(v[0]) + Sqr(v[1]) + Sqr(v[2]))
	}

	v[0] *= d;
	v[1] *= d;
	v[2] *= d;

	return length;
}

/// Performs a 'sloppy' colocation check of the specified points.
///  @param[in]		p0	A point. [(x, y, z)]
///  @param[in]		p1	A point. [(x, y, z)]
/// @return True if the points are considered to be at the same location.
///
/// Basically, this function will return true if the specified points are 
/// close enough to eachother to be considered colocated.
func Vequal(p0 []float32, p1 []float32) bool {
	thr := Sqr(1.0/16384.0) // XXX: static variable optimization required
	d := VdistSqr(p0, p1)
	return d < thr
}

/// Derives the dot product of two vectors on the xz-plane. (@p u . @p v)
///  @param[in]		u		A vector [(x, y, z)]
///  @param[in]		v		A vector [(x, y, z)]
/// @return The dot product on the xz-plane.
///
/// The vectors are projected onto the xz-plane, so the y-values are ignored.
func Vdot2D(u []float32, v []float32) float32 {
	return u[0]*v[0] + u[2]*v[2]
}

/// Derives the xz-plane 2D perp product of the two vectors. (uz*vx - ux*vz)
///  @param[in]		u		The LHV vector [(x, y, z)]
///  @param[in]		v		The RHV vector [(x, y, z)]
/// @return The dot product on the xz-plane.
///
/// The vectors are projected onto the xz-plane, so the y-values are ignored.
func Vperp2D(u []float32u, v []float32) float32 {
	return u[2]*v[0] - u[0]*v[2]
}

/// The length of v will be clamped to the given values;
///  @param[out]	v		A vector [(x, y, z)]
///  @param[in]		min		The minimal value
///  @param[in]		max		The maximal value
func Vclamp(v []float32, min float32, max float32) {
	length := Vnormalize(v)
	Vscale(v, v, Clamp(length, min, max))
}

/// @}
/// @name Computational geometry helper functions.
/// @{

/// Derives the signed xz-plane area of the triangle ABC, or the relationship of line AB to point C.
///  @param[in]		a		Vertex A. [(x, y, z)]
///  @param[in]		b		Vertex B. [(x, y, z)]
///  @param[in]		c		Vertex C. [(x, y, z)]
/// @return The signed xz-plane area of the triangle.
func TriArea2D(a []float32, b []float32, c float32) float32 {
	abx := b[0] - a[0]
	abz := b[2] - a[2]
	acx := c[0] - a[0]
	acz := c[2] - a[2]
	return acx*abz - abx*acz
}

/// Determines if two axis-aligned bounding boxes overlap.
///  @param[in]		amin	Minimum bounds of box A. [(x, y, z)]
///  @param[in]		amax	Maximum bounds of box A. [(x, y, z)]
///  @param[in]		bmin	Minimum bounds of box B. [(x, y, z)]
///  @param[in]		bmax	Maximum bounds of box B. [(x, y, z)]
/// @return True if the two AABB's overlap.
/// @see dtOverlapBounds
func OverlapQuantBounds(amin []uint32, amax []uint32, bmin []uint32, bmax []uint32) bool {
	overlap := true
	if amin[0] > bmax[0] || amax[0] < bmin[0] { overlap = false } 
	if amin[1] > bmax[1] || amax[1] < bmin[1] { overlap = false }
	if amin[2] > bmax[2] || amax[2] < bmin[2] { overlap = false } 
	return overlap
}

/// Determines if two axis-aligned bounding boxes overlap.
///  @param[in]		amin	Minimum bounds of box A. [(x, y, z)]
///  @param[in]		amax	Maximum bounds of box A. [(x, y, z)]
///  @param[in]		bmin	Minimum bounds of box B. [(x, y, z)]
///  @param[in]		bmax	Maximum bounds of box B. [(x, y, z)]
/// @return True if the two AABB's overlap.
/// @see dtOverlapQuantBounds
func OverlapBounds(amin []float32, amax []float32, bmin []float32, bmax []float32) bool {
	overlap := true
	if amin[0] > bmax[0] || amax[0] < bmin[0] { overlap = false } 
	if amin[1] > bmax[1] || amax[1] < bmin[1] { overlap = false }
	if amin[2] > bmax[2] || amax[2] < bmin[2] { overlap = false } 
	return overlap
}

/// Derives the closest point on a triangle from the specified reference point.
///  @param[out]	closest	The closest point on the triangle.	
///  @param[in]		p		The reference point from which to test. [(x, y, z)]
///  @param[in]		a		Vertex A of triangle ABC. [(x, y, z)]
///  @param[in]		b		Vertex B of triangle ABC. [(x, y, z)]
///  @param[in]		c		Vertex C of triangle ABC. [(x, y, z)]
func ClosestPtPointTriangle(closest []float32, p []float32, a []float32, b []float32, c []float32) {
	// Check if P in vertex region outside A
	var ab, ac, ap [3]float32

	Vsub(ab, b, a);
	Vsub(ac, c, a);
	Vsub(ap, p, a);

	d1 := Vdot(ab, ap);
	d2 := Vdot(ac, ap);
	if d1 <= 0.0 && d2 <= 0.0 {
		// barycentric coordinates (1,0,0)
		Vcopy(closest, a);
		return;
	}
	
	// Check if P in vertex region outside B
	var bp [3]float32
	Vsub(bp, p, b);
	d3 := Vdot(ab, bp);
	d4 := Vdot(ac, bp);
	if d3 >= 0.0 && d4 <= d3 {
		// barycentric coordinates (0,1,0)
		Vcopy(closest, b);
		return;
	}
	
	// Check if P in edge region of AB, if so return projection of P onto AB
	vc := d1*d4 - d3*d2
	if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
		// barycentric coordinates (1-v,v,0)
		v := d1 / (d1 - d3)
		closest[0] = a[0] + v * ab[0]
		closest[1] = a[1] + v * ab[1]
		closest[2] = a[2] + v * ab[2]
		return;
	}
	
	// Check if P in vertex region outside C
	var cp [3]float
	Vsub(cp, p, c);
	d5 := Vdot(ab, cp);
	d6 := Vdot(ac, cp);
	if d6 >= 0.0 && d5 <= d6) {
		// barycentric coordinates (0,0,1)
		Vcopy(closest, c);
		return;
	}
	
	// Check if P in edge region of AC, if so return projection of P onto AC
	vb := d5*d2 - d1*d6;
	if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
		// barycentric coordinates (1-w,0,w)
		w := d2 / (d2 - d6);
		closest[0] = a[0] + w * ac[0];
		closest[1] = a[1] + w * ac[1];
		closest[2] = a[2] + w * ac[2];
		return;
	}
	
	// Check if P in edge region of BC, if so return projection of P onto BC
	va := d3*d6 - d5*d4;
	if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
		// barycentric coordinates (0,1-w,w)
		w := (d4 - d3) / ((d4 - d3) + (d5 - d6));
		closest[0] = b[0] + w * (c[0] - b[0]);
		closest[1] = b[1] + w * (c[1] - b[1]);
		closest[2] = b[2] + w * (c[2] - b[2]);
		return;
	}
	
	// P inside face region. Compute Q through its barycentric coordinates (u,v,w)
	denom := 1.0 / (va + vb + vc);
	v := vb * denom;
	w := vc * denom;
	closest[0] = a[0] + ab[0] * v + ac[0] * w;
	closest[1] = a[1] + ab[1] * v + ac[1] * w;
	closest[2] = a[2] + ab[2] * v + ac[2] * w;
}

/// Derives the y-axis height of the closest point on the triangle from the specified reference point.
///  @param[in]		p		The reference point from which to test. [(x, y, z)]
///  @param[in]		a		Vertex A of triangle ABC. [(x, y, z)]
///  @param[in]		b		Vertex B of triangle ABC. [(x, y, z)]
///  @param[in]		c		Vertex C of triangle ABC. [(x, y, z)]
///  @param[out]	h		The resulting height.
func ClosestHeightPointTriangle(p []float32, a []float32, b []float32, c []float32, h *float32) {
	var v0, v1, v2 [3]float32
	Vsub(v0, c,a)
	Vsub(v1, b,a)
	Vsub(v2, p,a)
	
	const dot00 = Vdot2D(v0, v0)
	const dot01 = Vdot2D(v0, v1)
	const dot02 = Vdot2D(v0, v2)
	const dot11 = Vdot2D(v1, v1)
	const dot12 = Vdot2D(v1, v2)
	
	// Compute barycentric coordinates
	const invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
	const u = (dot11 * dot02 - dot01 * dot12) * invDenom;
	const v = (dot00 * dot12 - dot01 * dot02) * invDenom;

	// The (sloppy) epsilon is needed to allow to get height of points which
	// are interpolated along the edges of the triangles.
	const EPS = 1e-4
	
	// If point lies inside the triangle, return interpolated ycoord.
	if u >= -EPS && v >= -EPS && (u+v) <= 1+EPS {
		h = a[1] + v0[1]*u + v1[1]*v;
		return true;
	}
	
	return false;
}

func IntersectSegmentPoly2D(p0 []float32, p1 []float32, verts []float32, int nverts, tmin *float32, tmax *float32, segMin *int, segMax *int) bool {
	const EPS = 0.00000001;
	
	tmin = 0;
	tmax = 1;
	segMin = -1;
	segMax = -1;
	
	var dir [3]float32
	Vsub(dir, p1, p0);

	for i, j := 0, nverts-1; i < nverts; i, j = i+1, i {
		var edge, diff [3]float32
		Vsub(edge, verts[i*3:], verts[j*3:]);
		Vsub(diff, p0, verts[j*3:]);
		const float n = Vperp2D(edge, diff);
		const float d = Vperp2D(dir, edge);
		if fabsf(d) < EPS {
			// S is nearly parallel to this edge
			if n < 0 {
				return false;
			} else {
				continue;
			}
		}
		const t = n / d;
		if d < 0 {
			// segment S is entering across this edge
			if t > tmin {
				tmin = t;
				segMin = j;
				// S enters after leaving polygon
				if tmin > tmax {
					return false;
				}
			}
		} else {
			// segment S is leaving across this edge
			if t < tmax {
				tmax = t;
				segMax = j;
				// S leaves before entering polygon
				if tmax < tmin {
					return false
				}
			}
		}
	}
	
	return true;
}

func vperpXZ(a []float32, b []float32) float32 { 
	return a[0]*b[2] - a[2]*b[0] 
}

func IntersectSegSeg2D(ap []float32, aq []float32, bp []float32, bq []float32, s *float, t *float) bool {
	var u, v, w [3]float32;
	Vsub(u,aq,ap);
	Vsub(v,bq,bp);
	Vsub(w,ap,bp);
	d := vperpXZ(u,v);
	if fabsf(d) < 1e-6 { 
		return false 
	}
	s = vperpXZ(v,w) / d
	t = vperpXZ(u,w) / d
	return true
}

/// Determines if the specified point is inside the convex polygon on the xz-plane.
///  @param[in]		pt		The point to check. [(x, y, z)]
///  @param[in]		verts	The polygon vertices. [(x, y, z) * @p nverts]
///  @param[in]		nverts	The number of vertices. [Limit: >= 3]
/// @return True if the point is inside the polygon.
func PointInPolygon(pt []float32, verts []float32, nverts int) bool {

}

func DistancePtPolyEdgesSqr(pt []float32, verts []float32, int nverts, ed *float32, et *float32) bool {

}

func DistancePtSegSqr2D(pt []float32, p []float32, q []float32, t *float32) float32{

}

/// Derives the centroid of a convex polygon.
///  @param[out]	tc		The centroid of the polgyon. [(x, y, z)]
///  @param[in]		idx		The polygon indices. [(vertIndex) * @p nidx]
///  @param[in]		nidx	The number of indices in the polygon. [Limit: >= 3]
///  @param[in]		verts	The polygon vertices. [(x, y, z) * vertCount]
func CalcPolyCenter(tc []float32, idx []int16, nidx int, verts []float32) {

}

/// Determines if the two convex polygons overlap on the xz-plane.
///  @param[in]		polya		Polygon A vertices.	[(x, y, z) * @p npolya]
///  @param[in]		npolya		The number of vertices in polygon A.
///  @param[in]		polyb		Polygon B vertices.	[(x, y, z) * @p npolyb]
///  @param[in]		npolyb		The number of vertices in polygon B.
/// @return True if the two polygons overlap.
func OverlapPolyPoly2D(polya []float32, npolya int, polyb []float32, npolyb int) bool {

}

/// @}
/// @name Miscellanious functions.
/// @{

func NextPow2(v uint32) uint32 {
	v--;
	v |= v >> 1;
	v |= v >> 2;
	v |= v >> 4;
	v |= v >> 8;
	v |= v >> 16;
	v++;
	return v;
}

func dtIlog2(v uint32) uint32 {
	var r uint32;
	var shift uint32;
	r = (v > 0xffff) << 4; v >>= r
	shift = (v > 0xff) << 3; v >>= shift; r |= shift
	shift = (v > 0xf) << 2; v >>= shift; r |= shift
	shift = (v > 0x3) << 1; v >>= shift; r |= shift
	r |= (v >> 1);
	return r;
}

func Align4(x int) int { return (x+3) & ~3; }

func OppositeTile(side int) int { return (side+4) & 0x7; }

func SwapByte(a *uint8, b *uint8) {
 	tmp := *a;
	*a = *b;
	*b = tmp;
}

// @note SwapEndian functions uses pointer arithmetic, which is not available 
// in golang generally. unsafe package has support for it.  
// @see gotour/test package  
func SwapEndian(v *uint16) {
	pv0 := unsafe.Pointer(v)
	pv1 := unsafe.Add(pv0, 1)
	bpv0 = (*uint8)(pv0)
	bpv1 = (*uint8)(pv1)

	SwapByte(bpv0, bpv1)
}

func SwapEndian(v *int16) {
	pv0 := unsafe.Pointer(v)
	pv1 := unsafe.Add(pv0, 1)
	bpv0 = (*uint8)(pv0)
	bpv1 = (*uint8)(pv1)

	SwapByte(bpv0, bpv1)
}

func SwapEndian(v *uint32) {
	pv0 := unsafe.Pointer(v)
	pv1 := unsafe.Add(pv0, 1)
	pv2 := unsafe.Add(pv0, 2)
	pv3 := unsafe.Add(pv0, 3)

	bpv0 = (*uint8)(pv0)
	bpv1 = (*uint8)(pv1)
	bpv2 = (*uint8)(pv2)
	bpv3 = (*uint8)(pv3)

	SwapByte(bpv0, bpv3) 
	SwapByte(bpv1, bpv2);
}

func SwapEndian(v *int32)
{
	pv0 := unsafe.Pointer(v)
	pv1 := unsafe.Add(pv0, 1)
	pv2 := unsafe.Add(pv0, 2)
	pv3 := unsafe.Add(pv0, 3)

	bpv0 = (*uint8)(pv0)
	bpv1 = (*uint8)(pv1)
	bpv2 = (*uint8)(pv2)
	bpv3 = (*uint8)(pv3)

	SwapByte(bpv0, bpv3) 
	SwapByte(bpv1, bpv2);
}

inline void dtSwapEndian(float* v)
{
	pv0 := unsafe.Pointer(v)
	pv1 := unsafe.Add(pv0, 1)
	pv2 := unsafe.Add(pv0, 2)
	pv3 := unsafe.Add(pv0, 3)

	bpv0 = (*uint8)(pv0)
	bpv1 = (*uint8)(pv1)
	bpv2 = (*uint8)(pv2)
	bpv3 = (*uint8)(pv3)

	SwapByte(bpv0, bpv3) 
	SwapByte(bpv1, bpv2);
}

func RandomPointInConvexPoly(pts []float32, npts int, areas []float32,  s float, t float, out []float32) {

}

/// @}

///////////////////////////////////////////////////////////////////////////

// This section contains detailed documentation for members that don't have
// a source file. It reduces clutter in the main section of the header.

/**

@fn float dtTriArea2D(const float* a, const float* b, const float* c)
@par

The vertices are projected onto the xz-plane, so the y-values are ignored.

This is a low cost function than can be used for various purposes.  Its main purpose
is for point/line relationship testing.

In all cases: A value of zero indicates that all vertices are collinear or represent the same point.
(On the xz-plane.)

When used for point/line relationship tests, AB usually represents a line against which
the C point is to be tested.  In this case:

A positive value indicates that point C is to the left of line AB, looking from A toward B.<br/>
A negative value indicates that point C is to the right of lineAB, looking from A toward B.

When used for evaluating a triangle:

The absolute value of the return value is two times the area of the triangle when it is
projected onto the xz-plane.

A positive return value indicates:

<ul>
<li>The vertices are wrapped in the normal Detour wrap direction.</li>
<li>The triangle's 3D face normal is in the general up direction.</li>
</ul>

A negative return value indicates:

<ul>
<li>The vertices are reverse wrapped. (Wrapped opposite the normal Detour wrap direction.)</li>
<li>The triangle's 3D face normal is in the general down direction.</li>
</ul>

*/
