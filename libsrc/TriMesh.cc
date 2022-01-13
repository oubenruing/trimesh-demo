#include "TriMesh.h"
#include "emscripten/bind.h"
#include <emscripten/val.h>

using namespace emscripten;
// i+1 and i-1 modulo 3
#define NEXT(i) ((i) < 2 ? (i) + 1 : (i)-2)
#define PREV(i) ((i) > 0 ? (i)-1 : (i) + 2)

namespace trimesh
{
	std::vector<float> TriMesh::draw_boundaries(bool do_hidden)
	{
		std::vector<float>().swap(draw_points);
		need_faces();
		need_across_edge();
		for (size_t i = 0; i < faces.size(); i++)
		{
			for (int j = 0; j < 3; j++)
			{
				if (across_edge[i][j] >= 0)
					continue;
				int v1 = faces[i][(j + 1) % 3];
				int v2 = faces[i][(j + 2) % 3];
				draw_points.push_back(vertices[v1].x);
				draw_points.push_back(vertices[v1].y);
				draw_points.push_back(vertices[v1].z);
				draw_points.push_back(vertices[v2].x);
				draw_points.push_back(vertices[v2].y);
				draw_points.push_back(vertices[v2].z);
			}
		}
		return draw_points;
	}

	/****************************** 绘制特殊线条所需要的算法 start ******************************/

	// Compute gradient of (kr * sin^2 theta) at vertex i
	inline vec TriMesh::gradkr(int i)
	{
		vec viewdir = viewpos - vertices[i];
		float rlen_viewdir = 1.0f / len(viewdir);
		viewdir *= rlen_viewdir;

		float f_ndotv = viewdir DOT normals[i];
		float sintheta = sqrt(1.0f - sqr(f_ndotv));
		float csctheta = 1.0f / sintheta;
		float u = (viewdir DOT pdir1[i]) * csctheta;
		float v = (viewdir DOT pdir2[i]) * csctheta;
		float f_kr = curv1[i] * u * u + curv2[i] * v * v;
		float tr = u * v * (curv2[i] - curv1[i]);
		float kt = curv1[i] * (1.0f - u * u) +
							 curv2[i] * (1.0f - v * v);
		vec w = u * pdir1[i] + v * pdir2[i];
		vec wperp = u * pdir2[i] - v * pdir1[i];
		const Vec<4> &C = dcurv[i];

		vec g = pdir1[i] * (u * u * C[0] + 2.0f * u * v * C[1] + v * v * C[2]) +
						pdir2[i] * (u * u * C[1] + 2.0f * u * v * C[2] + v * v * C[3]) -
						2.0f * csctheta * tr * (rlen_viewdir * wperp + f_ndotv * (tr * w + kt * wperp));
		g *= (1.0f - sqr(f_ndotv));
		g -= 2.0f * f_kr * sintheta * f_ndotv * (f_kr * w + tr * wperp);
		return g;
	}

	// Find a zero crossing between val0 and val1 by linear interpolation
	// Returns 0 if zero crossing is at val0, 1 if at val1, etc.
	inline float TriMesh::find_zero_linear(float val0, float val1)
	{
		return val0 / (val0 - val1);
	}

	// Find a zero crossing using Hermite interpolation
	float TriMesh::find_zero_hermite(int v0, int v1, float val0, float val1,
													const vec &grad0, const vec &grad1)
	{
		if (unlikely(val0 == val1))
			return 0.5f;

		// Find derivatives along edge (of interpolation parameter in [0,1]
		// which means that e01 doesn't get normalized)
		vec e01 = vertices[v1] - vertices[v0];
		float d0 = e01 DOT grad0, d1 = e01 DOT grad1;

		// This next line would reduce val to linear interpolation
		//d0 = d1 = (val1 - val0);

		// Use hermite interpolation:
		//   val(s) = h1(s)*val0 + h2(s)*val1 + h3(s)*d0 + h4(s)*d1
		// where
		//  h1(s) = 2*s^3 - 3*s^2 + 1
		//  h2(s) = 3*s^2 - 2*s^3
		//  h3(s) = s^3 - 2*s^2 + s
		//  h4(s) = s^3 - s^2
		//
		//  val(s)  = [2(val0-val1) +d0+d1]*s^3 +
		//            [3(val1-val0)-2d0-d1]*s^2 + d0*s + val0
		// where
		//
		//  val(0) = val0; val(1) = val1; val'(0) = d0; val'(1) = d1
		//

		// Coeffs of cubic a*s^3 + b*s^2 + c*s + d
		float a = 2 * (val0 - val1) + d0 + d1;
		float b = 3 * (val1 - val0) - 2 * d0 - d1;
		float c = d0, d = val0;

		// -- Find a root by bisection
		// (as Newton can wander out of desired interval)

		// Start with entire [0,1] interval
		float sl = 0.0f, sr = 1.0f, valsl = val0, valsr = val1;

		// Check if we're in a (somewhat uncommon) 3-root situation, and pick
		// the middle root if it happens (given we aren't drawing curvy lines,
		// seems the best approach..)
		//
		// Find extrema of derivative (a -> 3a; b -> 2b, c -> c),
		// and check if they're both in [0,1] and have different signs
		float disc = 4 * b - 12 * a * c;
		if (disc > 0 && a != 0)
		{
			disc = sqrt(disc);
			float r1 = (-2 * b + disc) / (6 * a);
			float r2 = (-2 * b - disc) / (6 * a);
			if (r1 >= 0 && r1 <= 1 && r2 >= 0 && r2 <= 1)
			{
				float vr1 = (((a * r1 + b) * r1 + c) * r1) + d;
				float vr2 = (((a * r2 + b) * r2 + c) * r2) + d;
				// When extrema have different signs inside an
				// interval with endpoints with different signs,
				// the middle root is in between the two extrema
				if ((vr1 < 0.0f && vr2 >= 0.0f) ||
						(vr1 > 0.0f && vr2 <= 0.0f))
				{
					// 3 roots
					if (r1 < r2)
					{
						sl = r1;
						valsl = vr1;
						sr = r2;
						valsr = vr2;
					}
					else
					{
						sl = r2;
						valsl = vr2;
						sr = r1;
						valsr = vr1;
					}
				}
			}
		}

		// Bisection method (constant number of interations)
		for (int iter = 0; iter < 10; iter++)
		{
			float sbi = (sl + sr) / 2.0f;
			float valsbi = (((a * sbi + b) * sbi) + c) * sbi + d;

			// Keep the half which has different signs
			if ((valsl < 0.0f && valsbi >= 0.0f) ||
					(valsl > 0.0f && valsbi <= 0.0f))
			{
				sr = sbi;
				valsr = valsbi;
			}
			else
			{
				sl = sbi;
				valsl = valsbi;
			}
		}

		return 0.5f * (sl + sr);
	}

	// Draw part of a zero-crossing curve on one triangle face, but only if
	// "test_num/test_den" is positive.  v0,v1,v2 are the indices of the 3
	// vertices, "val" are the values of the scalar field whose zero
	// crossings we are finding, and "test_*" are the values we are testing
	// to make sure they are positive.  This function assumes that val0 has
	// opposite sign from val1 and val2 - the following function is the
	// general one that figures out which one actually has the different sign.
	void TriMesh::draw_face_isoline2(int v0, int v1, int v2,
													const std::vector<float> &val,
													const std::vector<float> &test_num,
													const std::vector<float> &test_den,
													bool do_hermite, bool do_test, float fade)
	{
		// How far along each edge?
		float w10 = do_hermite ? find_zero_hermite(v0, v1, val[v0], val[v1],
																							 gradkr(v0), gradkr(v1))
													 : find_zero_linear(val[v0], val[v1]);
		float w01 = 1.0f - w10;
		float w20 = do_hermite ? find_zero_hermite(v0, v2, val[v0], val[v2],
																							 gradkr(v0), gradkr(v2))
													 : find_zero_linear(val[v0], val[v2]);
		float w02 = 1.0f - w20;

		// Points along edges
		point p1 = w01 * vertices[v0] + w10 * vertices[v1];
		point p2 = w02 * vertices[v0] + w20 * vertices[v2];

		float test_num1 = 1.0f, test_num2 = 1.0f;
		float test_den1 = 1.0f, test_den2 = 1.0f;
		float z1 = 0.0f, z2 = 0.0f;
		bool valid1 = true;
		if (do_test)
		{
			// Interpolate to find value of test at p1, p2
			test_num1 = w01 * test_num[v0] + w10 * test_num[v1];
			test_num2 = w02 * test_num[v0] + w20 * test_num[v2];
			if (!test_den.empty())
			{
				test_den1 = w01 * test_den[v0] + w10 * test_den[v1];
				test_den2 = w02 * test_den[v0] + w20 * test_den[v2];
			}
			// First point is valid iff num1/den1 is positive,
			// i.e. the num and den have the same sign
			valid1 = ((test_num1 >= 0.0f) == (test_den1 >= 0.0f));
			// There are two possible zero crossings of the test,
			// corresponding to zeros of the num and den
			if ((test_num1 >= 0.0f) != (test_num2 >= 0.0f))
				z1 = test_num1 / (test_num1 - test_num2);
			if ((test_den1 >= 0.0f) != (test_den2 >= 0.0f))
				z2 = test_den1 / (test_den1 - test_den2);
			// Sort and order the zero crossings
			if (z1 == 0.0f)
				z1 = z2, z2 = 0.0f;
			else if (z2 < z1)
				std::swap(z1, z2);
		}

		// If the beginning of the segment was not valid, and
		// no zero crossings, then whole segment invalid
		if (!valid1 && !z1 && !z2)
			return;

		// Draw the valid piece(s)
		int npts = 0;
		if (valid1)
		{
			// glColor4f(currcolor[0], currcolor[1], currcolor[2],
			// 	  test_num1 / (test_den1 * fade + test_num1));
			// glVertex3fv(p1);
			draw_points.push_back(p1.x);
			draw_points.push_back(p1.y);
			draw_points.push_back(p1.z);
			npts++;
		}
		if (z1)
		{
			float num = (1.0f - z1) * test_num1 + z1 * test_num2;
			float den = (1.0f - z1) * test_den1 + z1 * test_den2;
			// glColor4f(currcolor[0], currcolor[1], currcolor[2],
			// 	  num / (den * fade + num));
			// glVertex3fv((1.0f - z1) * p1 + z1 * p2);
			const point p = (1.0f - z1) * p1 + z1 * p2;
			draw_points.push_back(p.x);
			draw_points.push_back(p.y);
			draw_points.push_back(p.z);
			npts++;
		}
		if (z2)
		{
			float num = (1.0f - z2) * test_num1 + z2 * test_num2;
			float den = (1.0f - z2) * test_den1 + z2 * test_den2;
			// glColor4f(currcolor[0], currcolor[1], currcolor[2],
			// 	  num / (den * fade + num));
			// glVertex3fv((1.0f - z2) * p1 + z2 * p2);
			const point p = (1.0f - z2) * p1 + z2 * p2;
			draw_points.push_back(p.x);
			draw_points.push_back(p.y);
			draw_points.push_back(p.z);
			npts++;
		}
		if (npts != 2)
		{
			// glColor4f(currcolor[0], currcolor[1], currcolor[2],
			// 	  test_num2 / (test_den2 * fade + test_num2));
			// glVertex3fv(p2);
			draw_points.push_back(p2.x);
			draw_points.push_back(p2.y);
			draw_points.push_back(p2.z);
		}
	}

	// See above.  This is the driver function that figures out which of
	// v0, v1, v2 has a different sign from the others.
	void TriMesh::draw_face_isoline(int v0, int v1, int v2,
												 const std::vector<float> &val,
												 const std::vector<float> &test_num,
												 const std::vector<float> &test_den,
												 bool do_bfcull, bool do_hermite,
												 bool do_test, float fade)
	{
		// Backface culling
		if (likely(do_bfcull && ndotv[v0] <= 0.0f &&
							 ndotv[v1] <= 0.0f && ndotv[v2] <= 0.0f))
			return;

		// Quick reject if derivs are negative
		if (do_test)
		{
			if (test_den.empty())
			{
				if (test_num[v0] <= 0.0f &&
						test_num[v1] <= 0.0f &&
						test_num[v2] <= 0.0f)
					return;
			}
			else
			{
				if (test_num[v0] <= 0.0f && test_den[v0] >= 0.0f &&
						test_num[v1] <= 0.0f && test_den[v1] >= 0.0f &&
						test_num[v2] <= 0.0f && test_den[v2] >= 0.0f)
					return;
				if (test_num[v0] >= 0.0f && test_den[v0] <= 0.0f &&
						test_num[v1] >= 0.0f && test_den[v1] <= 0.0f &&
						test_num[v2] >= 0.0f && test_den[v2] <= 0.0f)
					return;
			}
		}

		// Figure out which val has different sign, and draw
		if ((val[v0] < 0.0f && val[v1] >= 0.0f && val[v2] >= 0.0f) ||
				(val[v0] > 0.0f && val[v1] <= 0.0f && val[v2] <= 0.0f))
			draw_face_isoline2(v0, v1, v2,
												 val, test_num, test_den,
												 do_hermite, do_test, fade);
		else if ((val[v1] < 0.0f && val[v2] >= 0.0f && val[v0] >= 0.0f) ||
						 (val[v1] > 0.0f && val[v2] <= 0.0f && val[v0] <= 0.0f))
			draw_face_isoline2(v1, v2, v0,
												 val, test_num, test_den,
												 do_hermite, do_test, fade);
		else if ((val[v2] < 0.0f && val[v0] >= 0.0f && val[v1] >= 0.0f) ||
						 (val[v2] > 0.0f && val[v0] <= 0.0f && val[v1] <= 0.0f))
			draw_face_isoline2(v2, v0, v1,
												 val, test_num, test_den,
												 do_hermite, do_test, fade);
	}

	// Takes a scalar field and renders the zero crossings, but only where test_num/test_den is greater than 0.

	void TriMesh::draw_isolines(const std::vector<float> &val,
										 const std::vector<float> &test_num,
										 const std::vector<float> &test_den,
										 bool do_bfcull, bool do_hermite,
										 bool do_test, float fade)
	{
		const int *t = &tstrips[0];
		const int *stripend = t;
		const int *end = t + tstrips.size();

		// Walk through triangle strips
		while (1)
		{
			if (unlikely(t >= stripend))
			{
				if (unlikely(t >= end))
					return;
				// New strip: each strip is stored as
				// length followed by indices
				stripend = t + 1 + *t;
				// Skip over length plus first two indices of
				// first face
				t += 3;
			}
			// Draw a line if, among the values in this triangle,
			// at least one is positive and one is negative
			const float &v0 = val[*t], &v1 = val[*(t - 1)], &v2 = val[*(t - 2)];
			if (unlikely((v0 > 0.0f || v1 > 0.0f || v2 > 0.0f) &&
									 (v0 < 0.0f || v1 < 0.0f || v2 < 0.0f)))
				draw_face_isoline(*(t - 2), *(t - 1), *t,
													val, test_num, test_den,
													do_bfcull, do_hermite, do_test, fade);
			t++;
		}
	}

	/******************************* 绘制特殊线条所需要的算法 end *****************************/

	// Suggestive contours
	// 暗示线
	std::vector<float> TriMesh::draw_suggestive()
	{
		std::vector<float>().swap(draw_points);

		float fade = draw_faded ? 0.03f / sqr(f_feature_size) : 0.0f;
		draw_isolines(kr, sctest_num, sctest_den,
									true, use_hermite, true, fade);
		return draw_points;
	}

	// 轮廓线 contours
	std::vector<float> TriMesh::draw_contours()
	{
		std::vector<float>().swap(draw_points);
		draw_isolines(ndotv, kr, std::vector<float>(),
									false, false, true, 0.0f);
		return draw_points;
	}

	// Compute largest eigenvalue and associated eigenstd::vector of a
	// symmetric 2x2 matrix.  Solves characteristic equation.
	// Inputs: three elements of matrix (upper-left, diag, lower-right)
	// Outputs: largest (in magnitude) eigenstd::vector/value
	void TriMesh::largest_eig_2x2(float m1, float m12, float m2, vec2 &e1, float &l1)
	{
		l1 = 0.5f * (m1 + m2);
		// The result of the below sqrt is positive, so to get the largest
		// eigenvalue we add it if we were positive already, else subtract
		if (l1 > 0.0f)
			l1 += sqrt(sqr(m12) + 0.25f * sqr(m2 - m1));
		else
			l1 -= sqrt(sqr(m12) + 0.25f * sqr(m2 - m1));

		// Find corresponding eigenstd::vector
		e1 = vec2(m2 - l1, -m12);
		normalize(e1);
	}

	// Compute principal view-dependent curvatures and directions at vertex i.
	// ndotv = cosine of angle between normal and view direction
	// (u,v) = coordinates of w (projected view) in principal coordinates
	// Pass in u^2, u*v, and v^2, since those are readily available.
	// Fills in q1 and t1 (using the paper's notation).
	// Note that the latter is expressed in the (pdir1,pdir2) coordinate basis
	void TriMesh::compute_viewdep_curv(int i, float f_ndotv,
																		 float u2, float uv, float v2,
																		 float &f_q1, vec2 &vec2_t1)
	{
		// Find the entries in Q = S * P^-1
		//                       = S + (sec theta - 1) * S * w * w^T
		float sectheta_minus1 = 1.0f / fabs(f_ndotv) - 1.0f;
		float Q11 = curv1[i] * (1.0f + sectheta_minus1 * u2);
		float Q12 = curv1[i] * (sectheta_minus1 * uv);
		float Q21 = curv2[i] * (sectheta_minus1 * uv);
		float Q22 = curv2[i] * (1.0f + sectheta_minus1 * v2);

		// Find the three entries in the (symmetric) matrix Q^T Q
		float QTQ1 = Q11 * Q11 + Q21 * Q21;
		float QTQ12 = Q11 * Q12 + Q21 * Q22;
		float QTQ2 = Q12 * Q12 + Q22 * Q22;

		// Compute eigenstuff
		largest_eig_2x2(QTQ1, QTQ12, QTQ2, vec2_t1, f_q1);
	}

	// Compute D_{t_1} q_1 - the derivative of max view-dependent curvature
	// in the principal max view-dependent curvature direction.
	void TriMesh::compute_Dt1q1(int i, float f_ndotv, float &f_Dt1q1)
	{
		const point &v0 = vertices[i];
		float this_viewdep_curv = q1[i];
		vec world_t1 = t1[i][0] * pdir1[i] + t1[i][1] * pdir2[i];
		vec world_t2 = normals[i] CROSS world_t1;
		float v0_dot_t2 = v0 DOT world_t2;

		f_Dt1q1 = 0.0f;
		int n = 0;

		int naf = adjacentfaces[i].size();
		for (int j = 0; j < naf; j++)
		{
			// We're in a triangle adjacent to the vertex of interest.
			// The current vertex is v0 - let v1 and v2 be the other two
			int f = adjacentfaces[i][j];
			int ind = faces[f].indexof(i);
			int i1 = faces[f][NEXT(ind)];
			int i2 = faces[f][PREV(ind)];
			const point &v1 = vertices[i1];
			const point &v2 = vertices[i2];

			// Find the point p on the segment between v1 and v2 such that
			// its std::vector from v0 is along t1, i.e. perpendicular to t2.
			// Linear combination: p = w1*v1 + w2*v2, where w2 = 1-w1
			float v1_dot_t2 = v1 DOT world_t2;
			float v2_dot_t2 = v2 DOT world_t2;
			float w1 = (v2_dot_t2 - v0_dot_t2) / (v2_dot_t2 - v1_dot_t2);

			// If w1 is not in [0..1) then we're not interested.
			// Incidentally, the computation of w1 can result in infinity,
			// but the comparison should do the right thing...
			if (w1 < 0.0f || w1 >= 1.0f)
				continue;

			// Construct the opposite point
			float w2 = 1.0f - w1;
			point p = w1 * v1 + w2 * v2;

			// And interpolate to find the view-dependent curvature at that point
			float interp_viewdep_curv = w1 * q1[i1] + w2 * q1[i2];

			// Finally, take the *projected* view-dependent curvature derivative
			float proj_dist = (p - v0) DOT world_t1;
			proj_dist *= fabs(f_ndotv);
			f_Dt1q1 += (interp_viewdep_curv - this_viewdep_curv) / proj_dist;
			n++;

			// To save time, quit as soon as we have two estimates
			// (that's all we're going to get, anyway)
			if (n == 2)
			{
				f_Dt1q1 *= 0.5f;
				return;
			}
		}
	}

	void TriMesh::init()
	{
		need_tstrips();
		need_bsphere();
		need_normals();
		need_curvatures();
		need_dcurv();
		f_feature_size = feature_size();
	};

	void TriMesh::compute_perview(bool extra_sin2theta)
	{
		if (draw_apparent)
			need_adjacentfaces();

		int nv = vertices.size();

		float scthresh = sug_thresh / sqr(f_feature_size);
		float shthresh = sh_thresh / sqr(f_feature_size);
		bool need_DwKr = (draw_sc || draw_sh || draw_DwKr);

		ndotv.resize(nv);
		kr.resize(nv);
		if (draw_apparent)
		{
			q1.resize(nv);
			t1.resize(nv);
			Dt1q1.resize(nv);
		}
		if (need_DwKr)
		{
			sctest_num.resize(nv);
			sctest_den.resize(nv);
			if (draw_sh)
				shtest_num.resize(nv);
		}

		for (int i = 0; i < nv; i++)
		{
			// Compute n DOT v
			vec viewdir = viewpos - vertices[i];
			float rlv = 1.0f / len(viewdir);
			viewdir *= rlv;
			ndotv[i] = viewdir DOT normals[i];

			float u = viewdir DOT pdir1[i], u2 = u * u;
			float v = viewdir DOT pdir2[i], v2 = v * v;

			// Note:  this is actually Kr * sin^2 theta
			kr[i] = curv1[i] * u2 + curv2[i] * v2;
			if (draw_apparent)
			{
				float csc2theta = 1.0f / (u2 + v2);
				compute_viewdep_curv(i, ndotv[i],
														 u2 * csc2theta, u * v * csc2theta, v2 * csc2theta,
														 q1[i], t1[i]);
			}
			if (!need_DwKr)
				continue;

			// Use DwKr * sin(theta) / cos(theta) for cutoff test
			sctest_num[i] = u2 * (u * dcurv[i][0] +
														3.0f * v * dcurv[i][1]) +
											v2 * (3.0f * u * dcurv[i][2] +
														v * dcurv[i][3]);
			float csc2theta = 1.0f / (u2 + v2);
			sctest_num[i] *= csc2theta;
			float tr = (curv2[i] - curv1[i]) *
								 u * v * csc2theta;
			sctest_num[i] -= 2.0f * ndotv[i] * sqr(tr);
			if (extra_sin2theta)
				sctest_num[i] *= u2 + v2;

			sctest_den[i] = ndotv[i];

			if (draw_sh)
			{
				shtest_num[i] = -sctest_num[i];
				shtest_num[i] -= shthresh * sctest_den[i];
			}
			sctest_num[i] -= scthresh * sctest_den[i];
		}
		if (draw_apparent)
		{
			for (int i = 0; i < nv; i++)
				compute_Dt1q1(i, ndotv[i], Dt1q1[i]);
		}
	}

	// Draw the ridges (valleys) of the mesh
	std::vector<float> TriMesh::draw_mesh_ridges(bool do_ridge,
																							 bool do_bfcull)
	{
		std::vector<float>().swap(draw_points);
		bool do_test = test_rv;
		std::cout << rv_thresh << '\n';
		float thresh = rv_thresh / f_feature_size;
		const int *t = &tstrips[0];
		const int *stripend = t;
		const int *end = t + tstrips.size();

		// Walk through triangle strips
		while (1)
		{
			if (unlikely(t >= stripend))
			{
				if (unlikely(t >= end))
					return draw_points;
				// New strip: each strip is stored as
				// length followed by indices
				stripend = t + 1 + *t;
				// Skip over length plus first two indices of
				// first face
				t += 3;
			}

			draw_face_ridges(*(t - 2), *(t - 1), *t,
											 do_ridge, do_bfcull, do_test, thresh);
			t++;
		}
		return draw_points;
	}

	// Draw ridges or valleys (depending on do_ridge) in a triangle v0,v1,v2
	// - uses ndotv for backface culling (enabled with do_bfcull)
	// - do_test checks for curvature maxima/minina for ridges/valleys
	//   (when off, it draws positive minima and negative maxima)
	// Note: this computes ridges/valleys every time, instead of once at the
	//   start (given they aren't view dependent, this is wasteful)
	// Algorithm based on formulas of Ohtake et al., 2004.
	void TriMesh::draw_face_ridges(int v0, int v1, int v2,
																 bool do_ridge,
																 bool do_bfcull, bool do_test, float thresh)
	{
		// Backface culling
		if (likely(do_bfcull &&
							 ndotv[v0] <= 0.0f && ndotv[v1] <= 0.0f && ndotv[v2] <= 0.0f))
			return;

		// Check if ridge possible at vertices just based on curvatures
		if (do_ridge)
		{
			if ((curv1[v0] <= 0.0f) ||
					(curv1[v1] <= 0.0f) ||
					(curv1[v2] <= 0.0f))
				return;
		}
		else
		{
			if ((curv1[v0] >= 0.0f) ||
					(curv1[v1] >= 0.0f) ||
					(curv1[v2] >= 0.0f))
				return;
		}

		// Sign of curvature on ridge/valley
		float rv_sign = do_ridge ? 1.0f : -1.0f;

		// The "tmax" are the principal directions of maximal curvature,
		// flipped to point in the direction in which the curvature
		// is increasing (decreasing for valleys).  Note that this
		// is a bit different from the notation in Ohtake et al.,
		// but the tests below are equivalent.
		const float &emax0 = dcurv[v0][0];
		const float &emax1 = dcurv[v1][0];
		const float &emax2 = dcurv[v2][0];
		vec tmax0 = rv_sign * dcurv[v0][0] * pdir1[v0];
		vec tmax1 = rv_sign * dcurv[v1][0] * pdir1[v1];
		vec tmax2 = rv_sign * dcurv[v2][0] * pdir1[v2];

		// We have a "zero crossing" if the tmaxes along an edge
		// point in opposite directions
		bool z01 = ((tmax0 DOT tmax1) <= 0.0f);
		bool z12 = ((tmax1 DOT tmax2) <= 0.0f);
		bool z20 = ((tmax2 DOT tmax0) <= 0.0f);

		if (z01 + z12 + z20 < 2)
			return;

		if (do_test)
		{
			const point &p0 = vertices[v0],
									&p1 = vertices[v1],
									&p2 = vertices[v2];

			// Check whether we have the correct flavor of extremum:
			// Is the curvature increasing along the edge?
			z01 = z01 && ((tmax0 DOT(p1 - p0)) >= 0.0f ||
										(tmax1 DOT(p1 - p0)) <= 0.0f);
			z12 = z12 && ((tmax1 DOT(p2 - p1)) >= 0.0f ||
										(tmax2 DOT(p2 - p1)) <= 0.0f);
			z20 = z20 && ((tmax2 DOT(p0 - p2)) >= 0.0f ||
										(tmax0 DOT(p0 - p2)) <= 0.0f);

			if (z01 + z12 + z20 < 2)
				return;
		}

		// Draw line segment
		const float &kmax0 = curv1[v0];
		const float &kmax1 = curv1[v1];
		const float &kmax2 = curv1[v2];
		if (!z01)
		{
			draw_segment_ridge(v1, v2, v0,
												 emax1, emax2, emax0,
												 kmax1, kmax2, kmax0,
												 thresh, false);
		}
		else if (!z12)
		{
			draw_segment_ridge(v2, v0, v1,
												 emax2, emax0, emax1,
												 kmax2, kmax0, kmax1,
												 thresh, false);
		}
		else if (!z20)
		{
			draw_segment_ridge(v0, v1, v2,
												 emax0, emax1, emax2,
												 kmax0, kmax1, kmax2,
												 thresh, false);
		}
		else
		{
			// All three edges have crossings -- connect all to center
			draw_segment_ridge(v1, v2, v0,
												 emax1, emax2, emax0,
												 kmax1, kmax2, kmax0,
												 thresh, true);
			draw_segment_ridge(v2, v0, v1,
												 emax2, emax0, emax1,
												 kmax2, kmax0, kmax1,
												 thresh, true);
			draw_segment_ridge(v0, v1, v2,
												 emax0, emax1, emax2,
												 kmax0, kmax1, kmax2,
												 thresh, true);
		}
	}

	// Draw part of a ridge/valley curve on one triangle face.  v0,v1,v2
	// are the indices of the 3 vertices; this function assumes that the
	// curve connects points on the edges v0-v1 and v1-v2
	// (or connects point on v0-v1 to center if to_center is true)
	void TriMesh::draw_segment_ridge(int v0, int v1, int v2,
																	 float emax0, float emax1, float emax2,
																	 float kmax0, float kmax1, float kmax2,
																	 float thresh, bool to_center)
	{
		// Interpolate to find ridge/valley line segment endpoints
		// in this triangle and the curvatures there
		float w10 = fabs(emax0) / (fabs(emax0) + fabs(emax1));
		float w01 = 1.0f - w10;
		point p01 = w01 * vertices[v0] + w10 * vertices[v1];
		float k01 = fabs(w01 * kmax0 + w10 * kmax1);

		point p12;
		float k12;
		if (to_center)
		{
			// Connect first point to center of triangle
			p12 = (vertices[v0] +
						 vertices[v1] +
						 vertices[v2]) /
						3.0f;
			k12 = fabs(kmax0 + kmax1 + kmax2) / 3.0f;
		}
		else
		{
			// Connect first point to second one (on next edge)
			float w21 = fabs(emax1) / (fabs(emax1) + fabs(emax2));
			float w12 = 1.0f - w21;
			p12 = w12 * vertices[v1] + w21 * vertices[v2];
			k12 = fabs(w12 * kmax1 + w21 * kmax2);
		}

		// Don't draw below threshold
		k01 -= thresh;
		if (k01 < 0.0f)
			k01 = 0.0f;
		k12 -= thresh;
		if (k12 < 0.0f)
			k12 = 0.0f;

		// Skip lines that you can't see...
		if (k01 == 0.0f && k12 == 0.0f)
			return;

		// Fade lines
		if (draw_faded)
		{
			k01 /= (k01 + thresh);
			k12 /= (k12 + thresh);
		}
		else
		{
			k01 = k12 = 1.0f;
		}

		// Draw the line segment

		draw_points.push_back(p01.x);
		draw_points.push_back(p01.y);
		draw_points.push_back(p01.z);
		draw_points.push_back(p12.x);
		draw_points.push_back(p12.y);
		draw_points.push_back(p12.z);
	}

	EMSCRIPTEN_BINDINGS(stl_wrappers)
	{
		register_vector<float>("vector<float>");
		value_array<std::array<float, 4>>("array_float_4")
				.element(emscripten::index<0>())
				.element(emscripten::index<1>())
				.element(emscripten::index<2>())
				.element(emscripten::index<3>());
	}

	EMSCRIPTEN_BINDINGS(class)
	{
		class_<TriMesh>("TriMesh")
				.function("setViewpos", &TriMesh::setViewpos)
				.function("get_scaleParam", &TriMesh::get_scaleParam)
				// .function("need_faces", &TriMesh::need_faces)
				// .function("need_curvatures", &TriMesh::need_curvatures)
				// .function("need_tstrips", &TriMesh::need_tstrips)
				// .function("need_adjacentfaces", &TriMesh::need_adjacentfaces)
				.function("init", &TriMesh::init)
				.function("compute_perview", &TriMesh::compute_perview)
				// .function("draw_boundaries", &TriMesh::draw_boundaries)
				.function("draw_mesh_ridges", &TriMesh::draw_mesh_ridges)
				.function("draw_suggestive", &TriMesh::draw_suggestive)
				.function("draw_contours", &TriMesh::draw_contours)
				.class_function("read", select_overload<TriMesh *(const std::string &)>(&TriMesh::read), allow_raw_pointers());
	};
};
