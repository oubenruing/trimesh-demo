#ifndef TRIMESH_H
#define TRIMESH_H
/*
Szymon Rusinkiewicz
Princeton University

TriMesh.h
Class for triangle meshes.
*/

#include "Vec.h"
#include "Box.h"
#include "Color.h"
#include "strutil.h"
#include <vector>
#include <array>
#include <iostream>
#include <math.h>

namespace trimesh
{

	template <class T>
	static inline void clear_and_release(::std::vector<T> &v)
	{
		// Standard trick to release a vector's storage, since clear() doesn't
		::std::vector<T>().swap(v);
	}

	class TriMesh
	{
	public:
		//
		// Types
		//

		typedef Vec<3, int> Face;
		typedef Box<3, float> BBox;

		struct BSphere
		{
			point center;
			float r;
			bool valid;
			BSphere() : valid(false)
			{
			}
		};

		//
		// Enums
		//
		enum TstripRep
		{
			TSTRIP_LENGTH,
			TSTRIP_TERM
		};
		enum
		{
			GRID_INVALID = -1
		};
		enum StatOp
		{
			STAT_MIN,
			STAT_MINABS,
			STAT_MAX,
			STAT_MAXABS,
			STAT_SUM,
			STAT_SUMABS,
			STAT_SUMSQR,
			STAT_MEAN,
			STAT_MEANABS,
			STAT_RMS,
			STAT_MEDIAN,
			STAT_STDEV
		};
		enum StatVal
		{
			STAT_VALENCE,
			STAT_FACEAREA,
			STAT_ANGLE,
			STAT_DIHEDRAL,
			STAT_EDGELEN,
			STAT_X,
			STAT_Y,
			STAT_Z
		};

		//
		// Constructor
		//
		TriMesh() : grid_width(-1), grid_height(-1), flag_curr(0)
		{
			std::cout << "TriMesh"
								<< "\n";
		}

		//
		// Members
		//
		// The basics: vertices and faces
		::std::vector<point> vertices;
		::std::vector<Face> faces;

		// Triangle strips
		::std::vector<int> tstrips;

		// Grid, if present
		::std::vector<int> grid;
		int grid_width, grid_height;

		// Other per-vertex properties
		::std::vector<Color> colors;
		::std::vector<float> confidences;
		::std::vector<unsigned> flags;
		unsigned flag_curr;

		// Computed per-vertex properties
		::std::vector<vec> normals;
		::std::vector<vec> pdir1, pdir2;
		::std::vector<float> curv1, curv2;
		::std::vector<float> get_curv1()
		{
			return curv1;
		}
		::std::vector<Vec<4, float>> dcurv;
		::std::vector<vec> cornerareas;
		::std::vector<float> pointareas;

		// Bounding structures
		BBox bbox;
		BSphere bsphere;

		// Connectivity structures:
		//  For each vertex, all neighboring vertices
		::std::vector<::std::vector<int>> neighbors;
		//  For each vertex, all neighboring faces
		::std::vector<::std::vector<int>> adjacentfaces;
		//  For each face, the three faces attached to its edges
		//  (for example, across_edge[3][2] is the number of the face
		//   that's touching the edge opposite vertex 2 of face 3)
		::std::vector<Face> across_edge;

		::std::vector<float> ndotv, kr;
		::std::vector<float> sctest_num, sctest_den, shtest_num;
		::std::vector<float> q1, Dt1q1;
		::std::vector<vec2> t1;

		// compute members

		// Two cameras: the primary one, and an alternate one to fix the lines
		// and see them from a different direction
		float fov = 0.7f;
		double alt_projmatrix[16];
		point viewpos; // Current view position
		std::vector<float> draw_points;
		std::array<float, 4> scaleParam{};
		std::array<float, 4> get_scaleParam()
		{
			return scaleParam;
		}

		// Toggles for drawing various lines
		int draw_extsil = 0, draw_c = 1, draw_sc = 1;
		int draw_sh = 0, draw_phridges = 0, draw_phvalleys = 0;
		int draw_ridges = 0, draw_valleys = 0, draw_apparent = 0;
		int draw_K = 0, draw_H = 0, draw_DwKr = 0;
		int draw_bdy = 0, draw_isoph = 0, draw_topo = 0;
		int niso = 20, ntopo = 20;
		float topo_offset = 0.0f;

		// Toggles for tests we perform
		int draw_hidden = 0;
		int test_c = 1, test_sc = 1, test_sh = 1, test_ph = 1, test_rv = 1, test_ar = 1;
		float sug_thresh = 0.03, sh_thresh = 0.02, ph_thresh = 0.04;
		float rv_thresh = 0.03, ar_thresh = 0.1;

		// Toggles for style
		int use_texture = 0;
		int draw_faded = 1;
		int draw_colors = 0;
		int use_hermite = 0;
		int draw_edges = false;

		float f_feature_size; // Used to make thresholds dimensionless
		float currsmooth;			// Used in smoothing

		//
		// Compute all this stuff...
		//
		// void draw_line();

		void compute_viewdep_curv(int i, float f_ndotv,
															float u2, float uv, float v2,
															float &f_q1, trimesh::vec2 &vec2_t1);

		// Compute D_{t_1} q_1 - the derivative of max view-dependent curvature
		// in the principal max view-dependent curvature direction.
		void compute_Dt1q1(int i, float f_ndotv, float &f_Dt1q1);

		void compute_perview(bool extra_sin2theta = false);
		void largest_eig_2x2(float m1, float m12, float m2, vec2 &e1, float &l1);

		// Draw the ridges (valleys) of the mesh
		// 绘制峰谷图形
		std::vector<float> draw_mesh_ridges(bool do_ridge,
																				bool do_bfcull);

		void draw_face_ridges(int v0, int v1, int v2,
													bool do_ridge,
													bool do_bfcull, bool do_test, float thresh);

		void draw_segment_ridge(int v0, int v1, int v2,
														float emax0, float emax1, float emax2,
														float kmax0, float kmax1, float kmax2,
														float thresh, bool to_center);

		// 绘制边界
		std::vector<float> draw_boundaries(bool do_hidden);

		// 计算（kr * sin^2 theta）在顶点i的梯度
	 inline vec gradkr(int i);

		// 通过线性插值在val0和val1之间找到一个过零点。
		// 如果过零点在val0，则返回0，如果在val1，则返回1,之间差值。
	 inline float find_zero_linear(float val0, float val1);

		// 用Hermite插值法寻找零点交叉点
		float find_zero_hermite(int v0, int v1, float val0, float val1,
														const vec &grad0, const vec &grad1);

		// 绘制面上的分割线
		void draw_face_isoline(int v0, int v1, int v2,
													 const std::vector<float> &val,
													 const std::vector<float> &test_num,
													 const std::vector<float> &test_den,
													 bool do_bfcull, bool do_hermite,
													 bool do_test, float fade);

		void draw_face_isoline2(int v0, int v1, int v2,
														const std::vector<float> &val,
														const std::vector<float> &test_num,
														const std::vector<float> &test_den,
														bool do_hermite, bool do_test, float fade);

		// 绘制特殊分割线
		void draw_isolines(const std::vector<float> &val,
											 const std::vector<float> &test_num,
											 const std::vector<float> &test_den,
											 bool do_bfcull, bool do_hermite,
											 bool do_test, float fade);

		// 绘制暗示线
		std::vector<float> draw_suggestive();
		// 绘制轮廓线
		std::vector<float> draw_contours();

		void setViewpos(float x, float y, float z)
		{
			viewpos.x = x;
			viewpos.y = y;
			viewpos.z = z;
		};
		void need_faces()
		{
			if (!faces.empty())
				return;
			if (!tstrips.empty())
				unpack_tstrips();
			else if (!grid.empty())
				triangulate_grid();
		}
		void need_tstrips(TstripRep rep = TSTRIP_LENGTH);
		void convert_strips(TstripRep rep);
		void unpack_tstrips();
		void resize_grid(int width, int height)
		{
			grid_width = width;
			grid_height = height;
			grid.clear();
			grid.resize(grid_width * grid_height, GRID_INVALID);
		}
		void triangulate_grid(bool remove_slivers = true);
		void need_normals(bool simple_area_weighted = false);
		void need_curvatures();
		void need_dcurv();
		void need_pointareas();
		void need_bbox();
		void need_bsphere();
		void need_neighbors();
		void need_adjacentfaces();
		void need_across_edge();

		//
		// Delete everything and release storage
		//
		void clear_vertices() { clear_and_release(vertices); }
		void clear_faces() { clear_and_release(faces); }
		void clear_tstrips() { clear_and_release(tstrips); }
		void clear_grid()
		{
			clear_and_release(grid);
			grid_width = grid_height = -1;
		}
		void clear_colors() { clear_and_release(colors); }
		void clear_confidences() { clear_and_release(confidences); }
		void clear_flags()
		{
			clear_and_release(flags);
			flag_curr = 0;
		}
		void clear_normals() { clear_and_release(normals); }
		void clear_curvatures()
		{
			clear_and_release(pdir1);
			clear_and_release(pdir2);
			clear_and_release(curv1);
			clear_and_release(curv2);
		}
		void clear_dcurv() { clear_and_release(dcurv); }
		void clear_pointareas()
		{
			clear_and_release(pointareas);
			clear_and_release(cornerareas);
		}
		void clear_bbox() { bbox.clear(); }
		void clear_bsphere() { bsphere.valid = false; }
		void clear_neighbors() { clear_and_release(neighbors); }
		void clear_adjacentfaces() { clear_and_release(adjacentfaces); }
		void clear_across_edge() { clear_and_release(across_edge); }
		void clear()
		{
			clear_vertices();
			clear_faces();
			clear_tstrips();
			clear_grid();
			clear_colors();
			clear_confidences();
			clear_flags();
			clear_normals();
			clear_curvatures();
			clear_dcurv();
			clear_pointareas();
			clear_bbox();
			clear_bsphere();
			clear_neighbors();
			clear_adjacentfaces();
			clear_across_edge();
		}

		//
		// Input and output
		//
	protected:
		static bool read_helper(const char *filename, TriMesh *mesh);

	public:
		static TriMesh *read(const char *filename);
		static TriMesh *read(const ::std::string &filename);
		bool write(const char *filename);
		bool write(const ::std::string &filename);
		// static TriMesh *read_mem(const int & addr, const size_t & len);

		//
		// Useful queries
		//

		// Is vertex v on the mesh boundary?
		inline bool is_bdy(int v)
		{
			if (unlikely(neighbors.empty()))
				need_neighbors();
			if (unlikely(adjacentfaces.empty()))
				need_adjacentfaces();
			return neighbors[v].size() != adjacentfaces[v].size();
		}

		// Centroid of face f
		inline vec centroid(int f)
		{
			if (unlikely(faces.empty()))
				need_faces();
			return (1.0f / 3.0f) *
						 (vertices[faces[f][0]] +
							vertices[faces[f][1]] +
							vertices[faces[f][2]]);
		}

		// Normal of face f
		inline vec trinorm(int f)
		{
			if (unlikely(faces.empty()))
				need_faces();
			return trimesh::trinorm(vertices[faces[f][0]],
															vertices[faces[f][1]], vertices[faces[f][2]]);
		}

		// Angle of corner j in triangle i
		inline float cornerangle(int i, int j)
		{
			using namespace ::std;

			if (unlikely(faces.empty()))
				need_faces();
			const point &p0 = vertices[faces[i][j]];
			const point &p1 = vertices[faces[i][NEXT_MOD3(j)]];
			const point &p2 = vertices[faces[i][PREV_MOD3(j)]];
			return acos((p1 - p0) DOT(p2 - p0));
		}

		// Dihedral angle between face i and face across_edge[i][j]
		inline float dihedral(int i, int j)
		{
			if (unlikely(across_edge.empty()))
				need_across_edge();
			if (unlikely(across_edge[i][j] < 0))
				return 0.0f;
			vec mynorm = trinorm(i);
			vec othernorm = trinorm(across_edge[i][j]);
			float ang = angle(mynorm, othernorm);
			vec towards = 0.5f * (vertices[faces[i][NEXT_MOD3(j)]] +
														vertices[faces[i][PREV_MOD3(j)]]) -
										vertices[faces[i][j]];
			if ((towards DOT othernorm) < 0.0f)
				return M_PIf + ang;
			else
				return M_PIf - ang;
		}

		// Statistics
		float stat(StatOp op, StatVal val);
		void init();
		float feature_size();

		//
		// Debugging
		//

		// Debugging printout, controllable by a "verbose"ness parameter
		static int verbose;
		static void set_verbose(int);
		static void (*dprintf_hook)(const char *);
		static void set_dprintf_hook(void (*hook)(const char *));
		static void dprintf(const char *format, ...);

		// Same as above, but fatal-error printout
		static void (*eprintf_hook)(const char *);
		static void set_eprintf_hook(void (*hook)(const char *));
		static void eprintf(const char *format, ...);
	};

} // namespace trimesh

#endif
