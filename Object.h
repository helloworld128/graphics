#pragma once

#include "utils.h"
#include "bezier.h"

unsigned short mess[3] = { 1, 2, 3 };

class Object {
public:
	Texture texture;
	Object(Texture t) : texture(t) {}
	Object(Refl_t refl, Vector3 color, Vector3 emission, double brdf, std::string tname) :
		texture(tname, brdf, color, emission, refl) {}
	virtual std::pair<double, Vector3> intersect(Ray) { puts("virtual error in intersect!"); return std::make_pair(0.0, Vector3()); }
	// If no intersect, then return (INF, (0,0,0))
	virtual std::pair<Vector3, Vector3> aabb() { puts("virtual error in aabb!"); return std::make_pair(Vector3(), Vector3()); }
	virtual Vector3 norm(Vector3) { puts("virtual error in norm!"); return Vector3(); }
	// return norm vec out of obj
	virtual Vector3 change_for_bezier(Vector3) { puts("virtual error in bezier!"); return Vector3(); }
};

class BezierObject : public Object {
	// the curve will rotate line (x=pos.x and z=pos.z) as pivot
public:
	BezierCurve2D curve;
	Vector3 pos; // the buttom center point
	BezierObject(Vector3 pos_, BezierCurve2D c_, Texture t) :
		pos(pos_), curve(c_), Object(t) {}
	BezierObject(Vector3 pos_, BezierCurve2D c_, Refl_t refl, double brdf = 1.5, Vector3 color = Vector3(), Vector3 emission = Vector3(), std::string tname = "") :
		pos(pos_), curve(c_), Object(refl, color, emission, brdf, tname) {}
	double solve_t(double yc) { // solve y(t)=yc
						// assert(0 <= yc && yc <= curve.height);
		double t = .5, ft, dft;
		for (int i = 10; i--; )
		{
			if (t < 0) t = 0;
			else if (t > 1) t = 1;
			ft = curve.getpos(t).y - yc, dft = curve.getdir(t).y;
			if (std::abs(ft) < eps)
				return t;
			t -= ft / dft;
		}
		return -1;
	}
	virtual Vector3 change_for_bezier(Vector3 inter_p) {
		double t = solve_t(inter_p.y - pos.y);
		double u = atan2(inter_p.z - pos.z, inter_p.x - pos.x); // between -pi ~ pi
		if (u < 0)
			u += 2 * PI;
		return Vector3(u, t);
	}
	double get_sphere_intersect(Ray ray, Vector3 o, double r) {
		Vector3 ro = o - ray.o;
		double b = ray.d.dot(ro);
		double d = sqr(b) - ro.dot(ro) + sqr(r);
		if (d < 0) return -1;
		else d = sqrt(d);
		double t = b - d > eps ? b - d : b + d > eps ? b + d : -1;
		if (t < 0)
			return -1;
		return t;
	}
	virtual std::pair<double, Vector3> intersect(Ray ray) {
		double final_dis = INF;
		// check for |dy|<eps
		if (std::abs(ray.d.y) < 5e-4)
		{
			double dis_to_axis = (Vector3(pos.x, ray.o.y, pos.z) - ray.o).len();
			double hit = ray.get(dis_to_axis).y;
			if (hit < pos.y + eps || hit > pos.y + curve.height - eps)
				return std::make_pair(INF, Vector3());
			// solve function pos.y+y(t)=ray.o.y to get x(t)
			double t = solve_t(hit - pos.y);
			if (t < 0 || t > 1)
				return std::make_pair(INF, Vector3());
			Vector3 loc = curve.getpos(t);
			double ft = pos.y + loc.y - hit;
			if (std::abs(ft) > eps)
				return std::make_pair(INF, Vector3());
			// assume sphere (pos.x, pos.y + loc.y, pos.z) - loc.x
			final_dis = get_sphere_intersect(ray, Vector3(pos.x, pos.y + loc.y, pos.z), loc.x);
			if (final_dis < 0)
				return std::make_pair(INF, Vector3());
			Vector3 inter_p = ray.get(final_dis);
			// printf("y %f small!!!",std::abs((inter_p - Vector3(pos.x, inter_p.y, pos.z)).len2() - sqr(loc.x)));
			if (std::abs((inter_p - Vector3(pos.x, inter_p.y, pos.z)).len2() - sqr(loc.x)) > 1e-1)
				return std::make_pair(INF, Vector3());
			// second iteration, more accuracy
			hit = inter_p.y;
			if (hit < pos.y + eps || hit > pos.y + curve.height - eps)
				return std::make_pair(INF, Vector3());
			t = solve_t(hit - pos.y);
			loc = curve.getpos(t);
			ft = pos.y + loc.y - hit;
			if (std::abs(ft) > eps)
				return std::make_pair(INF, Vector3());
			final_dis = get_sphere_intersect(ray, Vector3(pos.x, pos.y + loc.y, pos.z), loc.x);
			if (final_dis < 0)
				return std::make_pair(INF, Vector3());
			inter_p = ray.get(final_dis);
			if (std::abs((inter_p - Vector3(pos.x, hit, pos.z)).len2() - sqr(loc.x)) > 1e-2)
				return std::make_pair(INF, Vector3());
			// printf("---y %f small!!!",std::abs((inter_p - Vector3(pos.x, inter_p.y, pos.z)).len2() - sqr(loc.x)));
			return std::make_pair(final_dis, inter_p);
		}
		// printf("y big\n");
		// check for top circle: the plane is y=pos.y + curve.height
		// TODO
		// check for buttom circle: the plane is y=pos.y
		// TODO
		// normal case
		// calc ay^2+by+c
		double a = 0, b = 0, c = 0, t1, t2;
		// (xo-x'+xd/yd*(y-yo))^2 -> (t1+t2*y)^2
		t1 = ray.o.x - pos.x - ray.d.x / ray.d.y * ray.o.y;
		t2 = ray.d.x / ray.d.y;
		a += t2 * t2;
		b += 2 * t1 * t2;
		c += t1 * t1;
		// (zo-z'+zd/yd*(y-yo))^2 -> (t1+t2*y)^2
		t1 = ray.o.z - pos.z - ray.d.z / ray.d.y * ray.o.y;
		t2 = ray.d.z / ray.d.y;
		a += sqr(t2);
		b += 2 * t1 * t2;
		c += sqr(t1);
		// ay^2+by+c -> a'(y-b')^2+c'
		c = c - b * b / 4 / a;
		b = -b / 2 / a - pos.y;
		// printf("%lf %lf %lf\n",a,b,c);
		if (0 <= b && b <= curve.height && c > curve.max2
			|| (b < 0 || b > curve.height) && std::min(sqr(b), sqr(curve.height - b)) * a + c > curve.max2) // no intersect
			return std::make_pair(INF, Vector3());
		// double pick[20] = {0, 0, 1}; int tot = 2;
		// for (double _ = 0; _ <= 1; _ += 0.1)
		// {
		// 	double t_pick = newton2(_, a, b, c);
		// 	if (0 <= t_pick && t_pick <= 1)
		// 	{
		// 		bool flag = 1;
		// 		for (int j = 1; j <= tot; ++j)
		// 			if (std::abs(t_pick - pick[j]) < eps)
		// 				flag = 0;
		// 		if (flag)
		// 			pick[++tot] = t_pick;
		// 	}
		// }
		// std::sort(pick + 1, pick + 1 + tot);
		// for (int j = 1; j < tot; ++j)
		// 	if (getft(pick[j], a, b, c) * getft(pick[j + 1], a, b, c) <= 0)
		// 		check(pick[j], pick[j+1], (pick[j] + pick[j + 1]) * .5, ray, a, b, c, final_dis);
		for (int ind = 0; ind <= curve.num; ++ind)
		{
			// y = curve.ckpt[ind] ~ curve.ckpt[ind+1]
			// calc min(a(y-b)^2+c)
			// double lower;
			// if (curve.data[ind].y0 <= b && b <= curve.data[ind].y1)
			// 	lower = c;
			// else
			// 	lower = a * std::min(sqr(curve.data[ind].y0 - b), sqr(curve.data[ind].y1 - b)) + c;
			double t0 = curve.data[ind].t0, t1 = curve.data[ind].t1;
			// if (t0 > eps) t0 += erand48(mess) * .01;
			// if (t1 < 1 - eps) t1 -= erand48(mess) * .01;
			// if (lower <= curve.data[ind].width2)
			{
				check(t0, t1, (t0 + t1 + t0) / 3, ray, a, b, c, final_dis);
				check(t0, t1, (t1 + t0 + t1) / 3, ray, a, b, c, final_dis);
			}
		}
		if (final_dis < INF / 2)
			return std::make_pair(final_dis, ray.get(final_dis));
		else
			return std::make_pair(INF, Vector3());
	}
	bool check(double low, double upp, double init, Ray ray, double a, double b, double c, double&final_dis)
	{
		double t = newton(init, a, b, c, low, upp);
		if (t <= 0 || t >= 1)
			return false;
		Vector3 loc = curve.getpos(t);
		double x = loc.x, y = loc.y;
		double ft = x - sqrt(a * sqr(y - b) + c);
		if (std::abs(ft) > eps)
			return false;
		// calc t for ray
		double dis = (pos.y + y - ray.o.y) / ray.d.y;
		if (dis < eps)
			return false;
		Vector3 inter_p = ray.get(dis);
		if (std::abs((Vector3(pos.x, pos.y + y, pos.z) - inter_p).len2() - x * x) > eps)
			return false;
		if (dis < final_dis)
		{
			final_dis = dis;
			// printf("%lf %lf %lf %lf\n",t,x , sqrt(a * sqr(y - b) + c), x - sqrt(a * sqr(y - b) + c));
			return true;
		}
		return false;
	}
	double getft(double t, double a, double b, double c)
	{
		if (t < 0) t = eps;
		if (t > 1) t = 1 - eps;
		Vector3 loc = curve.getpos(t);
		double x = loc.x, y = loc.y;
		return x - sqrt(a * sqr(y - b) + c);
	}
	double newton(double t, double a, double b, double c, double low = eps, double upp = 1 - eps)
	{
		// solve sqrt(a(y(t)+pos.y-b)^2+c)=x(t)
		// f(t) = x(t) - sqrt(a(y(t)+pos.y-b)^2+c)
		// f'(t) = x'(t) - a(y(t)+pos.y-b)*y'(t) / sqrt(...)
		// if t is not in [0, 1] then assume f(t) is a linear function
		double ft, dft, x, y, dx, dy, sq;
		Vector3 loc, dir;
		for (int i = 10; i--; )
		{
			if (t < 0) t = low;
			if (t > 1) t = upp;
			loc = curve.getpos(t), dir = curve.getdir(t);
			x = loc.x, dx = dir.x;
			y = loc.y, dy = dir.y;
			// printf("%lf %lf %lf\n",t,x,y);
			sq = sqrt(a * sqr(y - b) + c);
			ft = x - sq;
			dft = dx - a * (y - b) * dy / sq;
			if (std::abs(ft) < eps)
				return t;
			t -= ft / dft;
		}
		return -1;
	}
	double newton2(double t, double a, double b, double c)
	{
		double dft, ddft, y, dx, dy, ddx, ddy, sq;
		Vector3 loc, dir, dir2;
		for (int i = 5; i--; )
		{
			if (t < 0) t = eps;
			if (t > 1) t = 1 - eps;
			loc = curve.getpos(t), dir = curve.getdir(t), dir2 = curve.getdir2(t);
			y = loc.y, dx = dir.x, dy = dir.y;
			ddx = dir2.x, ddy = dir2.y;
			sq = sqrt(a * sqr(y - b) + c);
			dft = dx - a * (y - b) * dy / sq;
			ddft = ddx - a * ((y - b) * ddy + sqr(dy)) / sq + sqr(a * (y - b) * dy) / sq / sq / sq;
			if (std::abs(dft) < eps)
				return t;
			t -= dft / ddft;
		}
		return -1;
	}
	virtual std::pair<Vector3, Vector3> aabb() {
		return std::make_pair(Vector3(pos.x - curve.max, pos.y, pos.z - curve.max), Vector3(pos.x + curve.max, pos.y + curve.height, pos.z + curve.max));
	}
	virtual Vector3 norm(Vector3 p) {
		Vector3 tmp = change_for_bezier(p);
		Vector3 dir = curve.getdir(tmp.y);
		Vector3 d_surface = Vector3(cos(tmp.x), dir.y / dir.x, sin(tmp.x));
		Vector3 d_circ = Vector3(-sin(tmp.x), 0, cos(tmp.x));
		return d_circ.cross(d_surface).norm();
	}
};

class CubeObject : public Object {
	//store (x0, y0, z0) - (x1, y1, z1)
public:
	Vector3 m0, m1;
	CubeObject(Vector3 m0_, Vector3 m1_, Texture t) :
		m0(min(m0_, m1_)), m1(max(m0_, m1_)), Object(t) {}
	CubeObject(Vector3 m0_, Vector3 m1_, Refl_t refl, double brdf = 1.5, Vector3 color = Vector3(), Vector3 emission = Vector3(), std::string tname = "") :
		m0(min(m0_, m1_)), m1(max(m0_, m1_)), Object(refl, color, emission, brdf, tname) {}
	virtual Vector3 norm(Vector3 p) {
		if (std::abs(p.x - m0.x) < eps || std::abs(p.x - m1.x) < eps)
			return Vector3(std::abs(p.x - m1.x) < eps ? 1 : -1, 0, 0);
		if (std::abs(p.y - m0.y) < eps || std::abs(p.y - m1.y) < eps)
			return Vector3(0, std::abs(p.y - m1.y) < eps ? 1 : -1, 0);
		if (std::abs(p.z - m0.z) < eps || std::abs(p.z - m1.z) < eps)
			return Vector3(0, 0, std::abs(p.z - m1.z) < eps ? 1 : -1);
		assert(1 == 0);
	}
	virtual std::pair<double, Vector3> intersect(Ray ray) {
		double ft = INF, t;
		Vector3 fq = Vector3(), q;
		// x dir
		t = (m0.x - ray.o.x) / ray.d.x;
		if (0 < t && t < ft) {
			q = ray.get(t);
			if (m0.y <= q.y && q.y <= m1.y && m0.z <= q.z && q.z <= m1.z)
				ft = t, fq = q;
		}
		t = (m1.x - ray.o.x) / ray.d.x;
		if (0 < t && t < ft) {
			q = ray.get(t);
			if (m0.y <= q.y && q.y <= m1.y && m0.z <= q.z && q.z <= m1.z)
				ft = t, fq = q;
		}
		// y dir
		t = (m0.y - ray.o.y) / ray.d.y;
		if (0 < t && t < ft) {
			q = ray.get(t);
			if (m0.x <= q.x && q.x <= m1.x && m0.z <= q.z && q.z <= m1.z)
				ft = t, fq = q;
		}
		t = (m1.y - ray.o.y) / ray.d.y;
		if (0 < t && t < ft) {
			q = ray.get(t);
			if (m0.x <= q.x && q.x <= m1.x && m0.z <= q.z && q.z <= m1.z)
				ft = t, fq = q;
		}
		// z dir
		t = (m0.z - ray.o.z) / ray.d.z;
		if (0 < t && t < ft) {
			q = ray.get(t);
			if (m0.x <= q.x && q.x <= m1.x && m0.y <= q.y && q.y <= m1.y)
				ft = t, fq = q;
		}
		t = (m1.z - ray.o.z) / ray.d.z;
		if (0 < t && t < ft) {
			q = ray.get(t);
			if (m0.x <= q.x && q.x <= m1.x && m0.y <= q.y && q.y <= m1.y)
				ft = t, fq = q;
		}
		return std::make_pair(ft, fq);
	}
	virtual std::pair<Vector3, Vector3> aabb() {
		return std::make_pair(m0, m1);
	}
};

class SphereObject : public Object {
public:
	Vector3 o;
	double r;
	SphereObject(Vector3 o_, double r_, Texture t) :
		o(o_), r(r_), Object(t) {}
	SphereObject(Vector3 o_, double r_, Refl_t refl, double brdf = 1.5, Vector3 color = Vector3(), Vector3 emission = Vector3(), std::string tname = "") :
		o(o_), r(r_), Object(refl, color, emission, brdf, tname) {}
	virtual std::pair<double, Vector3> intersect(Ray ray) {
		Vector3 ro = o - ray.o;
		double b = ray.d.dot(ro);
		double d = sqr(b) - ro.dot(ro) + sqr(r);
		if (d < 0) return std::make_pair(INF, Vector3());
		else d = sqrt(d);
		double t = b - d > eps ? b - d : b + d > eps ? b + d : -1;
		if (t < 0)
			return std::make_pair(INF, Vector3());
		return std::make_pair(t, ray.get(t));
	}
	virtual std::pair<Vector3, Vector3> aabb() {
		return std::make_pair(o - r, o + r);
	}
	virtual Vector3 norm(Vector3 p) {
		double d = std::abs((p - o).len() - r);
		assert(d < eps);
		return (p - o).norm();
	}
};

class PlaneObject : public Object {
	// store ax+by+cz=1 n=(a,b,c)
public:
	Vector3 n, n0;
	PlaneObject(Vector3 n_, Texture t) :
		n(n_), n0(n_.norm()), Object(t) {}
	PlaneObject(Vector3 n_, Refl_t refl, double brdf = 1.5, Vector3 color = Vector3(), Vector3 emission = Vector3(), std::string tname = "") :
		n(n_), n0(n_.norm()), Object(refl, color, emission, brdf, tname) {}
	virtual std::pair<double, Vector3> intersect(Ray ray) {
		double t = (1 - ray.o.dot(n)) / ray.d.dot(n);
		if (t < eps)
			return std::make_pair(INF, Vector3());
		return std::make_pair(t, ray.get(t));
	}
	virtual std::pair<Vector3, Vector3> aabb() {
		Vector3 p0 = Vector3(min_p[0], min_p[1], min_p[2]);
		Vector3 p1 = Vector3(max_p[0], max_p[1], max_p[2]);
		if (std::abs(n.x) <= eps && std::abs(n.y) <= eps) { // horizontal plane
			p0.z = 1. / n.z - eps;
			p1.z = 1. / n.z + eps;
			return std::make_pair(p0, p1);
		}
		if (std::abs(n.y) <= eps && std::abs(n.z) <= eps) { // verticle plane
			p0.x = 1. / n.x - eps;
			p1.x = 1. / n.x + eps;
			return std::make_pair(p0, p1);
		}
		if (std::abs(n.x) <= eps && std::abs(n.z) <= eps) { // verticle plane
			p0.y = 1. / n.y - eps;
			p1.y = 1. / n.y + eps;
			return std::make_pair(p0, p1);
		}
		return std::make_pair(p0, p1);
	}
	virtual Vector3 norm(Vector3) {
		return n0;
	}
};
