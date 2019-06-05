#pragma once

#include "utils.h"
#include "ray.h"

class BezierCurve2D {
public:
	double *dx, *dy, max, height, max2, r, num;
	int n;
	struct D{
		double t0, t1, width, y0, y1, width2;
	}data[20];
	// x(t) = \sum_{i=0}^n dx_i * t^i
	// y(t) = \sum_{i=0}^n dy_i * t^i
	BezierCurve2D(double* px, double* py, int n_, int num_, double r_): num(num_), n(n_), r(r_) {
		dx = new double[n];
		dy = new double[n];
		assert(std::abs(py[0]) <= 1e-6);
		--n;
		// preproces
		for(int i = 0; i <= n; ++i)
		{
			dx[i] = px[0];
			dy[i] = py[0];
			for (int j = 0; j <= n - i; ++j)
			{
				px[j] = px[j + 1] - px[j];
				py[j] = py[j + 1] - py[j];
			}
		}
		double n_down = 1, fac = 1, nxt = n;
		for (int i = 0; i <= n; ++i, --nxt)
		{
			fac = fac * (i == 0 ? 1 : i);
			dx[i] = dx[i] * n_down / fac;
			dy[i] = dy[i] * n_down / fac;
			n_down *= nxt;
		}
		max = 0;
		double interval = 1. / (num - 1), c = 0;
		for (int cnt = 0; cnt <= num; c += interval, ++cnt)
		{
			data[cnt].width = 0;
			data[cnt].t0 = std::max(0., c - r);
			data[cnt].t1 = std::min(1., c + r);
			data[cnt].y0 = getpos(data[cnt].t0).y;
			data[cnt].y1 = getpos(data[cnt].t1).y;
			for (double t = data[cnt].t0; t <= data[cnt].t1; t += 0.00001)
			{
				Vector3 pos = getpos(t);
				if (data[cnt].width < pos.x)
					data[cnt].width = pos.x;
			}
			if (max < data[cnt].width)
				max = data[cnt].width;
			data[cnt].width += eps;
			data[cnt].width2 = sqr(data[cnt].width);
		}
		max += eps;
		max2 = max * max;
		height = getpos(1).y;
	}
	Vector3 getpos(double t)
	{
		double ans_x = 0, ans_y = 0, t_pow = 1;
		for (int i = 0; i <= n; ++i)
		{
			ans_x += dx[i] * t_pow;
			ans_y += dy[i] * t_pow;
			t_pow *= t;
		}
		return Vector3(ans_x, ans_y);
	}
	Vector3 getdir(double t)
	{
		double ans_x = 0, ans_y = 0, t_pow = 1;
		for(int i = 1; i <= n; ++i)
		{
			ans_x += dx[i] * i * t_pow;
			ans_y += dy[i] * i * t_pow;
			t_pow *= t;
		}
		return Vector3(ans_x, ans_y);
	}
	Vector3 getdir2(double t)
	{
		double ans_x = 0, ans_y = 0, t_pow = 1;
		for(int i = 2; i <= n; ++i)
		{
			ans_x += dx[i] * i * (i - 1) * t_pow;
			ans_y += dy[i] * i * (i - 1) * t_pow;
			t_pow *= t;
		}
		return Vector3(ans_x, ans_y);
	}
};