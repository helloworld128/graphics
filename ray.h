#pragma once

#include "utils.h"
#include "Vector3.h"

class Ray
{
public:
	//origin of the ray.
	Vector3 o;
	//direction of the ray.
	Vector3 d;
	Ray(Vector3 o_, Vector3 d_): o(o_), d(d_) {}
	Vector3 get(double t) { return o + d * t; }
	void print() { puts("Ray: "); o.print(); d.print(); }
};