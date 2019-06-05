#pragma once
#include "Vector3.h"

class Texture {
public:
	Vector3 color, emission;
	Refl_t refl;
	double brdf;
	std::string filename;
	unsigned char *buf;
	int w, h, c;
	Texture(const Texture&t): brdf(t.brdf), filename(t.filename), color(t.color), emission(t.emission), refl(t.refl) {
		if (t.filename != "")
			buf = stbi_load(filename.c_str(), &w, &h, &c, 0);
		else
			buf = NULL;
	}
	Texture(std::string _, double b, Vector3 col, Vector3 e, Refl_t r): brdf(b), filename(_), color(col), emission(e), refl(r) {
		if(_ != "")
			buf = stbi_load(filename.c_str(), &w, &h, &c, 0);
		else
			buf = NULL;
	}
	std::pair<Refl_t, Vector3> getcol(double a, double b) {
		if (buf == NULL)
			return std::make_pair(refl, color);
		int pw = (int(a * w) % w + w) % w, ph = (int(b * h) % h + h) % h;
		int idx = ph * w * c + pw * c;
		int x = buf[idx + 0], y = buf[idx + 1], z = buf[idx + 2];
		// printf("find point %d %d %lf %lf\n", ph, pw,a,b);
		if (x == 233 && y == 233 && z == 233) {
			return std::make_pair(SPEC, Vector3(1, 1, 1)*.999);
		}
		return std::make_pair(refl, Vector3(x, y, z) / 255.);
	}
};
