#pragma once
#include "texture.h"
#include "scene.h"
#include "kdtree.h"

std::pair<int, double> find_intersect_simple(Ray ray) {
	double t = INF;
	int id = -1;
	for (int i = 0; i < scene_num; ++i) {
		std::pair<double, Vector3> tmp = scene[i]->intersect(ray);
		if (tmp.first < INF / 2 && tmp.second.len2() > eps && tmp.first < t) {
			t = tmp.first;
			id = i;
		}
	}
	return std::make_pair(id, t);
}

Vector3 basic_render(Ray ray, int dep, unsigned short *X) {
	// printf("Dep %d\n",dep);
	// ray.print();
	int into = 0;
	std::pair<int, double> intersect_result = find_intersect_simple(ray);
	if (intersect_result.first == -1)
		return Vector3();
	Object* obj = scene[intersect_result.first];
	Texture& texture = obj->texture;
	Vector3 x = ray.get(intersect_result.second);
	std::pair<Refl_t, Vector3> feature = get_feature(obj, texture, x, X);
	Vector3 f = feature.second, n = obj->norm(x), nl = n.dot(ray.d) < 0 ? into = 1, n : -n;
	double p = f.max();
	if (f.max() < eps)
		return texture.emission;
	if (++dep > 5)
		if (erand48(X) < p) f /= p;
		else return texture.emission;
	if (feature.first == DIFF) {
		double r1 = 2 * PI * erand48(X), r2 = erand48(X), r2s = sqrt(r2);
		Vector3 w = nl, u=((fabs(w.x) > .1 ? Vector3(0, 1) : Vector3(1)) & w).norm(), v = w & u;
		Vector3 d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();
		return f.mult(basic_render(Ray(x, d), dep, X));
	}
	else {
		Ray reflray = Ray(x, ray.d.reflect(nl));
		if (feature.first == SPEC) {
			return f.mult(basic_render(reflray, dep, X));
		}
		else {
			Vector3 d = ray.d.refract(n, into ? 1 : texture.brdf, into ? texture.brdf : 1);
			if (d.len2() < eps) // Total internal reflection
				return f.mult(basic_render(reflray, dep, X));
			double a = texture.brdf - 1, b = texture.brdf + 1;
			double R0 = a * a / (b * b), c = 1 - (into ? -ray.d.dot(nl) : d.dot(n));
			double Re = R0 + (1 - R0) * c * c  * c * c * c, Tr = 1 - Re;
			double P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);
			return f.mult(dep > 2 ? (erand48(X) < P ?   // Russian roulette
				basic_render(reflray, dep, X) * RP : basic_render(Ray(x, d), dep, X) * TP)
			  : basic_render(reflray, dep, X) * Re + basic_render(Ray(x, d), dep, X) * Tr);
		}
	}
}

std::vector<SPPMnode> sppm_backtrace(Ray ray, int dep, int index, unsigned short* X, Vector3 pref = Vector3(1, 1, 1), double prob = 1.) {
// if index == -1 then the node is illegal
	std::vector<SPPMnode> result, tmp;
	if (pref.max() < eps || prob < eps) return result;
	int into = 0;
	std::pair<int, double> intersect_result = find_intersect_simple(ray);
	if (intersect_result.first == -1)
		return result;
	Object* obj = scene[intersect_result.first];
	Texture& texture = obj->texture;
	Vector3 x = ray.get(intersect_result.second);
	std::pair<Refl_t, Vector3> feature = get_feature(obj, texture, x, X);
	Vector3 f = feature.second, n = obj->norm(x), nl = n.dot(ray.d) < 0 ? into = 1, n : -n;
	double p = f.max();
	// if (debug)
	// 	printf("dep = %d\tf = (%.5f, %.5f, %.5f) %s col = (%.5f %.5f %.5f) hit = (%.5f %.5f %.5f)\n", dep, pref.x, pref.y, pref.z, 
	// 		feature.first == REFR ? "REFR" : feature.first == DIFF ? "DIFF" : "SPEC", f.x, f.y, f.z, x.x, x.y, x.z);
	if (f.max() < eps)
		return result;
	if (++dep > 5)
		if(erand48(X) < p) f /= p;
		else return result;
	Ray reflray = Ray(x, ray.d.reflect(nl));
	// result.push_back(SPPMnode(x, pref.mult(f), nl, 1, index, prob));
	// return result;
	if (feature.first == DIFF || texture.filename == "vase.png") { // vase: 0.8 prob
		result.push_back(SPPMnode(x, pref.mult(f), nl, 1, index, prob * (texture.filename == "vase.png" ? .9 : 1)));
	}
	if (feature.first == SPEC || texture.filename == "vase.png") { // vase: 0.2 prob
		tmp = sppm_backtrace(reflray, dep, index, X, pref.mult(f), prob * (texture.filename == "vase.png" ? .1 : 1.));
		result.insert(result.end(), tmp.begin(), tmp.end());
	}
	if (feature.first == REFR) {
		Vector3 d = ray.d.refract(n, into ? 1 : texture.brdf, into ? texture.brdf : 1);
		if (d.len2() < eps) // Total internal reflection
			return sppm_backtrace(reflray, dep, index, X, pref.mult(f), prob);
		double a = texture.brdf - 1, b = texture.brdf + 1;
		double R0 = a * a / (b * b), c = 1 - (into ? -ray.d.dot(nl) : d.dot(n));
		double Re = R0 + (1 - R0) * c * c  * c * c * c, Tr = 1 - Re;
		double P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);
		if (dep > 2)
			if (erand48(X) < P) {
				tmp = sppm_backtrace(reflray, dep, index, X, pref.mult(f), prob * RP);
				result.insert(result.end(), tmp.begin(), tmp.end());
			}
			else {
				tmp = sppm_backtrace(Ray(x, d), dep, index, X, pref.mult(f), prob * TP);
				result.insert(result.end(), tmp.begin(), tmp.end());
			}
		else {
			tmp = sppm_backtrace(reflray, dep, index, X, pref.mult(f), prob * Re);
			result.insert(result.end(), tmp.begin(), tmp.end());
			tmp = sppm_backtrace(Ray(x, d), dep, index, X, pref.mult(f), prob * Tr);
			result.insert(result.end(), tmp.begin(), tmp.end());
		}
	}
	return result;
}

void sppm_forward(Ray ray, int dep, Vector3 col, unsigned short *X, IMGbuf* c, KDTree* kdt, double prob = 1.) {
	if (col.max() < eps) return;
	int into = 0;
	std::pair<int, double> intersect_result = find_intersect_simple(ray);
	if (intersect_result.first == -1)
		return;
	Object* obj = scene[intersect_result.first];
	Texture& texture = obj->texture;
	Vector3 x = ray.get(intersect_result.second);
	std::pair<Refl_t, Vector3> feature = get_feature(obj, texture, x, X);
	Vector3 f = feature.second, n = obj->norm(x), nl = n.dot(ray.d) < 0 ? into = 1, n : -n;
	double p = f.max();
	if (f.max() < eps) {
		kdt->query(SPPMnode(x, col, nl), c);
		return;
	}
	if (++dep > 5)
		if (erand48(X) < p) f /= p;
		else {
			kdt->query(SPPMnode(x, col, nl), c);
			return;
		}
	if (feature.first == DIFF) {
		kdt->query(SPPMnode(x, col, nl), c); // query col
		double r1 = 2 * PI * erand48(X), r2 = erand48(X), r2s = sqrt(r2);
		Vector3 w = nl, u=((fabs(w.x) > .1 ? Vector3(0, 1) : Vector3(1)) & w).norm(), v = w & u;
		Vector3 d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();
		return sppm_forward(Ray(x, d), dep, col.mult(f), X, c, kdt, prob);
	}
	else {
		Ray reflray = Ray(x, ray.d.reflect(nl));
		if (feature.first == SPEC) {
			if (texture.filename == "vase.png")
				kdt->query(SPPMnode(x, col, nl), c); // query col
			return sppm_forward(reflray, dep, col.mult(f), X, c, kdt, prob);
		}
		else {
			Vector3 d = ray.d.refract(n, into ? 1 : texture.brdf, into ? texture.brdf : 1);
			if (d.len2() < eps) // Total internal reflection
				return sppm_forward(reflray, dep, col.mult(f), X, c, kdt, prob);
			double a = texture.brdf - 1, b = texture.brdf + 1;
			double R0 = a * a / (b * b), c0 = 1 - (into ? -ray.d.dot(nl) : d.dot(n));
			double Re = R0 + (1 - R0) * c0 * c0  * c0 * c0 * c0, Tr = 1 - Re;
			double P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);
			return dep > 2 ? (erand48(X) < P ?   // Russian roulette
				sppm_forward(reflray, dep, col.mult(f), X, c, kdt, prob * RP) : sppm_forward(Ray(x, d), dep, col.mult(f), X, c, kdt, prob * TP))
			 : (sppm_forward(reflray, dep, col.mult(f), X, c, kdt, prob * Re),  sppm_forward(Ray(x, d), dep, col.mult(f), X, c, kdt, prob * Tr));
		}
	}
}
