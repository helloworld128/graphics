#pragma once
#include "utils.h"

struct Vector3 {
	double x, y, z;
	Vector3(double x_ = 0, double y_ = 0, double z_ = 0) : x(x_), y(y_), z(z_) {}
	Vector3 operator-() const { return Vector3(-x, -y, -z); }
	Vector3 operator+(const Vector3&a) const { return Vector3(x + a.x, y + a.y, z + a.z); }
	Vector3 operator-(const Vector3&a) const { return Vector3(x - a.x, y - a.y, z - a.z); }
	Vector3 operator+(double p) const { return Vector3(x + p, y + p, z + p); }
	Vector3 operator-(double p) const { return Vector3(x - p, y - p, z - p); }
	Vector3 operator*(double p) const { return Vector3(x*p, y*p, z*p); }
	Vector3 operator/(double p) const { return Vector3(x / p, y / p, z / p); }
	bool operator==(const Vector3&a) const { return x == a.x && y == a.y && z == a.z; }
	bool operator!=(const Vector3&a) const { return x != a.x || y != a.y || z != a.z; }
	Vector3&operator+=(const Vector3&a) { return *this = *this + a; }
	Vector3&operator-=(const Vector3&a) { return *this = *this - a; }
	Vector3&operator+=(double p) { return *this = *this + p; }
	Vector3&operator-=(double p) { return *this = *this - p; }
	Vector3&operator*=(double p) { return *this = *this * p; }
	Vector3&operator/=(double p) { return *this = *this / p; }
	double operator|(const Vector3&a) const { return x*a.x + y*a.y + z*a.z; }
	double dot(const Vector3&a) const { return x*a.x + y*a.y + z*a.z; }
	double max() const { return x > y&&x > z ? x : y > z ? y : z; }
	Vector3 max(const Vector3&a) const { return Vector3(std::max(x, a.x), std::max(y, a.y), std::max(z, a.z)); }
	Vector3 min(const Vector3&a) const { return Vector3(std::min(x, a.x), std::min(y, a.y), std::min(z, a.z)); }
	double len() const { return sqrt(x*x + y*y + z*z); }
	double len2() const { return x*x + y*y + z*z; }
	Vector3 mult(const Vector3&a) const { return Vector3(x*a.x, y*a.y, z*a.z); }
	Vector3 operator&(const Vector3&a) const { return Vector3(y*a.z - z*a.y, z*a.x - x*a.z, x*a.y - y*a.x); }
	Vector3 cross(const Vector3&a) const { return Vector3(y*a.z - z*a.y, z*a.x - x*a.z, x*a.y - y*a.x); }
	Vector3 norm() const { return (*this) / len(); }
	Vector3 clip(double r0 = 0, double r1 = 1) const { return Vector3(x > r1 ? r1 : x<r0 ? r0 : x, y>r1 ? r1 : y<r0 ? r0 : y, z>r1 ? r1 : z < r0 ? r0 : z); }
	Vector3 reflect(const Vector3&n) const { return (*this) - n*2.*n.dot(*this); }
	Vector3 refract(const Vector3&n, double ni, double nr) const {
		double cosi = this->norm().dot(n);
		double nir = ni / nr;
		double cosr2 = 1. - nir*nir*(1 - cosi*cosi);
		if (cosr2 <= 0)
			return Vector3();
		double cosr = sqrt(cosr2);
		if (cosi > 0) // out
			cosr = -cosr;
		return ((*this)*nir - n*(nir*cosi + cosr)).norm();
	}
	void print() const { std::cout << x << " " << y << " " << z << std::endl; }
};

Vector3 min(Vector3 a, Vector3 b) { return Vector3(std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z)); }
Vector3 max(Vector3 a, Vector3 b) { return Vector3(std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z)); }