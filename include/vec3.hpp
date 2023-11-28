#ifndef vec3_hpp
#define vec3_hpp

#include <iostream>
#include "vec2.hpp"

class vec3;
vec3 V3_Normalised(const vec3& v);
vec3 V3_Cross(const vec3& v1, const vec3& v2);
vec3 V3_RandomUnit(const int &divisions=1000);
vec3 operator|(const vec2 &v, const float &f);

class vec3 {
public:
	friend vec3 operator+(const vec3& v1, const vec3& v2) {
		return {v1.x + v2.x, v1.y + v2.y, v1.z + v2.z};
	}
	friend vec3 operator-(const vec3& v1, const vec3& v2) {
		return {v1.x - v2.x, v1.y - v2.y, v1.z - v2.z};
	}
	friend float operator*(const vec3& v1, const vec3& v2){
		return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
	}
	friend vec3 operator*(const vec3& v, const float& a){
		return {a*v.x, a*v.y, a*v.z};
	}
	friend vec3 operator*(const float& a, const vec3& v){
		return {a*v.x, a*v.y, a*v.z};
	}
	friend vec3 operator/(const vec3& v, const float& a){
		return {v.x/a, v.y/a, v.z/a};
	}
	vec3& operator+=(const vec3& v){
		*this = *this + v;
		return *this;
	}
	vec3& operator-=(const vec3& v){
		*this = *this - v;
		return *this;
	}
	vec3& operator*=(const float& a){
		*this = *this * a;
		return *this;
	}
	vec3& operator/=(const float& a){
		*this = *this / a;
		return *this;
	}
	friend std::ostream& operator<<(std::ostream& stream, const vec3& v){
		stream << "(" << v.x << ", " << v.y << ", " << v.z << ")";
		return stream;
	}
	friend vec3 operator|(const vec2 &v, const float &f){
		return {v.x, v.y, f};
	}
	
	vec3& Normalise(){
		*this /= sqrt(SqMag());
		return *this;
	}
	
	friend vec3 V3_Normalised(const vec3& v){
		return v / sqrt(v.SqMag());
	}
	
	friend vec3 V3_Cross(const vec3& v1, const vec3& v2){
		return {v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x};
	}
	
	float SqMag() const {
		return x*x + y*y + z*z;
	}
	
	friend vec3 V3_RandomUnit(const int &divisions){
		const float phi = 2.0f*M_PI*(float)(rand() % divisions)/(float)divisions;
		const float costheta = 2.0f*(float)(rand() % divisions)/(float)divisions - 1.0f;
		const float theta = acosf(costheta);
		return {sinf(theta)*cosf(phi), sinf(theta)*sinf(phi), cosf(theta)};
	}
	
	float x, y, z;
};

#endif /* vec3_hpp */
