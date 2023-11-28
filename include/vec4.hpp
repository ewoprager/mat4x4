#ifndef vec4_hpp
#define vec4_hpp

#include "vec3.hpp"

class vec4;
vec4 V4_Normalised(const vec4& v);
vec4 operator|(const vec3 &v, const float &f);

class vec4 {
public:
	friend vec4 operator+(const vec4& v1, const vec4& v2) {
		return {v1.x + v2.x, v1.y + v2.y, v1.z + v2.z, v1.w + v2.w};
	}
	friend vec4 operator-(const vec4& v1, const vec4& v2) {
		return {v1.x - v2.x, v1.y - v2.y, v1.z - v2.z, v1.w - v2.w};
	}
	friend float operator*(const vec4& v1, const vec4& v2){
		return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z + v1.w*v2.w;
	}
	friend vec4 operator*(const vec4& v, const float& a){
		return {a*v.x, a*v.y, a*v.z, a*v.w};
	}
	friend vec4 operator*(const float& a, const vec4& v){
		return {a*v.x, a*v.y, a*v.z, a*v.w};
	}
	friend vec4 operator/(const vec4& v, const float& a){
		return {v.x/a, v.y/a, v.z/a, v.w/a};
	}
	vec4& operator+=(const vec4& v){
		*this = *this + v;
		return *this;
	}
	vec4& operator-=(const vec4& v){
		*this = *this - v;
		return *this;
	}
	vec4& operator*=(const float& a){
		*this = *this * a;
		return *this;
	}
	vec4& operator/=(const float& a){
		*this = *this / a;
		return *this;
	}
	friend std::ostream& operator<<(std::ostream& stream, const vec4& v){
		stream << "(" << v.x << ", " << v.y << ", " << v.z << ", " << v.w << ")";
		return stream;
	}
	float& operator[](const int& index){
		if(index == 0) return x;
		if(index == 1) return y;
		if(index == 2) return z;
		if(index == 3) return w;
		std::cout << "ERROR: vec4 index out of range." << std::endl;
		return x;
	}
	
	vec4& Normalise(){
		*this /= sqrt(SqMag());
		return *this;
	}
	friend vec4 V4_Normalised(const vec4& v){
		return v / sqrt(v.SqMag());
	}
	friend vec4 operator|(const vec3 &v, const float &f){
		return {v.x, v.y, v.z, f};
	}
	
	float SqMag() const {
		return x*x + y*y + z*z + w*w;
	}
	
	vec3 xyz() const {
		return {x, y, z};
	}
	
	float x, y, z, w;
};

#endif /* vec4_hpp */
