//============================================================
// FILE: vector3d_t.h
// Name: David Tu
// Program Description: 3D Vector Class Implementation. This class is able to construct 3D Vectors
// as well as be able to perfrom mathematical computations on 3D Vectors
//	Input: 3D Vectors, scalars
//	Output: 3D Vectors (After computation)
//============================================================

#pragma once
#ifndef __vector3d_T_H__
#define __vector3d_T_H__

#include <iostream>
#include <string>
#include <cmath>
#include <initializer_list>

template <typename T> class vector3d;
template <typename T> std::ostream&  operator<<(std::ostream& os, const vector3d<T>& v);

typedef vector3d<double> vector3dD;
typedef vector3d<float> vector3dF;
typedef vector3d<int> vector3dI;
typedef vector3d<long> vector3dL;

template <typename T>
class vector3d{
public:
	vector3d();
	vector3d(const std::string& name, int dims);
	vector3d(const std::string& name, int dims, const std::initializer_list<T>& li);
//----------------------------------------------------------------------
	T operator[](int i) const;
	T& operator[](int i);
//----------------------------------------------------------------------
	void name(const std::string& name);
	const std::string& name() const;
//----------------------------------------------------------------------
	vector3d<T>& operator+=(const vector3d<T>& v);
	vector3d<T>& operator-=(const vector3d<T>& v);
//----------------------------------------------------------------------
	vector3d<T>& operator+=(T k);
	vector3d<T>& operator-=(T k);
	vector3d<T>& operator*=(T k);
	vector3d<T>& operator/=(T k);
//----------------------------------------------------------------------
	vector3d<T> operator-();
	vector3d<T> operator+(const vector3d<T>& v);
	vector3d<T> operator-(const vector3d<T>& v);
//----------------------------------------------------------------------
	friend vector3d operator+(T k, const vector3d& v){
		return vector3d(std::to_string(k) + "+" + v.name_, v.dims_, {k + v[0], k + v[1], k + v[2], 0});
	}
	friend vector3d operator+(const vector3d& v, T k){
		return  k + v;
	}
	friend vector3d operator-(const vector3d& v, T k){
		return -k + v;
	}
	friend vector3d operator-(T k, const vector3d& v){
		return vector3d(std::to_string(k) + "-" + v.name_, v.dims_, {k - v[0], k - v[1], k - v[2], 0});
	}
	friend vector3d operator*(T k, const vector3d& v){
		return vector3d(std::to_string(k) + "*" + v.name_, v.dims_, {k * v[0], k * v[1], k * v[2], 0});
	}
	friend vector3d operator*(const vector3d& v, T k){
		return k * v;
	}
	friend vector3d operator/(const vector3d& v, T k){
		return vector3d(v.name_ + "/" + std::to_string(k), v.dims_, {v[0]/k, v[1]/k, v[2]/k, 0});
	}
//----------------------------------------------------------------------
	bool operator==(const vector3d<T>& v) const;
	bool operator!=(const vector3d<T>& v) const;
//----------------------------------------------------------------------
	T dot(const vector3d<T>& v) const;
	T magnitude() const;
	T angle(const vector3d<T>& v) const;
	vector3d<T> cross(const vector3d<T>& v) const;
//----------------------------------------------------------------------
	static vector3d<T> zero();
//----------------------------------------------------------------------
	friend std::ostream& operator<< <>(std::ostream& os, const vector3d<T>& v);
private:
	void check_equal_dims(const vector3d<T>& v) const;
	void check_bounds(int i) const;
private:
	constexpr static double EPSILON = 1.0e-10;
	std::string name_;
	int dims_;
	T data_[4];
};

//----------------------------------------------------------------------
	template <typename T> vector3d<T>::vector3d() : vector3d("", 3){
		// 3d default dims
	}
	template <typename T> vector3d<T>::vector3d(const std::string& name, int dims): name_(name), dims_(dims){
		std::memset(data_, 0, dims_ * sizeof(T));
		data_[3] = T();  // vectors have 0 at end, pts have 1
	}
	template <typename T> vector3d<T>::vector3d(const std::string& name, int dims, const std::initializer_list<T>& li): vector3d(name, dims){
		int i = 0;
		for (T value : li){
			if (i > dims_){ 
				break; 
			}        
			data_[i++] = value;    
		}    
		data_[3] = T(); 
	}
//----------------------------------------------------------------------
	template <typename T> T vector3d<T>::operator[](int i) const{// read-only index operator
		check_bounds(i);    
		return data_[i]; 
	}
	template <typename T> T& vector3d<T>::operator[](int i){// read-write index operator
		check_bounds(i);
		return data_[i];
	}
//----------------------------------------------------------------------
	template <typename T> void vector3d<T>::name(const std::string& name){ 
		name_ = name; 
	}
	template <typename T> const std::string& vector3d<T>::name() const{ 
		return name_; 
	}
//----------------------------------------------------------------------
	template <typename T> vector3d<T>& vector3d<T>::operator+=(const vector3d<T>& v){
		vector3d<T>& u = *this;
		for (int i = 0; i < 3; ++i){
			u[i] += v[i];
		}
		return *this;
	}
	template <typename T> vector3d<T>& vector3d<T>::operator-=(const vector3d<T>& v){
		vector3d<T>& u = *this;
		for (int i = 0; i < 3; ++i) {
			u[i] -= v[i];
		}
		return *this;
	}
//----------------------------------------------------------------------
	template <typename T> vector3d<T>& vector3d<T>::operator+=(T k){
		vector3d<T>& u = *this;
		for (int i = 0; i < 3; ++i) {
			u[i] += k;
		}
		return *this;
	}
	template <typename T> vector3d<T>& vector3d<T>::operator*=(T k){
		vector3d<T>& u = *this;
		for (int i = 0; i < 3; ++i) {
			u[i] *= k;
		}
		return *this;
	}
	template <typename T> vector3d<T>& vector3d<T>::operator-=(T k){
		vector3d<T>& u = *this;
		for (int i = 0; i < 3; ++i) {
			u[i] -= k;
		}
		return *this;
	}
	template <typename T> vector3d<T>& vector3d<T>::operator/=(T k){
		vector3d<T>& u = *this;
		for (int i = 0; i < 3; ++i) {
			u[i] /= k;
		}
		return *this;
	}
//----------------------------------------------------------------------
	template <typename T> vector3d<T> vector3d<T>::operator-(){
		return vector3d<T>("-" + name_, dims_, {-data_[0], -data_[1], -data_[2], 0});
	}
	template <typename T> vector3d<T> vector3d<T>::operator+(const vector3d& v){
		const vector3d<T>& u = *this;
		check_equal_dims(v);
		return vector3d<T>(u.name_ + "+" + v.name_, dims_, {u[0] + v[0], u[1] + v[1], u[2] + v[2], 0});
	}
	template <typename T> vector3d<T> vector3d<T>::operator-(const vector3d<T>& v){
		check_equal_dims(v);
		return vector3d<T>(name_ + "-" + v.name_, dims_, {data_[0] - v[0], data_[1] -v[1], data_[2] - v[2], 0});
	}
//----------------------------------------------------------------------
	template <typename T> bool vector3d<T>::operator==(const vector3d<T>& v) const{
		const vector3d<T>& u = *this;    
		check_equal_dims(v);    
		return  std::abs(u[0] - v[0]) < vector3d<T>::EPSILON && std::abs(u[1] - v[1]) < vector3d<T>::EPSILON && std::abs(u[2] - v[2]) < vector3d<T>::EPSILON;
	}
	template <typename T> bool vector3d<T>::operator!=(const vector3d<T>& v) const{ 
		return !(*this == v); 
	}
//----------------------------------------------------------------------
	template <typename T> T vector3d<T>::dot(const vector3d<T>& v) const{
		check_equal_dims(v);
		T dot = 0;
		for (int i = 0; i < dims_; ++i) {
			dot += data_[i] * v.data_[i];
		}
		return dot;
	}
	template <typename T> T vector3d<T>::magnitude() const{
		return sqrt(dot(*this));
	}
	template <typename T> T vector3d<T>::angle(const vector3d<T>& v) const{
		T dot = this->dot(v);
		T mag = this->magnitude();
		T othermag = v.magnitude();
		return acos(dot / (mag * othermag));
	}
	template <typename T> vector3d<T> vector3d<T>::cross(const vector3d<T>& v) const{
		const vector3d<T>& u = *this;
		check_equal_dims(v);
		if (v.dims_ != 3){
			throw new std::invalid_argument("cross_product only implemented for vector3d's");
		}
		return vector3d(name_ + " x " + v.name_, dims_, {u[1]*v[2] - u[2]*v[1], -(u[0]*v[2] - u[2]*v[0]), u[0]*v[1] - u[1]*v[0], 0 });
	}
//----------------------------------------------------------------------
	template <typename T> vector3d<T> vector3d<T>::zero(){ 
		return vector3d("zero", 3, { 0, 0, 0, 0 }); 
	}
//----------------------------------------------------------------------
	template <typename T> std::ostream& operator<<(std::ostream& os, const vector3d<T>& v){ 
		os << "<'" << v.name_ << "', ";    
		if (v.dims_ == 0){ 
			os << "empty>";
		} else{ 
			for (int i = 0; i < v.dims_ + 1; ++i){ 
				os << v[i];            
				if (i < v.dims_){ 
					os << " "; 
				} 
			}        
			os << ">"; 
		}    
		return os; 
	}
//----------------------------------------------------------------------
	template <typename T> void vector3d<T>::check_equal_dims(const vector3d<T>& v) const{
		if (dims_ != v.dims_){
			throw new std::invalid_argument("vector3d dims mismatch");
		}
	}
	template <typename T> void vector3d<T>::check_bounds(int i) const{
		if (i > dims_) { // 1 extra dimension for pts/vectors
			throw new std::invalid_argument("out of bounds");
		}
	}

#endif

//============================================================
// end of file: vector3d_t.h
//============================================================