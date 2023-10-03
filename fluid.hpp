#pragma once

#include <algorithm>
#include <cmath>

#include "parameters.hpp"

inline int INDEX(int x, int y)
{
    x = std::clamp(x, 0, parameters.N - 1);
    y = std::clamp(y, 0, parameters.N - 1);
    return x + y * parameters.N;
}

class FluidParticle
{
public:
    FluidParticle(int size, float dt, float diff, float visc);
    ~FluidParticle();

    void addDensity(int x, int y, float amount);
    void addVelocity(int x, int y, float amountX, float amountY);
    void fadeDensity();

    void set_bnd(int b, float *x);
    void lin_solve(int b, float *x, float *x0, float a, float c);
    void diffuse(int b, float *x, float *x0, float diff, float dt);
    void project(float* velocX, float* velocY, float* p, float* div);
    void advect(int b, float* d, float* d0, float* velocX, float* velocY, float dt);

    void step();

    inline int getSize() const { return size_; }
    inline float getDt() const { return dt_; }
    inline float getDiff() const { return diff_; }
    inline float getVisc() const { return visc_; }

    inline float* getDensity() const { return density_; }
    inline float* getVx() const { return Vx_; }
    inline float* getVy() const { return Vy_; }

    inline float* getVx0() const { return Vx0_; }
    inline float* getVy0() const { return Vy0_; }
private:
    int size_;
	float dt_;
	float diff_;
	float visc_;

	float *s_;
	float *density_;

	float *Vx_;
	float *Vy_;

	float *Vx0_;
	float *Vy0_;
};