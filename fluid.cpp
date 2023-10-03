#include "fluid.hpp"

FluidParticle::FluidParticle(int size, float dt, float diff, float visc)
    : size_(size), dt_(dt), diff_(diff), visc_(visc)
{
    int N = size_ * size_;

    s_ = new float[N];
    density_ = new float[N];

    Vx_ = new float[N];
    Vy_ = new float[N];

    Vx0_ = new float[N];
    Vy0_ = new float[N];

    std::fill(s_, s_ + N, 0.0f);
    std::fill(density_, density_ + N, 0.0f);

    std::fill(Vx_, Vx_ + N, 0.0f);
    std::fill(Vy_, Vy_ + N, 0.0f);

    std::fill(Vx0_, Vx0_ + N, 0.0f);
    std::fill(Vy0_, Vy0_ + N, 0.0f);
}

FluidParticle::~FluidParticle()
{
    delete[] s_;
    delete[] density_;

    delete[] Vx_;
    delete[] Vy_;

    delete[] Vx0_;
    delete[] Vy0_;
}

void FluidParticle::addDensity(int x, int y, float amount)
{
    density_[INDEX(x, y)] += amount;
}

void FluidParticle::addVelocity(int x, int y, float amountX, float amountY)
{
    int index = INDEX(x, y);

    Vx_[index] += amountX;
    Vy_[index] += amountY;
}

void FluidParticle::fadeDensity()
{
    for (int i = 0; i < size_ * size_; ++i)
    {
        density_[i] /= 1.02f;
    }
}

void FluidParticle::set_bnd(int b, float* x)
{
    for (int i = 1; i < size_ - 1; i++)
    {
        x[INDEX(i, 0)] = b == 2 ? -x[INDEX(i, 1)] : x[INDEX(i, 1)];
        x[INDEX(i, size_ - 1)] = b == 2 ? -x[INDEX(i, size_ - 2)] : x[INDEX(i, size_ - 2)];
    }

    for (int j = 1; j < size_ - 1; j++)
    {
        x[INDEX(0, j)] = b == 1 ? -x[INDEX(1, j)] : x[INDEX(1, j)];
        x[INDEX(size_ - 1, j)] = b == 1 ? -x[INDEX(size_ - 2, j)] : x[INDEX(size_ - 2, j)];
    }

    x[INDEX(0, 0)] = 0.5f * (x[INDEX(1, 0)] + x[INDEX(0, 1)]);
    x[INDEX(0, size_ - 1)] = 0.5f * (x[INDEX(1, size_ - 1)] + x[INDEX(0, size_ - 2)]);
    x[INDEX(size_ - 1, 0)] = 0.5f * (x[INDEX(size_ - 2, 0)] + x[INDEX(size_ - 1, 1)]);
    x[INDEX(size_ - 1, size_ - 1)] = 0.33f * (x[INDEX(size_ - 2, size_ - 1)] + x[INDEX(size_ - 1, size_ - 2)]);
}

void FluidParticle::lin_solve(int b, float* x, float* x0, float a, float c)
{
    float cRecip = 1.0f / c;
    for (int k = 0; k < parameters.iter; k++)
    {
        for (int j = 1; j < size_ - 1; j++)
        {
            for (int i = 1; i < size_ - 1; i++)
            {
                x[INDEX(i, j)] =
                    (x0[INDEX(i, j)]
                        + a * (x[INDEX(i + 1, j)]
                            + x[INDEX(i - 1, j)]
                            + x[INDEX(i, j + 1)]
                            + x[INDEX(i, j - 1)])) * cRecip;
            }
        }
        set_bnd(b, x);
    }
}

void FluidParticle::diffuse(int b, float* x, float* x0, float diff, float dt)
{
    float a = dt * diff * (size_ - 2) * (size_ - 2);
    lin_solve(b, x, x0, a, 1 + 6 * a);
}

void FluidParticle::project(float* velocX, float* velocY, float* p, float* div)
{
    for (int j = 1; j < size_ - 1; j++)
    {
        for (int i = 1; i < size_ - 1; i++)
        {
            div[INDEX(i, j)] = -0.5f * (
                velocX[INDEX(i + 1, j)]
                - velocX[INDEX(i - 1, j)]
                + velocY[INDEX(i, j + 1)]
                - velocY[INDEX(i, j - 1)]) / size_;
            p[INDEX(i, j)] = 0;
        }
    }

    set_bnd(0, div);
    set_bnd(0, p);
    lin_solve(0, p, div, 1, 6);

    for (int j = 1; j < size_ - 1; j++)
    {
        for (int i = 1; i < size_ - 1; i++)
        {
            velocX[INDEX(i, j)] -= 0.5f * (p[INDEX(i + 1, j)]
                - p[INDEX(i - 1, j)]) * size_;
            velocY[INDEX(i, j)] -= 0.5f * (p[INDEX(i, j + 1)]
                - p[INDEX(i, j - 1)]) * size_;
        }
    }

    set_bnd(1, velocX);
    set_bnd(2, velocY);
}

void FluidParticle::advect(int b, float* d, float* d0, float* velocX, float* velocY, float dt)
{
    float i0, i1, j0, j1;

    float dtx = dt * (size_ - 2);
    float dty = dt * (size_ - 2);

    float s0, s1, t0, t1;
    float tmp1, tmp2, x, y;

    float Nfloat = size_;
    float ifloat, jfloat;
    int i, j;

    for (j = 1, jfloat = 1; j < size_ - 1; j++, jfloat++)
    {
        for (i = 1, ifloat = 1; i < size_ - 1; i++, ifloat++)
        {
            tmp1 = dtx * velocX[INDEX(i, j)];
            tmp2 = dty * velocY[INDEX(i, j)];

            x = ifloat - tmp1;
            y = jfloat - tmp2;

            x = std::clamp(x, 0.5f, Nfloat + 0.5f);

            i0 = floorf(x);
            i1 = i0 + 1.0f;

            y = std::clamp(y, 0.5f, Nfloat + 0.5f);
            
            j0 = floorf(y);
            j1 = j0 + 1.0f;

            s1 = x - i0;
            s0 = 1.0f - s1;
            t1 = y - j0;
            t0 = 1.0f - t1;

            int i0i = i0;
            int i1i = i1;
            int j0i = j0;
            int j1i = j1;

            d[INDEX(i, j)] = s0
                    * (t0 * d0[INDEX(i0i, j0i)]
                       + t1 * d0[INDEX(i0i, j1i)])
                + s1
                    * (t0 * d0[INDEX(i1i, j0i)]
                       + t1 * d0[INDEX(i1i, j1i)]);
        }
    }

    set_bnd(b, d);
}

void FluidParticle::step()
{
    float visc = visc_;
    float diff = diff_;
    float dt = dt_;
    float* Vx = Vx_;
    float* Vy = Vy_;
    float* Vx0 = Vx0_;
    float* Vy0 = Vy0_;
    float* s = s_;
    float* density = density_;

    diffuse(1, Vx0, Vx, visc, dt);
    diffuse(2, Vy0, Vy, visc, dt);

    project(Vx0, Vy0, Vx, Vy);

    advect(1, Vx, Vx0, Vx0, Vy0, dt);
    advect(2, Vy, Vy0, Vx0, Vy0, dt);

    project(Vx, Vy, Vx0, Vy0);

    diffuse(0, s, density, diff, dt);
    advect(0, density, s, Vx, Vy, dt);
}