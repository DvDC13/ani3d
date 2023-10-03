#pragma once

struct Parameters
{
    int N = 128;
    int scale = 8;
    int iter = 16;

    int window_width = N * scale;
    int window_height = N * scale;

    float dt = 0.0001f;
    float diff = 0.2f;
    float visc = 0.0f;
};

static Parameters parameters;