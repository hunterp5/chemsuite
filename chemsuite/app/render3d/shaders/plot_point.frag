#version 440
// chemsuite 3D plots: a round marker with an antialiased edge (GPL-3.0-or-later).
layout(location = 0) in vec2 v_corner;
layout(location = 1) in vec4 v_color;
layout(location = 2) in float v_pixels;

layout(std140, binding = 0) uniform Frame {
    mat4 mvp;         // projection (with the backend's clip-space correction) x view
    mat4 view;        // world to eye, for lighting the surface
    vec4 viewport;    // x, y: size in device pixels; z: device pixels per point
    vec4 light;       // eye-space direction towards the light
};

layout(location = 0) out vec4 fragColor;

void main()
{
    float r = length(v_corner) * v_pixels * 0.5;   // pixels from the centre
    float edge = clamp(v_pixels * 0.5 - r, 0.0, 1.0);
    if (edge <= 0.0)
        discard;
    fragColor = vec4(v_color.rgb, v_color.a * edge);
}
