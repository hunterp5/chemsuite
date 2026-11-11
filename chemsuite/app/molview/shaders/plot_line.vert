#version 440
// chemsuite 3D plots: the cube, floor grid and drop lines (GPL-3.0-or-later).
layout(location = 0) in vec3 position;
layout(location = 1) in vec4 color;

layout(std140, binding = 0) uniform Frame {
    mat4 mvp;         // projection (with the backend's clip-space correction) x view
    mat4 view;        // world to eye, for lighting the surface
    vec4 viewport;    // x, y: size in device pixels; z: device pixels per point
    vec4 light;       // eye-space direction towards the light
};

layout(location = 0) out vec4 v_color;

void main()
{
    v_color = color;
    gl_Position = mvp * vec4(position, 1.0);
}
