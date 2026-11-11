#version 440
// chemsuite 3D plots: one point as a round marker of a fixed size on screen (GPL-3.0-or-later).
layout(location = 0) in vec2 corner;      // per vertex: (-1,-1) .. (1,1)
layout(location = 1) in vec4 point;       // per instance: position xyz, diameter in points
layout(location = 2) in vec4 color;       // per instance: rgba

layout(std140, binding = 0) uniform Frame {
    mat4 mvp;         // projection (with the backend's clip-space correction) x view
    mat4 view;        // world to eye, for lighting the surface
    vec4 viewport;    // x, y: size in device pixels; z: device pixels per point
    vec4 light;       // eye-space direction towards the light
};

layout(location = 0) out vec2 v_corner;
layout(location = 1) out vec4 v_color;
layout(location = 2) out float v_pixels;

void main()
{
    vec4 clip = mvp * vec4(point.xyz, 1.0);
    float pixels = point.w * viewport.z + 1.0;  // one pixel more, for the soft edge
    clip.xy += corner * pixels / viewport.xy * clip.w;
    v_corner = corner;
    v_color = color;
    v_pixels = pixels;
    gl_Position = clip;
}
