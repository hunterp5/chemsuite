#version 440
// chemsuite 3D plots: the surface lit from the eye's side, both faces alike (GPL-3.0-or-later).
layout(location = 0) in vec3 v_normal;
layout(location = 1) in vec4 v_color;

layout(std140, binding = 0) uniform Frame {
    mat4 mvp;         // projection (with the backend's clip-space correction) x view
    mat4 view;        // world to eye, for lighting the surface
    vec4 viewport;    // x, y: size in device pixels; z: device pixels per point
    vec4 light;       // eye-space direction towards the light
};

layout(location = 0) out vec4 fragColor;

void main()
{
    vec3 n = normalize(v_normal);
    float diffuse = abs(dot(n, normalize(light.xyz)));
    fragColor = vec4(v_color.rgb * (0.45 + 0.55 * diffuse), v_color.a);
}
