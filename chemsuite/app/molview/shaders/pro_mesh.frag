#version 440
// chemsuite protein view: mesh shading with slab, selection and hover (GPL-3.0-or-later).
layout(location = 0) in vec3 v_pos;
layout(location = 1) in vec3 v_normal;
layout(location = 2) in vec4 v_color;
layout(location = 3) in float v_flags;

layout(std140, binding = 0) uniform Frame {
    mat4 view;
    mat4 proj;
    vec4 light;
    vec4 params;
    vec4 background;
    vec4 slab;
    vec4 highlight;
};

layout(location = 0) out vec4 fragColor;

void main()
{
    if (-v_pos.z < slab.x || -v_pos.z > slab.y)
        discard;
    vec3 dir = params.w > 0.5 ? vec3(0.0, 0.0, -1.0) : normalize(v_pos);
    vec3 n = normalize(v_normal);
    if (dot(n, dir) > 0.0)
        n = -n;                                // two-sided: the inside of a cut surface
    vec3 l = normalize(light.xyz);
    float diffuse = max(dot(n, l), 0.0);
    float spec = pow(max(dot(reflect(-l, n), -dir), 0.0), 32.0);
    float rim = pow(1.0 - max(dot(n, -dir), 0.0), 2.0);
    float flags = floor(v_flags + 0.5);       // whole: Direct3D interpolates 4.0 as 3.9999
    float selected = mod(flags, 2.0);
    float hover = mod(floor(flags / 2.0), 2.0);
    float ghost = mod(floor(flags / 4.0), 2.0);
    vec3 colour = mix(v_color.rgb, highlight.rgb, 0.45 * selected);
    vec3 shade = colour * (0.35 + 0.65 * diffuse) + vec3(0.22 * spec);
    shade += highlight.rgb * rim * 0.6 * selected + vec3(rim * 0.6 * hover);
    float fog = clamp((-v_pos.z - params.y) / max(params.z - params.y, 1e-4), 0.0, 1.0);
    shade = mix(shade, background.rgb, fog * 0.65);
    shade = mix(shade, background.rgb, slab.z * ghost);
    fragColor = vec4(shade, v_color.a);
}
