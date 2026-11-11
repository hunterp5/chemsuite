#version 440
// chemsuite protein view: ray-sphere hit, slab, selection and hover (GPL-3.0-or-later).
layout(location = 0) in vec3 v_pos;
layout(location = 1) in vec4 v_sphere;
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

vec3 look(vec3 base, vec3 n, vec3 dir, float raw, float depth)
{
    float flags = floor(raw + 0.5);           // whole: Direct3D interpolates 4.0 as 3.9999
    vec3 l = normalize(light.xyz);
    float diffuse = max(dot(n, l), 0.0);
    float spec = pow(max(dot(reflect(-l, n), -dir), 0.0), 48.0);
    float rim = pow(1.0 - max(dot(n, -dir), 0.0), 2.0);
    float selected = mod(flags, 2.0);
    float hover = mod(floor(flags / 2.0), 2.0);
    float ghost = mod(floor(flags / 4.0), 2.0);
    vec3 colour = mix(base, highlight.rgb, 0.45 * selected);
    vec3 shade = colour * (0.32 + 0.68 * diffuse) + vec3(0.28 * spec);
    shade += highlight.rgb * rim * 0.9 * selected + vec3(rim * 0.8 * hover);
    float fog = clamp((depth - params.y) / max(params.z - params.y, 1e-4), 0.0, 1.0);
    shade = mix(shade, background.rgb, fog * 0.65);
    return mix(shade, background.rgb, slab.z * ghost);
}

void main()
{
    bool ortho = params.w > 0.5;               // rays straight down -z from the eye's plane
    vec3 dir = ortho ? vec3(0.0, 0.0, -1.0) : normalize(v_pos);
    vec3 o = ortho ? vec3(v_pos.xy, 0.0) : vec3(0.0);
    vec3 c = v_sphere.xyz - o;
    float r = v_sphere.w;
    float b = dot(dir, c);
    float h = b * b - dot(c, c) + r * r;
    if (h < 0.0)
        discard;
    vec3 p = o + dir * (b - sqrt(h));
    c = v_sphere.xyz;
    if (-p.z < slab.x || -p.z > slab.y)
        discard;
    vec3 n = normalize(p - c);
    vec4 clip = proj * vec4(p, 1.0);
    float z = clip.z / clip.w;
    gl_FragDepth = params.x > 0.5 ? z : z * 0.5 + 0.5;
    fragColor = vec4(look(v_color.rgb, n, dir, v_flags, -p.z), v_color.a);
}
