#version 440
// chemsuite protein view: ray-cylinder hit, two-colour bonds, slab, selection and hover
// (GPL-3.0-or-later).
layout(location = 0) in vec3 v_pos;
layout(location = 1) in vec4 v_start;
layout(location = 2) in vec4 v_axis;
layout(location = 3) in vec4 v_color0;
layout(location = 4) in vec4 v_color1;
layout(location = 5) in vec2 v_flags;

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
    vec3 dir = normalize(v_pos);
    vec3 a = v_start.xyz;
    float r = v_start.w;
    vec3 ax = v_axis.xyz;
    float len = v_axis.w;
    vec3 oc = -a;
    float dd = dot(dir, ax);
    float od = dot(oc, ax);
    vec3 dp = dir - ax * dd;
    vec3 op = oc - ax * od;
    float qa = dot(dp, dp);
    float qb = 2.0 * dot(dp, op);
    float qc = dot(op, op) - r * r;
    float disc = qb * qb - 4.0 * qa * qc;
    if (qa < 1e-10 || disc < 0.0)
        discard;
    float t = (-qb - sqrt(disc)) / (2.0 * qa);
    float s = od + t * dd;
    if (s < 0.0 || s > len)
        discard;
    vec3 p = dir * t;
    if (-p.z < slab.x || -p.z > slab.y)
        discard;
    vec3 n = normalize(p - (a + ax * s));
    vec4 clip = proj * vec4(p, 1.0);
    float z = clip.z / clip.w;
    gl_FragDepth = params.x > 0.5 ? z : z * 0.5 + 0.5;
    bool first = s < 0.5 * len;
    vec4 base = first ? v_color0 : v_color1;
    float flags = first ? v_flags.x : v_flags.y;
    fragColor = vec4(look(base.rgb, n, dir, flags, -p.z), base.a);
}
